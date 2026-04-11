import os
import sys
import time
import argparse
import psutil
from rich.console import Console
from rich.layout import Layout
from rich.panel import Panel
from rich.table import Table
from rich.live import Live
from rich.columns import Columns
from rich.progress_bar import ProgressBar
from rich.text import Text
from datetime import datetime

def read_config(root_dir):
    config_path = os.path.join(root_dir, "config")
    config = {}
    if not os.path.exists(config_path):
        return config
    with open(config_path, "r") as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                config[parts[0]] = parts[1]
    return config

def get_cluster_status(root_dir, num_clusters, round_idx):
    status = {}
    cluster_dir = os.path.join(root_dir, "CLUSTER")
    exists_any = False
    now = time.time()
    for i in range(1, int(num_clusters) + 1):
        d = os.path.join(cluster_dir, f"{round_idx}.{i}.cluster")
        if not os.path.exists(d):
            status[i] = "waiting"
        else:
            exists_any = True
            # Check for activity (modified in last 2 minutes)
            is_active = False
            for root_sub, dirs, files in os.walk(d):
                for f in files:
                    if f.startswith("task") and f.endswith(".started"):
                        f_path = os.path.join(root_sub, f)
                        if now - os.path.getmtime(f_path) < 120:
                            is_active = True
                            break
                if is_active: break

            if os.path.exists(os.path.join(d, "expand.DONE")):
                status[i] = "done"
            elif os.path.exists(os.path.join(d, "cmsearch.DONE")):
                status[i] = "active_searching" if is_active else "searching"
            elif os.path.exists(os.path.join(d, "model.DONE")):
                status[i] = "active_modeled" if is_active else "modeled"
            elif os.path.exists(os.path.join(d, "pp.DONE")):
                status[i] = "active_aligned" if is_active else "aligned"
            else:
                status[i] = "active_init" if is_active else "init"
    return status, exists_any

def get_sys_info():
    try:
        cpu = psutil.cpu_percent(interval=None)
    except:
        cpu = 0
    workers = 0
    active_procs = []
    try:
        for proc in psutil.process_iter(['name']):
            try:
                name = proc.info.get('name') if proc.info else None
                if name is None:
                    continue
                if name in ['cmsearch', 'mlocarna', 'cmcalibrate', 'cmbuild']:
                    workers += 1
                    try:
                        cmd = " ".join(proc.cmdline())
                        if "CLUSTER/" in cmd:
                            import re
                            m = re.search(r"CLUSTER/(\d+\.\d+)\.cluster", cmd)
                            if m:
                                active_procs.append(m.group(1))
                    except (psutil.AccessDenied, SystemError, PermissionError):
                        pass
            except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess,
                    SystemError, PermissionError):
                pass
    except (SystemError, PermissionError, OSError):
        pass
    return cpu, workers, sorted(list(set(active_procs)))

def calculate_eta(status, start_time):
    done_count = sum(1 for s in status.values() if s == "done")
    if done_count == 0:
        return "Calculating..."
    elapsed = time.time() - start_time
    avg_per_cluster = elapsed / done_count
    remaining = (len(status) - done_count) * avg_per_cluster
    return time.strftime("%H:%M:%S", time.gmtime(remaining))

def get_stage_status(root_dir):
    """Check progress of all pipeline stages."""
    result = []

    # Stage 0: Init — done when FASTA/ directory exists
    fasta_dir = os.path.join(root_dir, "FASTA")
    stage0_done = os.path.exists(fasta_dir)
    result.append({
        "label": "Stage 0", "desc": "Initialisation",
        "total": 1, "finished": 1 if stage0_done else 0,
        "errors": 0, "done": stage0_done,
    })

    # Stage 1: Preprocessing — done when FASTA/data.fasta exists
    data_fasta = os.path.join(fasta_dir, "data.fasta")
    stage1_done = os.path.exists(data_fasta)
    result.append({
        "label": "Stage 1", "desc": "Preprocessing",
        "total": 1, "finished": 1 if stage1_done else 0,
        "errors": 0, "done": stage1_done,
    })

    # Stages with sge_log-based progress tracking
    batch_stages = [
        ("Stage 2", "Graph Generation", os.path.join(root_dir, "GSPAN"), "gspan.DONE"),
        ("Stage 3", "NSPDK Vectors", os.path.join(root_dir, "SVECTOR"), "svector.groups.DONE"),
    ]
    for label, desc, stage_dir, done_file in batch_stages:
        if not os.path.exists(stage_dir):
            result.append({
                "label": label, "desc": desc,
                "total": 0, "finished": 0, "errors": 0, "done": False,
            })
            continue
        is_done = os.path.exists(os.path.join(stage_dir, done_file))
        sge_log = os.path.join(stage_dir, "sge_log")
        total_jobs = 0
        finished_jobs = 0
        error_jobs = 0
        if os.path.exists(sge_log):
            task_jobs_file = os.path.join(sge_log, "task.jobs")
            if os.path.exists(task_jobs_file):
                try:
                    with open(task_jobs_file) as f:
                        total_jobs = int(f.read().strip())
                except (ValueError, IOError):
                    pass
            if total_jobs > 0:
                for i in range(1, total_jobs + 1):
                    if os.path.exists(os.path.join(sge_log, f"task-{i}.finished")):
                        finished_jobs += 1
                    elif os.path.exists(os.path.join(sge_log, f"task-{i}.error")):
                        error_jobs += 1
        result.append({
            "label": label, "desc": desc,
            "total": total_jobs, "finished": finished_jobs,
            "errors": error_jobs, "done": is_done,
        })
    # Stage 10: Final Results — done when RESULTS/partitions/final_partition.soft exists
    results_done = os.path.exists(os.path.join(root_dir, "RESULTS", "partitions", "final_partition.soft"))
    result.append({
        "label": "Stage 10", "desc": "Final Results",
        "total": 1, "finished": 1 if results_done else 0,
        "errors": 0, "done": results_done,
    })

    return result


def get_stage_panel(stages):
    """Render a panel showing stage progress bars."""
    if not stages:
        return Panel("[dim]No stage data yet[/dim]", title="Pipeline Stages", border_style="magenta")
    table = Table.grid(expand=True, padding=(0, 1))
    table.add_column(width=28)  # label
    table.add_column(ratio=1)   # bar
    table.add_column(width=16, justify="right")  # counts
    for s in stages:
        total = s["total"]
        finished = s["finished"]
        errors = s["errors"]
        if s["done"]:
            status_str = f"[bold green]{s['label']}: {s['desc']}[/]"
            bar_text = Text("█" * 30, style="green")
            count_str = "[green]✓ DONE[/green]"
        elif total == 0:
            status_str = f"[dim]{s['label']}: {s['desc']}[/]"
            bar_text = Text("░" * 30, style="grey37")
            count_str = "[dim]waiting[/dim]"
        else:
            pct = finished / total if total > 0 else 0
            filled = int(pct * 30)
            status_str = f"[bold cyan]{s['label']}: {s['desc']}[/]"
            bar_chars = "█" * filled + "░" * (30 - filled)
            color = "yellow" if pct < 1.0 else "green"
            bar_text = Text(bar_chars, style=color)
            err_str = f" [red]({errors} err)[/red]" if errors > 0 else ""
            count_str = f"[white]{finished}/{total}[/white] ({pct:.0%}){err_str}"
        table.add_row(status_str, bar_text, count_str)
    return Panel(table, title="Pipeline Stages", border_style="magenta")


def make_layout():
    layout = Layout()
    layout.split(
        Layout(name="header", size=6),
        Layout(name="stages", size=7),
        Layout(name="main"),
        Layout(name="footer", size=10)
    )
    layout["main"].split_row(
        Layout(name="body", ratio=3),
        Layout(name="cumulative", ratio=1, minimum_size=30),
    )
    return layout

def get_header_panel(cpu, workers, current_round, eta, active_clusters):
    grid = Table.grid(expand=True)
    grid.add_column(justify="left")
    grid.add_column(justify="right")
    
    # Legend string with larger dots
    legend = (
        "[grey37]●[/] Waiting "
        "[blue]●[/] Init "
        "[cyan]●[/] Aligned "
        "[yellow]●[/] Modeled "
        "[orange3]●[/] Searching "
        "[green]●[/] Done "
        "[bold white blink]◎[/] ACTIVE"
    )

    grid.add_row(
        f"[bold cyan]Pipeline:[/bold cyan] Round {current_round} [dim]|[/dim] {legend}",
        f"[bold yellow]CPU:[/bold yellow] {cpu}% | [bold yellow]Workers:[/bold yellow] {workers}"
    )
    grid.add_row(
        f"[bold cyan]ETA (Curr Round):[/bold cyan] {eta}",
        f"[bold magenta]Active Clusters:[/bold magenta] {', '.join(active_clusters) if active_clusters else 'None'}"
    )
    return Panel(grid, title="System Status & Legend", border_style="bright_blue")

def get_round_stage5_info(root_dir, round_idx):
    """Get Stage 5 NSPDK Clustering progress for a specific round."""
    import glob as globmod
    import re
    stage5_out = os.path.join(root_dir, "SGE_LOG", f"stage5-1.out.{round_idx}")
    if not os.path.exists(stage5_out):
        # try without round suffix (Round 1 format)
        stage5_out = os.path.join(root_dir, "SGE_LOG", "stage5-1.out")
    done_file = os.path.join(root_dir, "SVECTOR", f"data.svector.fast_cluster.{round_idx}.DONE")
    is_done = os.path.exists(done_file)
    total = 0
    finished = 0
    if os.path.exists(stage5_out):
        try:
            with open(stage5_out, "r") as f:
                content = f.read()
            m = re.search(r"(\d+)\s+instances", content)
            if m:
                total = int(m.group(1))
            last_k = list(re.finditer(r"(\d+)K", content))
            if last_k:
                finished = int(last_k[-1].group(1)) * 1000
                finished = min(finished, total)
            if is_done:
                finished = total
        except (IOError, OSError):
            pass
    return {"total": total, "finished": finished, "done": is_done}

def _add_stage_row(table, label, desc, info, use_k=False):
    """Add a stage progress row to a table."""
    t, f_val, done = info["total"], info["finished"], info["done"]
    if done:
        table.add_row(f"[bold green]{label}: {desc}[/]", Text("█" * 30, style="green"), "[green]✓ DONE[/green]")
    elif t == 0:
        table.add_row(f"[dim]{label}: {desc}[/]", Text("░" * 30, style="grey37"), "[dim]waiting[/dim]")
    else:
        pct = f_val / t if t > 0 else 0
        filled = int(pct * 30)
        bar = Text("█" * filled + "░" * (30 - filled), style="yellow")
        if use_k and t >= 1000:
            c_str = f"[white]{f_val//1000}K/{t//1000}K[/white] ({pct:.0%})"
        else:
            c_str = f"[white]{f_val}/{t}[/white] ({pct:.0%})"
        table.add_row(f"[bold cyan]{label}: {desc}[/]", bar, c_str)

def get_matrix_panel(status, round_idx, root_dir=None):
    round_summary = None
    # Side-by-side layout: grid left, stages right
    side_by_side = Table.grid(expand=True)
    side_by_side.add_column(width=32)  # grid
    side_by_side.add_column(ratio=1)   # stages

    # Cluster grid (left side)
    grid = Table.grid(padding=(0, 1))
    num_cols = 10
    for _ in range(num_cols):
        grid.add_column(justify="center")

    current_row = []
    for i in sorted(status.keys()):
        s = status[i]
        is_active = s.startswith("active_")
        base_status = s.replace("active_", "")

        color = "grey37"
        if base_status == "done": color = "green"
        elif base_status == "searching": color = "orange3"
        elif base_status == "modeled": color = "yellow"
        elif base_status == "aligned": color = "cyan"
        elif base_status == "init": color = "blue"

        char = "●" if not is_active else "◎"
        style = f"bold {color}"
        if is_active: style += " blink"

        current_row.append(f"[{style}]{char}[/]")

        if len(current_row) == num_cols:
            grid.add_row(*current_row)
            current_row = []

    if current_row:
        current_row.extend([""] * (num_cols - len(current_row)))
        grid.add_row(*current_row)

    # Stages (right side)
    stages_table = Table.grid(expand=True, padding=(0, 1))
    stages_table.add_column(width=24)
    stages_table.add_column(ratio=1)
    stages_table.add_column(width=14, justify="right")

    if root_dir:
        # Stage 5: NSPDK Clustering
        stage5 = get_round_stage5_info(root_dir, round_idx)
        _add_stage_row(stages_table, "Stage 5", "NSPDK", stage5, use_k=True)

        # Stage 6: LocARNA Align
        total_clusters = len(status)
        s6_total = total_clusters if stage5["done"] else 0
        s6_done = sum(1 for s in status.values() if s.replace("active_", "") in ("aligned", "modeled", "searching", "done"))
        _add_stage_row(stages_table, "Stage 6", "Align",
                       {"total": s6_total, "finished": s6_done, "done": s6_done == s6_total and s6_total > 0})

        # Stage 7: Build Model
        s7_total = s6_total
        s7_done = sum(1 for s in status.values() if s.replace("active_", "") in ("modeled", "searching", "done"))
        _add_stage_row(stages_table, "Stage 7", "Model",
                       {"total": s7_total, "finished": s7_done, "done": s7_done == s7_total and s7_total > 0})

        # Stage 8: CM Search
        s8_total = s7_total
        s8_done = sum(1 for s in status.values() if s.replace("active_", "") == "done")
        _add_stage_row(stages_table, "Stage 8", "CM Search",
                       {"total": s8_total, "finished": s8_done, "done": s8_done == s8_total and s8_total > 0})

        # Stage 9: Collect Results
        round_done = os.path.exists(os.path.join(root_dir, f"{round_idx}.round.DONE"))
        s9_all_done = s8_done == s8_total and s8_total > 0
        _add_stage_row(stages_table, "Stage 9", "Results",
                       {"total": 1 if s9_all_done else 0, "finished": 1 if round_done else 0, "done": round_done})

        # Round summary info — only this round's contribution
        round_summary = None
        if round_done:
            summary_parts = []

            # Calculate this round's clusters/hits by subtracting previous round
            merged_file = os.path.join(root_dir, "EVAL", "partitions", f"{round_idx}.hard.merged")
            prev_merged = os.path.join(root_dir, "EVAL", "partitions", f"{round_idx - 1}.hard.merged") if round_idx > 1 else None

            if os.path.exists(merged_file):
                try:
                    curr_clusters = set()
                    curr_hits = 0
                    with open(merged_file) as f:
                        for line in f:
                            curr_hits += 1
                            parts = line.strip().split()
                            if parts:
                                curr_clusters.add(parts[-1])

                    prev_clusters = set()
                    prev_hits = 0
                    if prev_merged and os.path.exists(prev_merged):
                        with open(prev_merged) as f:
                            for line in f:
                                prev_hits += 1
                                parts = line.strip().split()
                                if parts:
                                    prev_clusters.add(parts[-1])

                    new_clusters = len(curr_clusters) - len(prev_clusters)
                    new_hits = curr_hits - prev_hits
                    summary_parts.append(f"[green]{new_clusters}[/] clusters")
                    summary_parts.append(f"[green]{new_hits}[/] hits")
                except (IOError, OSError):
                    pass

            # Errors (clusters without cmsearch.DONE)
            n_errors = sum(1 for s in status.values() if s.replace("active_", "") != "done" and s != "waiting")
            if n_errors > 0:
                summary_parts.append(f"[red]{n_errors}[/] errors")

            if summary_parts:
                round_summary = " [dim]·[/] ".join(summary_parts)

    side_by_side.add_row(grid, stages_table)

    outer = Table.grid(expand=True)
    outer.add_column()
    outer.add_row(side_by_side)
    if round_summary:
        outer.add_row(Text.from_markup(f"  {round_summary}"))

    return Panel(outer, title=f"Round {round_idx} Progress", border_style="cyan", padding=(0, 1))

def get_cumulative_panel(root_dir, max_rounds):
    """Show cumulative results across all completed rounds."""
    # Find the latest completed round's merged partition
    latest_merged = None
    latest_round = 0
    for r in range(max_rounds, 0, -1):
        f = os.path.join(root_dir, "EVAL", "partitions", f"{r}.hard.merged")
        if os.path.exists(f):
            latest_merged = f
            latest_round = r
            break

    if not latest_merged:
        return Panel("[dim]No results yet[/dim]", title="Cumulative Results", border_style="green", padding=(0, 1))

    try:
        clusters = set()
        hits = 0
        with open(latest_merged) as f:
            for line in f:
                hits += 1
                parts = line.strip().split()
                if parts:
                    clusters.add(parts[-1])

        # Blacklist for info
        bl_file = os.path.join(root_dir, "SVECTOR", f"data.svector.blacklist.{latest_round + 1}")
        bl_size = 0
        if os.path.exists(bl_file):
            with open(bl_file) as f:
                bl_size = sum(1 for _ in f)

        # Total sequences
        data_fasta = os.path.join(root_dir, "FASTA", "data.fasta")
        total_seqs = 0
        if os.path.exists(data_fasta):
            with open(data_fasta) as f:
                total_seqs = sum(1 for line in f if line.startswith(">"))

        info = Table.grid(expand=True, padding=(0, 1))
        info.add_column(justify="right", width=8)
        info.add_column()

        info.add_row(f"[bold green]{len(clusters)}[/]", "clusters")
        info.add_row(f"[bold green]{hits}[/]", "hits")
        info.add_row(f"[bold yellow]{bl_size}[/]", "blacklisted")
        remaining = total_seqs - bl_size
        info.add_row(f"[dim]{remaining}[/]", f"[dim]/ {total_seqs} remaining[/]")

        return Panel(info, title=f"Total (R1-{latest_round})", border_style="green", padding=(0, 1))

    except (IOError, OSError):
        return Panel("[dim]Error reading results[/dim]", title="Cumulative Results", border_style="green", padding=(0, 1))

def get_log_panel(root_dir):
    log_path = os.path.join(root_dir, "pipeline.log")
    lines = []
    if os.path.exists(log_path):
        try:
            with open(log_path, "r") as f:
                lines = f.readlines()[-8:]
        except Exception:
            lines = ["Error reading pipeline.log"]
    return Panel("".join(lines), title="Recent Logs", border_style="dim")

def main():
    parser = argparse.ArgumentParser(description="GraphClust Live Monitor")
    parser.add_argument("--root-dir", default=".", help="GraphClust root directory")
    parser.add_argument("--interval", type=int, default=2, help="Update interval in seconds")
    args = parser.parse_args()

    root = args.root_dir
    config = read_config(root)
    num_clusters = int(config.get("GLOBAL_num_clusters", 80))
    max_rounds = int(config.get("GLOBAL_iterations", 3))
    
    # We track start time of the CURRENT round for better ETA
    last_round_seen = 0
    round_start_time = time.time()
    
    layout = make_layout()
    with Live(layout, refresh_per_second=1, screen=True):
        try:
            while True:
                cpu, workers, active_clusters = get_sys_info()
                
                round_panels = []
                current_round = 1
                current_status = {}
                
                for r in range(1, max_rounds + 1):
                    status, exists = get_cluster_status(root, num_clusters, r)
                    if exists:
                        round_panels.append(get_matrix_panel(status, r, root))
                        current_round = r
                        current_status = status
                
                # Update round start time if we just entered a new round
                if current_round > last_round_seen:
                    round_start_time = time.time()
                    last_round_seen = current_round
                
                eta = calculate_eta(current_status, round_start_time)
                
                stages = get_stage_status(root)

                layout["header"].update(get_header_panel(cpu, workers, current_round, eta, active_clusters))
                layout["stages"].update(get_stage_panel(stages))
                layout["body"].update(Columns(round_panels, equal=True, expand=True))
                layout["cumulative"].update(get_cumulative_panel(root, max_rounds))
                layout["footer"].update(get_log_panel(root))
                
                time.sleep(args.interval)
        except KeyboardInterrupt:
            pass

if __name__ == "__main__":
    main()
