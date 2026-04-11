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
from rich.progress import ProgressBar

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

def get_cluster_status(root_dir, num_clusters):
    status = {}
    cluster_dir = os.path.join(root_dir, "CLUSTER")
    for i in range(1, int(num_clusters) + 1):
        # We assume Round 1 for now, can be extended for CI loop
        d = os.path.join(cluster_dir, f"1.{i}.cluster")
        if not os.path.exists(d):
            status[i] = "waiting"
        elif os.path.exists(os.path.join(d, "expand.DONE")):
            status[i] = "done"
        elif os.path.exists(os.path.join(d, "cmsearch.DONE")):
            status[i] = "expanding"
        elif os.path.exists(os.path.join(d, "cluster.DONE")):
            status[i] = "searching"
        else:
            status[i] = "init"
    return status

def get_sys_info():
    cpu = psutil.cpu_percent(interval=None)
    workers = 0
    for proc in psutil.process_iter(['name']):
        try:
            if proc.info['name'] in ['cmsearch', 'mlocarna']:
                workers += 1
        except (psutil.NoSuchProcess, psutil.AccessDenied, psutil.ZombieProcess):
            pass
    return cpu, workers

def calculate_eta(status, start_time):
    done_count = sum(1 for s in status.values() if s == "done")
    if done_count == 0:
        return "Calculating..."
    elapsed = time.time() - start_time
    avg_per_cluster = elapsed / done_count
    remaining = (len(status) - done_count) * avg_per_cluster
    return time.strftime("%H:%M:%S", time.gmtime(remaining))

def make_layout():
    layout = Layout()
    layout.split(
        Layout(name="header", size=5),
        Layout(name="body"),
        Layout(name="footer", size=10)
    )
    return layout

def get_header_panel(cpu, workers, round_num, eta):
    grid = Table.grid(expand=True)
    grid.add_column(justify="left")
    grid.add_column(justify="right")
    grid.add_row(
        f"[bold cyan]Pipeline:[/bold cyan] Round {round_num}",
        f"[bold yellow]CPU:[/bold yellow] {cpu}% | [bold yellow]Workers:[/bold yellow] {workers}"
    )
    grid.add_row(
        f"[bold cyan]ETA:[/bold cyan] {eta}",
        ""
    )
    return Panel(grid, title="System Status")

def get_matrix_panel(status):
    table = Table.grid(padding=1)
    # Automatically determine columns based on status size
    num_cols = 10
    for _ in range(num_cols):
        table.add_column()
    
    current_row = []
    # Sort status keys to ensure correct order
    for i in sorted(status.keys()):
        s = status[i]
        color = "white"
        if s == "done": color = "green"
        elif s == "expanding": color = "orange3"
        elif s == "searching": color = "yellow"
        elif s == "init": color = "blue"
        
        current_row.append(f"[{color}]●[/{color}]")
        if len(current_row) == num_cols:
            table.add_row(*current_row)
            current_row = []
    # Add any remaining icons
    if current_row:
        while len(current_row) < num_cols:
            current_row.append("")
        table.add_row(*current_row)
        
    return Panel(table, title=f"Cluster Matrix (1-{len(status)})")

def get_log_panel(root_dir):
    log_path = os.path.join(root_dir, "pipeline.log")
    lines = []
    if os.path.exists(log_path):
        try:
            with open(log_path, "r") as f:
                # Read last 8 lines
                lines = f.readlines()[-8:]
        except Exception:
            lines = ["Error reading pipeline.log"]
    return Panel("".join(lines), title="Recent Logs")

def main():
    parser = argparse.ArgumentParser(description="GraphClust Live Monitor")
    parser.add_argument("--root-dir", default=".", help="GraphClust root directory")
    parser.add_argument("--interval", type=int, default=2, help="Update interval in seconds")
    args = parser.parse_args()

    root = args.root_dir
    config = read_config(root)
    num_clusters = int(config.get("GLOBAL_num_clusters", 80))
    start_time = time.time()
    
    layout = make_layout()
    with Live(layout, refresh_per_second=1, screen=True):
        try:
            while True:
                cpu, workers = get_sys_info()
                status = get_cluster_status(root, num_clusters)
                eta = calculate_eta(status, start_time)
                
                # Determine current round (simple heuristic for now)
                round_num = 1
                if os.path.exists(os.path.join(root, "1.round.DONE")):
                    round_num = 2
                
                layout["header"].update(get_header_panel(cpu, workers, round_num, eta))
                layout["body"].update(get_matrix_panel(status))
                layout["footer"].update(get_log_panel(root))
                time.sleep(args.interval)
        except KeyboardInterrupt:
            pass

if __name__ == "__main__":
    main()
