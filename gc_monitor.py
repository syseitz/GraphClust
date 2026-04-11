import os
import sys
import time
import argparse
import psutil
from rich.console import Console
from rich.layout import Layout
from rich.panel import Panel

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

def main():
    parser = argparse.ArgumentParser(description="GraphClust Live Monitor")
    parser.add_argument("--root-dir", default=".", help="GraphClust root directory")
    parser.add_argument("--interval", type=int, default=2, help="Update interval in seconds")
    args = parser.parse_args()

    console = Console()
    console.print(f"[bold green]Starting GraphClust Monitor in {args.root_dir}...[/bold green]")

if __name__ == "__main__":
    main()
