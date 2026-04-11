import os
import sys
import time
import argparse
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

def main():
    parser = argparse.ArgumentParser(description="GraphClust Live Monitor")
    parser.add_argument("--root-dir", default=".", help="GraphClust root directory")
    parser.add_argument("--interval", type=int, default=2, help="Update interval in seconds")
    args = parser.parse_args()

    console = Console()
    console.print(f"[bold green]Starting GraphClust Monitor in {args.root_dir}...[/bold green]")

if __name__ == "__main__":
    main()
