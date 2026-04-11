import os
import sys
import time
import argparse
from rich.console import Console
from rich.layout import Layout
from rich.panel import Panel

def main():
    parser = argparse.ArgumentParser(description="GraphClust Live Monitor")
    parser.add_argument("--root-dir", default=".", help="GraphClust root directory")
    parser.add_argument("--interval", type=int, default=2, help="Update interval in seconds")
    args = parser.parse_args()

    console = Console()
    console.print(f"[bold green]Starting GraphClust Monitor in {args.root_dir}...[/bold green]")

if __name__ == "__main__":
    main()
