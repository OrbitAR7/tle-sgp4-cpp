#!/usr/bin/env python3
"""
Quick plotting script for the SGP4 examples.

Runs the C++ examples, grabs the output tables, and generates some plots.
Just needs matplotlib - install with: python3 -m pip install matplotlib
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import List, Tuple

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT / "build"
ISS_BIN = BUILD_DIR / "iss_tracking"
CONST_BIN = BUILD_DIR / "constellation_demo"


def run_command(cmd: List[str]) -> str:
    """Run command and return stdout; raise on non-zero exit."""
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed ({proc.returncode}): {' '.join(cmd)}\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
    return proc.stdout


def parse_iss_ground_track(text: str) -> Tuple[List[float], List[float], List[float], List[float]]:
    """Pull out the ground track data from ISS output."""
    rows: List[Tuple[float, float, float, float]] = []
    grab = False
    for line in text.splitlines():
        if line.startswith("Time [min]"):
            grab = True
            continue
        if grab:
            if not line.strip():
                break
            parts = line.split()
            if len(parts) >= 4:
                try:
                    t, lat, lon, alt = map(float, parts[:4])
                except ValueError:
                    # Skip separator or malformed rows
                    continue
                rows.append((t, lat, lon, alt))
    if not rows:
        raise ValueError("No ground-track rows found in iss_tracking output.")
    t, lat, lon, alt = zip(*rows)
    return list(t), list(lat), list(lon), list(alt)


def parse_constellation_altitudes(text: str) -> Tuple[List[str], List[float], List[List[float]]]:
    """Parse the altitude table from constellation demo output."""
    names: List[str] = []
    rows: List[Tuple[float, List[float]]] = []
    for line in text.splitlines():
        if line.startswith("#ALT_NAMES"):
            parts = line.split(",")
            names = [p.replace("_", " ") for p in parts[1:] if p]
            continue
        if line.startswith("T [min]"):
            # Skip the formatted header; we already have names from #ALT_NAMES
            continue
        if names and line[:1].isdigit():
            parts = line.split()
            t = float(parts[0])
            vals = list(map(float, parts[1:]))
            rows.append((t, vals))
    if not rows:
        raise ValueError("No altitude table found in constellation_demo output.")
    times = [r[0] for r in rows]

    # Pad/truncate each row to match number of names to keep alignment
    padded_rows: List[List[float]] = []
    for _, vals in rows:
        row = list(vals)
        if len(row) < len(names):
            row.extend([float('nan')] * (len(names) - len(row)))
        elif len(row) > len(names):
            row = row[:len(names)]
        padded_rows.append(row)

    # Build series per satellite index (aligned with names)
    series: List[List[float]] = []
    for i in range(len(names)):
        series.append([row[i] for row in padded_rows])

    return names, times, series


def plot_iss(t: List[float], lat: List[float], lon: List[float], alt: List[float], outdir: Path) -> None:
    (outdir / "plots").mkdir(parents=True, exist_ok=True)
    # Ground track
    plt.figure(figsize=(6, 3))
    plt.plot(lon, lat, marker="o")
    plt.xlabel("Lon [deg]")
    plt.ylabel("Lat [deg]")
    plt.title("ISS Ground Track (one orbit)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outdir / "plots" / "iss_ground_track.png", dpi=200)
    plt.close()

    # Altitude vs time
    plt.figure(figsize=(6, 3))
    plt.plot(t, alt, marker="o")
    plt.xlabel("Time [min]")
    plt.ylabel("Alt [km]")
    plt.title("ISS Altitude (sampled one orbit)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(outdir / "plots" / "iss_altitude.png", dpi=200)
    plt.close()


def plot_constellation(names: List[str], times: List[float], series: List[List[float]], outdir: Path) -> None:
    (outdir / "plots").mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(7, 4))
    # Make Molniya visually distinct (red) if present
    for name, vals in zip(names, series):
        color = "red" if "MOLNIYA" in name.upper() else None
        plt.plot(times, vals, marker="o", label=name, color=color)
    plt.xlabel("Time [min]")
    plt.ylabel("Altitude [km]")
    plt.title("Altitude Evolution (24 hours)")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outdir / "plots" / "constellation_altitude.png", dpi=200)
    plt.close()


def main() -> int:
    if not ISS_BIN.exists() or not CONST_BIN.exists():
        print("Build outputs not found. Please run 'cmake -S . -B build && cmake --build build' first.", file=sys.stderr)
        return 1

    print("Running iss_tracking ...", flush=True)
    iss_out = run_command([str(ISS_BIN)])
    t, lat, lon, alt = parse_iss_ground_track(iss_out)
    plot_iss(t, lat, lon, alt, ROOT)
    print("Created plots/iss_ground_track.png and plots/iss_altitude.png")

    print("Running constellation_demo ...", flush=True)
    const_out = run_command([str(CONST_BIN)])
    names, times, series = parse_constellation_altitudes(const_out)
    plot_constellation(names, times, series, ROOT)
    print("Created plots/constellation_altitude.png")

    print("Done. Plots saved under", ROOT / "plots")
    return 0


if __name__ == "__main__":
    sys.exit(main())
