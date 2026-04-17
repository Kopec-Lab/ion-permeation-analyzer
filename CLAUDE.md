# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment

All work runs in the `mdanalysis` conda environment:

```bash
conda activate mdanalysis
# or, to run without activating:
conda run -n mdanalysis python ion-permeation-analyzer.py ...
```

When piping interactive stdin through `conda run`, wrap in bash to preserve the pipe:

```bash
conda run -n mdanalysis bash -c "printf '25\n7\n' | python ion-permeation-analyzer.py ..."
```

## Running the script

```bash
python ion-permeation-analyzer.py -s structure.gro -f traj.xtc -n index.ndx -o out.dat
```

The script is interactive: it prints the NDX groups and asks the user to select (by number) the phosphate group, the ion group, and optionally a protein/cylinder group (`-cyl` flag). Test files are in `TEST-TRAJ/`. The reference run uses group 25 (P) for phosphates and group 7 (K) for ions.

## Architecture

`count_ions4.py` is the original pmx-based version (kept for reference). `ion-permeation-analyzer.py` is the active script — a rewrite using MDAnalysis as the trajectory engine.

**Data flow in `main()`:**

1. `mda.Universe(strFile, trjFile)` loads structure + trajectory. Positions are in Ångströms.
2. `read_ndx()` parses the GROMACS `.ndx` file into 0-based atom indices (subtracts 1 from GROMACS 1-based indices). `select_ndx()` handles interactive group selection.
3. Trajectory loop: per frame, `get_limits(phosCrd)` divides the z-axis into 4 regions (plow/pmid/phigh from phosphate atom positions). Each ion's region assignment is appended to `ions[idx]` and its raw z-coordinate to `ion_zpos[idx]`.
4. `track_permeations()` runs a state machine over each ion's region sequence. An upward permeation requires the uninterrupted sequence 1→2→3→4; downward requires 4→3→2→1. The state machine handles PBC wrapping (4→1 and 1→4 steps).
5. Two plots are saved: `_permeations.png` (cumulative counts + linear fit) and `_positions.png` (z vs time per ion species).

**Key design decisions to be aware of:**

- `ion_resname` (keyed by 0-based atom index) is used to split ions into subplots by residue name. Non-permeating ions are plotted as gray dots; permeating ions get unique colors from `gist_rainbow` sampled as `i / n_perm` (not `i / (n_perm-1)`) to avoid the colormap's red wrap-around making the first and last ion look identical.
- The cylinder mode (`-cyl`) tracks a secondary state machine (`CYLstate0/1`) that requires the ion to also be within the cylinder radius in regions 2 and 3.
- `count_ions4.py` uses pmx atom coordinates in nm; `ion-permeation-analyzer.py` uses MDAnalysis positions in Å. The analysis functions (`get_limits`, `get_cylinder`, etc.) are unit-agnostic since boundaries are computed from the same coordinate system each frame.
