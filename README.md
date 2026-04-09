# ion permeation counter

Ion permeation counter and density / PMF analysis tool for GROMACS molecular
dynamics trajectories. Built on [MDAnalysis](https://www.mdanalysis.org/).

The tool counts how many times each ion of a chosen species crosses a lipid
bilayer (or any membrane defined by a phosphate-like reference group), and
produces ready-to-use plots and numerical output for the permeation events,
ion z-trajectories, 1D z-density / PMF profiles, and 2D XZ / YZ density / PMF
profiles useful for wider channels.

## What it does

For every frame the script:

1. Computes three z-axis reference levels from the user-selected phosphate
   group: `plow` (lower headgroups), `pmid` (bilayer center) and `phigh`
   (upper headgroups). The z-axis is divided into four regions:

   ```
   (bulk)  4
   ----  phigh  ----   upper leaflet headgroups
            3
   ----  pmid   ----   bilayer centre
            2
   ----  plow   ----   lower leaflet headgroups
   (bulk)  1
   ```

2. Assigns every ion to one of those four regions and records its (x, y, z)
   position.

3. Runs a state machine over each ion's region sequence. An **upward**
   permeation requires the uninterrupted sequence `1 → 2 → 3 → 4`, a
   **downward** permeation requires `4 → 3 → 2 → 1`. Periodic-boundary
   wrap-arounds (`4 → 1` and `1 → 4`) are handled. Partial crossings or
   reversals reset the state machine — they are *not* counted as events.

4. Optionally (`-cyl`) restricts counting to ions that also stay inside a
   cylinder defined by a protein/pore index group while traversing regions
   2 and 3.

5. From the per-ion event lists, computes ionic currents (pA), per-ion transit
   times (entry into the membrane region until completion of the crossing),
   and 1D / 2D density and PMF profiles.

## Installation

The script requires Python ≥ 3.10 and the following packages:

- MDAnalysis
- numpy
- scipy
- matplotlib

A conda environment named `mdanalysis` is the supported workflow:

```bash
conda create -n mdanalysis python=3.11
conda activate mdanalysis
pip install MDAnalysis numpy scipy matplotlib
```

## Usage

```bash
python count-ions-mdanalysis.py -s <structure> -f <trajectory> -n <index> -o <output> [-cyl]
```

| Argument | Description |
|----------|-------------|
| `-s`     | Input structure file (`.gro` / `.pdb`) |
| `-f`     | Input trajectory file (`.xtc` / `.trr`) |
| `-n`     | GROMACS index file (`.ndx`) — default `index.ndx` |
| `-o`     | Output base name for the `.dat` summary and all plots — default `out.dat` |
| `-cyl`   | Enable cylinder-restricted counting (asks for an extra index group) |

The script is **interactive**: it prints the groups found in the index file
and asks you (by number) to pick:

1. the **phosphate** group (defines the membrane reference)
2. the **ion** group (only one ion species per run — re-run the script for
   each species you want to count)
3. *(only with `-cyl`)* the **protein / pore** group used to construct the
   cylinder

### Piping the interactive answers

When running through `conda run`, wrap the pipe in `bash -c` so stdin is
preserved:

```bash
conda run -n mdanalysis bash -c \
  "printf '25\n7\n' | python count-ions-mdanalysis.py \
     -s md.gro -f md.xtc -n index.ndx -o results/ions.dat"
```

## Output

All output uses the base name from `-o` (without the `.dat` extension).
Times are reported in **nanoseconds**, distances in **Ångströms**.

### `<base>.dat`
Plain-text summary:

- Total simulation time (ns)
- Total upward / downward permeation counts
- Ionic current estimates (pA): `I = N · e / t`
- Per-ion permeation times and transit times (ns)
- The same quantities restricted to the cylinder (if `-cyl` was used)

### `<base>_permeations.png`
Cumulative permeation count vs simulation time, with a least-squares linear
fit overlaid for both directions. The slope of the fit is reported as the
permeation rate in events / ns.

### `<base>_positions.png`
Z-position of every ion vs simulation time, one subplot per ion species
(residue name).

- **Gray dots** — ions that did not complete a full membrane traversal.
  Some gray traces may dip into the headgroup / interior region and still
  remain gray: these ions entered the membrane but reversed before completing
  the full `1→2→3→4` (or `4→3→2→1`) sequence required for an event. This is
  the physically correct behaviour, not a counting error.
- **Coloured dots** — ions that permeated at least once. Each permeating ion
  gets a unique colour from an evenly-spaced rainbow palette so all traces
  remain visually distinct.
- **▼ / ▲** — markers for upward / downward permeation events.
- Horizontal lines show the time-averaged `plow / pmid / phigh` boundaries.

### `<base>_<resname>_density_pmf.png` / `.xvg`
Per ion species — 1D z-density profile and PMF.

- **Top panel** — z-density. All ion z-positions from every frame are wrapped
  to `[0, ⟨L_z⟩)` (to fold PBC artefacts), histogrammed into 1 Å bins, divided
  by `n_frames · bin_width`, and normalised by the 95th percentile of the
  Gaussian-smoothed (σ = 2 Å) density. **No artificial pinning to 1** is
  applied, so asymmetric profiles (e.g. under applied voltage) are preserved.
- **Bottom panel** — `PMF(z) = −ln(ρ(z) / ρ_bulk)` in kBT, computed from the
  smoothed density. Positive values inside the membrane, near zero in bulk.
- The `.xvg` file contains the numerical data with columns
  `z (Å)  ρ_raw  ρ_smooth  PMF (kBT)` and standard xmgrace headers.

### `<base>_<resname>_density_pmf_XZ.png` / `.dat` and `..._YZ.png` / `.dat`
Per ion species — 2D density and PMF heatmaps in the XZ and YZ planes,
useful for wider channels where lateral structure matters.

- The 2D histogram uses 1 Å square bins, is normalised by `n_frames · bin_area`,
  and rescaled by the 95th percentile of the Gaussian-smoothed 2D density
  (same logic as the 1D version, no pinning).
- Two panels per file: normalised density and PMF (kBT). Membrane reference
  lines (`plow`, `phigh`) are overlaid.
- The `.dat` file contains a flat grid table with columns
  `H_center (Å)  Z_center (Å)  ρ_raw  ρ_smooth  PMF (kBT)`,
  blank lines between H rows for compatibility with `gnuplot`'s `splot`.


## Example

```bash
conda activate mdanalysis
python count-ions-mdanalysis.py \
    -s structure.gro \
    -f trajectory.xtc \
    -n index.ndx \
    -o out.dat
# select group 25 (P) when asked for the phosphate group
# select group  7 (K) when asked for the ion group
```

This produces

```
out.dat
out_permeations.png
out_positions.png
out_K_density_pmf.png
out_K_density_pmf.xvg
out_K_density_pmf_XZ.png
out_K_density_pmf_XZ.dat
out_K_density_pmf_YZ.png
out_K_density_pmf_YZ.dat
```
