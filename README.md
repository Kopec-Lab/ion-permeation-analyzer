# ion permeation analyzer

Ion permeation analyzer and density / PMF analysis tool for GROMACS molecular
dynamics trajectories. Built on [MDAnalysis](https://www.mdanalysis.org/).

The tool counts how many times each ion of a chosen species crosses a lipid
bilayer (or any membrane defined by a phosphate-like reference group), and
produces ready-to-use plots and numerical output for the permeation events,
ion z-trajectories, 1D / 2D density and PMF profiles, and automatic
binding-site detection from the free energy landscape. Based on the earlier code written by Vytautas Gapsys.

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
   and 1D / 2D density and PMF profiles with automatic binding-site detection.

## Installation

The script requires Python >= 3.10 and the following packages:

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
python ion-permeation-analyzer.py -s <structure> -f <trajectory> -n <index> -o <output> [-cyl] [-coord]
```

| Argument | Description |
|----------|-------------|
| `-s`     | Input structure file (`.gro` / `.pdb`) |
| `-f`     | Input trajectory file (`.xtc` / `.trr`) |
| `-n`     | GROMACS index file (`.ndx`) — default `index.ndx` |
| `-o`     | Output base name for the `.dat` summary and all plots — default `out.dat` |
| `-cyl`   | Enable cylinder-restricted counting and pore PMF (asks for an extra index group) |
| `-coord` | Compute coordination analysis for permeating ions (slow — requires second trajectory pass) |

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
  "printf '19\n13\n1\n' | python ion-permeation-analyzer.py \
     -s md.pdb -f md.xtc -n index.ndx -o out.dat -cyl"
```

### When to use `-cyl`

Use `-cyl` whenever you have a protein-defined pore and want to:
- Count only permeation events that pass through the pore (not lipid defects).
- Resolve **binding sites** inside the pore. Without `-cyl`, bulk ions (~200)
  swamp the ~4 filter ions in the histogram, and binding sites inside narrow
  channels (e.g. K+ channel selectivity filter) cannot be resolved. With
  `-cyl`, only ion positions inside the cylinder are histogrammed, giving
  clean pore-only PMF profiles where individual sites become visible.

## Output

All output uses the base name from `-o` (without the `.dat` extension).
Times are reported in **nanoseconds**, distances in **Angstroms**.

### Permeation counting

#### `<base>.dat`
Plain-text summary:

- Total simulation time (ns)
- Total upward / downward permeation counts
- Ionic current estimates (pA): `I = N * e / t`
- Per-ion permeation times and transit times (ns)
- The same quantities restricted to the cylinder (if `-cyl` was used)

#### `<base>_permeations.png`
Cumulative permeation count vs simulation time, with a least-squares linear
fit overlaid for both directions. The slope of the fit is reported as the
permeation rate in events / ns.

#### `<base>_positions.png`
Z-position of every ion vs simulation time, one subplot per ion species
(residue name).

- **Gray dots** — ions that did not complete a full membrane traversal.
  Some gray traces may dip into the headgroup / interior region and still
  remain gray: these ions entered the membrane but reversed before completing
  the full `1->2->3->4` (or `4->3->2->1`) sequence required for an event.
  This is the physically correct behaviour, not a counting error.
- **Coloured dots** — ions that permeated at least once. Each permeating ion
  gets a unique colour from an evenly-spaced rainbow palette so all traces
  remain visually distinct.
- **down-triangle / up-triangle** — markers for upward / downward permeation
  events.
- Horizontal lines show the time-averaged `plow / pmid / phigh` boundaries.

### Density and PMF profiles

The script produces density and PMF output at three levels of detail,
from coarse (full box) to fine (pore-only). All three run automatically;
the pore mode additionally requires `-cyl`.

#### Global: `<base>_<resname>_density_pmf.png` / `.xvg`
Full-box 1D z-density and PMF normalised to bulk.

- **Top panel** — z-density. All ion z-positions from every frame are wrapped
  to `[0, <L_z>)` (to fold PBC artefacts), histogrammed into 1 A bins, and
  normalised by the 95th percentile of the Gaussian-smoothed (sigma = 2 A)
  density. **No artificial pinning** is applied, so asymmetric profiles
  (e.g. under applied voltage) are preserved.
- **Bottom panel** — `PMF(z) = -ln(rho / rho_bulk)` in kBT.
- The `.xvg` file contains columns `z  rho_raw  rho_smooth  PMF`.
- Binding-site detection: min depth 0.5 kBT, search window `[plow-5, phigh+5]`.
  Results in `<base>_<resname>_binding_sites.dat`.

#### Local: `<base>_<resname>_density_pmf_local.png` / `.xvg`
Zoomed to the membrane region `[plow-10, phigh+10]` with lighter smoothing
and detrending to expose local features.

- Uses sigma = 1 A smoothing (preserves finer structure than the global mode).
- Density normalised to the local maximum within the window (not bulk).
- **Detrended PMF**: the PMF is decomposed into a heavy baseline
  (sigma = 10 A Gaussian) capturing the large-scale barrier shape, and a
  residual that exposes local wells. A third panel shows the detrended PMF
  with detected binding sites annotated.
- Binding-site detection: min depth 0.2 kBT on the detrended profile, search
  restricted to `[plow, phigh]` to avoid edge artefacts.
  Results in `<base>_<resname>_binding_sites_local.dat`.
- The `.xvg` file contains columns `z  rho_raw  rho_smooth  PMF  baseline  detrended`.

#### Pore (requires `-cyl`): `<base>_<resname>_density_pmf_pore.png` / `.xvg`
Cylinder-filtered density and PMF using **only** ion positions that were
inside the pore at that frame. This is the most sensitive mode for resolving
binding sites in narrow channels.

- Uses 0.5 A bins (finer than the other modes) to resolve sites spaced ~3.5 A
  apart (e.g. K+ channel selectivity filter S1-S4).
- sigma = 0.5 A smoothing.
- Detrended with a sigma = 5 A baseline.
- Binding-site detection: min depth 0.1 kBT, contour 0.15 kBT.
  Results in `<base>_<resname>_binding_sites_pore.dat`.
- The `.xvg` file contains columns `z  count  density_smooth  PMF  baseline  detrended`.

**Why three levels?** The global view shows the overall barrier shape and
bulk reference. The local view reveals membrane-interior structure without
requiring a pore definition. The pore view gives the cleanest signal for
channels but needs the `-cyl` flag and a suitable protein index group.

#### Baseline, detrending, and the sigma parameter

Inside the membrane, the PMF has two components superimposed:

1. A **large-scale barrier** — the overall energetic cost of moving an ion
   from bulk water into the membrane interior. This is a smooth, broad shape
   spanning ~40 A, typically 2-5 kBT high.
2. **Local wells** — individual binding sites (e.g. selectivity filter sites
   in a K+ channel), spaced ~3-6 A apart, with depths of only 0.2-0.7 kBT.

When you plot the raw PMF, the local wells appear as tiny ripples on a
massive slope — invisible to both the eye and to peak-finding algorithms.

The **baseline** is a heavily smoothed version of the PMF (Gaussian filter
with a large sigma). It captures only the broad barrier shape and is blind
to features narrower than ~2 x sigma. In the plots, the baseline is shown
as a dashed gray line alongside the raw PMF, so you can visually see what
the detrending subtracts out.

The **detrended PMF** = PMF - baseline. It removes the large-scale slope
and leaves only the local landscape: wells become clearly visible minima
around zero, barriers become positive peaks. This is what binding-site
detection runs on.

**Sigma** is the width of the Gaussian smoothing kernel in Angstroms. It
controls the scale separation — features narrower than ~2*sigma end up in
the detrended curve, features broader than ~2*sigma stay in the baseline:

| Mode  | Fine smoothing | Baseline sigma | Rationale |
|-------|---------------|----------------|-----------|
| Local | 1 A           | 10 A           | The membrane is ~30-40 A across. sigma=10 A captures the overall barrier shape while anything at the ~3-6 A binding-site scale passes through to the detrended profile. |
| Pore  | 0.5 A         | 5 A            | The pore data is already spatially restricted (~15-30 A), so a smaller baseline sigma is appropriate. sigma=5 A still captures the vestibule-to-filter slope but preserves ~3.5 A selectivity filter features. |

The choice of sigma is a trade-off: too small and the baseline follows the
wells (detrending removes what you are looking for), too large and the
baseline is too flat to remove the slope effectively. The current values
were tuned on K+ channel trajectories where the known S1-S4 selectivity
filter sites at ~3.5 A spacing are correctly resolved.

In short: **baseline** = what the PMF would look like without binding sites;
**detrended** = what binding sites look like without the overall barrier.

#### Why detrending is needed (comparison with specialised tools)

Specialised tools for specific channel types (e.g. xtck for K+ channels)
can avoid the need for detrending entirely. They do this by selecting only
ions that are already inside the selectivity filter — a ~15 A segment —
and reporting positions relative to a known binding site (e.g. S4). In
that narrow window, the density comes almost entirely from ions hopping
between adjacent filter sites. There is no vestibule-to-filter gradient,
so the PMF is nearly flat with clear wells and no detrending is needed.

This tool takes a different, general-purpose approach. The cylinder
(`-cyl`) is defined by the entire protein, encompassing the full
cross-section and the full z-range between plow and phigh (~40 A). The
histogram therefore includes a mix of tightly-bound filter ions and
loosely-associated vestibule or cavity ions. The ion density naturally
drops as you move deeper into the membrane, creating a large-scale
energetic slope (typically 2-5 kBT) on top of which the local binding
wells (0.2-0.7 kBT) are superimposed. Without detrending, these wells
are invisible as tiny ripples on a massive gradient.

In principle, defining a much tighter cylinder (e.g. 5 A radius,
restricted to just the filter z-range) would reproduce the specialised-
tool result without detrending. But that requires knowing where the filter
is ahead of time, which defeats the purpose of an automated analysis.
The broad-cylinder + detrending approach used here works for any channel
geometry — narrow selectivity filters, wide pores, ligand-gated channels
— without prior structural knowledge of the binding-site locations.

### Binding-site detection

Binding sites are identified as local minima of the PMF profile using
`scipy.signal.find_peaks`. Each site is characterised by:

| Field | Description |
|-------|-------------|
| `z_min` | Z-position of the PMF minimum |
| `z_lo`, `z_hi` | Extent of the site (contour above the minimum) |
| `depth_kBT` | Prominence of the minimum in kBT |
| `occ_frac` | Fraction of frames with at least one ion in `[z_lo, z_hi]` |
| `mean_count` | Average number of ions per frame in the site |
| `nearby_residues` | Protein residues near the site (pore mode only, see below) |

Sites are annotated on the PMF panels (orange shading + labels) and written
to `_binding_sites.dat` / `_binding_sites_local.dat` / `_binding_sites_pore.dat`.

#### Nearby protein residues (pore mode)

When `-cyl` is used, the pore binding-sites output additionally identifies
which protein residues line each site. For each site, protein atoms from
the cylinder group are selected if their z-position falls within the site
extent (+/- 2 A margin) and their xy-distance from the pore axis is within
10 A. The unique residues (name + number) are listed in the
`nearby_residues` column of `_binding_sites_pore.dat` and annotated on the
detrended PMF panel (up to 6 residues shown in the plot; full list in the
data file).

This is useful for direct comparison with structurally proposed binding
sites — e.g. verifying that the detected selectivity filter sites in a K+
channel correspond to the expected threonine/glycine-rich signature
residues.

#### Binding-site visualisation PDB (pore mode)

When `-cyl` is used and binding sites are detected, the script writes a PDB
file `<base>_<resname>_binding_sites_pore.pdb` containing:

- All protein atoms from the cylinder group (ATOM records).
- One dummy atom per binding site (HETATM records, residue name `DUM`,
  chain `Z`), placed at the pore axis center (x, y) and the site's
  z-position. The B-factor column stores the site depth in kBT, which
  can be used for colour-coding.

This allows quick visualisation of binding-site locations relative to the
protein structure in VMD, PyMOL, or ChimeraX. To select the site markers:
- VMD: `resname DUM`
- PyMOL: `select resn DUM`
- ChimeraX: `select :DUM`

### 2D density and PMF

#### `<base>_<resname>_density_pmf_XZ.png` / `.dat` and `..._YZ.png` / `.dat`
Per ion species — 2D density and PMF heatmaps in the XZ and YZ planes,
useful for wider channels where lateral structure matters.

- 1 A square bins, normalised by `n_frames * bin_area`, rescaled by the 95th
  percentile of the Gaussian-smoothed 2D density (no pinning).
- Two panels per file: normalised density and PMF (kBT). Membrane reference
  lines (`plow`, `phigh`) are overlaid.
- The `.dat` file contains a flat grid table with columns
  `H_center  Z_center  rho_raw  rho_smooth  PMF`,
  blank lines between H rows for compatibility with gnuplot's `splot`.

### Pore occupancy and ion-ion distances (requires `-cyl`)

#### `<base>_<resname>_pore_occupancy.png` / `.dat`
Number of ions inside the pore (cylinder AND membrane z-range `[plow, phigh]`)
at each frame. Two panels:

- **Left** — occupancy time series with a running average overlaid.
- **Right** — histogram of occupancy states (fraction of frames with 0, 1,
  2, ... ions in the pore), with the mean occupancy marked.

This immediately reveals the dominant occupancy state of the pore. For
example, a K+ channel selectivity filter typically has 2-3 K+ ions at all
times (plus additional ions in the vestibules within the cylinder), while a
narrow single-ion channel would show predominantly 0-1 ions.

The `.dat` file contains a header with the occupancy distribution and mean,
followed by per-frame columns `time_ns  n_ions`.

#### `<base>_<resname>_ion_ion_distances.png` / `.dat`
Histogram of all pairwise distances between ions simultaneously inside the
pore. Only frames with >= 2 ions contribute.

- Sharp peaks at short distances (e.g. ~3.5 A and ~7 A) indicate ions
  occupying adjacent and next-nearest binding sites — the signature of
  multi-ion knock-on permeation as seen in K+ channel selectivity filters.
- A smooth, broad distribution peaking at large distances indicates ions
  permeate largely independently of each other.

The `.dat` file contains a header with summary statistics (total pairs,
median, mean), followed by the histogram as columns
`bin_center_A  count  probability_density`.

### Coordination analysis (requires `-coord`)

#### `<base>_<resname>_coordination.png` / `.xvg`
Per-ion coordination number vs z-position for all permeating ions, averaged
and plotted as three curves:

- **Water coordination** — number of water oxygen atoms within the first
  coordination shell (default cutoff: 3.5 A).
- **Protein coordination** — number of protein heavy atoms within the same
  cutoff. Only computed when `-cyl` is also used (protein group required);
  otherwise this curve is zero.
- **Total coordination** — sum of water + protein.

This reveals how the solvation shell changes during permeation: whether ions
cross the membrane fully hydrated, partially desolvated, or with water
replaced by protein contacts (e.g. carbonyl oxygens in a K+ channel
selectivity filter).

The analysis requires a **second pass** through the entire trajectory, which
can be slow for long simulations. It is therefore gated behind the `-coord`
flag and not run by default.

- The `.xvg` file contains columns `z  water_coord  protein_coord  total_coord`.
- The plot shows the three coordination curves with membrane reference lines
  (`plow`, `pmid`, `phigh`) overlaid.

> **Note — planned development:** The coordination analysis is functional but
> will be extended in future versions. Planned improvements include:
> configurable cutoff radius, detection of coordinating atom types and their
> physicochemical properties (e.g. distinguishing backbone carbonyls from
> side-chain hydroxyls), and per-site coordination breakdown aligned with
> the detected binding sites.

## Example

### Basic run (no cylinder)

```bash
conda activate mdanalysis
python ion-permeation-analyzer.py \
    -s structure.gro \
    -f trajectory.xtc \
    -n index.ndx \
    -o out.dat
# select phosphate group, then ion group
```

Output files:
```
out.dat                          # permeation summary
out_permeations.png              # cumulative count plot
out_positions.png                # z vs time plot
out_K_density_pmf.png / .xvg    # global density + PMF
out_K_density_pmf_local.png / .xvg  # local density + PMF (membrane zoom)
out_K_binding_sites.dat          # binding sites (global)
out_K_binding_sites_local.dat    # binding sites (local)
out_K_density_pmf_XZ.png / .dat  # 2D XZ density + PMF
out_K_density_pmf_YZ.png / .dat  # 2D YZ density + PMF
```

### With cylinder (recommended for channels with a protein pore)

```bash
python ion-permeation-analyzer.py \
    -s channel.pdb -f channel.xtc -n index.ndx -o out.dat -cyl
# select phosphate group, ion group, then protein group
```

Additional output:
```
out_POT_density_pmf_pore.png / .xvg    # pore-only density + PMF
out_POT_binding_sites_pore.dat          # binding sites (pore, finest resolution)
out_POT_binding_sites_pore.pdb          # protein + dummy atoms at site centers
out_POT_pore_occupancy.png / .dat      # pore ion occupancy vs time + histogram
out_POT_ion_ion_distances.png / .dat   # pairwise ion-ion distance histogram
```

### With coordination analysis

```bash
python ion-permeation-analyzer.py \
    -s channel.pdb -f channel.xtc -n index.ndx -o out.dat -cyl -coord
# select phosphate group, ion group, then protein group
```

Additional output:
```
out_POT_coordination.png / .xvg        # coordination number vs z for permeating ions
```
