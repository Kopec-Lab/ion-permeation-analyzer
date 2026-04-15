
__doc__ = """
ion-permeation-counter.py  —  Ion permeation counter for GROMACS MD trajectories
=================================================================================

Usage
-----
    python ion-permeation-counter.py -s <structure> -f <trajectory> -n <index> -o <output> [-cyl] [-coord]

Required arguments
------------------
    -s <file>   Input structure file (.gro or .pdb).  Used to initialise the
                MDAnalysis Universe and read atom/residue metadata.
    -f <file>   Input trajectory file (.xtc or .trr).
    -n <file>   GROMACS index file (.ndx).  The script will interactively ask
                you to select three groups from this file:
                  1. Phosphate group  — atoms used to define the membrane
                     boundaries (plow, pmid, phigh) along the z-axis.
                  2. Ion group        — atoms whose permeation events are counted.
                  3. Protein/cylinder group (only with -cyl) — atoms used to
                     define the cylinder through which permeation is tracked.
    -o <file>   Output data file (default: out.dat).  Base name is also used
                for the output plots (see below).

Optional arguments
------------------
    -cyl        Enable cylinder-restricted counting.  Ions must pass through a
                cylinder (defined by a protein/pore group) to be counted as
                permeation events.  The cylinder axis is parallel to z.
    -coord      Compute coordination numbers (water + protein) for permeating
                ions vs z-position.  Requires a second trajectory pass (slow).

How it works
------------
Each frame the z-axis is divided into four regions using the mean positions of
the phosphate atoms:

    (bulk)  4
    ----  phigh  ----   (upper leaflet headgroups)
            3
    ----  pmid   ----   (bilayer centre)
            2
    ----  plow   ----   (lower leaflet headgroups)
    (bulk)  1

A permeation event (upward) is registered when an ion traverses the full
sequence 1→2→3→4 without reversing.  A downward event requires 4→3→2→1.
Partial crossings and backward steps reset the state machine.

With -cyl, only ions that additionally pass through the cylinder interior
(xy-plane distance from centre ≤ cylinder radius) in regions 2 and 3 are
counted as pore permeation events.

Output
------
  <output>.dat        Plain-text summary containing:
                        • Total simulation time (ns)
                        • Total permeation counts (up / down)
                        • Ionic current estimates (pA) = N × e / t
                        • Per-ion permeation times and transit times (ns)
                        • Same quantities restricted to the cylinder (if -cyl)

  <output>_permeations.png
                      Cumulative permeation count vs simulation time for
                      upward and downward events, with a least-squares linear
                      fit overlaid.  The slope of the fit gives the permeation
                      rate in events/ns.

  <output>_positions.png
                      Z-position of every ion vs simulation time, one subplot
                      per ion species (residue name).

                      Colour coding:
                        • Gray dots  — ions that did NOT complete a full
                          membrane traversal.  Note that some gray traces may
                          dip into the headgroup/interior region (between plow
                          and phigh) and still remain gray: these ions entered
                          the membrane but reversed before completing the full
                          1→2→3→4 (or 4→3→2→1) sequence required to register
                          a permeation event.  This is physically correct
                          behaviour and is not a counting error.
                        • Colored dots — ions that permeated at least once.
                          Each permeating ion gets a unique color drawn from an
                          evenly-spaced rainbow palette so all traces are
                          visually distinct regardless of how many there are.

                      Event markers on colored traces:
                        ▼ downward triangle — upward permeation event
                          (ion completed 1→2→3→4)
                        ▲ upward triangle   — downward permeation event
                          (ion completed 4→3→2→1)

                      Horizontal lines show the time-averaged plow / pmid /
                      phigh boundaries derived from the phosphate group.

  <output>_<resname>_density_pmf.png / .xvg
                      Per ion species — two-panel figure and plain-text data.
                      The .xvg columns are: z (Å), raw density, smoothed density,
                      PMF (kBT).  Can be plotted directly with xmgrace or numpy.

                      Figure: two-panel figure per ion species:

                      Top panel — z-density profile.  All ion z-positions
                      from every frame are histogrammed into 1-Å bins and
                      normalised so that the mean value in the bulk regions
                      (z < plow and z > phigh) equals 1.  Both the raw
                      histogram and a lightly smoothed curve (Gaussian,
                      σ = 2 Å) are shown.

                      Bottom panel — potential of mean force (PMF).
                      PMF(z) = −ln(ρ(z)/ρ_bulk) in units of kBT, computed
                      from the smoothed density.  The PMF is 0 by definition
                      in the bulk and positive inside the membrane where ions
                      are excluded.  Negative values indicate ion accumulation
                      relative to bulk (e.g. at the headgroup region).

Dependencies
------------
    MDAnalysis, numpy, matplotlib
    Install:  pip install MDAnalysis  (numpy and matplotlib are pulled in automatically)

Example
-------
    python ion-permeation-counter.py -s md.gro -f md.xtc -n index.ndx -o results/ions.dat
"""

import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import MDAnalysis as mda


def track_permeations(ions, trjTimes, bCyl=False, ionsCyl=None):
    if ionsCyl is None:
        ionsCyl = {}
    transitionsUp = {}
    transitionsDown = {}
    transitionsUpTimes = {}
    transitionsDownTimes = {}
    transitionsUpTransit = {}    # transit times (ns) for each upward permeation
    transitionsDownTransit = {}  # transit times (ns) for each downward permeation
    totalUp = 0
    totalDown = 0
    # counting within cylinder (optional)
    CYLtransitionsUp = {}
    CYLtransitionsDown = {}
    CYLtransitionsUpTimes = {}
    CYLtransitionsDownTimes = {}
    CYLtotalUp = 0
    CYLtotalDown = 0

    for key in ions.keys():
        transitionsUp[key] = 0
        transitionsDown[key] = 0
        transitionsUpTimes[key] = []
        transitionsDownTimes[key] = []
        transitionsUpTransit[key] = []
        transitionsDownTransit[key] = []
        # counting within cylinder (optional)
        CYLtransitionsUp[key] = 0
        CYLtransitionsDown[key] = 0
        CYLtransitionsUpTimes[key] = []
        CYLtransitionsDownTimes[key] = []

        state0 = 0  # for tracking 1->2->3->4
        state1 = 5  # for tracking 4->3->2->1
        CYLstate0 = 0  # for tracking 2->3 in the cylinder
        CYLstate1 = 0  # for tracking 3->2 in the cylinder
        entry_time_up = None    # time when ion entered lower headgroup (region 1->2)
        entry_time_down = None  # time when ion entered upper headgroup (region 4->3)

        counter = 0
        for i in ions[key]:
            # Record membrane entry times before state updates.
            # Upward: ion steps from region 1 into region 2 (lower headgroup entry).
            # entry_time_up is reset each time the ion re-enters from region 1,
            # so transit is always measured from the most recent membrane entry.
            if state0 == 1 and i == 2:
                entry_time_up = trjTimes[counter]
            # Downward: ion steps from region 4 into region 3 (upper headgroup entry).
            if state1 == 4 and i == 3:
                entry_time_down = trjTimes[counter]
            if i == state0 + 1:  # step forward
                state0 += 1
            elif state0 == 4 and i == 1:  # periodic step forward
                state0 = 1
            elif i == state0 - 1:  # step backward
                state0 -= 1
            elif state0 == 1 and i == 4:  # periodic step backward
                state0 = 0
            # WARNING: jump over two positions
            elif state0 == 2 and i == 4:
                print('WARNING: large ion jump, your trj output frequency is too low')
                state0 = 0
            elif state0 == 4 and i == 2:
                print('WARNING: large ion jump, your trj output frequency is too low')
                state0 = 2
            elif state0 == 3 and i == 1:
                print('WARNING: large ion jump, your trj output frequency is too low')
                state0 = 4
            elif state0 == 1 and i == 3:
                print('WARNING: large ion jump, your trj output frequency is too low')
                state0 = 0
            ######## cylinder #######
            if bCyl:
                if state0 == 2:
                    if ionsCyl[key][counter] == 1:
                        CYLstate0 = 1
                    else:
                        CYLstate0 = 0
                if state0 == 3:
                    if ionsCyl[key][counter] == 1 and (CYLstate0 == 1 or CYLstate0 == 2):
                        CYLstate0 = 2
                    else:
                        CYLstate0 = 0
            #########################

            if i == state1 - 1:  # step forward
                state1 -= 1
            elif state1 == 1 and i == 4:  # periodic step forward
                state1 = 4
            elif i == state1 + 1:  # step backward
                state1 += 1
            elif state1 == 4 and i == 1:  # periodic step backward
                state1 = 5
            elif state1 == 2 and i == 4:
                print('WARNING: large ion jump, your trj output frequency is too low')
                state1 = 1
            elif state1 == 4 and i == 2:
                print('WARNING: large ion jump, your trj output frequency is too low')
                state1 = 5
            elif state1 == 3 and i == 1:
                print('WARNING: large ion jump, your trj output frequency is too low')
                state1 = 5
            elif state1 == 1 and i == 3:
                print('WARNING: large ion jump, your trj output frequency is too low')
                state1 = 3
            ######## cylinder #######
            if bCyl:
                if state1 == 3:
                    if ionsCyl[key][counter] == 1:
                        CYLstate1 = 1
                    else:
                        CYLstate1 = 0
                if state1 == 2:
                    if ionsCyl[key][counter] == 1 and (CYLstate1 == 1 or CYLstate1 == 2):
                        CYLstate1 = 2
                    else:
                        CYLstate1 = 0
            #########################

            if state0 == 4:  # transition up
                state0 = 0
                transitionsUp[key] += 1
                transitionsUpTimes[key].append(trjTimes[counter])
                if entry_time_up is not None:
                    transitionsUpTransit[key].append(trjTimes[counter] - entry_time_up)
                    entry_time_up = None
                totalUp += 1
                if bCyl == True and CYLstate0 == 2:
                    CYLstate0 = 0
                    CYLtransitionsUp[key] += 1
                    CYLtransitionsUpTimes[key].append(trjTimes[counter])
                    CYLtotalUp += 1
            if state1 == 1:  # transition down
                state1 = 5
                transitionsDown[key] += 1
                transitionsDownTimes[key].append(trjTimes[counter])
                if entry_time_down is not None:
                    transitionsDownTransit[key].append(trjTimes[counter] - entry_time_down)
                    entry_time_down = None
                totalDown += 1
                if bCyl == True and CYLstate1 == 2:
                    CYLstate1 = 0
                    CYLtransitionsDown[key] += 1
                    CYLtransitionsDownTimes[key].append(trjTimes[counter])
                    CYLtotalDown += 1

            counter += 1

    return (transitionsUp, transitionsDown, transitionsUpTimes, transitionsDownTimes,
            transitionsUpTransit, transitionsDownTransit,
            totalUp, totalDown,
            CYLtransitionsUp, CYLtransitionsDown, CYLtransitionsUpTimes, CYLtransitionsDownTimes, CYLtotalUp, CYLtotalDown)


def identify_range(z, plow, pmid, phigh):
    if z < plow:
        return 1
    elif z < pmid:
        return 2
    elif z < phigh:
        return 3
    else:
        return 4


def identify_if_in_cylinder(x, y, cylCenter, cylRadSq):
    """Check if (x,y) is within the cylinder. cylRadSq is radius squared."""
    d2 = (x - cylCenter[0]) ** 2 + (y - cylCenter[1]) ** 2
    return 1 if d2 <= cylRadSq else 0


def get_cylinder(crd):
    """Return the [x,y] center and radius (max xy-distance) of the cylinder."""
    center = crd[:, :2].mean(axis=0)
    rad = float(np.max(np.sqrt(np.sum((crd[:, :2] - center) ** 2, axis=1))))
    return center, rad


def get_limits(crd, z=2):
    pmid = np.mean(crd, axis=0)[z]
    indhigh = np.where(crd[:, z] > pmid)
    indlow = np.where(crd[:, z] < pmid)
    phigh = np.mean(crd[indhigh, z])
    plow = np.mean(crd[indlow, z])
    return plow, pmid, phigh


def read_ndx(fname):
    """Parse a GROMACS NDX file. Returns dict of {group_name: list of 0-based atom indices}."""
    groups = {}
    names = []
    current_group = None
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if line.startswith('['):
                current_group = line[1:line.index(']')].strip()
                groups[current_group] = []
                names.append(current_group)
            elif current_group is not None and line:
                groups[current_group].extend(int(x) - 1 for x in line.split())  # convert to 0-based
    return groups, names


def select_ndx(fname="index.ndx", message=False):
    if not os.path.isfile(fname):
        return False
    groups, names = read_ndx(fname)
    for i, name in enumerate(names):
        sys.stdout.write('%d %s: %d atoms\n' % (i, name, len(groups[name])))
    sys.stdout.write('\n')

    if not message:
        sys.stdout.write('Select a group for analysis:\n')
    else:
        sys.stdout.write(message + '\n')

    ndxNum = -1
    while ndxNum == -1:
        ndxNum = input()
        if not ndxNum.isdigit():
            sys.stdout.write('Wrong index group number selected (use an integer number)\n')
            ndxNum = -1
            continue

        ndxNum = int(ndxNum)
        if (ndxNum >= len(names)) or (ndxNum < 0):
            sys.stdout.write('Wrong index group number selected\n')
            ndxNum = -1
    sys.stdout.write('Selected group %d\n\n' % ndxNum)

    res = np.array(groups[names[ndxNum]])
    return res, ndxNum


def plot_permeations(transitionsUpTimes, transitionsDownTimes, trjTimes, outbase):
    """Plot cumulative permeation counts vs time with linear fits."""
    fig, ax = plt.subplots(figsize=(8, 5))
    t_end = trjTimes[-1]

    for times_dict, color, label in [
        (transitionsUpTimes,   'tab:blue', 'Up'),
        (transitionsDownTimes, 'tab:red',  'Down'),
    ]:
        all_times = sorted(t for k in times_dict for t in times_dict[k])
        if not all_times:
            continue
        n = len(all_times)

        # step function: starts at 0, increments at each event
        ax.step([0.0] + all_times, list(range(n + 1)),
                where='post', color=color, label=label, linewidth=1.5)

        # linear fit through the event times vs cumulative count
        if n >= 2:
            coeffs = np.polyfit(all_times, np.arange(1, n + 1, dtype=float), 1)
            t_fit = np.array([0.0, t_end])
            rate_per_ns = coeffs[0]  # slope is already events/ns since time is in ns
            ax.plot(t_fit, np.polyval(coeffs, t_fit), '--', color=color, alpha=0.7,
                    label=f'{label} fit: {rate_per_ns:.4f} events/ns')

    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('Cumulative permeation count')
    ax.set_title('Ion permeation events')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    outpath = outbase + '_permeations.png'
    plt.savefig(outpath, dpi=150)
    plt.close(fig)
    print('Saved: %s' % outpath)


def plot_ion_positions(ion_zpos, ion_resname, trjTimes, plows, pmids, phighs,
                       transitionsUpTimes, transitionsDownTimes, outbase):
    """Plot z-position of each ion vs time, one subplot per ion type.

    Non-permeating ions are drawn as thin gray lines.
    Permeating ions are drawn in distinct colours on top, with triangle
    markers at each permeation event (up = ▼ reaching top, down = ▲ reaching bottom).
    """
    resnames = sorted(set(ion_resname.values()))
    n_types = len(resnames)
    t_arr = np.array(trjTimes)
    mean_plow  = np.mean(plows)
    mean_pmid  = np.mean(pmids)
    mean_phigh = np.mean(phighs)

    # ions that permeated at least once in either direction
    permeating = set()
    for key, times in transitionsUpTimes.items():
        if times:
            permeating.add(key)
    for key, times in transitionsDownTimes.items():
        if times:
            permeating.add(key)

    fig, axes = plt.subplots(n_types, 1, figsize=(10, 4 * n_types), squeeze=False)

    for ax_idx, resname in enumerate(resnames):
        ax = axes[ax_idx, 0]
        ion_keys = [k for k, v in ion_resname.items() if v == resname]
        perm_keys = [k for k in ion_keys if k in permeating]

        # evenly spaced colors across the full spectrum — avoids the paired
        # similar shades that tab20 produces for small numbers of ions.
        # Divide by n_perm (not n_perm-1) so we never reach 1.0: gist_rainbow
        # wraps (start and end are both red/magenta), so 0/(n-1)…(n-1)/(n-1)
        # would make the first and last ion look identical.
        n_perm = len(perm_keys)
        perm_colors = [plt.cm.gist_rainbow(i / max(n_perm, 1))
                       for i in range(n_perm)]

        # --- non-permeating ions: small gray dots as background ---
        for key in ion_keys:
            if key not in permeating:
                ax.plot(t_arr, ion_zpos[key],
                        color='lightgray', linestyle='none', marker='.', markersize=1,
                        alpha=0.6, zorder=1)

        # --- permeating ions: colored dots, drawn on top ---
        for ci, key in enumerate(perm_keys):
            color = perm_colors[ci]
            z_arr = np.array(ion_zpos[key])
            ax.plot(t_arr, z_arr,
                    color=color, linestyle='none', marker='.', markersize=2,
                    alpha=0.85, zorder=2, label='ion %d' % int(key))

            # mark permeation events: downward triangle = crossed upward,
            #                         upward triangle   = crossed downward
            for t_event in transitionsUpTimes.get(key, []):
                fi = np.argmin(np.abs(t_arr - t_event))
                ax.plot(t_event, z_arr[fi], 'v',
                        color=color, markersize=7, markeredgecolor='black',
                        markeredgewidth=0.5, zorder=3)
            for t_event in transitionsDownTimes.get(key, []):
                fi = np.argmin(np.abs(t_arr - t_event))
                ax.plot(t_event, z_arr[fi], '^',
                        color=color, markersize=7, markeredgecolor='black',
                        markeredgewidth=0.5, zorder=3)

        # membrane boundaries
        ax.axhline(mean_plow,  color='steelblue', linestyle='--', linewidth=1.0,
                   label='plow  %.1f Å' % mean_plow)
        ax.axhline(mean_pmid,  color='steelblue', linestyle='-',  linewidth=1.0,
                   label='pmid  %.1f Å' % mean_pmid)
        ax.axhline(mean_phigh, color='steelblue', linestyle=':',  linewidth=1.0,
                   label='phigh %.1f Å' % mean_phigh)

        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Z position (Å)')
        ax.set_title('%s ions  (n=%d,  %d permeating)' % (resname, len(ion_keys), len(perm_keys)))
        ax.grid(True, alpha=0.2)

    plt.tight_layout()
    outpath = outbase + '_positions.png'
    plt.savefig(outpath, dpi=150)
    plt.close(fig)
    print('Saved: %s' % outpath)


def find_binding_sites(bin_centers, pmf, mean_plow, mean_phigh,
                       min_depth=0.5, min_width=1.0, margin=5.0,
                       contour_kBT=1.0):
    """Identify binding sites as local minima of the 1D PMF profile.

    Parameters
    ----------
    bin_centers : 1D array
        Z-positions of the PMF bins (Å).
    pmf : 1D array
        PMF values in kBT (NaN allowed).
    mean_plow, mean_phigh : float
        Membrane reference levels; sites are searched in
        [plow - margin, phigh + margin].
    min_depth : float
        Minimum prominence of the minimum, in kBT (default 0.5).
    min_width : float
        Minimum width of the minimum at half-prominence, in Å (default 1).
    margin : float
        How far outside the headgroup region to still allow sites, in Å.
    contour_kBT : float
        Energy contour above the minimum used to define each site's extent
        (default 1 kBT).

    Returns
    -------
    list of dicts, each with keys:
        z_min, depth_kBT, z_lo, z_hi
    """
    from scipy.signal import find_peaks

    z_lo_search = mean_plow  - margin
    z_hi_search = mean_phigh + margin

    # Mask: only consider bins inside the search window with finite PMF
    mask = (bin_centers >= z_lo_search) & (bin_centers <= z_hi_search) & np.isfinite(pmf)
    if not np.any(mask):
        return []

    idx_window = np.where(mask)[0]
    pmf_window = pmf[idx_window]
    z_window   = bin_centers[idx_window]

    # find_peaks works on maxima → invert PMF
    bin_width = float(bin_centers[1] - bin_centers[0])
    peaks, props = find_peaks(-pmf_window,
                              prominence=min_depth,
                              width=max(min_width / bin_width, 1.0))

    sites = []
    for pi in peaks:
        z_min      = float(z_window[pi])
        pmf_min    = float(pmf_window[pi])
        prominence = float(props['prominences'][np.where(peaks == pi)[0][0]])

        # Extent: walk left and right while PMF stays below pmf_min + contour_kBT
        threshold = pmf_min + contour_kBT
        lo = pi
        while lo > 0 and pmf_window[lo - 1] <= threshold:
            lo -= 1
        hi = pi
        while hi < len(pmf_window) - 1 and pmf_window[hi + 1] <= threshold:
            hi += 1
        sites.append({
            'z_min':      z_min,
            'pmf_min':    pmf_min,
            'depth_kBT':  prominence,
            'z_lo':       float(z_window[lo]),
            'z_hi':       float(z_window[hi]),
        })

    sites.sort(key=lambda s: s['z_min'])
    return sites


def compute_site_occupancy(ion_zpos, ion_keys, mean_box_z, sites):
    """Fraction of frames in which at least one ion of the species lies inside the site.

    Also returns mean occupancy = average number of ions per frame in the site.
    """
    if not sites or not ion_keys:
        return [(0.0, 0.0) for _ in sites]
    n_frames = len(ion_zpos[ion_keys[0]])
    z_arr = np.stack([np.array(ion_zpos[k]) % mean_box_z for k in ion_keys], axis=0)
    results = []
    for s in sites:
        in_site = (z_arr >= s['z_lo']) & (z_arr <= s['z_hi'])  # (n_ions, n_frames)
        per_frame_count = in_site.sum(axis=0)
        frac_occupied   = float(np.mean(per_frame_count > 0))
        mean_count      = float(np.mean(per_frame_count))
        results.append((frac_occupied, mean_count))
    return results


def plot_density_pmf(ion_zpos, ion_resname, plows, pmids, phighs, mean_box_z, outbase):
    """Plot z-density profile and PMF for each ion species.

    Steps:
      1. Wrap all ion z-positions into [0, box_z) to eliminate PBC artefacts.
      2. Histogram into 1-Å bins over the full box length.
      3. Divide by (n_frames * bin_width) → average ion count per Å per frame.
      4. Normalise: bulk reference = 95th percentile of the Gaussian-smoothed
         density.  For ions excluded from the membrane, the bulk plateau is the
         dominant density level so the 95th percentile lands squarely on it
         regardless of box geometry.  No additional pinning is applied so that
         asymmetric profiles (e.g. under applied voltage) are preserved.
      5. PMF = -ln(normalised density) in kBT.

    Outputs per ion species:
      <outbase>_<resname>_density_pmf.png  — two-panel figure
      <outbase>_<resname>_density_pmf.xvg  — plain-text data (z, raw, smooth, PMF)
    """
    from scipy.ndimage import gaussian_filter1d

    resnames   = sorted(set(ion_resname.values()))
    mean_plow  = np.mean(plows)
    mean_pmid  = np.mean(pmids)
    mean_phigh = np.mean(phighs)
    bin_width  = 1.0  # Å

    for resname in resnames:
        ion_keys = [k for k, v in ion_resname.items() if v == resname]
        n_frames = len(ion_zpos[ion_keys[0]])

        # Wrap z-coordinates into [0, box_z) to fold PBC artefact positions
        # (ions with z slightly > L_z) back into the correct bin.
        all_z = np.concatenate([np.array(ion_zpos[k]) for k in ion_keys])
        all_z = all_z % mean_box_z

        # Histogram over the full box length with uniform 1-Å bins
        bins        = np.arange(0.0, mean_box_z + bin_width, bin_width)
        hist, edges = np.histogram(all_z, bins=bins)
        bin_centers = 0.5 * (edges[:-1] + edges[1:])

        # Average number of ions per Å per frame
        density = hist / (n_frames * bin_width)

        # Smooth first — used for both the reference estimate and the PMF curve
        density_smooth = gaussian_filter1d(density.astype(float), sigma=2)

        # Bulk reference: 95th percentile of the smoothed density.
        # No secondary pinning is applied so that asymmetric systems
        # (e.g. applied voltage) retain their natural bulk asymmetry.
        bulk_ref     = np.percentile(density_smooth, 95)
        norm_density = density        / bulk_ref
        norm_smooth  = density_smooth / bulk_ref

        # PMF = -ln(rho / rho_bulk) in kBT; NaN where density is zero
        with np.errstate(divide='ignore', invalid='ignore'):
            pmf = np.where(norm_smooth > 0, -np.log(norm_smooth), np.nan)

        # ---- binding-site detection (local minima of the PMF) ----
        sites = find_binding_sites(bin_centers, pmf, mean_plow, mean_phigh)
        occ   = compute_site_occupancy(ion_zpos, ion_keys, mean_box_z, sites)

        # ---- plot ----
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True,
                                        gridspec_kw={'hspace': 0.05})

        # -- density panel --
        ax1.plot(bin_centers, norm_density, color='tab:blue', linewidth=0.7,
                 alpha=0.35, label='raw (1 Å bins)')
        ax1.plot(bin_centers, norm_smooth,  color='tab:blue', linewidth=2.0,
                 label='smoothed (σ = 2 Å)')
        ax1.axhline(1.0, color='gray', linestyle='--', linewidth=0.8, label='bulk = 1')
        ax1.axvline(mean_plow,  color='steelblue', linestyle='--', linewidth=1.0,
                    label='plow  %.1f Å' % mean_plow)
        ax1.axvline(mean_pmid,  color='steelblue', linestyle='-',  linewidth=1.0,
                    label='pmid  %.1f Å' % mean_pmid)
        ax1.axvline(mean_phigh, color='steelblue', linestyle=':',  linewidth=1.0,
                    label='phigh %.1f Å' % mean_phigh)
        ax1.set_ylabel('Density (normalised to bulk)')
        ax1.set_title('%s  —  z-density and potential of mean force' % resname)
        ax1.legend(fontsize=8, loc='upper right', ncol=2)
        ax1.grid(True, alpha=0.3)

        # -- PMF panel --
        ax2.plot(bin_centers, pmf, color='tab:red', linewidth=2.0)
        ax2.axhline(0.0, color='gray', linestyle='--', linewidth=0.8, label='0 kBT (bulk)')
        ax2.axvline(mean_plow,  color='steelblue', linestyle='--', linewidth=1.0,
                    label='plow')
        ax2.axvline(mean_pmid,  color='steelblue', linestyle='-',  linewidth=1.0,
                    label='pmid')
        ax2.axvline(mean_phigh, color='steelblue', linestyle=':',  linewidth=1.0,
                    label='phigh')
        # binding-site annotations on both panels
        for si, (s, (frac, mean_n)) in enumerate(zip(sites, occ), start=1):
            label = 'S%d' % si
            # shaded extent on density panel
            ax1.axvspan(s['z_lo'], s['z_hi'], color='orange', alpha=0.15)
            # marker + label on PMF panel
            ax2.axvspan(s['z_lo'], s['z_hi'], color='orange', alpha=0.15)
            ax2.plot(s['z_min'], s['pmf_min'], 'o', color='orange',
                     markeredgecolor='black', markersize=7, zorder=5)
            ax2.annotate('%s\n%.1f kBT\nocc %.2f' % (label, s['depth_kBT'], frac),
                         xy=(s['z_min'], s['pmf_min']),
                         xytext=(0, -22), textcoords='offset points',
                         ha='center', va='top', fontsize=7,
                         color='black')

        ax2.set_xlabel('Z position (Å)')
        ax2.set_ylabel('PMF (kBT)')
        ax2.legend(fontsize=8, loc='upper right')
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        outpath = '%s_%s_density_pmf.png' % (outbase, resname)
        plt.savefig(outpath, dpi=150)
        plt.close(fig)
        print('Saved: %s' % outpath)

        # ---- XVG output ----
        xvg_path = '%s_%s_density_pmf.xvg' % (outbase, resname)
        with open(xvg_path, 'w') as fxvg:
            fxvg.write('# Ion z-density and PMF for species: %s\n' % resname)
            fxvg.write('# Generated by count-ions-mdanalysis.py\n')
            fxvg.write('@    title "z-density and PMF: %s"\n' % resname)
            fxvg.write('@    xaxis label "Z position (A)"\n')
            fxvg.write('@    yaxis label "Density (norm.) / PMF (kBT)"\n')
            fxvg.write('@TYPE xy\n')
            fxvg.write('@ s0 legend "density raw"\n')
            fxvg.write('@ s1 legend "density smooth"\n')
            fxvg.write('@ s2 legend "PMF (kBT)"\n')
            for j in range(len(bin_centers)):
                pmf_val = pmf[j] if not np.isnan(pmf[j]) else 0.0
                fxvg.write('%10.4f  %12.6f  %12.6f  %12.6f\n' % (
                    bin_centers[j], norm_density[j], norm_smooth[j], pmf_val))
        print('Saved: %s' % xvg_path)

        # ---- binding sites table ----
        sites_path = '%s_%s_binding_sites.dat' % (outbase, resname)
        with open(sites_path, 'w') as fs:
            fs.write('# Binding sites for %s identified as local minima of the 1D PMF\n' % resname)
            fs.write('# Search window: [plow - 5, phigh + 5] Å, min depth 0.5 kBT, contour 1 kBT\n')
            fs.write('# %-4s %10s %10s %10s %12s %12s %12s\n' %
                     ('id', 'z_min(A)', 'z_lo(A)', 'z_hi(A)', 'depth(kBT)',
                      'occ_frac', 'mean_count'))
            for si, (s, (frac, mean_n)) in enumerate(zip(sites, occ), start=1):
                fs.write('  %-4d %10.3f %10.3f %10.3f %12.3f %12.4f %12.4f\n' %
                         (si, s['z_min'], s['z_lo'], s['z_hi'],
                          s['depth_kBT'], frac, mean_n))
        print('Saved: %s  (%d site%s)' % (sites_path, len(sites),
                                          '' if len(sites) == 1 else 's'))


def plot_density_pmf_local(ion_zpos, ion_resname, plows, pmids, phighs,
                           mean_box_z, outbase, margin=10.0):
    """Local-mode density and PMF zoomed to the membrane region.

    Complements the global plot by:
      • Restricting the view to [plow - margin, phigh + margin].
      • Using lighter Gaussian smoothing (σ = 1 Å) to preserve fine structure
        inside the pore / selectivity filter.
      • Normalising density to its local maximum within the window (not bulk).
      • Computing a detrended PMF: the PMF from σ=1 smoothing minus a heavy
        baseline (σ=10 Å).  The baseline captures the large-scale barrier
        shape; subtracting it exposes local wells (binding sites) that would
        otherwise be invisible as tiny ripples on a multi-kBT slope.
      • Running binding-site detection on the detrended PMF with a 0.3 kBT
        depth threshold.

    Outputs per ion species:
      <outbase>_<resname>_density_pmf_local.png   — three-panel figure
      <outbase>_<resname>_density_pmf_local.xvg   — numerical data
      <outbase>_<resname>_binding_sites_local.dat  — detected sites
    """
    from scipy.ndimage import gaussian_filter1d

    resnames   = sorted(set(ion_resname.values()))
    mean_plow  = np.mean(plows)
    mean_pmid  = np.mean(pmids)
    mean_phigh = np.mean(phighs)
    bin_width  = 1.0  # Å
    z_lo_view  = mean_plow  - margin
    z_hi_view  = mean_phigh + margin

    for resname in resnames:
        ion_keys = [k for k, v in ion_resname.items() if v == resname]
        n_frames = len(ion_zpos[ion_keys[0]])

        # Wrap and restrict to the membrane window
        all_z = np.concatenate([np.array(ion_zpos[k]) for k in ion_keys])
        all_z = all_z % mean_box_z

        # Histogram over the full box …
        bins_full   = np.arange(0.0, mean_box_z + bin_width, bin_width)
        hist, edges = np.histogram(all_z, bins=bins_full)
        bc_full     = 0.5 * (edges[:-1] + edges[1:])
        density_full = hist / (n_frames * bin_width)

        # … then crop to the window
        mask = (bc_full >= z_lo_view) & (bc_full <= z_hi_view)
        bin_centers = bc_full[mask]
        density     = density_full[mask]

        # Light smoothing (σ = 1 Å) to preserve local structure
        density_smooth = gaussian_filter1d(density.astype(float), sigma=1)

        # Local normalisation to peak density in window
        local_ref = density_smooth.max()
        if local_ref <= 0:
            local_ref = 1.0
        norm_density = density        / local_ref
        norm_smooth  = density_smooth / local_ref

        # PMF (locally normalised)
        with np.errstate(divide='ignore', invalid='ignore'):
            pmf = np.where(norm_smooth > 0, -np.log(norm_smooth), np.nan)

        # Detrended PMF: subtract a heavy baseline (σ = 10 Å) to remove the
        # large-scale barrier shape and expose local wells.
        pmf_for_baseline = np.where(np.isfinite(pmf), pmf, 0.0)
        pmf_baseline = gaussian_filter1d(pmf_for_baseline, sigma=10)
        pmf_detrended = pmf - pmf_baseline

        # ---- binding-site detection on the detrended profile ----
        # margin=0: search strictly inside [plow, phigh] to avoid spurious
        # minima from baseline edge artefacts in bulk.
        sites = find_binding_sites(bin_centers, pmf_detrended, mean_plow, mean_phigh,
                                   min_depth=0.2, min_width=1.0, margin=0.0,
                                   contour_kBT=0.2)
        occ   = compute_site_occupancy(ion_zpos, ion_keys, mean_box_z, sites)

        # ---- plot (3 panels) ----
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9, 10), sharex=True,
                                             gridspec_kw={'hspace': 0.06})

        # -- density panel --
        ax1.plot(bin_centers, norm_density, color='tab:blue', linewidth=0.7,
                 alpha=0.35, label='raw (1 Å bins)')
        ax1.plot(bin_centers, norm_smooth,  color='tab:blue', linewidth=2.0,
                 label='smoothed (σ = 1 Å)')
        ax1.axvline(mean_plow,  color='steelblue', linestyle='--', linewidth=1.0,
                    label='plow  %.1f Å' % mean_plow)
        ax1.axvline(mean_pmid,  color='steelblue', linestyle='-',  linewidth=1.0,
                    label='pmid  %.1f Å' % mean_pmid)
        ax1.axvline(mean_phigh, color='steelblue', linestyle=':',  linewidth=1.0,
                    label='phigh %.1f Å' % mean_phigh)
        for si, (s, (frac, mean_n)) in enumerate(zip(sites, occ), start=1):
            ax1.axvspan(s['z_lo'], s['z_hi'], color='orange', alpha=0.15)
        ax1.set_ylabel('Density (normalised to local max)')
        ax1.set_title('%s  —  local z-density and PMF (membrane region)' % resname)
        ax1.legend(fontsize=8, loc='upper right', ncol=2)
        ax1.grid(True, alpha=0.3)

        # -- PMF panel --
        ax2.plot(bin_centers, pmf, color='tab:red', linewidth=2.0, label='PMF (σ = 1 Å)')
        ax2.plot(bin_centers, pmf_baseline, color='gray', linewidth=1.5,
                 linestyle='--', label='baseline (σ = 10 Å)')
        ax2.axhline(0.0, color='gray', linestyle=':', linewidth=0.6)
        ax2.axvline(mean_plow,  color='steelblue', linestyle='--', linewidth=1.0)
        ax2.axvline(mean_pmid,  color='steelblue', linestyle='-',  linewidth=1.0)
        ax2.axvline(mean_phigh, color='steelblue', linestyle=':',  linewidth=1.0)
        for si, (s, (frac, mean_n)) in enumerate(zip(sites, occ), start=1):
            ax2.axvspan(s['z_lo'], s['z_hi'], color='orange', alpha=0.15)
        ax2.set_ylabel('PMF (kBT, relative to local max)')
        ax2.legend(fontsize=8, loc='upper right')
        ax2.grid(True, alpha=0.3)

        # -- detrended PMF panel --
        ax3.plot(bin_centers, pmf_detrended, color='tab:purple', linewidth=2.0)
        ax3.axhline(0.0, color='gray', linestyle='--', linewidth=0.8, label='baseline')
        ax3.axvline(mean_plow,  color='steelblue', linestyle='--', linewidth=1.0,
                    label='plow')
        ax3.axvline(mean_pmid,  color='steelblue', linestyle='-',  linewidth=1.0,
                    label='pmid')
        ax3.axvline(mean_phigh, color='steelblue', linestyle=':',  linewidth=1.0,
                    label='phigh')
        for si, (s, (frac, mean_n)) in enumerate(zip(sites, occ), start=1):
            slabel = 'S%d' % si
            ax3.axvspan(s['z_lo'], s['z_hi'], color='orange', alpha=0.15)
            ax3.plot(s['z_min'], s['pmf_min'], 'o', color='orange',
                     markeredgecolor='black', markersize=7, zorder=5)
            ax3.annotate('%s\n%.2f kBT\nocc %.2f' % (slabel, s['depth_kBT'], frac),
                         xy=(s['z_min'], s['pmf_min']),
                         xytext=(0, -22), textcoords='offset points',
                         ha='center', va='top', fontsize=7,
                         color='black')
        ax3.set_xlabel('Z position (Å)')
        ax3.set_ylabel('Detrended PMF (kBT)')
        ax3.legend(fontsize=8, loc='upper right')
        ax3.grid(True, alpha=0.3)

        plt.tight_layout()
        outpath = '%s_%s_density_pmf_local.png' % (outbase, resname)
        plt.savefig(outpath, dpi=150)
        plt.close(fig)
        print('Saved: %s' % outpath)

        # ---- XVG output ----
        xvg_path = '%s_%s_density_pmf_local.xvg' % (outbase, resname)
        with open(xvg_path, 'w') as fxvg:
            fxvg.write('# Local density and PMF for %s (membrane region)\n' % resname)
            fxvg.write('# Smoothing: sigma = 1 A (fine), sigma = 10 A (baseline)\n')
            fxvg.write('# Detrended PMF = PMF - baseline\n')
            fxvg.write('@    title "local density and PMF: %s"\n' % resname)
            fxvg.write('@    xaxis label "Z position (A)"\n')
            fxvg.write('@TYPE xy\n')
            fxvg.write('@ s0 legend "density raw"\n')
            fxvg.write('@ s1 legend "density smooth"\n')
            fxvg.write('@ s2 legend "PMF (kBT)"\n')
            fxvg.write('@ s3 legend "PMF baseline (kBT)"\n')
            fxvg.write('@ s4 legend "PMF detrended (kBT)"\n')
            for j in range(len(bin_centers)):
                pmf_val = pmf[j] if not np.isnan(pmf[j]) else 0.0
                det_val = pmf_detrended[j] if not np.isnan(pmf_detrended[j]) else 0.0
                fxvg.write('%10.4f  %12.6f  %12.6f  %12.6f  %12.6f  %12.6f\n' % (
                    bin_centers[j], norm_density[j], norm_smooth[j],
                    pmf_val, pmf_baseline[j], det_val))
        print('Saved: %s' % xvg_path)

        # ---- binding sites table ----
        sites_path = '%s_%s_binding_sites_local.dat' % (outbase, resname)
        with open(sites_path, 'w') as fs:
            fs.write('# Binding sites (detrended local PMF) for %s\n' % resname)
            fs.write('# Membrane window: [%.1f, %.1f] Å, min depth 0.2 kBT, '
                     'contour 0.2 kBT\n' % (z_lo_view, z_hi_view))
            fs.write('# Fine smooth σ = 1 Å, baseline σ = 10 Å\n')
            fs.write('# %-4s %10s %10s %10s %12s %12s %12s\n' %
                     ('id', 'z_min(A)', 'z_lo(A)', 'z_hi(A)', 'depth(kBT)',
                      'occ_frac', 'mean_count'))
            for si, (s, (frac, mean_n)) in enumerate(zip(sites, occ), start=1):
                fs.write('  %-4d %10.3f %10.3f %10.3f %12.3f %12.4f %12.4f\n' %
                         (si, s['z_min'], s['z_lo'], s['z_hi'],
                          s['depth_kBT'], frac, mean_n))
        print('Saved: %s  (%d site%s)' % (sites_path, len(sites),
                                          '' if len(sites) == 1 else 's'))


def identify_site_residues(u, ndxCyl, sites, z_margin=2.0, r_cutoff=10.0):
    """Find protein residues near each binding site.

    For each site, selects protein atoms (from the cylinder group) whose
    mean z-position falls within [z_lo - z_margin, z_hi + z_margin].
    Returns a list (one per site) of sorted unique residue strings
    like 'THR75', 'VAL76'.

    Parameters
    ----------
    u : mda.Universe
        The MDAnalysis Universe (positioned at any frame — uses current
        positions as a representative snapshot).
    ndxCyl : array
        0-based atom indices for the cylinder/protein group.
    sites : list of dicts
        Binding sites with 'z_lo' and 'z_hi' keys.
    z_margin : float
        Extra Å above/below the site extent to search (default 2).
    r_cutoff : float
        Maximum xy-distance (Å) from the protein group centroid to include
        an atom (default 10). Filters out lipid-facing residues in large
        protein complexes.
    """
    if not sites or len(ndxCyl) == 0:
        return [[] for _ in sites]

    positions = u.atoms.positions
    prot_pos = positions[ndxCyl]
    # Cylinder center in xy
    cx = np.mean(prot_pos[:, 0])
    cy = np.mean(prot_pos[:, 1])

    results = []
    for s in sites:
        z_lo = s['z_lo'] - z_margin
        z_hi = s['z_hi'] + z_margin
        residues = set()
        for i, idx in enumerate(ndxCyl):
            az = prot_pos[i, 2]
            if az < z_lo or az > z_hi:
                continue
            # xy distance from pore center
            dx = prot_pos[i, 0] - cx
            dy = prot_pos[i, 1] - cy
            if np.sqrt(dx*dx + dy*dy) > r_cutoff:
                continue
            atom = u.atoms[int(idx)]
            residues.add((atom.resname, atom.resid))
        # Sort by resid, format as "RESNAME RESID"
        sorted_res = sorted(residues, key=lambda x: x[1])
        results.append(['%s%d' % (rn, ri) for rn, ri in sorted_res])
    return results


def write_binding_sites_pdb(u, ndxCyl, sites, outpath):
    """Write a PDB file with the protein and dummy atoms at binding-site centers.

    The protein atoms are taken from the cylinder group (current Universe
    frame).  For each binding site, a dummy atom (element X, residue name
    DUM, chain Z) is placed at (pore_center_x, pore_center_y, z_min) so
    the sites can be visualised in VMD / PyMOL / ChimeraX overlaid on the
    protein structure.
    """
    positions = u.atoms.positions
    prot_pos = positions[ndxCyl]
    cx = np.mean(prot_pos[:, 0])
    cy = np.mean(prot_pos[:, 1])

    with open(outpath, 'w') as fp:
        fp.write('REMARK  Binding sites visualisation\n')
        fp.write('REMARK  Protein atoms from cylinder group + DUM atoms at site centers\n')

        serial = 1
        # Write protein atoms
        for i, idx in enumerate(ndxCyl):
            atom = u.atoms[int(idx)]
            x, y, z = prot_pos[i]
            aname = atom.name[:4].ljust(4)
            resname = atom.resname[:3].rjust(3)
            chain = atom.chainID if hasattr(atom, 'chainID') and atom.chainID else 'A'
            chain = chain[0] if len(chain) > 0 else 'A'
            resid = atom.resid % 10000  # PDB resid field is 4 digits
            fp.write('ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00\n' %
                     (serial % 100000, aname, resname, chain, resid, x, y, z))
            serial += 1

        fp.write('TER\n')

        # Write dummy atoms at binding-site centers
        for si, s in enumerate(sites, start=1):
            resid = si
            fp.write('HETATM%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00%6.2f\n' %
                     (serial % 100000, 'DU', 'DUM', 'Z', resid,
                      cx, cy, s['z_min'], s['depth_kBT']))
            serial += 1

        fp.write('END\n')
    print('Saved: %s  (%d dummy atoms for binding sites)' % (outpath, len(sites)))


def plot_density_pmf_pore(ion_zpos, ionsCyl, ion_resname, plows, pmids, phighs,
                          mean_box_z, outbase, u=None, ndxCyl=None, margin=5.0):
    """Pore-only density and PMF using only ion positions inside the cylinder.

    When -cyl is used, ionsCyl stores 0/1 per frame per ion indicating whether
    the ion was inside the pore cylinder.  By restricting the histogram to
    those positions, bulk ions are excluded and the fine structure of the
    selectivity filter (binding sites) becomes clearly resolved.

    Uses 0.5 Å bins, σ = 0.5 Å smoothing, and no detrending — because bulk
    ions are excluded, the PMF directly reflects the pore landscape without
    needing baseline subtraction.

    Outputs per ion species:
      <outbase>_<resname>_density_pmf_pore.png
      <outbase>_<resname>_density_pmf_pore.xvg
      <outbase>_<resname>_binding_sites_pore.dat
    """
    from scipy.ndimage import gaussian_filter1d

    resnames   = sorted(set(ion_resname.values()))
    mean_plow  = np.mean(plows)
    mean_pmid  = np.mean(pmids)
    mean_phigh = np.mean(phighs)
    bin_width  = 0.5  # Å — fine bins for resolving ~3.5 Å-spaced sites
    z_lo_view  = mean_plow  - margin
    z_hi_view  = mean_phigh + margin

    for resname in resnames:
        ion_keys = [k for k, v in ion_resname.items() if v == resname]

        # Collect only z-positions where the ion was inside the cylinder
        pore_z = []
        for k in ion_keys:
            zs = np.array(ion_zpos[k])
            cyl = np.array(ionsCyl[k])
            pore_z.append(zs[cyl == 1])
        all_z = np.concatenate(pore_z)
        all_z = all_z % mean_box_z

        if len(all_z) == 0:
            print('No in-cylinder data for %s — skipping pore PMF' % resname)
            continue

        n_frames = len(ion_zpos[ion_keys[0]])

        # Histogram restricted to the membrane window
        bins = np.arange(z_lo_view, z_hi_view + bin_width, bin_width)
        hist, edges = np.histogram(all_z, bins=bins)
        bin_centers = 0.5 * (edges[:-1] + edges[1:])

        density = hist / (n_frames * bin_width)
        density_smooth = gaussian_filter1d(density.astype(float), sigma=1)  # σ = 0.5 Å in bins

        # PMF = -ln(density) + const, shifted so minimum = 0
        with np.errstate(divide='ignore', invalid='ignore'):
            pmf = np.where(density_smooth > 0, -np.log(density_smooth), np.nan)
        pmf -= np.nanmin(pmf)

        # Detrended PMF: subtract a heavy baseline (σ = 5 Å = 10 bins) to
        # remove the large-scale slope from vestibule→filter and expose
        # individual binding sites.
        pmf_for_baseline = np.where(np.isfinite(pmf), pmf, 0.0)
        pmf_baseline = gaussian_filter1d(pmf_for_baseline, sigma=10)  # σ=5 Å
        pmf_detrended = pmf - pmf_baseline

        # ---- binding-site detection on detrended profile ----
        sites = find_binding_sites(bin_centers, pmf_detrended, mean_plow, mean_phigh,
                                   min_depth=0.1, min_width=0.5, margin=2.0,
                                   contour_kBT=0.15)
        occ = compute_site_occupancy(ion_zpos, ion_keys, mean_box_z, sites)
        if u is not None and ndxCyl is not None:
            site_residues = identify_site_residues(u, ndxCyl, sites)
        else:
            site_residues = [[] for _ in sites]

        # ---- plot (3 panels) ----
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9, 10), sharex=True,
                                             gridspec_kw={'hspace': 0.06})

        ax1.bar(bin_centers, hist, width=bin_width, color='tab:blue', alpha=0.3,
                label='raw (0.5 Å bins)')
        ax1.plot(bin_centers, density_smooth * n_frames * bin_width,
                 color='tab:blue', linewidth=2, label='smoothed (σ = 0.5 Å)')
        ax1.axvline(mean_plow,  color='steelblue', linestyle='--', linewidth=1.0,
                    label='plow  %.1f Å' % mean_plow)
        ax1.axvline(mean_pmid,  color='steelblue', linestyle='-',  linewidth=1.0,
                    label='pmid  %.1f Å' % mean_pmid)
        ax1.axvline(mean_phigh, color='steelblue', linestyle=':',  linewidth=1.0,
                    label='phigh %.1f Å' % mean_phigh)
        for si, (s, _) in enumerate(zip(sites, occ), start=1):
            ax1.axvspan(s['z_lo'], s['z_hi'], color='orange', alpha=0.15)
        ax1.set_ylabel('Count (in-cylinder ions)')
        ax1.set_title('%s  —  pore density and PMF (cylinder-filtered)' % resname)
        ax1.legend(fontsize=8, loc='upper right', ncol=2)
        ax1.grid(True, alpha=0.3)

        ax2.plot(bin_centers, pmf, color='tab:red', linewidth=2.0, label='PMF')
        ax2.plot(bin_centers, pmf_baseline, color='gray', linewidth=1.5,
                 linestyle='--', label='baseline (σ = 5 Å)')
        ax2.axhline(0.0, color='gray', linestyle=':', linewidth=0.6)
        ax2.axvline(mean_plow,  color='steelblue', linestyle='--', linewidth=1.0)
        ax2.axvline(mean_pmid,  color='steelblue', linestyle='-',  linewidth=1.0)
        ax2.axvline(mean_phigh, color='steelblue', linestyle=':',  linewidth=1.0)
        for si, (s, _) in enumerate(zip(sites, occ), start=1):
            ax2.axvspan(s['z_lo'], s['z_hi'], color='orange', alpha=0.15)
        ax2.set_ylabel('PMF (kBT)')
        ax2.legend(fontsize=8, loc='upper right')
        ax2.grid(True, alpha=0.3)

        ax3.plot(bin_centers, pmf_detrended, color='tab:purple', linewidth=2.0)
        ax3.axhline(0.0, color='gray', linestyle='--', linewidth=0.8)
        ax3.axvline(mean_plow,  color='steelblue', linestyle='--', linewidth=1.0)
        ax3.axvline(mean_pmid,  color='steelblue', linestyle='-',  linewidth=1.0)
        ax3.axvline(mean_phigh, color='steelblue', linestyle=':',  linewidth=1.0)
        for si, (s, (frac, mean_n), reslist) in enumerate(
                zip(sites, occ, site_residues), start=1):
            slabel = 'S%d' % si
            ax3.axvspan(s['z_lo'], s['z_hi'], color='orange', alpha=0.15)
            ax3.plot(s['z_min'], s['pmf_min'], 'o', color='orange',
                     markeredgecolor='black', markersize=7, zorder=5)
            # Build annotation: site label + depth + occupancy + residues
            ann_text = '%s\n%.2f kBT\nocc %.2f' % (slabel, s['depth_kBT'], frac)
            if reslist:
                # Show up to 6 residues in annotation; full list in .dat file
                res_str = ', '.join(reslist[:6])
                if len(reslist) > 6:
                    res_str += ', ...'
                ann_text += '\n' + res_str
            ax3.annotate(ann_text,
                         xy=(s['z_min'], s['pmf_min']),
                         xytext=(0, -18), textcoords='offset points',
                         ha='center', va='top', fontsize=6,
                         color='black')
        ax3.set_xlabel('Z position (Å)')
        ax3.set_ylabel('Detrended PMF (kBT)')
        ax3.grid(True, alpha=0.3)

        plt.tight_layout()
        outpath = '%s_%s_density_pmf_pore.png' % (outbase, resname)
        plt.savefig(outpath, dpi=150)
        plt.close(fig)
        print('Saved: %s' % outpath)

        # ---- XVG output ----
        xvg_path = '%s_%s_density_pmf_pore.xvg' % (outbase, resname)
        with open(xvg_path, 'w') as fxvg:
            fxvg.write('# Pore density and PMF for %s (cylinder-filtered)\n' % resname)
            fxvg.write('# Only positions where ion was inside the cylinder\n')
            fxvg.write('# Bin width = 0.5 A, smoothing sigma = 0.5 A\n')
            fxvg.write('@    title "pore density and PMF: %s"\n' % resname)
            fxvg.write('@    xaxis label "Z position (A)"\n')
            fxvg.write('@TYPE xy\n')
            fxvg.write('@ s0 legend "count raw"\n')
            fxvg.write('@ s1 legend "density smooth"\n')
            fxvg.write('@ s2 legend "PMF (kBT)"\n')
            fxvg.write('@ s3 legend "PMF baseline (kBT)"\n')
            fxvg.write('@ s4 legend "PMF detrended (kBT)"\n')
            for j in range(len(bin_centers)):
                pmf_val = pmf[j] if not np.isnan(pmf[j]) else 0.0
                det_val = pmf_detrended[j] if not np.isnan(pmf_detrended[j]) else 0.0
                fxvg.write('%10.4f  %12.0f  %12.6f  %12.6f  %12.6f  %12.6f\n' % (
                    bin_centers[j], hist[j], density_smooth[j],
                    pmf_val, pmf_baseline[j], det_val))
        print('Saved: %s' % xvg_path)

        # ---- binding sites table ----
        sites_path = '%s_%s_binding_sites_pore.dat' % (outbase, resname)
        with open(sites_path, 'w') as fs:
            fs.write('# Binding sites (pore PMF) for %s\n' % resname)
            fs.write('# Cylinder-filtered, 0.5 Å bins, σ = 0.5 Å\n')
            fs.write('# %-4s %10s %10s %10s %12s %12s %12s  %s\n' %
                     ('id', 'z_min(A)', 'z_lo(A)', 'z_hi(A)', 'depth(kBT)',
                      'occ_frac', 'mean_count', 'nearby_residues'))
            for si, (s, (frac, mean_n), reslist) in enumerate(
                    zip(sites, occ, site_residues), start=1):
                res_str = ', '.join(reslist) if reslist else '-'
                fs.write('  %-4d %10.3f %10.3f %10.3f %12.3f %12.4f %12.4f  %s\n' %
                         (si, s['z_min'], s['z_lo'], s['z_hi'],
                          s['depth_kBT'], frac, mean_n, res_str))
        print('Saved: %s  (%d site%s)' % (sites_path, len(sites),
                                          '' if len(sites) == 1 else 's'))

        # ---- PDB with protein + dummy atoms at binding sites ----
        if u is not None and ndxCyl is not None and sites:
            pdb_path = '%s_%s_binding_sites_pore.pdb' % (outbase, resname)
            write_binding_sites_pdb(u, ndxCyl, sites, pdb_path)


def plot_density_pmf_2d(ion_xpos, ion_ypos, ion_zpos, ion_resname,
                        mean_box_x, mean_box_y, mean_box_z,
                        plows, phighs, outbase):
    """Plot 2D density profiles and PMFs in the XZ and YZ planes per ion species.

    Same logic as the 1D version:
      • Wrap coordinates into [0, box) to fold PBC artefacts.
      • 2D histogram with 1-Å square bins.
      • Normalise by (n_frames * bin_area) -> ions/Å² per frame.
      • Bulk reference = 95th percentile of the smoothed 2D density (no pinning).
      • PMF = -ln(rho/rho_bulk) in kBT.

    Outputs per ion species and plane:
      <outbase>_<resname>_density_pmf_<plane>.png  (4-panel: density + PMF)
      <outbase>_<resname>_density_pmf_<plane>.dat  (grid data; columns: h, z, rho_raw, rho_smooth, PMF)
    """
    from scipy.ndimage import gaussian_filter

    resnames   = sorted(set(ion_resname.values()))
    mean_plow  = np.mean(plows)
    mean_phigh = np.mean(phighs)
    bin_width  = 1.0  # Å

    planes = [
        ('XZ', ion_xpos, mean_box_x, 'X (Å)'),
        ('YZ', ion_ypos, mean_box_y, 'Y (Å)'),
    ]

    for resname in resnames:
        ion_keys = [k for k, v in ion_resname.items() if v == resname]
        n_frames = len(ion_zpos[ion_keys[0]])

        for plane_name, hpos_dict, box_h, hlabel in planes:
            all_h = np.concatenate([np.array(hpos_dict[k]) for k in ion_keys]) % box_h
            all_z = np.concatenate([np.array(ion_zpos[k]) for k in ion_keys]) % mean_box_z

            h_bins = np.arange(0.0, box_h + bin_width, bin_width)
            z_bins = np.arange(0.0, mean_box_z + bin_width, bin_width)

            hist, h_edges, z_edges = np.histogram2d(all_h, all_z, bins=[h_bins, z_bins])
            h_centers = 0.5 * (h_edges[:-1] + h_edges[1:])
            z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])

            # ions per Å² per frame
            density = hist / (n_frames * bin_width * bin_width)
            density_smooth = gaussian_filter(density.astype(float), sigma=2)

            bulk_ref = np.percentile(density_smooth, 95)
            if bulk_ref <= 0:
                bulk_ref = density_smooth.max() if density_smooth.max() > 0 else 1.0
            norm_density = density        / bulk_ref
            norm_smooth  = density_smooth / bulk_ref

            with np.errstate(divide='ignore', invalid='ignore'):
                pmf = np.where(norm_smooth > 0, -np.log(norm_smooth), np.nan)

            # ---- plot ----
            fig, axes = plt.subplots(1, 2, figsize=(12, 5))
            extent = [h_centers[0], h_centers[-1], z_centers[0], z_centers[-1]]

            # density panel
            im1 = axes[0].imshow(norm_smooth.T, origin='lower', extent=extent,
                                 aspect='auto', cmap='viridis',
                                 vmin=0.0, vmax=max(1.5, np.nanpercentile(norm_smooth, 99)))
            axes[0].axhline(mean_plow,  color='white', linestyle='--', linewidth=0.8)
            axes[0].axhline(mean_phigh, color='white', linestyle=':',  linewidth=0.8)
            axes[0].set_xlabel(hlabel)
            axes[0].set_ylabel('Z (Å)')
            axes[0].set_title('%s — density (norm.) [%s plane]' % (resname, plane_name))
            plt.colorbar(im1, ax=axes[0], label='ρ / ρ_bulk')

            # PMF panel — clip for display
            pmf_disp = np.where(np.isnan(pmf), np.nanmax(pmf[np.isfinite(pmf)]) if np.any(np.isfinite(pmf)) else 0.0, pmf)
            vmax = np.nanpercentile(pmf_disp, 99) if np.any(np.isfinite(pmf_disp)) else 5.0
            vmin = np.nanpercentile(pmf_disp, 1)  if np.any(np.isfinite(pmf_disp)) else -1.0
            im2 = axes[1].imshow(pmf_disp.T, origin='lower', extent=extent,
                                 aspect='auto', cmap='RdBu_r', vmin=vmin, vmax=vmax)
            axes[1].axhline(mean_plow,  color='black', linestyle='--', linewidth=0.8)
            axes[1].axhline(mean_phigh, color='black', linestyle=':',  linewidth=0.8)
            axes[1].set_xlabel(hlabel)
            axes[1].set_ylabel('Z (Å)')
            axes[1].set_title('%s — PMF (kBT) [%s plane]' % (resname, plane_name))
            plt.colorbar(im2, ax=axes[1], label='PMF (kBT)')

            plt.tight_layout()
            outpath = '%s_%s_density_pmf_%s.png' % (outbase, resname, plane_name)
            plt.savefig(outpath, dpi=150)
            plt.close(fig)
            print('Saved: %s' % outpath)

            # ---- data file ----
            dat_path = '%s_%s_density_pmf_%s.dat' % (outbase, resname, plane_name)
            with open(dat_path, 'w') as fdat:
                fdat.write('# 2D density and PMF for %s in %s plane\n' % (resname, plane_name))
                fdat.write('# Generated by count-ions-mdanalysis.py\n')
                fdat.write('# Columns: %s_center(A)  Z_center(A)  rho_raw_norm  rho_smooth_norm  PMF(kBT)\n'
                           % plane_name[0])
                for ih in range(len(h_centers)):
                    for iz in range(len(z_centers)):
                        pmf_val = pmf[ih, iz] if not np.isnan(pmf[ih, iz]) else 0.0
                        fdat.write('%10.4f  %10.4f  %12.6f  %12.6f  %12.6f\n' % (
                            h_centers[ih], z_centers[iz],
                            norm_density[ih, iz], norm_smooth[ih, iz], pmf_val))
                    fdat.write('\n')
            print('Saved: %s' % dat_path)


def compute_coordination(u, permeating_keys, ion_resname, ndxCyl=None,
                          cutoff=3.5):
    """Second trajectory pass: compute coordination numbers for permeating ions.

    For each permeating ion and each frame, counts the number of water
    oxygens (name OW) and protein heavy atoms within `cutoff` Å.

    Parameters
    ----------
    u : mda.Universe
    permeating_keys : set of int
        0-based atom indices of ions that permeated at least once.
    ion_resname : dict
        Mapping ion index → residue name.
    ndxCyl : array or None
        Protein group indices (for protein coordination). If None, only
        water coordination is computed.
    cutoff : float
        First coordination shell cutoff in Å (default 3.5).

    Returns
    -------
    coord_data : dict
        {ion_idx: {'z': [...], 'n_water': [...], 'n_prot': [...]}}
    """
    from MDAnalysis.lib.distances import distance_array

    water_ox = u.select_atoms('name OW')
    if len(water_ox) == 0:
        # Try alternative water oxygen names
        water_ox = u.select_atoms('name OH2 or name O and resname HOH WAT TIP3')
    n_water_atoms = len(water_ox)

    if ndxCyl is not None and len(ndxCyl) > 0:
        prot_ag = u.atoms[ndxCyl]
        # Heavy atoms only (mass > 2 excludes H regardless of naming)
        prot_heavy = prot_ag.select_atoms('not name H*')
        if len(prot_heavy) == 0:
            prot_heavy = prot_ag.atoms[[i for i in range(len(prot_ag))
                                        if prot_ag[i].mass > 2.0]]
        n_prot_atoms = len(prot_heavy)
    else:
        prot_heavy = None
        n_prot_atoms = 0

    print('Coordination analysis: %d permeating ions, %d water O, %d protein heavy atoms, cutoff %.1f Å'
          % (len(permeating_keys), n_water_atoms, n_prot_atoms, cutoff))

    perm_list = sorted(permeating_keys)
    coord_data = {k: {'z': [], 'n_water': [], 'n_prot': []} for k in perm_list}

    n_frames = len(u.trajectory)
    for fi, ts in enumerate(u.trajectory):
        if fi % 500 == 0:
            sys.stdout.write('Coordination: frame %d / %d\r' % (fi, n_frames))
            sys.stdout.flush()

        box = ts.dimensions
        positions = u.atoms.positions
        wat_pos = water_ox.positions
        prot_pos = prot_heavy.positions if prot_heavy is not None else None

        for key in perm_list:
            ion_pos = positions[key].reshape(1, 3)
            z = positions[key][2]

            # Water coordination
            dists_w = distance_array(ion_pos, wat_pos, box=box)[0]
            n_w = int(np.sum(dists_w <= cutoff))

            # Protein coordination
            if prot_pos is not None:
                dists_p = distance_array(ion_pos, prot_pos, box=box)[0]
                n_p = int(np.sum(dists_p <= cutoff))
            else:
                n_p = 0

            coord_data[key]['z'].append(z)
            coord_data[key]['n_water'].append(n_w)
            coord_data[key]['n_prot'].append(n_p)

    sys.stdout.write('Coordination: done (%d frames)          \n' % n_frames)
    return coord_data


def plot_coordination(coord_data, ion_resname, plows, pmids, phighs,
                      mean_box_z, outbase, bin_width=1.0):
    """Plot average coordination number vs z-position for permeating ions.

    Three curves: water coordination, protein coordination, total.
    """
    from scipy.ndimage import gaussian_filter1d

    mean_plow  = np.mean(plows)
    mean_pmid  = np.mean(pmids)
    mean_phigh = np.mean(phighs)

    resnames = sorted(set(ion_resname[k] for k in coord_data))

    for resname in resnames:
        keys = [k for k in coord_data if ion_resname[k] == resname]
        if not keys:
            continue

        # Pool all (z, n_water, n_prot) data from all permeating ions
        all_z = []
        all_nw = []
        all_np = []
        for k in keys:
            all_z.extend(coord_data[k]['z'])
            all_nw.extend(coord_data[k]['n_water'])
            all_np.extend(coord_data[k]['n_prot'])
        all_z  = np.array(all_z) % mean_box_z
        all_nw = np.array(all_nw, dtype=float)
        all_np = np.array(all_np, dtype=float)
        all_nt = all_nw + all_np

        # Bin by z and compute mean in each bin
        bins = np.arange(0.0, mean_box_z + bin_width, bin_width)
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        bin_idx = np.digitize(all_z, bins) - 1
        n_bins = len(bin_centers)

        mean_nw = np.full(n_bins, np.nan)
        mean_np = np.full(n_bins, np.nan)
        mean_nt = np.full(n_bins, np.nan)
        for bi in range(n_bins):
            mask = bin_idx == bi
            if np.sum(mask) > 0:
                mean_nw[bi] = np.mean(all_nw[mask])
                mean_np[bi] = np.mean(all_np[mask])
                mean_nt[bi] = np.mean(all_nt[mask])

        # Light smoothing
        valid = np.isfinite(mean_nw)
        nw_smooth = np.full(n_bins, np.nan)
        np_smooth = np.full(n_bins, np.nan)
        nt_smooth = np.full(n_bins, np.nan)
        if np.sum(valid) > 3:
            nw_smooth[valid] = gaussian_filter1d(mean_nw[valid], sigma=1)
            np_smooth[valid] = gaussian_filter1d(mean_np[valid], sigma=1)
            nt_smooth[valid] = gaussian_filter1d(mean_nt[valid], sigma=1)

        # ---- plot ----
        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(bin_centers, nw_smooth, color='tab:blue', linewidth=2.0,
                label='Water O')
        ax.plot(bin_centers, np_smooth, color='tab:orange', linewidth=2.0,
                label='Protein')
        ax.plot(bin_centers, nt_smooth, color='black', linewidth=2.0,
                linestyle='--', label='Total')
        # Raw data as faint lines
        ax.plot(bin_centers, mean_nw, color='tab:blue', linewidth=0.5, alpha=0.3)
        ax.plot(bin_centers, mean_np, color='tab:orange', linewidth=0.5, alpha=0.3)

        ax.axvline(mean_plow,  color='steelblue', linestyle='--', linewidth=1.0,
                   label='plow  %.1f Å' % mean_plow)
        ax.axvline(mean_pmid,  color='steelblue', linestyle='-',  linewidth=1.0,
                   label='pmid  %.1f Å' % mean_pmid)
        ax.axvline(mean_phigh, color='steelblue', linestyle=':',  linewidth=1.0,
                   label='phigh %.1f Å' % mean_phigh)

        ax.set_xlabel('Z position (Å)')
        ax.set_ylabel('Coordination number')
        ax.set_title('%s  —  coordination of permeating ions (n=%d)' % (resname, len(keys)))
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(mean_plow - 15, mean_phigh + 15)

        plt.tight_layout()
        outpath = '%s_%s_coordination.png' % (outbase, resname)
        plt.savefig(outpath, dpi=150)
        plt.close(fig)
        print('Saved: %s' % outpath)

        # ---- XVG output ----
        xvg_path = '%s_%s_coordination.xvg' % (outbase, resname)
        with open(xvg_path, 'w') as fxvg:
            fxvg.write('# Coordination numbers for permeating %s ions\n' % resname)
            fxvg.write('# Averaged over %d permeating ions\n' % len(keys))
            fxvg.write('@    title "Coordination: %s"\n' % resname)
            fxvg.write('@    xaxis label "Z position (A)"\n')
            fxvg.write('@    yaxis label "Coordination number"\n')
            fxvg.write('@TYPE xy\n')
            fxvg.write('@ s0 legend "water O (raw)"\n')
            fxvg.write('@ s1 legend "protein (raw)"\n')
            fxvg.write('@ s2 legend "total (raw)"\n')
            fxvg.write('@ s3 legend "water O (smooth)"\n')
            fxvg.write('@ s4 legend "protein (smooth)"\n')
            fxvg.write('@ s5 legend "total (smooth)"\n')
            for j in range(n_bins):
                nw_r = mean_nw[j] if np.isfinite(mean_nw[j]) else 0.0
                np_r = mean_np[j] if np.isfinite(mean_np[j]) else 0.0
                nt_r = mean_nt[j] if np.isfinite(mean_nt[j]) else 0.0
                nw_s = nw_smooth[j] if np.isfinite(nw_smooth[j]) else 0.0
                np_s = np_smooth[j] if np.isfinite(np_smooth[j]) else 0.0
                nt_s = nt_smooth[j] if np.isfinite(nt_smooth[j]) else 0.0
                fxvg.write('%10.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f  %8.4f\n' %
                           (bin_centers[j], nw_r, np_r, nt_r, nw_s, np_s, nt_s))
        print('Saved: %s' % xvg_path)


def analyze_pore_occupancy(ionsCyl, ion_zpos, ion_xpos, ion_ypos,
                           ion_resname, trjTimes, plows, phighs,
                           cyl_zlows, cyl_zhighs, outbase):
    """Analyze ion occupancy of the pore cylinder vs time.

    Only ions that are inside the cylinder AND within the effective z-range
    are counted.  The effective z-range per frame is the intersection of the
    membrane boundaries (plow, phigh) and the cylinder group z-extent
    (cyl_zlow, cyl_zhigh) — whichever is tighter.

    Produces:
      - <base>_<resname>_pore_occupancy.png  : occupancy time series + histogram
      - <base>_<resname>_pore_occupancy.dat  : frame-by-frame occupancy data
      - <base>_<resname>_ion_ion_distances.png : histogram of pairwise distances
      - <base>_<resname>_ion_ion_distances.dat : raw distance data
    """
    from itertools import combinations

    # --- determine number of frames and ion keys ---
    keys = sorted(ionsCyl.keys())
    n_frames = len(ionsCyl[keys[0]])
    times = np.array(trjTimes[:n_frames])

    # --- per-frame occupancy (cylinder + effective z-range) ---
    occupancy = np.zeros(n_frames, dtype=int)
    for fi in range(n_frames):
        zlo = max(plows[fi], cyl_zlows[fi])
        zhi = min(phighs[fi], cyl_zhighs[fi])
        for k in keys:
            if ionsCyl[k][fi] == 1:
                z = ion_zpos[k][fi]
                if zlo <= z <= zhi:
                    occupancy[fi] += 1

    # --- per-frame pairwise distances for frames with >= 2 ions ---
    all_distances = []
    for fi in range(n_frames):
        zlo = max(plows[fi], cyl_zlows[fi])
        zhi = min(phighs[fi], cyl_zhighs[fi])
        # collect positions of ions inside cylinder AND effective z-range
        in_cyl = []
        for k in keys:
            if ionsCyl[k][fi] == 1:
                z = ion_zpos[k][fi]
                if zlo <= z <= zhi:
                    in_cyl.append((ion_xpos[k][fi], ion_ypos[k][fi], z))
        if len(in_cyl) >= 2:
            for (a, b) in combinations(range(len(in_cyl)), 2):
                dx = in_cyl[a][0] - in_cyl[b][0]
                dy = in_cyl[a][1] - in_cyl[b][1]
                dz = in_cyl[a][2] - in_cyl[b][2]
                all_distances.append(np.sqrt(dx*dx + dy*dy + dz*dz))

    all_distances = np.array(all_distances)

    # --- determine resname (use first ion key) ---
    resname = ion_resname[keys[0]]

    # ===== PLOT 1: occupancy time series + histogram =====
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4),
                                    gridspec_kw={'width_ratios': [3, 1]})

    ax1.plot(times, occupancy, linewidth=0.3, color='steelblue', alpha=0.7)
    # running average for clarity
    window = max(1, n_frames // 200)
    if window > 1:
        occ_smooth = np.convolve(occupancy, np.ones(window)/window, mode='same')
        ax1.plot(times, occ_smooth, linewidth=1.5, color='darkblue',
                 label='running avg (w=%d)' % window)
        ax1.legend(fontsize=8)
    ax1.set_xlabel('Time (ns)')
    ax1.set_ylabel('Number of %s in pore' % resname)
    ax1.set_title('Pore occupancy vs time')
    ax1.set_ylim(bottom=-0.2)

    # histogram
    max_occ = int(occupancy.max())
    bins_hist = np.arange(-0.5, max_occ + 1.5, 1)
    counts, _ = np.histogram(occupancy, bins=bins_hist)
    fractions = counts / n_frames
    ax2.bar(range(max_occ + 1), fractions, color='steelblue', edgecolor='black')
    ax2.set_xlabel('Number of %s in pore' % resname)
    ax2.set_ylabel('Fraction of frames')
    ax2.set_title('Occupancy distribution')
    if max_occ <= 15:
        ax2.set_xticks(range(max_occ + 1))
    else:
        step = max(1, (max_occ + 1) // 10)
        ax2.set_xticks(range(0, max_occ + 1, step))

    # annotate mean
    mean_occ = occupancy.mean()
    ax2.axvline(mean_occ, color='red', linestyle='--', linewidth=1.5)
    ax2.text(mean_occ + 0.1, fractions.max() * 0.9, 'mean=%.2f' % mean_occ,
             color='red', fontsize=9)

    plt.tight_layout()
    occ_png = '%s_%s_pore_occupancy.png' % (outbase, resname)
    plt.savefig(occ_png, dpi=150)
    plt.close()
    print('Saved: %s' % occ_png)

    # write occupancy dat file
    occ_dat = '%s_%s_pore_occupancy.dat' % (outbase, resname)
    with open(occ_dat, 'w') as f:
        f.write('# Pore occupancy: number of %s ions inside cylinder per frame\n' % resname)
        f.write('# z-range: max(plow, cyl_zmin) to min(phigh, cyl_zmax)\n')
        f.write('# Mean effective z-range: %.1f - %.1f A\n'
                % (np.mean([max(plows[i], cyl_zlows[i]) for i in range(n_frames)]),
                   np.mean([min(phighs[i], cyl_zhighs[i]) for i in range(n_frames)])))
        f.write('# Mean occupancy: %.4f\n' % mean_occ)
        f.write('# Occupancy distribution (n_ions  fraction_of_frames):\n')
        for n_ions in range(max_occ + 1):
            f.write('#   %d  %.6f\n' % (n_ions, fractions[n_ions] if n_ions < len(fractions) else 0.0))
        f.write('#\n')
        f.write('# %10s  %6s\n' % ('time_ns', 'n_ions'))
        for fi in range(n_frames):
            f.write('%12.4f  %6d\n' % (times[fi], occupancy[fi]))
    print('Saved: %s' % occ_dat)

    # ===== PLOT 2: ion-ion distance histogram =====
    if len(all_distances) > 0:
        fig, ax = plt.subplots(figsize=(7, 4))
        bin_width = 0.5  # Angstroms
        d_max = min(all_distances.max(), 50.0)  # cap at 50 A for clarity
        bins_d = np.arange(0, d_max + bin_width, bin_width)
        ax.hist(all_distances[all_distances <= d_max], bins=bins_d,
                color='coral', edgecolor='black', linewidth=0.3, density=True)
        ax.set_xlabel('Pairwise ion-ion distance (A)')
        ax.set_ylabel('Probability density')
        ax.set_title('Ion-ion distances inside pore (%s, %d pairs from %d frames with >= 2 ions)'
                      % (resname, len(all_distances),
                         int(np.sum(occupancy >= 2))))
        # annotate median
        median_d = np.median(all_distances)
        ax.axvline(median_d, color='red', linestyle='--', linewidth=1.5)
        ax.text(median_d + 0.3, ax.get_ylim()[1] * 0.9,
                'median=%.1f A' % median_d, color='red', fontsize=9)
        plt.tight_layout()
        dist_png = '%s_%s_ion_ion_distances.png' % (outbase, resname)
        plt.savefig(dist_png, dpi=150)
        plt.close()
        print('Saved: %s' % dist_png)

        # write distance dat file
        dist_dat = '%s_%s_ion_ion_distances.dat' % (outbase, resname)
        with open(dist_dat, 'w') as f:
            f.write('# Pairwise ion-ion distances (A) for %s ions inside the pore cylinder\n' % resname)
            f.write('# Total pairs: %d  (from %d frames with >= 2 ions)\n'
                    % (len(all_distances), int(np.sum(occupancy >= 2))))
            f.write('# Median distance: %.4f A\n' % median_d)
            f.write('# Mean distance: %.4f A\n' % np.mean(all_distances))
            f.write('#\n')
            # write histogram
            counts_d, edges_d = np.histogram(all_distances, bins=np.arange(0, all_distances.max() + bin_width, bin_width))
            f.write('# Histogram (bin_center_A  count  probability_density)\n')
            total = float(len(all_distances))
            for i in range(len(counts_d)):
                center = (edges_d[i] + edges_d[i+1]) / 2.0
                density = counts_d[i] / (total * bin_width) if total > 0 else 0.0
                f.write('%10.4f  %8d  %12.6f\n' % (center, counts_d[i], density))
        print('Saved: %s' % dist_dat)
    else:
        print('No frames with >= 2 ions in pore — skipping ion-ion distance analysis')

    return occupancy, all_distances


def main(argv):

    ################################
    # scheme for region counting
    ################################
    ################################
    #
    # (bulk) 4
    # PPPP  phigh  PPPPPPPPPPPPPPPPP
    #  3
    #------  middle  -------------
    #  2
    # PPPP  plow   PPPPPPPPPPPPPPPP
    # (bulk) 1
    #
    ##################################
    ##################################
    # when counting within the cylinder #
    # only regions 2 and 3 need to be within its radius #
    #####################################################

    parser = argparse.ArgumentParser(
        description='Count ions permeating through a membrane.\n'
                    'Takes into account all ion crossings from one side of the bilayer to the other.\n'
                    'Considers only one ion species at a time.\n'
                    'To count several ion species run script several times separately.\n\n'
                    'Optionally, can consider ions passing through a cylinder.\n'
                    'The cylinder is constructed from the atoms in a separate index group, e.g. protein atoms.\n'
                    'These atoms are used to define height and radius of the cylinder parallel to z-axis.\n'
                    'For the -cyl option to work properly, the system needs to be superimposed accordingly.'
    )
    parser.add_argument('-s', required=True, metavar='init.pdb',
                        help='Input structure file (pdb/gro)')
    parser.add_argument('-f', required=True, metavar='traj.xtc',
                        help='Input trajectory file (xtc/trr)')
    parser.add_argument('-n', default='index.ndx', metavar='index.ndx',
                        help='Input index file')
    parser.add_argument('-o', default='out.dat', metavar='out.dat',
                        help='Output counts file')
    parser.add_argument('-cyl', action='store_true',
                        help='Count ions passing through a cylinder defined by a protein group')
    parser.add_argument('-coord', action='store_true',
                        help='Compute coordination numbers for permeating ions (second trajectory pass, slow)')

    args = parser.parse_args(argv[1:])

    strFile = args.s
    trjFile = args.f
    ndxFile = args.n
    predFile = args.o
    bCyl = args.cyl
    bCoord = args.coord

    # Load structure and trajectory with MDAnalysis
    u = mda.Universe(strFile, trjFile)

    #### index ####
    ndxPhos = []
    ndxIons = []
    ndxCyl = []
    (ndxPhos, ndxPhosNum) = select_ndx(ndxFile, message='Select phosphate group\n')
    (ndxIons, ndxIonsNum) = select_ndx(ndxFile, message='Select ion group\n')
    if bCyl:
        (ndxCyl, ndxCylNum) = select_ndx(ndxFile, message='Select protein group for cylinder\n')
    ionNum = len(ndxIons)

    ### initialize ion structures ###
    ions = {}
    ionsCyl = {}  # track whether within a cylinder (optional): 0 - not in cylinder, 1 - in cylinder
    ion_xpos = {}  # x-coordinate time series for each ion
    ion_ypos = {}  # y-coordinate time series for each ion
    ion_zpos = {}  # z-coordinate time series for each ion
    ion_resname = {}  # residue name for each ion atom
    for i in range(ionNum):
        idx = ndxIons[i]
        ions[idx] = []
        ionsCyl[idx] = []
        ion_xpos[idx] = []
        ion_ypos[idx] = []
        ion_zpos[idx] = []
        ion_resname[idx] = u.atoms[int(idx)].resname

    ##############################
    #### reading trajectory ####
    frnum = 0
    trjTimes = []
    plows = []
    pmids = []
    phighs = []
    box_xs = []
    box_ys = []
    box_zs = []
    cyl_zlows = []
    cyl_zhighs = []

    for ts in u.trajectory:
        # u.atoms.positions is shape (N, 3) in Angstroms
        positions = u.atoms.positions

        # identify 1-2-3-4 regions and, if needed, cylinder (for each frame)
        phosCrd = positions[ndxPhos]
        plow, pmid, phigh = get_limits(phosCrd)
        plows.append(plow)
        pmids.append(pmid)
        phighs.append(phigh)
        if bCyl:
            cylCrd = positions[ndxCyl]
            cylCenter, cylRad = get_cylinder(cylCrd)
            cyl_zlows.append(float(cylCrd[:, 2].min()))
            cyl_zhighs.append(float(cylCrd[:, 2].max()))

        # analyze ions
        for i in range(ionNum):
            pos = positions[ndxIons[i]]
            x = pos[0]
            y = pos[1]
            z = pos[2]
            # identify the range
            r = identify_range(z, plow, pmid, phigh)
            ions[ndxIons[i]].append(r)
            ion_xpos[ndxIons[i]].append(x)
            ion_ypos[ndxIons[i]].append(y)
            ion_zpos[ndxIons[i]].append(z)
            # track if in cylinder
            if bCyl:
                c = identify_if_in_cylinder(x, y, cylCenter, cylRad * cylRad)
                ionsCyl[ndxIons[i]].append(c)

        # store trj time (converted ps -> ns) and box dimensions
        trjTimes.append(ts.time / 1000.0)
        box_xs.append(ts.dimensions[0])
        box_ys.append(ts.dimensions[1])
        box_zs.append(ts.dimensions[2])

        sys.stdout.write('Step: %d       Time: %f                  Frame: %d\r' % (ts.frame, ts.time, frnum))
        sys.stdout.flush()
        frnum += 1

    sys.stdout.write('\n')
    mean_box_x = float(np.mean(box_xs))
    mean_box_y = float(np.mean(box_ys))
    mean_box_z = float(np.mean(box_zs))

    # track permeations
    (transitionsUp, transitionsDown, transitionsUpTimes, transitionsDownTimes,
     transitionsUpTransit, transitionsDownTransit,
     totalUp, totalDown,
     CYLtransitionsUp, CYLtransitionsDown, CYLtransitionsUpTimes, CYLtransitionsDownTimes,
     CYLtotalUp, CYLtotalDown) = track_permeations(ions, trjTimes, bCyl, ionsCyl)

    ######################################################
    ######################################################
    ###################### output ########################
    ######################################################
    ######################################################
    with open(predFile, 'w') as fp:
        # summary
        fp.write('-----------------------------------\n')
        fp.write('-----------------------------------\n')
        fp.write('Total simulation time: {0} ns\n'.format(trjTimes[-1]))
        fp.write('Total permeations up: {0}\n'.format(totalUp))
        fp.write('Total permeations down: {0}\n'.format(totalDown))
        fp.write('-----------------------------------\n')
        currentUp = totalUp * 1.602176634 / trjTimes[-1] * 100.0  # pA
        currentDown = totalDown * 1.602176634 / trjTimes[-1] * 100.0  # pA
        fp.write('Current up: {0} pA\n'.format(round(currentUp, 5)))
        fp.write('Current down: {0} pA\n'.format(round(currentDown, 5)))
        fp.write('-----------------------------------\n')
        fp.write('-----------------------------------\n')

        if bCyl:
            fp.write('\n**************************\n')
            fp.write('******** Cylinder ********\n')
            fp.write('Cylinder permeations up: {0}\n'.format(CYLtotalUp))
            fp.write('Cylinder permeations down: {0}\n'.format(CYLtotalDown))
            fp.write('---------\n')
            currentUp = CYLtotalUp * 1.602176634 / trjTimes[-1] * 100.0  # pA
            currentDown = CYLtotalDown * 1.602176634 / trjTimes[-1] * 100.0  # pA
            fp.write('Cylinder current up: {0} pA\n'.format(round(currentUp, 5)))
            fp.write('Cylinder current down: {0} pA\n'.format(round(currentDown, 5)))
            fp.write('******** Cylinder ********\n')
            fp.write('**************************\n\n')

        fp.write('-----------\n')
        fp.write('--Details--\n')
        fp.write('-----------\n')
        for direction, trans, trans_times, trans_transit, total, prefix in [
            ('Up',   transitionsUp,   transitionsUpTimes,   transitionsUpTransit,   totalUp,   'Up'),
            ('Down', transitionsDown, transitionsDownTimes, transitionsDownTransit, totalDown, 'Down'),
        ]:
            if total > 0:
                for key in trans:
                    if trans[key] > 0:
                        fp.write('{0}: ion {1} had {2} permeations at times (ns):'.format(
                            prefix, key, trans[key]))
                        for t in trans_times[key]:
                            fp.write(' {0}'.format(t))
                        fp.write('\n')
                        fp.write('{0}: ion {1} transit times (ns):'.format(prefix, key))
                        for tt in trans_transit[key]:
                            fp.write(' {0}'.format(round(tt, 1)))
                        fp.write('\n')

        if bCyl:
            fp.write('\n***********\n')
            fp.write('**Details**\n')
            fp.write('***********\n')
            for direction, trans, trans_times, total, prefix in [
                ('Up',   CYLtransitionsUp,   CYLtransitionsUpTimes,   CYLtotalUp,   'Cylinder up'),
                ('Down', CYLtransitionsDown, CYLtransitionsDownTimes, CYLtotalDown, 'Cylinder down'),
            ]:
                if total > 0:
                    for key in trans:
                        if trans[key] > 0:
                            fp.write('{0}: ion {1} had {2} permeations at times (ns):'.format(
                                prefix, key, trans[key]))
                            for t in trans_times[key]:
                                fp.write(' {0}'.format(t))
                            fp.write('\n')

    ######################################################
    ###################### plots ########################
    ######################################################
    outbase = os.path.splitext(predFile)[0]
    plot_permeations(transitionsUpTimes, transitionsDownTimes, trjTimes, outbase)
    plot_ion_positions(ion_zpos, ion_resname, trjTimes, plows, pmids, phighs,
                       transitionsUpTimes, transitionsDownTimes, outbase)
    plot_density_pmf(ion_zpos, ion_resname, plows, pmids, phighs, mean_box_z, outbase)
    plot_density_pmf_local(ion_zpos, ion_resname, plows, pmids, phighs, mean_box_z, outbase)
    if bCyl:
        plot_density_pmf_pore(ion_zpos, ionsCyl, ion_resname, plows, pmids, phighs,
                              mean_box_z, outbase, u=u, ndxCyl=ndxCyl)
        analyze_pore_occupancy(ionsCyl, ion_zpos, ion_xpos, ion_ypos,
                               ion_resname, trjTimes, plows, phighs,
                               cyl_zlows, cyl_zhighs, outbase)
    plot_density_pmf_2d(ion_xpos, ion_ypos, ion_zpos, ion_resname,
                        mean_box_x, mean_box_y, mean_box_z,
                        plows, phighs, outbase)

    # ---- coordination analysis for permeating ions (optional, slow) ----
    if bCoord:
        permeating = set()
        for key, times in transitionsUpTimes.items():
            if times:
                permeating.add(key)
        for key, times in transitionsDownTimes.items():
            if times:
                permeating.add(key)
        if permeating:
            coord_data = compute_coordination(u, permeating, ion_resname,
                                              ndxCyl=ndxCyl if bCyl else None)
            plot_coordination(coord_data, ion_resname, plows, pmids, phighs,
                              mean_box_z, outbase)
        else:
            print('No permeating ions — skipping coordination analysis')


if __name__ == '__main__':
    main(sys.argv)
