
__doc__ = """
count-ions-mdanalysis.py  —  Ion permeation counter for GROMACS MD trajectories
=================================================================================

Usage
-----
    python count-ions-mdanalysis.py -s <structure> -f <trajectory> -n <index> -o <output> [-cyl]

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
    python count-ions-mdanalysis.py -s md.gro -f md.xtc -n index.ndx -o results/ions.dat
"""

import sys
import os
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import MDAnalysis as mda


mainchain = ['N', 'NH', 'C', 'O', 'CA', 'OC1', 'OC2', 'OX', 'H', 'HA', 'HA1', 'HA2']
bckb = ['N', 'C', 'CA', 'NH']


def track_permeations(ions, trjTimes, bCyl=False, ionsCyl={}):
    transitionsUp = {}
    transitionsDown = {}
    transitionsUpTimes = {}
    transitionsDownTimes = {}
    transitionsUpTransit = {}    # transit times (ps) for each upward permeation
    transitionsDownTransit = {}  # transit times (ps) for each downward permeation
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
            if bCyl == True:
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

            # WARNING: jump over two positions
            elif (state1 == 2 and i == 4) or \
                 (state1 == 4 and i == 2) or \
                 (state1 == 3 and i == 1) or \
                 (state1 == 1 and i == 3):
                print('WARNING: large ion jump, your trj output frequency is too low')
                state1 = 5
            ######## cylinder #######
            if bCyl == True:
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
    r = 0
    if z < plow:
        r = 1
    elif z >= plow and z < pmid:
        r = 2
    elif z >= pmid and z < phigh:
        r = 3
    elif z >= phigh:
        r = 4
    return r


def identify_if_in_cylinder(x, y, cylCenter, cylRad):
    d = np.sqrt(np.power(x - cylCenter[0], 2.0) + np.power(y - cylCenter[1], 2.0))
    if d <= cylRad:
        return 1
    return 0


def get_cylinder(crd):
    # return the [x,y] coordinates of the center point and radius
    center = [0.0, 0.0]
    rad = 0.0

    # 1. get the center
    count = 0
    for c in crd:
        center[0] += c[0]
        center[1] += c[1]
        count += 1
    center[0] /= float(count)
    center[1] /= float(count)

    # 2. get the radius (as a maximal distance from the center)
    for c in crd:
        d = np.sqrt(np.power(c[0] - center[0], 2.0) + np.power(c[1] - center[1], 2.0))
        if d > rad:
            rad = d

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

    if message == False:
        sys.stdout.write('Select a group for analysis:\n')
    else:
        sys.stdout.write(message + '\n')

    ndxNum = -1
    while ndxNum == -1:
        ndxNum = input()
        if ndxNum.isdigit() == False:
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
        #ax.legend(loc='upper right', fontsize=7, ncol=max(1, len(perm_keys) // 10 + 1))
        ax.grid(True, alpha=0.2)

    plt.tight_layout()
    outpath = outbase + '_positions.png'
    plt.savefig(outpath, dpi=150)
    plt.close(fig)
    print('Saved: %s' % outpath)


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

    args = parser.parse_args(argv[1:])

    strFile = args.s
    trjFile = args.f
    ndxFile = args.n
    predFile = args.o
    bCyl = args.cyl

    # Load structure and trajectory with MDAnalysis
    u = mda.Universe(strFile, trjFile)

    #### index ####
    ndxPhos = []
    ndxIons = []
    ndxCyl = []
    (ndxPhos, ndxPhosNum) = select_ndx(ndxFile, message='Select phosphate group\n')
    (ndxIons, ndxIonsNum) = select_ndx(ndxFile, message='Select ion group\n')
    if bCyl == True:
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
    plow = 0.0
    pmid = 0.0
    phigh = 0.0
    plows = []
    pmids = []
    phighs = []
    box_xs = []
    box_ys = []
    box_zs = []

    for ts in u.trajectory:
        # u.atoms.positions is shape (N, 3) in Angstroms
        positions = u.atoms.positions

        # identify 1-2-3-4 regions and, if needed, cylinder (for each frame)
        phosCrd = positions[ndxPhos]
        plow, pmid, phigh = get_limits(phosCrd)
        plows.append(plow)
        pmids.append(pmid)
        phighs.append(phigh)
        if bCyl == True:
            cylCrd = positions[ndxCyl]
            cylCenter, cylRad = get_cylinder(cylCrd)

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
            if bCyl == True:
                c = identify_if_in_cylinder(x, y, cylCenter, cylRad)
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
    fp = open(predFile, 'w')
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

    #########################
    ##### cylinder ##########
    #########################
    if bCyl == True:
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
    if totalUp > 0:
        for key in transitionsUp.keys():
            if transitionsUp[key] > 0:
                fp.write('Up: ion {0} had {1} permeations at times (ns):'.format(key, transitionsUp[key]))
                for i in range(transitionsUp[key]):
                    fp.write(' {0}'.format(transitionsUpTimes[key][i]))
                fp.write('\n')
                fp.write('Up: ion {0} transit times (ns):'.format(key))
                for i in range(len(transitionsUpTransit[key])):
                    fp.write(' {0}'.format(round(transitionsUpTransit[key][i], 1)))
                fp.write('\n')
    if totalDown > 0:
        for key in transitionsDown.keys():
            if transitionsDown[key] > 0:
                fp.write('Down: ion {0} had {1} permeations at times (ns):'.format(key, transitionsDown[key]))
                for i in range(transitionsDown[key]):
                    fp.write(' {0}'.format(transitionsDownTimes[key][i]))
                fp.write('\n')
                fp.write('Down: ion {0} transit times (ns):'.format(key))
                for i in range(len(transitionsDownTransit[key])):
                    fp.write(' {0}'.format(round(transitionsDownTransit[key][i], 1)))
                fp.write('\n')

    #########################
    ##### cylinder ##########
    #########################
    if bCyl == True:
        fp.write('\n***********\n')
        fp.write('**Details**\n')
        fp.write('***********\n')
        if CYLtotalUp > 0:
            for key in CYLtransitionsUp.keys():
                if CYLtransitionsUp[key] > 0:
                    fp.write('Cylinder up: ion {0} had {1} permeations at times (ns):'.format(key, CYLtransitionsUp[key]))
                    for i in range(CYLtransitionsUp[key]):
                        fp.write(' {0}'.format(CYLtransitionsUpTimes[key][i]))
                    fp.write('\n')
        if CYLtotalDown > 0:
            for key in CYLtransitionsDown.keys():
                if CYLtransitionsDown[key] > 0:
                    fp.write('Cylinder down: ion {0} had {1} permeations at times (ns):'.format(key, CYLtransitionsDown[key]))
                    for i in range(CYLtransitionsDown[key]):
                        fp.write(' {0}'.format(CYLtransitionsDownTimes[key][i]))
                    fp.write('\n')
    fp.close()

    ######################################################
    ###################### plots ########################
    ######################################################
    outbase = os.path.splitext(predFile)[0]
    plot_permeations(transitionsUpTimes, transitionsDownTimes, trjTimes, outbase)
    plot_ion_positions(ion_zpos, ion_resname, trjTimes, plows, pmids, phighs,
                       transitionsUpTimes, transitionsDownTimes, outbase)
    plot_density_pmf(ion_zpos, ion_resname, plows, pmids, phighs, mean_box_z, outbase)
    plot_density_pmf_2d(ion_xpos, ion_ypos, ion_zpos, ion_resname,
                        mean_box_x, mean_box_y, mean_box_z,
                        plows, phighs, outbase)


if __name__ == '__main__':
    main(sys.argv)
