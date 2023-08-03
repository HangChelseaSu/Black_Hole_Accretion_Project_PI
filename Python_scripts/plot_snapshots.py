from pathlib import Path
import sys
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('ggplot')
from discopy import util

def R_centroid(rm, rp):
    
    rc = 2*(rm*rm + rm*rp + rp*rp) / (3*(rm + rp))

    return rc

def plot_snap(snap_path, Rmd_max=1.0, Nr_md=30):

    # Get just the filename
    name = snap_path.stem


    # Load the snapshot data ("_" data are redundant, already returned by Qrz)
    t, r, z, Qrz, rf, zf, planetDat = util.loadSnapshotRZ(snap_path)
    t, Qarr, _, _, _ = util.loadSnapshotArr(snap_path)

    # Reshape the raw Qarr into a more usable shape
    # [Num_planets, Num_R_minidisk, Num_quantities+1]
    Qmd = Qarr.reshape(2, Nr_md, 20)

    # Grab the planet Masses
    M0 = planetDat[0, 0]
    M1 = planetDat[1, 0]

    # These are the cell "centers" (really, centroid) in radius.
    R = R_centroid(rf[:-1], rf[1:])

    # These are the bin edges for the radial bins for the minidisks
    Rf_md = np.linspace(0.0, Rmd_max, Nr_md+1)

    # And the minidisk radial centroids
    R_md = R_centroid(Rf_md[:-1], Rf_md[1:])

    # Circumbinary density, radial velocity, angular velocity
    sig_cb = Qrz[0, :, 0]
    vr_cb = Qrz[0, :, 1] / Qrz[0, :, 0]
    om_cb = Qrz[0, :, 2] / (R**2 * Qrz[0, :, 0])

    # Minidisk density, velocities for BH 0
    sig_0 = Qmd[0, :, 1] / Qmd[0, :, 0]
    vr_0 = Qmd[0, :, 2] / Qmd[0, :, 1]
    om_0 = Qmd[0, :, 3] / (R_md**2 * Qmd[0, :, 1])

    # Minidisk density, velocities for BH 1
    sig_1 = Qmd[1, :, 1] / Qmd[1, :, 0]
    vr_1 = Qmd[1, :, 2] / Qmd[1, :, 1]
    om_1 = Qmd[1, :, 3] / (R_md**2 * Qmd[1, :, 1])

    M = M0 + M1
    Om_kep_cb = np.sqrt(M / R**3)
    Om_kep_md0 = np.sqrt(M0 / R_md**3)
    Om_kep_md1 = np.sqrt(M1 / R_md**3)

    # Make a plot
    fig, ax = plt.subplots(3, 3, figsize=(9, 6))
    
    # Plot the densities
    ax[0, 0].plot(R, sig_cb)
    ax[0, 1].plot(R_md, sig_0)
    ax[0, 2].plot(R_md, sig_1)
    
    # Plot the radial velocities densities
    ax[1, 0].plot(R, vr_cb)
    ax[1, 1].plot(R_md, vr_0)
    ax[1, 2].plot(R_md, vr_1)

    # Plot the keplerian angular velocities
    ax[2, 0].plot(R, Om_kep_cb, alpha=0.5, lw=2, color='grey')
    ax[2, 1].plot(R_md, Om_kep_md0, alpha=0.5, lw=2, color='grey')
    ax[2, 2].plot(R_md, Om_kep_md1, alpha=0.5, lw=2, color='grey')

    # Plot the actual angular velocities
    ax[2, 0].plot(R, om_cb)
    ax[2, 1].plot(R_md, om_0)
    ax[2, 2].plot(R_md, om_1)

    # Column titles!
    ax[0, 0].set_title("CB Disk")
    ax[0, 1].set_title("Minidisk 0")
    ax[0, 2].set_title("Minidisk 1")

    # Make densty plots nice 
    ax[0, 0].set(ylabel=r'$\Sigma$', ylim=(0, 500))
    ax[0, 1].set(ylim=(0, 100))
    ax[0, 2].set(ylim=(0, 100))

    # Make radial velocity plots nice 
    ax[1, 0].set(ylabel=r'$v_r$', ylim=(-1, 1))
    ax[1, 1].set(ylim=(-1, 1))
    ax[1, 2].set(ylim=(-1, 1))
    
    # Make angular velocity plots nice 
    ax[2, 0].set(ylabel=r'$\Omega$', xlabel=r'$r$',
                 ylim=(0.01, 1000), yscale='log')
    ax[2, 1].set(xlabel=r'$r$',
                 ylim=(0.1, 1000), yscale='log')
    ax[2, 2].set(xlabel=r'$r$',
                 ylim=(0.1, 1000), yscale='log')

    # Make sure everything fits nicely
    fig.tight_layout()

    # Save and close!
    figname = name + ".png"
    print("Saving", figname)
    fig.savefig(figname)
    plt.close(fig)



if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Need some snapshot files to nom")

    filenames = [Path(x) for x in sys.argv[1:]]

    for f in filenames:
        plot_snap(f)