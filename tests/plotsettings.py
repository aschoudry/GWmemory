from matplotlib import rc

rc('text', usetex=True)
rc('font', family='serif')
rc('font', serif='times')
#rc('font', weight='bolder')
rc('mathtext', default='sf')
rc("lines", markeredgewidth=0.6)
rc("lines", linewidth=2)
rc('axes', labelsize=24) #24
rc("axes", linewidth=0.5) #2)
rc('xtick', labelsize=15)
rc('ytick', labelsize=15)
rc('legend', fontsize=20) #16
rc('xtick.major', pad=8) #8)
rc('ytick.major', pad=8) #8)
rc('xtick.major', size=25) #8)
rc('ytick.major', size=25) #8)
rc('xtick.minor', size=25) #8)
rc('ytick.minor', size=25) #8)

def set_tick_sizes(ax, major, minor):
    for l in ax.get_xticklines() + ax.get_yticklines():
        l.set_markersize(major)
    for tick in ax.xaxis.get_minor_ticks() + ax.yaxis.get_minor_ticks():
        tick.tick1line.set_markersize(minor)
        tick.tick2line.set_markersize(minor)
    ax.xaxis.LABELPAD=10.
    ax.xaxis.OFFSETTEXTPAD=10.

