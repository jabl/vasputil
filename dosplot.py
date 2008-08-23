#!/usr/bin/pythoon

# Plot DOS from vasp
# First use split_dos to create DOSx files.

from pylab import *
#from matplotlib.numerix import sin, cos, exp, pi, arange

def makedosplot(dos, s=True, p=True, d=True):
    t = arange(0.0, 2.0, 0.01)
    s1 = sin(2*pi*t)
    s2 = exp(-t)
    s3 = sin(2*pi*t)*exp(-t)
    s4 = sin(2*pi*t)*cos(4*pi*t)
     
    fig = figure()
    t = arange(0.0, 2.0, 0.01)
    
    yprops = dict(rotation=0,
               horizontalalignment='right',
               verticalalignment='center',
               x=-0.01)
    
    axprops = dict(yticks=[])
    
    ax1 =fig.add_axes([0.1, 0.7, 0.8, 0.2], **axprops)
    ax1.plot(t, s1)
    ax1.set_ylabel('S1', **yprops)
    
    axprops['sharex'] = ax1
    axprops['sharey'] = ax1
    # force x axes to remain in register, even with toolbar navigation
    ax2 = fig.add_axes([0.1, 0.5, 0.8, 0.2], **axprops)
    
    ax2.plot(t, s2)
    ax2.set_ylabel('S2', **yprops)
    
    ax3 = fig.add_axes([0.1, 0.3, 0.8, 0.2], **axprops)
    ax3.plot(t, s4)
    ax3.set_ylabel('S3', **yprops)
    
    ax4 = fig.add_axes([0.1, 0.1, 0.8, 0.2], **axprops)
    ax4.plot(t, s4)
    ax4.set_ylabel('S4', **yprops)
    
    # turn off x ticklabels for all but the lower axes
    for ax in ax1, ax2, ax3:
        setp(ax.get_xticklabels(), visible=False)
    
    show()

def readdos(fname):
    dos = load(fname)

def showdp():
    xlabel("E-E$_\mathrm{f}$ (eV)", size=16)
    #ylabel("LDOS", size=16)
    #f = figure(0)
    #text(0., 0.4, "LDOS", rotation='vertical', transform = f.axes.transAxes)
    #figtext(0., 0.4, "LDOS", rotation='vertical', figure=f)
    figtext(0.03, 0.45, "LDOS", rotation='vertical', size=16)
    loc, lab = xticks()
#    set(l, size = 16)
    lab.set_size = 16
    loc, lab = yticks()
    lab.set_size=16
    legend()
    subplots_adjust(hspace=0.0)
    show()
