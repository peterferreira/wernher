import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt


def plot_orbit_3d(o):
    km = 1000
    π = np.pi
    T = o.period
    npoints = 100
    tt = np.linspace(0,T,npoints)
    f,g,df,fg = o.lagrange_coefficients_at_time(tt)
    xx,yy,zz = o.position_at_lagrange_coefficients(f,g)/km

    mpl.rcParams['legend.fontsize'] = 10

    def axisEqual3D(ax):
        extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
        sz = extents[:,1] - extents[:,0]
        centers = np.mean(extents, axis=1)
        maxsize = max(abs(sz))
        r = maxsize/2
        for ctr, dim in zip(centers, 'xyz'):
            getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

    def plot_sphere(ax, position, radius, *args, **kwargs):
        u = np.linspace(0, 2 * np.pi, 360)
        v = np.linspace(0, np.pi, 180)
        x0,y0,z0 = position
        x = x0 + radius * np.outer(np.cos(u), np.sin(v))
        y = y0 + radius * np.outer(np.sin(u), np.sin(v))
        z = z0 + radius * np.outer(np.ones(np.size(u)), np.cos(v))
        ax.plot_surface(x, y, z, *args, **kwargs)


    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(xx,yy,zz,label='orbit',lw=3)
    plot_sphere(ax,(0,0,0),6378, alpha=0.7, color='lightblue')

    x,y,z = o.position_at_epoch/km
    ax.plot([x],[y],[z], marker='o', markersize=10, color='green')

    x,y,z = o.position_at_ascending_node/km
    ax.plot([x],[y],[z], marker='o', markersize=10, color='cyan')

    x,y,z = o.position_at_descending_node/km
    ax.plot([x],[y],[z], marker='o', markersize=10, color='magenta')

    ax.legend()
    axisEqual3D(ax)

    return fig, ax
