import matplotlib.pyplot as plt
import numpy as np
import matplotlib.collections as mcoll
import matplotlib.path as mpath

def colorline(ax, x, y, z, **kwargs):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, **kwargs)

    ax.add_collection(lc)
    if ax.get_autoscale_on():
        ax.autoscale_view()

    return lc


def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


if __name__ == '__main__':
    from matplotlib.colors import Normalize

    N = 10
    np.random.seed(101)
    x = np.random.rand(N)
    y = np.random.rand(N)
    fig, ax = plt.subplots()

    path = mpath.Path(np.column_stack([x, y]))
    verts = path.interpolated(steps=3).vertices
    x, y = verts[:, 0], verts[:, 1]
    z = np.linspace(0, 100, len(x))
    colorline(ax, x, y, z, cmap=plt.get_cmap('jet'), linewidth=2)

    fig,ax = plt.subplots()

    x = np.linspace(0,1,50)
    y = 2 - x**2
    z = np.linspace(0,25,len(x))
    colorline(ax,x,y,z,norm=Normalize(vmin=0,vmax=50))

    x = np.linspace(0,1,50)
    y = 1 - x**2
    z = np.linspace(0,50,len(x))
    colorline(ax,x,y,z,norm=Normalize(vmin=0,vmax=50))


    plt.show()
