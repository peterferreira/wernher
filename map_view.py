import os
import numpy as np
from matplotlib import pyplot, ticker

''' example using basemap:

from mpl_toolkits.basemap import Basemap

bmap = Basemap(resolution=None,rsphere=body.equatorial_radius,ax=ax)
im = bmap.imshow(image)

major_parallels = np.linspace(-90,90,2*3+1)
major_meridians = np.linspace(-180,180,2*4+1)
line_opts = dict(labels=[1,0,0,1], dashes=[3,9])
bmap.drawparallels(major_parallels, **line_opts)
bmap.drawmeridians(major_meridians, **line_opts)
'''

class MapView(object):
    def __init__(self,body,**kwargs):
        self.body = body
        self.maptype = kwargs.pop('maptype','sat')
        self.zoomlevel = kwargs.pop('zoomlevel',3)

    @staticmethod
    def set_ticks(ax):
        opts = dict(nbins=9, steps=[1, 2, 3, 6, 15, 18])
        ax.xaxis.set_major_locator(ticker.MaxNLocator(**opts))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(**opts))
        def _lon_fmt(x,pos):
            def _lon_dir(x):
                if np.isclose(x,0) or np.isclose(abs(x),180):
                    return ''
                elif x < 0:
                    return ' W'
                else:
                    return ' E'
            s = '{x:g}˚{dir}'.format(x=abs(x),dir=_lon_dir(x))
            return s
        def _lat_fmt(x,pos):
            def _lat_dir(x):
                if np.isclose(x,0):
                    return ''
                elif x < 0:
                    return ' S'
                else:
                    return ' N'
            s ='{x:g}˚{dir}'.format(x=abs(x),dir=_lat_dir(x))
            return s
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(_lon_fmt))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(_lat_fmt))

    @staticmethod
    def split_tracks(lat,lon):
        '''assumes eastward motion'''
        tracks = []
        lt,ln = [lat[0]],[lon[0]]
        for i in range(1,len(lon)):
            lt.append(lat[i])
            d1 = abs(lon[i] - lon[i-1])
            d2 = abs((lon[i-1] + 360) - lon[i])
            d3 = abs(lon[i-1] - (lon[i] + 360))
            if d2 < d1:
                ln.append(lon[i]-360)
                tracks.append((np.array(lt),np.array(ln)))
                lt = [lat[i-1],lat[i]]
                ln = [lon[i-1]+360,lon[i]]
            elif d3 < d1:
                ln.append(lon[i]+360)
                tracks.append((np.array(lt),np.array(ln)))
                lt = [lat[i-1],lat[i]]
                ln = [lon[i-1]-360,lon[i]]
            else:
                ln.append(lon[i])
        if len(lt):
            tracks.append((np.array(lt),np.array(ln)))
        return tracks

    def plot_basemap(self,ax):
        ax.imshow(self.image,
            extent=[-180,180,-90,90],
            aspect='equal',
            origin='lower')
        ax.grid(True)
        MapView.set_ticks(ax)
        ax.autoscale(False)

    @staticmethod
    def plot_marker(ax,lat,lon,**kwargs):
        kw = dict(
            color = 'red',
            marker = 'o',
            markersize = 10)
        kw.update(**kwargs)
        return ax.plot(lon,lat,**kw)

    @staticmethod
    def plot_track(ax,lat,lon,**kwargs):
        kw = dict(
            color = 'cyan',
            alpha = 0.5,
            lw = 4)
        kw.update(**kwargs)
        tracks = MapView.split_tracks(lat,lon)
        pts = []
        for lt,ln in tracks:
            pts.append(ax.plot(ln,lt,**kw))
        return pts

    @property
    def image(self):
        bodyname = self.body.name.lower()
        maptype = self.maptype
        zoomlevel = self.zoomlevel

        key = (bodyname,maptype,zoomlevel)
        if hasattr(self,'_images'):
            if key in self._images:
                return self._images[key]
        else:
            self._images = {}

        curdir = os.path.dirname(os.path.realpath(__file__))
        fpath = os.path.join(curdir,
            'map_images',bodyname,maptype,str(zoomlevel))
        ffmt = '{col}_{row}.png'

        images = []

        for col in range(2**(zoomlevel+1)):
            imrow = []
            for row in range(2**zoomlevel):
                fname = os.path.join(fpath,ffmt.format(col=col,row=row))
                data = pyplot.imread(fname)
                data_u8 = (data[::-1,...,...] * np.iinfo(np.uint8).max)\
                            .astype(np.uint8,casting='unsafe')
                imrow.append(data_u8)
                del data
            images.append(np.vstack(imrow))
        self._images[key] = np.hstack(images)
        return self._images[key]

if __name__ == '__main__':
    class CelestialBody(object):
        def __init__(self,name='kerbin',radius=600000):
            self.name = name
            self.equatorial_radius = radius

    body = CelestialBody() #'duna',320000)
    mview = MapView(body)
    fig,ax = pyplot.subplots(figsize=(16,8))
    mview.plot_basemap(ax)
    pyplot.show()
