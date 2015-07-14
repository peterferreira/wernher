import os
from urllib import request

msgfmt = '{body} {maptype} {zoom} {col} {row}'
dfmt = '{body}/{maptype}/{zoom}'
ffmt = '{body}/{maptype}/{zoom}/{col}_{row}.png'
urlfmt = 'http://tiles.kerbalmaps.com/{body}/{maptype}/{zoom}/{col}/{row}.png'

bodies = '''\
    moho
    eve
    gilly
    kerbin
    mun
    minmus
    duna
    ike
    dres
    laythe
    vall
    tylo
    bop
    pol
    eeloo
'''.split()

maptypes = '''\
    sat
    color
    slope
    biome
'''.split()


# change to range(6) to get *all* maps (about 2GB of data total)
zooms = range(3)


def cols(zoom):
    return range(2**(zoom+1))

def rows(zoom):
    return range(2**zoom)


fmt = {}
for maptype in maptypes:
    fmt['maptype'] = maptype
    for body in bodies:
        fmt['body'] = body
        for zoom in zooms:
            fmt['zoom'] = str(zoom)

            dname = dfmt.format(**fmt)
            if not os.path.exists(dname):
                os.makedirs(dname)

            for col in cols(zoom):
                fmt['col'] = str(int(col))
                for row in rows(zoom):
                    fmt['row'] = str(int(row))

                    fname = ffmt.format(**fmt)

                    if not os.path.exists(fname):

                        print(msgfmt.format(**fmt))

                        try:
                            url = urlfmt.format(**fmt)
                            req = request.urlopen(url)

                            with open(fname,'wb') as fout:
                                fout.write(req.read())
                        except:
                            print('    failed!')
