import cf
import h5py
import cfplot as cfp

h5_infile = 'adios_test.h5'
nc_infile = '../data/u_plev_192x145x17.nc'
ilev = 8

with h5py.File(h5_infile, 'r') as f:
    lon = f.attrs['longitude']
    lat = f.attrs['latitude']
    lev = f.attrs['p']
    print (lon)
    print (lat)
    print (lev)
    data = f["Step5/u_copy_2"]
    print (data.shape)
    cfp.con(f=data[ilev], x=lon, y=lat, ptype=1)

g = cf.read(nc_infile, 'r')[0]
g.dump()
cfp.con(g[0,ilev])
