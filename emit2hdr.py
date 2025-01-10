import netCDF4 as nc
from osgeo import gdal, osr, ogr
import numpy as np
import sys, string, os, re
from pyproj import Proj
from pyresample import kd_tree, geometry
def emit2hdr(pathname):
    dataset = nc.Dataset(pathname)
    lon=dataset['location/lon'][:].data
    lat=dataset['location/lat'][:].data
    spectrum=dataset['reflectance'][:].data
    w=dataset['sensor_band_parameters']['wavelengths'][:].data/1000
    bp=dataset['sensor_band_parameters']['fwhm'][:].data/1000
    spectrum[spectrum<0]=0
    spectrum[spectrum>1]=1
    spectrum[spectrum==-9999]=0
    name=os.path.splitext(os.path.split(pathname)[1])[0]
    folder=os.path.split(pathname)[0]
    p=Proj('+proj=utm +zone='+int(np.ceil(lon.mean()/6)+30).__str__() +' datum=WGS84 +units=m +no_defs')
    extent=[*p(lon.min(),lat.min()-0.05),*p(lon.max(),lat.max())]
    area_def = geometry.AreaDefinition('areaD', 'custom', 'areaD',
                               '+proj=utm +zone='+int(np.ceil(lon.mean()/6)+30).__str__() +' datum=WGS84 +units=m +no_defs',
                               int((extent[2]-extent[0])/60),int((extent[3]-extent[1])/60),
                               extent)
    swath_def = geometry.SwathDefinition(lons=lon, lats=lat)
    result = kd_tree.resample_nearest(swath_def,     (spectrum*10000).astype('uint16')  ,
                                area_def, radius_of_influence=90)
    rows,columns,samples=result.shape
    driver = gdal.GetDriverByName('ENVI')
    raster = driver.Create(folder+'/'+name+'_', columns,rows,samples,gdal.GDT_UInt16)
    raster.SetMetadataItem('AREA_OR_POINT', 'Point')
    raster.SetGeoTransform((extent[0],60,0,extent[-1],0,-60))
    raster.SetProjection('+proj=utm +zone='+int(np.ceil(lon.mean()/6)+30).__str__() +' datum=WGS84 +units=m +no_defs')
    for i in range(samples):
        raster.GetRasterBand(i+1).WriteArray(result[:,:,i])
        raster.FlushCache()
    raster=None
    with open(folder+'/'+name+'_.hdr','a') as f:
        f.write('\nwavelength = {\n'+  ', '.join([str(i) for i in w])     +'}')
        f.write('\nfwhm = {\n'+  ', '.join([str(i) for i in bp])    +'}')
        f.write('\nwavelength units = Micrometers')
        f.write('\nreflectance scale factor = 10000.000000')
        f.close()
