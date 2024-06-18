import os
import zipfile
import shutil
import numpy as np
import xarray as xr
import rioxarray as rxr
import geopandas as gpd
from typing import Union 
from shapely.geometry import box

def create_bbox_kml(file_path:Union[str,os.PathLike], 
                    buffer_dist:float=0, 
                    file_type:str='kml',
                    out_type:str='kml', 
                    out_crs:int=4326):
    '''Create a bounding box from a kml file.

    input:
        file_path: str, path to the kml file
        buffer_dist: float, buffer distance around the bounding box
        type: str, output file type, either 'kml' or 'geojson'
        out_crs: int, output crs, default is 4326
    
    returns:
        nothing, saves the file to the same directory as the input file.
    '''
    file_path = os.path.abspath(file_path)

    out_crs = 'EPSG:' + str(out_crs)
    if file_type == 'kml':
        try:
            file = gpd.read_file(file_path, driver='LIBKML')
        except:
            gpd.io.file.fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'
            file = gpd.read_file(file_path, driver='LIBKML')

    elif file_type == 'geojson' or file_type == 'shapefile':
        file = gpd.read_file(file_path)

    if out_crs[5:] != file.crs.to_epsg():
        file.to_crs(out_crs, inplace=True)
        
    bounds = file.total_bounds
    bbox = box(*bounds).buffer(buffer_dist)    
    bbox_gdf = gpd.GeoDataFrame(gpd.GeoSeries(bbox), columns=['geometry'])

    if out_type == 'shapefile':
        if not os.path.exists(os.path.join(os.path.dirname(file_path), 'shp')):
            os.mkdir(os.path.join(os.path.dirname(file_path), 'shp'))
        outpath = os.path.join(os.path.dirname(file_path), 'shp', os.path.basename(file_path).split('.')[0])
    else:
        outpath = os.path.join(os.path.dirname(file_path), os.path.basename(file_path).split('.')[0])

    if out_type == 'kml':
        try:
            bbox_gdf.to_file(outpath + '_bbox.kml', driver='LIBKML', crs=out_crs)
        except:
            gpd.io.file.fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'
            bbox_gdf.to_file(outpath + '_bbox.kml', driver='LIBKML', crs=out_crs)
    if out_type == 'shapefile':
        bbox_gdf.to_file(outpath + '_bbox.shp', driver='ESRI Shapefile', crs=out_crs)
    elif out_type == 'geojson':
        bbox_gdf.to_file(outpath + '_bbox.geojson', driver='GeoJSON', crs=out_crs) 

def zip_shapefiles(file_path:Union[str,os.PathLike], delete_dir:bool=False):
    '''Zip shapefiles in the same directory as the input file.

    input:
        file_path: str, path to the shapefile components
    output:
        zipped shapefile
    '''
    path = os.path.join(os.path.dirname(file_path))
    fname = os.path.basename(file_path).split('.')[0]

    ext = ('.shp', '.dbf', '.prj', '.shx', '.cpg')
    
    shp_list = []
    with zipfile.ZipFile(f'{fname}_bbox.zip', 'w') as zipMe:   
        for path, dirc, files in os.walk(path):
            for name in files:
                if name.endswith(ext):
                    shp_list.append(name)
                    zipMe.write(os.path.join(path, name), compress_type=zipfile.ZIP_DEFLATED)
    
    if delete_dir == True:
        shutil.rmtree(path)
    
def make_a_raster(path:Union[str, os.PathLike], 
                  res:float, 
                  left:float, 
                  bottom:float, 
                  right:float, 
                  top:float, 
                  crs:str = '4326'):
    '''Makes an empty raster given bounding box extents.

    inputs:
    res: resolution of the raster
    left: left extent of the raster
    bottom: bottom extent of the raster
    right: right extent of the raster
    top: top extent of the raster
    path: path to save the raster

    returns:
    da: xarray dataset with the raster
    '''

    x = np.arange(left, right, res) 
    y = np.arange(bottom, top, res) 
    y = y[::-1]

    X, _ = np.meshgrid(x, y)
    da = xr.DataArray(X, coords={"y":y, "x":x}, dims=['y', 'x'])
    da.rio.write_crs(crs, inplace=True)
    da.rio.to_raster(path)

def get_centroid(path:Union[os.PathLike, str], 
                 file_type:str):
    '''Get the centroid of a shapefile or kml file.

    inputs:
        path: str, path to the file
        file_type: str, type of file: 'shp', 'geojson', or 'kml'
    returns:
        centroid: tuple, centroid of the file
    '''
    
    if file_type == 'kml':
		
        file_type = 'LIBKML'
        try:
            df = gpd.read_file(path, driver=file_type)
        except:
            gpd.io.file.fiona.drvsupport.supported_drivers['LIBKML'] = 'rw'
            df = gpd.read_file(path, driver=file_type)
    else:
        df = gpd.read_file(path, driver=file_type)
	
    centroid = df['geometry'].centroid
    centroid = (centroid.x.item(), centroid.y.item())
    return centroid

def nearest_species(
        detections: gpd.GeoDataFrame, 
        field: gpd.GeoDataFrame
    ):
    '''Match attribute (species) in a point layer with nearest field measurement.

    inputs:
        detections: (geopandas.GeoDataFrame) point layer with species attribute
        field: (geopandas.GeoDataFrame) point layer with field measurements
    returns:
        detections_spp: (geopandas.GeoDataFrame) point layer with nearest field measurement (species)
    '''

    detections_spp = gpd.GeoDataFrame()
    field_alt = field.copy()

    for index, _ in detections.iterrows():

        df = gpd.GeoDataFrame(detections.iloc[index:index+1])
        joined = gpd.sjoin_nearest(df, field_alt)

        #remove the joined row from the field df
        field_alt = field_alt[~field_alt.index.isin(joined['index_right'])]
        detections_spp = pd.concat([detections_spp, joined])

        if len(field_alt) == 0:
            break
                
    leftovers = detections.iloc[index+1:]
    leftover_species = gpd.sjoin_nearest(leftovers, field)
    detections_spp = pd.concat([detections_spp, leftover_species])

    return detections_spp