import os
import zipfile
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
        if not os.path.exists(os.path.join(os.path.dirname(file_path) + 'shp')):
            os.mkdir(os.path.join(os.path.dirname(file_path) + 'shp'))
        outpath = os.path.join(os.path.dirname(file_path) + 'shp', os.path.basename(file_path).split('.')[0])
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

def zip_shapefiles(file_path:Union[str,os.PathLike]):
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

    