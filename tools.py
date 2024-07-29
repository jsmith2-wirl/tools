import os
import zipfile
import shutil
import numpy as np
import pandas as pd
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

def generate_table4(input_df: Union[gpd.GeoDataFrame, pd.DataFrame], total_area, n_sapling):
    """Generates Table 4 for the report template, carbon summary by tree species.
    input:
        biomass (gpd.GeoDataFrame): Geopandas dataframe containing the final biomass output
        total_area (float): Total area of the plot in hectares
        n_sapling (int): Number of saplings in the project (Dbh < 7 cm).
    returns:
        Nothing. Writes the output to a CSV file.
    """
    # Extract necessary data
    total_tree_count = len(input_df)
    total_tCO2e = input_df["co2e_total"].sum().round(2)
    avg_height = input_df["height"].mean().round(1)
    std_height = input_df["height"].std().round(1)
    min_height = input_df["height"].min().round(1)
    max_height = input_df["height"].max().round(1)
    avg_dbh = input_df["predicted_dbh"].mean().round(1)
    std_dbh = input_df["predicted_dbh"].std().round(1)
    min_dbh = input_df["predicted_dbh"].min().round(1)
    max_dbh = input_df["predicted_dbh"].max().round(1)
    avg_tCO2e = input_df["co2e_total"].mean().round(2)

    # Start building the CSV content
    csv_content = f"Total Area (ha),,{total_area},,,,\n"
    csv_content += f"Sapling count,,{n_sapling},,,,\n"
    csv_content += ",,Tree count,tCO2e,Mean Tree Height (m),Mean Tree Dbh (cm),Mean tCO2e\n"
    csv_content += (
        f"Total,,{total_tree_count},{total_tCO2e},{avg_height} ± "
        f"{std_height} [{min_height}-{max_height}],{avg_dbh} ± {std_dbh} [{min_dbh}-{max_dbh}],{avg_tCO2e} \n"
    )

    # Add species-specific rows
    spps = input_df["species"].value_counts().keys().tolist()

    for i, species in enumerate(spps):
        species_df = input_df[input_df["species"] == species]
        species_tree_count = len(species_df)
        species_tCO2e = species_df["co2e_total"].sum().round(2)
        species_avg_height = species_df["height"].mean().round(1)
        species_std_height = species_df["height"].std().round(1)
        species_min_height = species_df["height"].min().round(1)
        species_max_height = species_df["height"].max().round(1)
        species_avg_dbh = species_df["predicted_dbh"].mean().round(1)
        species_std_dbh = species_df["predicted_dbh"].std().round(1)
        species_min_dbh = species_df["predicted_dbh"].min().round(1)
        species_max_dbh = species_df["predicted_dbh"].max().round(1)
        species_avg_tCO2e = species_df["co2e_total"].mean().round(2)

        if i == 0:
            csv_content += (
                f"By Species,{species},{species_tree_count},{species_tCO2e},{species_avg_height} ± "
                f"{species_std_height} [{species_min_height} - {species_max_height}],{species_avg_dbh} ± "
                f"{species_std_dbh} [{species_min_dbh}-{species_max_dbh}], {species_avg_tCO2e} \n"
            )
        else:
            csv_content += (
                f",{species},{species_tree_count},{species_tCO2e},{species_avg_height} ± "
                f"{species_std_height} [{species_min_height} - {species_max_height}],{species_avg_dbh} ± "
                f"{species_std_dbh} [{species_min_dbh}-{species_max_dbh}], {species_avg_tCO2e} \n"
            )

    # Write the content to a CSV file
    out_dir = os.path.join("./treeconomy_utils/report/tables/output")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    with open(os.path.join(out_dir, "table4.csv"), "w") as file:
        file.write(csv_content)


def compartmentalizer(input_df: gpd.GeoDataFrame):
    """Looping function to iterate through each project compartment. To be used for generate_appendix6.
    input:
        biomass (gpd.GeoDataFrame): Geopandas dataframe containing the final biomass output
    returns:
        Nothing. Writes the output to a CSV file.
    """
    comps = input_df["compartment"].unique()
    for comp in comps:
        comp_data = input_df[input_df["compartment"] == comp]
        generate_appendix6(comp_data, comp)


def generate_appendix3(equations_df: pd.DataFrame) -> pd.DataFrame:
    """Generates appendix 3, the table indicating equations used for Dbh predictions.
    input:
        equations_df (pd.DataFrame): DataFrame containing the equations and RMSE values, generated from height_to_dbh.py script.
    returns:
        appendix_3 (pd.DataFrame): DataFrame containing the formatted appendix 3 table.
    """
    appendix_3 = equations_df.copy()

    for _, row in appendix_3.iterrows():
        ast = "a="
        aen = ","
        bst = "b="
        a = f"{float(re.search(f'{ast}(.*){aen}', row['Parameters']).group(1)):.5f}"
        b = f"{float(re.search(f'{bst}(.*)', row['Parameters']).group(1)):.5f}"

        if row["Fitting_Function"] == "power_law_function":
            appendix_3.loc[_, "Parameters"] = f"Dbh = {a} * height^{b}"

        elif row["Fitting_Function"] == "exponential_function":
            appendix_3.loc[_, "Parameters"] = f"Dbh = {a} * {b} ^ height"

    appendix_3["RMSE"] = appendix_3["RMSE"].round(1)

    appendix_3.rename(columns={"Parameters": "Equation", "Common_Name": "Species"}, inplace=True)
    appendix_3.drop(columns=["Fitting_Function", "Tallo_Used", "Species_Latin"], inplace=True)
    appendix_3.loc[0, "Source"] = "Models fit by Treeconomy from the Tallo Global database (Jucker et al, 2022)."

    out_dir = os.path.join("./treeconomy_utils/report/tables/output")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    appendix_3.to_csv(os.path.join(out_dir, "appendix3.csv"), index=False)

    return appendix_3


def generate_appendix6(input_df: Union[gpd.GeoDataFrame, pd.DataFrame], compartment):
    """Generate Appendix 6, compartment-based summaries of per-species carbon.
    input:
        input_df (gpd.GeoDataFrame): Geopandas dataframe containing the final biomass output.
        compartment (str): Compartment number/ID of the target compartment.
    returns:
        Nothing. Writes the output to a CSV file.
    """
    compartment_name = (
        f"Compartment {input_df['compartment'].unique()[0]}"  # Assuming single compartment for simplicity
    )
    compartment_name = f"Compartment {compartment}"

    total_area = input_df["area"].iloc[0]  # Assuming same area for all rows for simplicity
    total_tree_count = len(input_df)
    avg_density = input_df["density"].iloc[0]
    avg_spacing = input_df["spacing"].iloc[0]
    total_tCO2e = input_df["co2e_total"].sum().round(2)
    avg_height = input_df["height"].mean().round(1)
    avg_dbh = input_df["predicted_dbh"].mean().round(1)
    avg_tCO2e = input_df["co2e_total"].mean().round(2)

    # Start building the CSV content
    csv_content = f"{compartment_name},,,,,,\n"
    csv_content += f"Total Area (ha),,{total_area},,,,\n"
    csv_content += (
        ",,Tree count,Density (ha),Spacing (m),tCO2e,Mean Tree Height (m),Mean Tree Dbh (cm),Mean tree tCO2e\n"
    )
    csv_content += (
        f"Total,,{total_tree_count},{avg_density},{avg_spacing},{total_tCO2e},{avg_height},{avg_dbh},{avg_tCO2e}\n"
    )

    input_df = input_df[input_df['compartment'] == compartment]
    spps = input_df["species"].value_counts().keys().tolist()

    # Add species-specific rows
    for i, species in enumerate(spps):
        species_df = input_df[input_df["species"] == species]
        species_tree_count = len(species_df)
        species_tCO2e = species_df["co2e_total"].sum().round(2)
        species_avg_height = species_df["height"].mean().round(1)
        species_avg_dbh = species_df["predicted_dbh"].mean().round(1)
        species_avg_tCO2e = species_df["co2e_total"].mean().round(2)

        if i == 0:
            csv_content += f"By Species,{species},{species_tree_count},,,{species_tCO2e},{species_avg_height},{species_avg_dbh},{species_avg_tCO2e}\n"
        else:
            csv_content += f",{species},{species_tree_count},,,{species_tCO2e},{species_avg_height},{species_avg_dbh}, {species_avg_tCO2e}\n"

    out_dir = os.path.join("./treeconomy_utils/report/tables/output")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Write the content to a CSV file
    with open(os.path.join(out_dir, f"{compartment}_appendix6.csv"), "w") as file:
        file.write(csv_content)


def generate_appendix7(species_dict: dict):
    """Generate Appendix 7, table showing species mapped to WCC-compatible species.
    input:
        species_dict (dict): Dictionary containing the original species as keys and the re-mapped species as values.
    returns:
        Nothing. Writes the output to a CSV file.
    """

    df = pd.DataFrame(species_dict.items(), columns=["Original species", "Re-mapped species"])

    out_dir = os.path.join("./treeconomy_utils/report/tables/output")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    df.to_csv(os.path.join(out_dir, "appendix7.csv"), index=False)

    return df