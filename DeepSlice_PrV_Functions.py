#!/usr/bin/env python
# coding: utf-8

import numpy as np
import scipy.stats as stats
import math
import pandas as pd
import glob2 as glob
from statistics import mean
import numpy as np
from scipy.stats import norm
import math

##### from DeepSlice plane_alignment.py and depth_estimation.py 
# https://github.com/PolarBean/DeepSlice/DeepSlice/DeepSlice/coord_post_processing/depth_estimation.py 
# https://github.com/PolarBean/DeepSlice/DeepSlice/coord_post_processing/plane_alignment_functions/plane_alignment.py

def find_plane_equation(plane):
    """
    Finds the plane equation of a plane
    :param plane: the plane to find the equation of
    :type plane: :any:`numpy.ndarray`
    :returns: the normal vector of the plane and the constant k
    :rtype: :any:`numpy.ndarray`, float
    """
    a, b, c = (
        np.array(plane[0:3], dtype=np.float64),
        np.array(plane[3:6], dtype=np.float64),
        np.array(plane[6:9], dtype=np.float64),
    )
    cross = np.cross(b, c)
    cross /= 9
    k = -((a[0] * cross[0]) + (a[1] * cross[1]) + (a[2] * cross[2]))
    return (cross, k)

def calculate_brain_center_depth(section):
    """
    Calculates the depth of the brain center for a given section

    :param section: the section coordinates as an array consisting of Oxyz,Uxyz,Vxyz 
    :type section: np.array
    :return: the depth of the brain center
    :rtype: float
    """
    cross, k = find_plane_equation(section)
    translated_volume = np.array((456, 0, 320))
    linear_point = (
        ((translated_volume[0] / 2) * cross[0])
        + ((translated_volume[2] / 2) * cross[2])
    ) + k
    depth = -(linear_point / cross[1])
    return depth

def calculate_brain_center_depths(predictions):
    """
    Calculates the depths of the brain center for a series of predictions
    
    :param predictions: dataframe of predictions
    :type predictions: pandas.DataFrame
    :return: a list of depths
    :rtype: list[float]
    """
    depths = []
    for prediction in predictions[
        ["ox", "oy", "oz", "ux", "uy", "uz", "vx", "vy", "vz"]
    ].values:
        depths.append(calculate_brain_center_depth(prediction))
    return depths


# Custom functions

def interpolate_depths(reference_df, sample_df):
    # DeepSlice_wDepth = interpolate_depths(DeepSlice_PrV, DeepSlice_rawOutput)
    
    lookup_df = pd.DataFrame(data= {'depths': reference_df["depths"][1:397], 'vole_depths':reference_df["vole_depths"][1:397]})
    lookup_df = lookup_df.sort_values(by='depths', ascending=True)
    min_vals = [np.min(lookup_df["depths"]),np.min(lookup_df["vole_depths"])]
    max_vals = [np.max(lookup_df["depths"]),np.max(lookup_df["vole_depths"])]
    
    input_df = sample_df.copy()
    depths_list = []
    vole_depths_list = []
    
    for alignment in input_df.iterrows():
        m = alignment[1][["ox", "oy", "oz", "ux", "uy", "uz", "vx", "vy", "vz"]].values.astype(np.float64)
        depth = calculate_brain_center_depth(m)
        
        if depth < min_vals[0]: 
            vole_depth = min_vals[1] # values more posterior than end of look-up table are matched with last value of look-up table
            print('depth = ' + str(depth) + ', which is outside of range. Matched with nearest reference value.')
        elif depth > max_vals[0]: 
            vole_depth = max_vals[1] # values more anterior than start of look-up table are matched with first value of look-up table
            print('depth = ' + str(depth) + ', which is outside of range. Matched with nearest reference value.')
        else: vole_depth = np.interp(depth, lookup_df["depths"], lookup_df["vole_depths"])

        depths_list.append(depth)
        vole_depths_list.append(vole_depth)
    
    input_df.insert(12, "depths",depths_list,False)
    input_df.insert(13, "vole_depths",vole_depths_list,False)
   
    return input_df


def adjust_size_center(reference_df, section):
    # o_adj, u_adj, v_adj = adjust_size_center(DeepSlice_PrV, alignment[1][["ox", "oy", "oz", "ux", "uy", "uz", "vx", "vy", "vz","vole_depths"]])
    
    ######################################## Match section with relevant reference section ##################################################################

    lookup_df = pd.DataFrame(data= {'ox':reference_df["ox"][1:397],'oy':reference_df["oy"][1:397],'oz':reference_df["oz"][1:397], 
                                'ux':reference_df["ux"][1:397],'uy':reference_df["uy"][1:397],'uz':reference_df["uz"][1:397],
                                'vx':reference_df["vx"][1:397],'vy':reference_df["vy"][1:397],'vz':reference_df["vz"][1:397],
                                'vole_depths':reference_df["vole_depths"][1:397]})
    lookup_df = lookup_df.sort_values(by='vole_depths', ascending=True)

    sample_depth = section["vole_depths"]

    o_ref = [np.interp(sample_depth, lookup_df["vole_depths"], lookup_df["ox"]),
                   np.interp(sample_depth, lookup_df["vole_depths"], lookup_df["oy"]),
                   np.interp(sample_depth, lookup_df["vole_depths"], lookup_df["oz"])]
    u_ref = [np.interp(sample_depth, lookup_df["vole_depths"], lookup_df["ux"]),
                   np.interp(sample_depth, lookup_df["vole_depths"], lookup_df["uy"]),
                   np.interp(sample_depth, lookup_df["vole_depths"], lookup_df["uz"])]
    v_ref = [np.interp(sample_depth, lookup_df["vole_depths"], lookup_df["vx"]),
             np.interp(sample_depth, lookup_df["vole_depths"], lookup_df["vy"]),
             np.interp(sample_depth, lookup_df["vole_depths"], lookup_df["vz"])]

    o_sample = section[["ox","oy","oz"]].values.astype(np.float64)
    u_sample = section[["ux","uy","uz"]].values.astype(np.float64)
    v_sample = section[["vx","vy","vz"]].values.astype(np.float64)

    ####################################################### Move centroid and adjust scale ##################################################################

    centroid_ref = centroid_from_vectors(o_ref,u_ref,v_ref) # finds rough center of reference section (in mouse vx)
    centroid_sample = centroid_from_vectors(o_sample,u_sample,v_sample) # finds rough center of sample section (in mouse vx)

    width_ref = magnitude_3d(u_ref) # finds width of reference (in mouse vx)
    width_sample = magnitude_3d(u_sample) # finds width of sample (in mouse vx)

    height_ref = magnitude_3d(v_ref) # finds height of reference (in mouse vx)
    height_sample = magnitude_3d(v_sample) # finds height of sample (in mouse vx)

    centroid_dif = [s-r for s,r in zip(centroid_sample,centroid_ref)] # finds difference between sample and reference centroid (in mouse vx)
    adj_centr_dif = [(centroid_dif[0]*417/width_ref),0,(centroid_dif[2]*350/height_ref)] # scales the x and z components to vole dimensions
    centroid_ref_vole = [417/2, 0, 350/2] # center of reference section in vole space
    centroid_sample_vole = [r+dif for r,dif in zip(centroid_ref_vole,adj_centr_dif)] # uses scaled difference vector to set centroid of section

    adjusted_width = width_sample*417/width_ref # scales the sample width to vole dimensions
    adjusted_height = height_sample*350/height_ref # scales the sample height to vole dimensions

    o_adj = [(centroid_sample_vole[0]+(adjusted_width/2)),sample_depth,(centroid_sample_vole[2]+(adjusted_height/2))] # new o based on w/h and centroid (in vole vx)
    u_adj = [-adjusted_width,0,0] # new u from adjusted width (vole vx)
    v_adj = [0,0,-adjusted_height] # new v from adjusted height (vole vx)

    ##################################################### Adjust angles of U and V, find new O ##############################################################

    u_sample_flat = [u_sample[0],u_sample[2]] # flattened vector u from sample section
    v_sample_flat = [v_sample[0],v_sample[2]] # flattened vector v from sample section
    u_ref_flat = [u_ref[0],u_ref[2]] # flattened vector u from reference section
    v_ref_flat = [v_ref[0],v_ref[2]] # flattened vector v from reference section

    ### find new ux/z and vx/z
    u_rot,angle_u_dif = uv_rotation(u_sample_flat,u_ref_flat,adjusted_width,'u')
    v_rot,angle_v_dif = uv_rotation(v_sample_flat,v_ref_flat,adjusted_height,'v')
    
    ### find new ox, oz
    o_rot = o_rotation(o_adj, angle_v_dif, centroid_ref_vole)

    
    return o_rot,u_rot,v_rot

def uv_rotation(sample_flat,ref_flat,adj_length,u_or_v):
    # u_rot = uv_rotation(u_sample_flat,u_ref_flat,adjusted_width,'u')
    
    angle_sample = angle_xz(sample_flat,u_or_v,'rad') # calculate original angle of sample u in radians
    angle_ref = angle_xz(ref_flat,u_or_v,'rad') # calculate original angle of reference u in radians
    angle_dif = angle_ref-angle_sample
    
    if 'u' in u_or_v.lower():
        x_rot = adj_length * np.cos(angle_dif) # find rotated u (adjusted u angle with scaled width)
        z_rot = adj_length * np.sin(angle_dif) 
        vector_rot = [-x_rot,0,z_rot]
    elif 'v' in u_or_v.lower():     
        x_rot = adj_length * np.sin(angle_dif) # find rotated v (adjusted v angle with scaled height)
        z_rot = adj_length * np.cos(angle_dif)
        vector_rot = [x_rot,0,-z_rot]
    return vector_rot, angle_dif

def o_rotation(o_adj, angle_v_dif, centroid):
    cO = [O-c for O,c in zip(o_adj,centroid)]
    cO_flat = [cO[0], cO[2]]
    len_cO = magnitude_3d([cO[0], 0, cO[2]])
    cO_angle = np.arctan(cO_flat[1]/cO_flat[0])*180/np.pi

    new_opposite_angle = (90-cO_angle) - (angle_v_dif * 180/np.pi)

    delta_x_rot = np.sin((new_opposite_angle * np.pi/180)) * len_cO
    delta_z_rot = np.cos((new_opposite_angle * np.pi/180)) * len_cO
    cO_rot = [delta_x_rot, 0, delta_z_rot]

    o_rot = [cO_rot[0] + centroid[0], cO_rot[1] + o_adj[1], cO_rot[2]+centroid[2]]
    
    return o_rot

def angle_xz(vector,u_or_v,units):
    #e.g., angle_u = angle_xz(u_sample_flat,'u','degrees')
    
    if len(vector) == 2:
        x_val = vector[0]
        z_val = vector[1]
    elif len(vector) == 3:
        x_val = vector[0]
        z_val = vector[2]
    else: print('Error: vector should either be formatted as [x,y,z] or [x,z]')
    
    if u_or_v == 'u': angle_rad = np.arctan(z_val/x_val)
    elif u_or_v == 'v': angle_rad = np.arctan(x_val/z_val)
    else: print("Error: must specify 'u' or 'v'")
    
    if 'deg' in units.lower(): angle = angle_rad*180/np.pi
    elif 'rad' in units.lower(): angle = angle_rad
    else: print("Error: units must contain 'rad' or 'deg'")

    return angle

def substitute_PrVRef_values(sample_df, o_list,u_list,v_list):
    # DeepSlice_processed = substitute_PrVRef_values(DeepSlice_wDepth,o_adj,u_adj,v_adj)
    
    input_df = sample_df.copy()
    output_df = pd.DataFrame().reindex_like(input_df.filter(['Filenames', 'ox','oy', 'oz', 'ux', 'uy', 'uz', 'vx', 'vy', 'vz', 'width', 'height','vole_depths'], axis=1))
    
    nr_list = []
    for i in range(len(output_df)): # Copy filenames, width, and height from the raw output
        output_df.loc[i,"Filenames"] = input_df.loc[i,"Filenames"]
        output_df.loc[i,"width"] = input_df.loc[i,"width"]
        output_df.loc[i,"height"] = input_df.loc[i,"height"]
        output_df.loc[i,"vole_depths"] = input_df.loc[i,"vole_depths"]
        
        fname = input_df.loc[i,"Filenames"]
        if '_s' in fname: nr_list.append(int(fname[fname.find('_s')+2:].split('_')[0]))

        o_adj = o_list[i]
        u_adj = u_list[i]
        v_adj = v_list[i]
        
        output_df.loc[i,"ox"] = o_adj[0]
        output_df.loc[i,"oy"] = o_adj[1]
        output_df.loc[i,"oz"] = o_adj[2]

        output_df.loc[i,"ux"] = u_adj[0]
        output_df.loc[i,"uy"] = 0
        output_df.loc[i,"uz"] = u_adj[2]

        output_df.loc[i,"vx"] = v_adj[0]
        output_df.loc[i,"vy"] = 0
        output_df.loc[i,"vz"] = v_adj[2]
        
    if 'nr' not in output_df.columns:
        if len(nr_list)==len(output_df): output_df.insert(len(output_df.columns),'nr',nr_list,allow_duplicates = False)
    
    if 'nr' in output_df.columns: output_df = output_df.sort_values(by=['nr'],ignore_index=True)
    else:
        if 'vole_depths' in output_df.columns: output_df = output_df.sort_values(by=['vole_depths'],ascending=False,ignore_index=True)
    
    return output_df

def centroid_from_vectors(o, u, v):
    # eg. centroid = centroid_from_vectors([ox, oy, oz], [ux, uy, uz], [vx, vy, vz])
    
    top_left_corner = o
    top_right_corner = list(map(sum, zip(o,u)))
    bottom_left_corner = list(map(sum, zip(o,v)))
    bottom_right_corner = list(map(sum, zip(o,u,v)))
    
    centroid = [(top_right_corner[0]+bottom_left_corner[0])/2, (top_right_corner[1]+bottom_left_corner[1])/2,(top_right_corner[2]+bottom_left_corner[2])/2]
    
    return centroid #depth is centroid y component

def magnitude_3d(a):
    
    magnitude = np.sqrt(np.sum([a[0]**2, a[1]**2, a[2]**2]))
    
    return magnitude


##### from DeepSlice QuickNII_functions.py (https://github.com/PolarBean/DeepSlice/DeepSlice/read_and_write/)
import pandas as pd
import xml.etree.ElementTree as ET
import json
import numpy as np


def write_QuickNII_XML(df: pd.DataFrame, filename: str, aligner: str) -> None: # Note: aligner is DeepSlice version (JS)
    """
    Converts a pandas DataFrame to a quickNII compatible XML
    """
    df_temp = df.copy()
    if "nr" not in df_temp.columns:
        df_temp["nr"] = np.arange(len(df_temp)) + 1
    df_temp[["ox", "oy", "oz", "ux", "uy", "uz", "vx", "vy", "vz", "nr"]] = df[
        ["ox", "oy", "oz", "ux", "uy", "uz", "vx", "vy", "vz", "nr"]
    ].astype(str)
    out_df = pd.DataFrame(
        {
            "anchoring": "ox="
            + (df_temp.ox)
            + "&oy="
            + (df_temp.oy)
            + "&oz="
            + (df_temp.oz)
            + "&ux="
            + (df_temp.ux)
            + "&uy="
            + (df_temp.uy)
            + "&uz="
            + (df_temp.uz)
            + "&vx="
            + (df_temp.vx)
            + "&vy="
            + (df_temp.vy)
            + "&vz="
            + (df_temp.vz),
            "filename": df_temp.Filenames,
            "height": df_temp.height,
            "width": df_temp.width,
            "nr": df_temp.nr,
        }
    )
    print(f"saving to {filename}.xml")

    out_df.to_xml(
        filename + ".xml",
        index=False,
        root_name="series",
        row_name="slice",
        attr_cols=list(out_df.columns),
        namespaces={
            "first": df_temp.nr.values[0],
            "last": df_temp.nr.values[-1],
            "name": filename,
            "aligner": aligner,
            "": "",
        },
    )


def read_QuickNII_XML(filename: str) -> pd.DataFrame:
    """
    Converts a QuickNII XML to a pandas dataframe

    :param xml: The path to the QuickNII XML
    :type xml: str
    :return: A pandas dataframe
    :rtype: pd.DataFrame
    """
    df = pd.read_xml(filename)
    # split the anchoring string into separate columns
    anchoring = df.anchoring.str.split("&", expand=True).values
    # lambda function to remove non_numeric characters besides '.', we need this as all the 'ox=' etc is still in the strings
    strip = lambda x: "".join(
        c for c in x if c.isdigit() or c == "." or c == "-" or c == "e"
    )
    ##vectorise the lambda function and apply it to all elements
    anchoring = np.vectorize(strip)(anchoring)
    anchoring = anchoring.astype(np.float64)
    out_df = pd.DataFrame({"Filenames": df.filename})
    out_df[["ox", "oy", "oz", "ux", "uy", "uz", "vx", "vy", "vz"]] = anchoring
    return out_df


def write_QUINT_JSON(
    df: pd.DataFrame, filename: str, aligner: str, target: str 
) -> None: # Note: target is name of cutlas file (JS)
    """
    Converts a pandas DataFrame to a QUINT (QuickNII, Visualign, & Nutil) compatible JSON
    """
    if "nr" not in df.columns:
        df["nr"] = np.arange(len(df)) + 1
    alignments = df[["ox", "oy", "oz", "ux", "uy", "uz", "vx", "vy", "vz"]].values
    if "markers" in df.columns:
        markers = df.markers.values
    else:
        markers = [[]] * len(df)
    #print(len(markers)) commented out by JS
    alignment_metadata = [
        {
            "filename": fn,
            "anchoring": list(alignment),
            "height": h,
            "width": w,
            "nr": nr,
            "markers": marker[0] if len(marker)> 0 else [],
        }
        for fn, alignment, nr, marker, h, w in zip(
            df.Filenames, alignments, df.nr, markers, df.height, df.width
        )
    ]
    QUINT_json = {
        "name": "",
        "target": target,
        "aligner": aligner,
        "slices": alignment_metadata,
    }
    print(f"saving to {filename}.json")
    with open(filename + ".json", "w") as f:
        json.dump(QUINT_json, f)
    with open(filename + ".json", "w") as outfile:
        json.dump(QUINT_json, outfile)


def read_QUINT_JSON(filename: str) -> pd.DataFrame:
    """
    Converts a QUINT JSON to a pandas dataframe
    
    :param json: The path to the QUINT JSON
    :type json: str
    :return: A pandas dataframe
    :rtype: pd.DataFrame
    """
    with open(filename, "r") as f:
        data = json.load(f)
    sections = data["slices"]
    target_volume = data["target"]
    alignments = [
        row["anchoring"] if "anchoring" in row else 9 * [np.nan] for row in sections
    ]
    height = [row["height"] if "height" in row else [] for row in sections]
    width = [row["width"] if "width" in row else [] for row in sections]
    filenames = [row["filename"] if "filename" in row else [] for row in sections]
    section_numbers = [row["nr"] if "nr" in row else [] for row in sections]
    markers = [row["markers"] if "markers" in row else [] for row in sections]
    df = pd.DataFrame({"Filenames": filenames, "nr": section_numbers})
    df[["ox", "oy", "oz", "ux", "uy", "uz", "vx", "vy", "vz"]] = alignments
    df["markers"] = markers
    df["height"] = height
    df["width"] = width
    return df, target_volume
