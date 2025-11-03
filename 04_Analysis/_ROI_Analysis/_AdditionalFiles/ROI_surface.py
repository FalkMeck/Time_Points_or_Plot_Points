# -*- coding: utf-8 -*-
"""
Created on Sat Oct 25 14:22:53 2025

@author: f_meck01
"""

import numpy as np
from nilearn import surface, datasets
import warnings

import pyvista as pv
from matplotlib import cm
from os.path import join
from matplotlib.colors import LinearSegmentedColormap
import glob


warnings.filterwarnings('ignore')


def sulc_to_gray(sulc_values, vmin=None, vmax=None, gray_min=0.3, gray_max=0.7):
    if vmin is None:
        vmin = np.percentile(sulc_values, 2)
    if vmax is None:
        vmax = np.percentile(sulc_values, 98)

    sulc_norm = np.clip((sulc_values - vmin) / (vmax - vmin), 0, 1)
    sulc_rescaled = gray_min + (gray_max - gray_min) * sulc_norm
    rgb_gray = np.tile(sulc_rescaled[:, None], (1, 3))
    return rgb_gray


def normalize(data):
    return (data - data.min()) / (data.max() - data.min())

def rescale_alpha(data, alpha_min=0.2, alpha_max=1.0):
    norm = normalize(np.abs(data))  # range 0 to 1
    return alpha_min + (alpha_max - alpha_min) * norm

def scaled_alpha(data, lower=0.2, upper=1.0, power=0.5):
    """Apply nonlinear scaling to emphasize low values."""
    normed = normalize(np.abs(data))
    curved = np.power(normed, power)  # e.g. sqrt if power=0.5
    return lower + (upper - lower) * curved

def symmetric_normalize(data):
    max_abs = np.nanmax(np.abs(data))
    return 0.5 + 0.5 * data / max_abs  # Ensures zero maps to 0.5

def make_afni_style_colormap():
    return LinearSegmentedColormap.from_list(
        'afni_hotcold',
        [(0.0,"#00ffff"),  # cyan
         (0.4,"#0000ff"),  # blue
         (0.5,"#ffffff"),  # neutral (at midpoint = 0)
         (0.6,"#ff0000"),  # red
         (1.0,"#ffff00")  # yellow
         ],
        N=256
    )

def rgba_from_scalar(scalars, colormap="afni_hotcold", alpha=None):
    if colormap == 'afni_hotcold':
        cmap = make_afni_style_colormap()
    else:
        cmap = cm.get_cmap(colormap)

    normed = symmetric_normalize(scalars)
    rgba = cmap(normed)
    if alpha is not None:
        rgba[:, 3] = alpha
        return rgba[:, :3], rgba[:, 3]
    else:
        return rgba, None
    
def hex_to_rgba(hex_code, opacity=1.0):
    h = hex_code.lstrip('#')
    rgb = np.array([int(h[i:i+2], 16) for i in (0, 2, 4)]) / 255.0
    return np.hstack([rgb, opacity])


fsaverage = datasets.fetch_surf_fsaverage(mesh = 'fsaverage7')
fsavg_maps = {'pial_left': fsaverage.pial_left,
 'pial_right': fsaverage.pial_right,
 'infl_left': fsaverage.infl_left, 
 'infl_right': fsaverage.infl_right,
 'sulc_left': fsaverage.sulc_left,
 'sulc_right': fsaverage.sulc_right}

  
hemi = 'left'  # switch to right
surf_key = f'pial_{hemi}'
inflated = True

outDir = '.../Material_for_Figure/ROI'

corr_img_path = '.../NIFTI/_AFNI_Analysis/INT02_36sSceneShot_04sSceneShot/Con_INT02_36sSceneShot_04sSceneShot_corr.nii.gz'
stat_img_path = '.../NIFTI/_AFNI_Analysis/INT02_36sSceneShot_04sSceneShot/Con_INT02_36sSceneShot_04sSceneShot_tstat.nii.gz'

corr_data = surface.vol_to_surf(corr_img_path, fsavg_maps[surf_key], radius=3, interpolation='linear', kind='auto')
t_data = surface.vol_to_surf(stat_img_path, fsavg_maps[surf_key], radius=3, interpolation='linear', kind='auto')

#baldassanoROI = '.../Baldassano_2017_ROIs'
#additionalROI = '.../Additional_ROIs'
STSROI = '.../STS_freeroi'

roi_img_paths = sorted(glob.glob(join(STSROI, "l_*.nii"))) # switch to r

# select appropriate, order of colors must match order of ROI files
# hexcolors = ['#404788','#FDE725','#73D055','#440154','#1F968B']# Baldassano order 
#hexcolors = ['#FDE725','#404788','#1F968B','#440154','#73D055']# Addtional
hexcolors = ['#440154','#1F968B','#FDE725'] #STS

# init plotter
plotter = pv.Plotter(off_screen = True)
plotter.enable_anti_aliasing('ssaa')

if inflated:
    surf_mesh = fsavg_maps[f'infl_{hemi}']    
    bg_map = fsavg_maps[f'sulc_{hemi}']
    sulc_map =  surface.load_surf_data(bg_map)
    sulc_grey = sulc_to_gray(sulc_map)
else:
    surf_mesh = fsavg_maps[f'pial_{hemi}']    
    bg_map = None

# Load surface coordinates and faces
coords_surf, faces_surf = surface.load_surf_mesh(surf_mesh)


mesh = pv.PolyData(coords_surf, np.hstack([
        np.full((faces_surf.shape[0], 1), 3),  # Number of points per face (3)
        faces_surf
    ]))


# --- Mesh 1: Sulcal depth background on inflated surface ---
mesh_bg = mesh.copy()
if inflated:
    mesh_bg.point_data['sulc_grey'] = sulc_grey
    plotter.add_mesh(mesh_bg, scalars='sulc_grey', rgb=True)
else:
    plotter.add_mesh(mesh_bg, color='lightgrey')

# --- Mesh 2: Alpha-modulated statistical map ---
colors_t, alphas_t = rgba_from_scalar(
    corr_data,
    colormap="afni_hotcold",
    alpha=rescale_alpha(t_data)
)
rgba = np.concatenate([colors_t, alphas_t[:, None]], axis=1)


# Compute normals
normals = mesh.point_normals

# Offset overlay mesh slightly
offset_distance = 0.1  # in mm; tweak if needed
coords_offset = coords_surf + offset_distance * normals

mesh_overlay = pv.PolyData(coords_surf, mesh.faces.copy())
mesh_overlay.point_data['trans'] = rgba
#mesh_overlay.point_data['trans'] = np.concatenate([rgba[:,0:3], np.ones((colors_t.shape[0], 1))], axis = 1)
plotter.add_mesh(mesh_overlay, scalars='trans', rgba=True, show_scalar_bar=True)

# Mesh that only contans the significant clusters

for path, color in zip(roi_img_paths, hexcolors):

    roi_data = surface.vol_to_surf(path, fsavg_maps[surf_key], radius=3, interpolation='linear', kind='auto')

    mask = roi_data > 0
    cluster_faces = faces_surf[np.all(mask[faces_surf], axis=1)]

    if cluster_faces.size > 0:
        cluster_faces_pv = np.hstack([np.full((cluster_faces.shape[0], 1), 3), cluster_faces])
        mesh_cluster = pv.PolyData(coords_offset, cluster_faces_pv)
        # Choose your hex color
        rgba_color = hex_to_rgba(color, opacity=0.5)  # fully opaque
    
        # Duplicate this color for every point in the mesh
        rgba_opac = np.tile(rgba_color, (mesh_cluster.n_points, 1))
        mesh_cluster.point_data['colors'] = rgba_opac  # full array OK
       
        
        outline_edges = mesh_cluster.extract_feature_edges(
           boundary_edges=True,
           non_manifold_edges=False,
           feature_edges=False,
           manifold_edges=False )
        
        plotter.add_mesh(outline_edges, color=color, style='wireframe', line_width=15)
        plotter.add_mesh(mesh_cluster, scalars='colors', rgba=True)


views = {
    'right':  [[1, 0, 0], [0, 0, 0], [0, 0, 1]],
    'left':   [[-1, 0, 0], [0, 0, 0], [0, 0, 1]],
    'back':   [[0, -1, 0], [0, 0, 0], [0, 0, 1]],
    'front':  [[0, 1, 0], [0, 0, 0], [0, 0, 1]],
    'top':    [[0, 0, 1], [0, 0, 0], [0, 1, 0]],
    'bottom': [[0, 0, -1], [0, 0, 0], [0, 1, 0]],
}

   
    
for name, cpos in views.items():
     plotter.camera_position = cpos
     plotter.reset_camera()  # Ensures full view in frame
     plotter.screenshot(join(outDir, f"STS_LH_{name}.png")) # switch to RH
plotter.close()