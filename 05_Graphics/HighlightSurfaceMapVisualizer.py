# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 15:45:28 2025

@author: f_meck01
"""

import numpy as np
import warnings
from nilearn import surface, datasets
import pyvista as pv
from matplotlib import cm
from os.path import join, exists
from os import makedirs
from matplotlib.colors import LinearSegmentedColormap

warnings.filterwarnings('ignore')


class HighlightSurfaceMapVisualizer:
    def __init__(
        self,
        stat_img_path,
        corr_img_path,
        clust_img_path,
        out_dir,
        hemi='left',
        fsaverage_res='fsaverage5',
        img_prefix='brain',
        inflated=True,
        cmap='afni_hotcold',
        views=None,
        alpha_range=(0.2, 1.0),
        sulc_gray_range=(0.3, 0.7),
        offset=0.1,
    ):
        self.stat_img_path = stat_img_path
        self.corr_img_path = corr_img_path
        self.clust_img_path = clust_img_path
        self.out_dir = out_dir
        self.hemi = hemi
        self.inflated = inflated
        self.cmap = cmap
        self.img_prefix = img_prefix
        self.alpha_range = alpha_range
        self.sulc_gray_range = sulc_gray_range
        self.offset = offset
        self.views = views or {
            'right':  [[1, 0, 0], [0, 0, 0], [0, 0, 1]],
            'left':   [[-1, 0, 0], [0, 0, 0], [0, 0, 1]],
            'back':   [[0, -1, 0], [0, 0, 0], [0, 0, 1]],
            'front':  [[0, 1, 0], [0, 0, 0], [0, 0, 1]],
            'top':    [[0, 0, 1], [0, 0, 0], [0, 1, 0]],
            'bottom': [[0, 0, -1], [0, 0, 0], [0, 1, 0]],
        }

        if not exists(out_dir):
            makedirs(out_dir, exist_ok=True)

        self.fsavg = datasets.fetch_surf_fsaverage(mesh=fsaverage_res)
        self._prepare_surface_data()

    def _sulc_to_gray(self, sulc_values):
        vmin = np.percentile(sulc_values, 2)
        vmax = np.percentile(sulc_values, 98)
        norm = np.clip((sulc_values - vmin) / (vmax - vmin), 0, 1)
        gmin, gmax = self.sulc_gray_range
        rescaled = gmin + (gmax - gmin) * norm
        return np.tile(rescaled[:, None], (1, 3))

    def _normalize(self, data):
        data_range = data.max() - data.min()
        if data_range == 0:
            return np.zeros_like(data)
        return (data - data.min()) / data_range

    def _scaled_alpha(self, data, power=0.5):
        norm = self._normalize(np.abs(data))
        curved = np.power(norm, power)
        amin, amax = self.alpha_range
        return amin + (amax - amin) * curved

    def _symmetric_normalize(self, data):
        max_abs = np.max(np.abs(data))
        if max_abs == 0:
            return np.full_like(data, 0.5)
        return 0.5 + 0.5 * data / max_abs
    
    
    def _make_afni_style_colormap(self):
        return LinearSegmentedColormap.from_list(
            'afni_hotcold',
            [(0.0,"#00ffff"),  # cyan
             (0.4,"#0000ff"),  # blue
             (0.5,"#ffffff"),  # neutral (at midpoint = 0)
             (0.6,"#ff0000"),  # red
             (1.0,"#ffff00")  # yellow
             ], N=256)

    def _rgba_from_scalar(self, scalars, alpha=None):
        if self.cmap == 'afni_hotcold':
            cmap = self._make_afni_style_colormap()
        else:
            cmap = cm.get_cmap(self.cmap)
        normed = self._symmetric_normalize(scalars)
        rgba = cmap(normed)
        if alpha is not None:
            rgba[:, 3] = alpha
        return rgba[:, :3], rgba[:, 3]

    def _prepare_surface_data(self):
        
        other_hemi = 'right' if self.hemi == 'left' else 'left'
        
        # Projection surfaces (pial is usually recommended for data mapping)
        surf_key_for_projection = f'pial_{self.hemi}'
        surf_key_for_projection_other = f'pial_{other_hemi}'
        
        projection_surface = self.fsavg[surf_key_for_projection]
        projection_surface_other = self.fsavg[surf_key_for_projection_other]
        
        # Visualization mesh
        self.surface_mesh_path = self.fsavg[f'infl_{self.hemi}' if self.inflated else surf_key_for_projection]
        self.sulc_map = surface.load_surf_data(self.fsavg[f'sulc_{self.hemi}']) if self.inflated else None
        self.coords, self.faces = surface.load_surf_mesh(self.surface_mesh_path)
    
        # --- Load both hemispheres' data for normalization ---
        stat_data_self  = surface.vol_to_surf(self.stat_img_path, projection_surface)
        stat_data_other = surface.vol_to_surf(self.stat_img_path, projection_surface_other)
        
        corr_data_self  = surface.vol_to_surf(self.corr_img_path, projection_surface)
        corr_data_other = surface.vol_to_surf(self.corr_img_path, projection_surface_other)
        
        clust_data_self  = surface.vol_to_surf(self.clust_img_path, projection_surface)
        clust_data_other = surface.vol_to_surf(self.clust_img_path, projection_surface_other)
        
        # Store the full datasets for normalization
        self.stat_data_full  = np.concatenate([stat_data_self,  stat_data_other])
        self.corr_data_full  = np.concatenate([corr_data_self,  corr_data_other])
        self.clust_data_full = np.concatenate([clust_data_self, clust_data_other])
        
        # Store the hemisphere-specific data for plotting
        self.stat_data  = stat_data_self
        self.corr_data  = corr_data_self
        self.clust_data = clust_data_self
    
        

    def render(self):
        plotter = pv.Plotter(off_screen=True)
        plotter.enable_anti_aliasing('ssaa')
        # Base mesh
        mesh = pv.PolyData(self.coords, np.hstack([
            np.full((self.faces.shape[0], 1), 3),
            self.faces
        ]))
        mesh_bg = mesh.copy()
        
        if self.inflated:
            sulc_grey = self._sulc_to_gray(self.sulc_map)
            mesh_bg.point_data['sulc_grey'] = sulc_grey
            plotter.add_mesh(mesh_bg, scalars='sulc_grey', rgb=True)
            label = 'inflated'
        else:
            plotter.add_mesh(mesh_bg, color='lightgrey')
            label = 'pial'

        # Calculate colors and alpha once for both meshes
        colors, alphas = self._rgba_from_scalar(
            self.corr_data_full, alpha=self._scaled_alpha(self.stat_data_full)
        )
        rgba = np.concatenate([colors[:len(self.stat_data),], alphas[:len(self.stat_data), None]], axis=1)
        
        # Compute normals and offset coordinates
        normals = mesh.point_normals
        coords_offset = self.coords + self.offset * normals

        # Overlay (with transparency) - entire surface
        mesh_overlay = pv.PolyData(coords_offset, mesh.faces.copy())
        mesh_overlay.point_data['trans'] = rgba
        plotter.add_mesh(mesh_overlay, scalars='trans', rgba=True, show_scalar_bar=True)

        # Cluster overlay - only significant clusters with full opacity
        mask = np.abs(self.clust_data) > 0
        cluster_faces = self.faces[np.all(mask[self.faces], axis=1)]

        if cluster_faces.size > 0:
            cluster_faces_pv = np.hstack([np.full((cluster_faces.shape[0], 1), 3), cluster_faces])
            mesh_cluster = pv.PolyData(self.coords, cluster_faces_pv)
            
            # Same colors but with full opacity for clusters
            rgba_opaque = np.concatenate([colors[:len(self.stat_data)], mask.astype(np.float64).reshape(-1, 1)], axis=1)
            mesh_cluster.point_data['colors'] = rgba_opaque
            plotter.add_mesh(mesh_cluster, scalars='colors', rgba=True)
            
            # Add outline edges
            outline_edges = mesh_cluster.extract_feature_edges(
                boundary_edges=True,
                non_manifold_edges=False,
                feature_edges=False,
                manifold_edges=False
            )
            plotter.add_mesh(outline_edges, color='black', style='wireframe', line_width=5)

        # Render each view
        for name, cpos in self.views.items():
            plotter.camera_position = cpos
            plotter.reset_camera()
            out_path = join(self.out_dir, f"{self.img_prefix}_{label}_{self.hemi}_{name}.png")
            plotter.screenshot(out_path)
            print(f"Saved: {out_path}")
        print(f"Min: {self.corr_data_full.min()}, Max: {self.corr_data_full.max()}")

        plotter.close()