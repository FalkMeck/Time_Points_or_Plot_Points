# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 15:49:35 2025

@author: f_meck01
"""

### Create the plots necessary!

import os
import HighlightSurfaceMapVisualizer_bothHemi as hsmv

os.chdir('...05_Graphics')
study_dir = '.../_AFNI_Analysis/'
         
contrasts= ["01_Scene_Shot", "02_04s_12s","03_04s_36s","04_12s_36s",
           "INT01_12sSceneShot_04sSceneShot","INT02_36sSceneShot_04sSceneShot"]


for i, con in enumerate(contrasts):
    
    corr_img = study_dir +  con + '/Con_' + con + '_corr.nii.gz'
    t_img = study_dir +  con + '/Con_' + con + '_tstat.nii.gz'
    sig_clus_img = study_dir +  con + '/ClusterizeOutput_report' + '/Con_' + con + '_EffEst.nii.gz'
    outDir = study_dir +  con +  '/Result_Images'
    
    for h in ('left', 'right'):
        for inflation in (True, False):
            visualizer = hsmv.HighlightSurfaceMapVisualizer(stat_img_path = t_img,
                              corr_img_path = corr_img, 
                              clust_img_path = sig_clus_img,
                              img_prefix = con,
                              out_dir = outDir,
                              hemi = h,
                              fsaverage_res='fsaverage7',
                              inflated=inflation,
                              cmap='afni_hotcold',  # or 'bwr'
                              views=None, # or leave as None for all 6    
                              offset = 0.1
                              )
            visualizer.render()