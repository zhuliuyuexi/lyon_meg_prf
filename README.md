# pRF analysis for Lyon MEG project

this is a bit cludgy, but this collection of scripts runs post-fmriprep preprocessing, as well as both grid and iterative prf fitting.
It isn't using all the newest software, but we're aiming for a set of simple scripts that will run the required analysis. The settings for the analysis are in `settings.yml`, where any and all settings should end up. All files that run analyses will read this `yml` file and internalize these settings as a `settings` dictionary.

# To prepare

- run `python update_cortex.py`, which will update the default locations of the pycortex database and colormaps. 
- convert all your data to BIDS

# Pipeline Description

1.  One runs fmriprep and mriqc on the BIDS data. 
2. freesurfer anatomies need to be checked and edited by hand.
3. freesurfer surfaces need to be cut and flattened following [the instructions on our wiki](https://github.com/VU-Cog-Sci/wiki/wiki/Cutting-inflated-FreeSurfer-surfaces-to-obtain-a-flattened-cortical-surface), making sure that you save the flattened files as instructed.
4. run `python prepost.py --subject $SJNR`, where `$SJNR` is the subject number.
   - This will perform temporal filtering and percent-signal change conversion of the fmriprep output `.nii.gz` files.
   - It will also average (median) across all repetitions of the runs, and
   - Copy over files like example epi file and brainmask file.
5. run `python run_surf.py --subject $SJNR --import_subject --import_session`, which will import your subject's anatomy into pycortex's folders.
6. run `python run_fit.py --subject $SJNR --grid`, which will perform grid fitting using pre-calculated regressors, and save out the results to `grid_rsq.nii.gz` and `grid_params.nii.gz`
7. run `python run_surf.py --subject $SJNR --flatmap_pRFs --webgl_pRFs 1`, which will perform visualizations based on `grid_rsq.nii.gz` and `grid_params.nii.gz`

# Other stuff

There are some additional files that aid analysis. First, some shell scripts to run fmriprep and visualize results. Then, an ipython notebook to create the design matrix used in the analysis. This latter design matrix should not have to change after we've created it correctly.
 
