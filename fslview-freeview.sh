subject_nr=01
subjects_folder=/Users/knapen/Downloads/prf_lyon/derivatives/out/freesurfer
freeview -v ${subjects_folder}/sub-${subject_nr}/mri/T1.mgz \
${subjects_folder}/sub-${subject_nr}/mri/wm.mgz:visible=0 \
${subjects_folder}/sub-${subject_nr}/mri/brainmask.mgz:visible=0 \
${subjects_folder}/sub-${subject_nr}/mri/T2.mgz:colormap=jet:colorscale=100,600 \
-f ${subjects_folder}/sub-${subject_nr}/surf/lh.white:edgecolor=white \
${subjects_folder}/sub-${subject_nr}/surf/lh.pial:edgecolor=yellow \
${subjects_folder}/sub-${subject_nr}/surf/lh.woT2.pial:edgecolor=red \
${subjects_folder}/sub-${subject_nr}/surf/rh.white:edgecolor=white \
${subjects_folder}/sub-${subject_nr}/surf/rh.pial:edgecolor=yellow \
${subjects_folder}/sub-${subject_nr}/surf/rh.woT2.pial:edgecolor=red



pp_dir=/Users/knapen/Downloads/prf_lyon/derivatives/out/pp/sub-${subject_nr}/ses-01/
anat_dir=/Users/knapen/Downloads/prf_lyon/derivatives/out/fmriprep/sub-${subject_nr}/anat/

fsleyes -b ${pp_dir}sub-${subject_nr}_task-prf_acq-median_T1w_desc-preproc_bold.nii.gz\
{anat_dir}sub-${subject_nr}_desc-preproc_T1w.nii.gz