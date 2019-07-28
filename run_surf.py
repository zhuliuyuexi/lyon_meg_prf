from prf_fit import *
import argparse
import cortex
import numpy as np
import nibabel as nb
import shutil
import glob
import os
import sys
import yaml

import matplotlib.pyplot as plt
import matplotlib.colors as colors


###################################################################
###
# start by parsing arguments
###
###################################################################

parser = argparse.ArgumentParser(
    description='Run the pRF fitting on given subject')
parser.add_argument('--subject', type=int, default=1,
                    help='BIDS integer for this subject')

parser.add_argument('--import_subject', type=int, default=0,
                    help='whether to import subject into pycortex')
parser.add_argument('--import_session', type=int, default=0,
                    help='whether to import session into pycortex subject')

parser.add_argument('--flatmap_pRFs', type=int, default=0,
                    help='visualize pRF results in pdfs')
parser.add_argument('--webgl_pRFs', type=int, default=0,
                    help='visualize pRF results in web browser')
args = parser.parse_args()
subject = str(args.subject).zfill(2)

with open('settings.yml', 'r') as f:
    settings = yaml.load(f)


###################################################################
###
# check the system we're at from the settings
###
###################################################################

base_dir = None
uname = str(os.uname())
possible_systems = list(settings['systems'].keys())
possible_unames = [settings['systems'][s]['uname'] for s in possible_systems]
for ps, pu in zip(possible_systems, possible_unames):
    if pu in uname:
        base_dir = settings['systems'][ps]['base_dir']
        system = settings['systems'][ps]['uname']
        N_PROCS = settings['systems'][ps]['threads']

if base_dir == None:
    sys.exit()

###################################################################
###
# load data, mask, and example epi, along with FS T1 file
###
###################################################################

fs_dir = os.path.join(base_dir, settings['result_dir'], 'freesurfer')
# cortex.config.update(
# {'default_filestore': os.path.join(os.path.split(fs_dir)[0], 'cortex')})

fs_T1_mgz_file = os.path.join(
    base_dir, settings['result_dir'], 'freesurfer/sub-{subject}/mri/T1.mgz'.format(subject=subject, session=settings['session']))
fs_T1_nii_file = fs_T1_mgz_file.replace('.mgz', '.nii.gz')
if not os.path.isfile(fs_T1_nii_file):
    os.system('mri_convert {mgz} {nii}'.format(
        mgz=fs_T1_mgz_file, nii=fs_T1_nii_file))

example_epi_file = os.path.join(
    base_dir, settings['result_dir'], settings['example_epi_filename'].format(subject=subject, session=settings['session']))

###################################################################
###
# import the subject's anatomy.
###
###################################################################

if args.import_subject == 1:
    cortex.freesurfer.import_subj(
        subject='sub-{subject}'.format(subject=subject),
        sname='sub-{subject}'.format(subject=subject),
        freesurfer_subject_dir=fs_dir)

    if not os.path.isfile(os.path.join(fs_dir, 'sub-{subject}'.format(subject=subject), 'surf', 'rh.full.flat.patch.3d')):
        print("""No flattened patches with name 'rh-lh.full.flat.patch.3d' found. Will not import flat, but we will not be able to draw ROIs without.
        Please see https://github.com/VU-Cog-Sci/wiki/wiki/Cutting-inflated-FreeSurfer-surfaces-to-obtain-a-flattened-cortical-surface
        for a tutorial on how to cut and flatten inflated cortex. Please remember to run this code again once you've created the flattened surfaces""")
    else:
        cortex.freesurfer.import_flat(subject='sub-{subject}'.format(subject=subject),
                                      patch='full',
                                      sname='sub-{subject}'.format(
                                          subject=subject),
                                      freesurfer_subject_dir=fs_dir)

###################################################################
###
# import the functional data as session.
###
###################################################################

if args.import_session == 1:

    mtx_file = os.path.join(base_dir, settings['result_dir'], 'pp/sub-{subject}/ses-{session}/flirt.mtx'.format(
        subject=subject, session=settings['session']))

    # first, we register the subject using bbregister
    os.environ.update({'SUBJECTS_DIR': fs_dir})
    bbr_cmd = '''bbregister --s sub-{subject} --mov {example_epi} \
    --reg {reg_file} --fslmat {mtx_file} --init-fsl --bold'''
    bbr_cmd = bbr_cmd.format(
        example_epi=example_epi_file,
        reg_file=os.path.join(base_dir, settings['result_dir'], 'pp/sub-{subject}/ses-{session}/register.dat'.format(
            subject=subject, session=settings['session'])),
        mtx_file=mtx_file,
        subject=subject
    )

    # retain the existing matrix
    if not os.path.isfile(mtx_file):
        print("Running for registration: " + bbr_cmd)
        os.system(bbr_cmd)
    else:
        print("Working with existing coregistration matrix file {mtx_file}".format(
            mtx_file=mtx_file))

    # then, we use the generated flirt matrix to import the subject into the session.
    xfm_data = np.loadtxt(mtx_file)
    xfm = cortex.xfm.Transform.from_fsl(xfm=xfm_data,
                                        func_nii=example_epi_file,
                                        anat_nii=fs_T1_nii_file)
    # # Save as pycortex 'coord' transform
    xfm.save('sub-{subject}'.format(subject=subject), 'fmriprep_T1', 'coord')

    ###################################################################
    ###
    # and, while we're at it, create a cortical mask and save it.
    ###
    ###################################################################

    mask = cortex.get_cortical_mask(
        'sub-{subject}'.format(subject=subject), 'fmriprep_T1', type='thick')

    epi_nii = nb.load(example_epi_file)
    dims = epi_nii.shape

    mask_img = nb.Nifti1Image(dataobj=mask.transpose(
        (2, 1, 0)), affine=epi_nii.affine, header=epi_nii.header)

    cortical_mask_file = os.path.join(base_dir, settings['result_dir'], settings['cortical_mask_filename'].format(
        subject=subject,
        session=settings['session']))
    mask_img.to_filename(cortical_mask_file)


###################################################################
###
# create webgl representation of prf results
###
###################################################################

if args.webgl_pRFs == 1 or args.webgl_pRFs == 2 or args.flatmap_pRFs == 1:
    params = nb.load(os.path.join(base_dir, 'derivatives', 'out', 'pp', 'sub-{subject}'.format(
        subject=subject), 'ses-{ses}'.format(ses=settings['session']), 'grid_params.nii.gz'))
    p_data = params.get_data()
    rsq = nb.load(os.path.join(base_dir, 'derivatives', 'out', 'pp', 'sub-{subject}'.format(
        subject=subject), 'ses-{ses}'.format(ses=settings['session']), 'grid_rsq.nii.gz')).get_data()

    # now construct polar angle and eccentricity values
    complex_location = p_data[..., 0] + p_data[..., 1] * 1j
    polar_angle = np.angle(complex_location)
    eccentricity = np.abs(complex_location)

    size = p_data[..., 2]
    # baseline and beta were swapped when running the first fits.
    baseline = p_data[..., 3]
    beta = p_data[..., 4]

    polar_angle_n = (polar_angle + np.pi) / (np.pi * 2.0)

    # make discrete angles for clarity
    angle_offset = 0.1
    polar_angle_n = np.fmod(polar_angle_n+angle_offset, 1.0)

    # convert angles to colors, using correlations as weights
    hsv = np.zeros(list(polar_angle_n.shape) + [3])
    hsv[..., 0] = polar_angle_n  # angs_discrete  # angs_n
    # np.sqrt(rsq) #np.ones_like(rsq)  # np.sqrt(rsq)
    hsv[..., 1] = np.ones_like(rsq)
    hsv[..., 2] = np.ones_like(rsq)  # np.sqrt(rsq)# np.ones_like(rsq)

    alpha_mask = (rsq <= settings['rsq_threshold']).T
    alpha = np.sqrt(rsq).T * 5  # for graded rsq viz
    # alpha[alpha_mask] = 0
    alpha = np.ones(alpha.shape)
    alpha[alpha_mask] = 0

    rgb = colors.hsv_to_rgb(hsv)

    # use the above calculations to create pycortex representations
    vrgba = cortex.VolumeRGB(red=rgb[..., 0].T,
                             green=rgb[..., 1].T,
                             blue=rgb[..., 2].T,
                             subject='sub-{subject}'.format(
        subject=subject),
        alpha=alpha,
        xfmname='fmriprep_T1')
    vecc = cortex.Volume2D(eccentricity.T, rsq.T, 'sub-{subject}'.format(
        subject=subject), 'fmriprep_T1',
        vmin=0, vmax=10,
        vmin2=settings['rsq_threshold'], vmax2=1.0, cmap='BROYG_2D')
    vsize = cortex.Volume2D(size.T, rsq.T, 'sub-{subject}'.format(
        subject=subject), 'fmriprep_T1',
        vmin=0, vmax=10,
        vmin2=settings['rsq_threshold'], vmax2=1.0, cmap='BROYG_2D')
    vbeta = cortex.Volume2D(beta.T, rsq.T, 'sub-{subject}'.format(
        subject=subject), 'fmriprep_T1',
        vmin=-2.5, vmax=2.5,
        vmin2=settings['rsq_threshold'], vmax2=1.0, cmap='RdBu_r_alpha')
    vbaseline = cortex.Volume2D(baseline.T, rsq.T, 'sub-{subject}'.format(
        subject=subject), 'fmriprep_T1',
        vmin=-1, vmax=1,
        vmin2=settings['rsq_threshold'], vmax2=1.0, cmap='RdBu_r_alpha')
    vrsq = cortex.Volume2D(rsq.T, rsq.T, 'sub-{subject}'.format(
        subject=subject), 'fmriprep_T1',
        vmin=0, vmax=0.8,
        vmin2=settings['rsq_threshold'], vmax2=1.0, cmap='fire_alpha')

    ds = cortex.Dataset(polar=vrgba, ecc=vecc, size=vsize,
                        amplitude=vbeta, baseline=vbaseline, rsq=vrsq)
    ds.save(os.path.join(base_dir, 'derivatives', 'out', 'pp', 'sub-{subject}'.format(
        subject=subject), 'ses-{ses}'.format(ses=settings['session']), 'pycortex_ds.h5'))

if args.webgl_pRFs == 1:
    cortex.webgl.make_static(outpath=os.path.join(base_dir, 'derivatives', 'out', 'pp', 'sub-{subject}'.format(
        subject=subject), 'ses-{ses}'.format(ses=settings['session']), 'webgl'), data=ds, recache=True)
elif args.webgl_pRFs == 2:
    cortex.webgl.show(data=ds, recache=True, port=12001)

if args.flatmap_pRFs == 1:
    for vol, name in zip([vrgba, vecc, vsize, vbeta, vbaseline, vrsq],
                         ['polar', 'ecc', 'size', 'amplitude', 'baseline', 'rsq']):
        fig = cortex.quickflat.make_figure(vol, with_dropout=True)
        plt.savefig(outpath=os.path.join(base_dir, 'derivatives', 'out', 'pp', 'sub-{subject}'.format(
            subject=subject), 'ses-{ses}'.format(ses=settings['session']), name+'.pdf'))
        plt.close('all')
