from prf_fit import *
import argparse
# import cortex
import numpy as np
import nibabel as nb
import shutil
import glob
import os
import sys
import yaml
import multiprocessing

from IPython.core.debugger import set_trace
# call as `set_trace()` for interactive debugging

from utils import *


###################################################################
###
# start by parsing arguments
###
###################################################################

parser = argparse.ArgumentParser(
    description='Run the pRF fitting on given subject')
parser.add_argument('--subject', type=int, default=1,
                    help='BIDS integer for this subject')
parser.add_argument('--grid', type=str2bool, nargs='?',
                    const=True, default=False,
                    help='whether to perform grid fitting')
parser.add_argument('--iter', type=str2bool, nargs='?',
                    const=True, default=False,
                    help='whether to perform iterative fitting')
parser.set_defaults(grid=False, iter=False)
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
    print('cannot find base_dir in "{uname}", exiting'.format(uname=uname))
    sys.exit()


###################################################################
###
# load data, mask, and example epi, along with FS T1 file
###
###################################################################


input_nii = nb.load(os.path.join(
    base_dir, settings['result_dir'], settings['median_psc_filename'].format(subject=subject, session=settings['session'])))
input_data = input_nii.get_data()
data = input_data.reshape((-1, input_data.shape[-1]))
mask_nii = nb.load(os.path.join(
    base_dir, settings['result_dir'], settings['cortical_mask_filename'].format(subject=subject, session=settings['session'])))
mask = mask_nii.get_data().astype(bool)
example_epi_file = os.path.join(
    base_dir, settings['result_dir'], settings['example_epi_filename'].format(subject=subject, session=settings['session']))

fs_T1_mgz_file = os.path.join(
    base_dir, settings['result_dir'], 'freesurfer/sub-{subject}/mri/T1.mgz'.format(subject=subject))
fs_T1_nii_file = fs_T1_mgz_file.replace('.mgz', '.nii.gz')
if not os.path.isfile(fs_T1_nii_file):
    os.system('mri_convert {mgz} {nii}'.format(
        mgz=fs_T1_mgz_file, nii=fs_T1_nii_file))

# load design matrix
visual_dm = np.load(os.path.join(
    base_dir, 'code', settings['dm_file'])).T


###################################################################
###
# set up fitting parameters from settings
###
###################################################################


# Fit: define search grids
x_grid_bound = (-settings["max_eccen"], settings["max_eccen"])
y_grid_bound = (-settings["max_eccen"], settings["max_eccen"])
sigma_grid_bound = (settings["min_size"], settings["max_size"])
n_grid_bound = (settings["min_n"], settings["max_n"])
grid_steps = settings["grid_steps"]

# Fit: define search bounds, these are more lenient than the grid search versions
x_fit_bound = (-settings["max_eccen"]*2, settings["max_eccen"]*2)
y_fit_bound = (-settings["max_eccen"]*2, settings["max_eccen"]*2)
sigma_fit_bound = (1e-6, 1e2)
n_fit_bound = (1e-6, 2)
beta_fit_bound = (-1e6, 1e6)
baseline_fit_bound = (-1e6, 1e6)

if settings["fit_model"] == 'gauss' or settings["fit_model"] == 'gauss_sg':
    bound_grids = (x_grid_bound, y_grid_bound, sigma_grid_bound)
    bound_fits = (x_fit_bound, y_fit_bound, sigma_fit_bound,
                  beta_fit_bound, baseline_fit_bound)
elif settings["fit_model"] == 'css' or settings["fit_model"] == 'css_sg':
    bound_grids = (x_grid_bound, y_grid_bound, sigma_grid_bound, n_grid_bound)
    bound_fits = (x_fit_bound, y_fit_bound, sigma_fit_bound,
                  n_fit_bound, beta_fit_bound, baseline_fit_bound)

###################################################################
###
# intitialize prf analysis object
###
###################################################################

prf = PRF_fit(data=data[mask.ravel()].astype(np.float32),
              fit_model=settings["fit_model"],
              visual_design=visual_dm,
              screen_distance=settings["screen_distance"],
              screen_width=settings["screen_width"],
              scale_factor=1/4.0,
              tr=settings["TR"],
              bound_grids=bound_grids,
              grid_steps=grid_steps,
              bound_fits=bound_fits,
              n_jobs=N_PROCS,
              sg_filter_window_length=settings["sg_filt_window_length"],
              sg_filter_polyorder=settings["sg_filt_polyorder"],
              sg_filter_deriv=settings["sg_filt_deriv"],
              )


###################################################################
###
# grid fit stage
###
###################################################################


if args.grid == 1:
    # will need to move/delete this file for new predictions
    prediction_file = os.path.join(
        base_dir, 'derivatives', 'out', 'pp', 'predictions.npy')
    if os.path.isfile(prediction_file):
        prf.load_grid_predictions(prediction_file=prediction_file)
    else:
        prf.make_predictions(out_file=prediction_file)

    print('starting grid prf fitting on {nvox} voxels using {nreg} regressors'.format(
        nvox=np.prod(mask_nii.shape), nreg=prf.predictions.shape[-1]))

    # actual fit
    prf.grid_fit()

    rsq_output = np.zeros(mask_nii.shape)
    rsq_output[mask] = prf.gridsearch_r2
    rsq_out_nii = nb.Nifti1Image(
        rsq_output, affine=mask_nii.affine, header=mask_nii.header)
    rsq_out_nii.to_filename(
        os.path.join(base_dir, 'derivatives', 'out', 'pp', 'sub-{subject}'.format(subject=subject), 'ses-{ses}'.format(ses=settings['session']), 'grid_rsq.nii.gz'))

    params_output = np.zeros(list(mask_nii.shape) +
                             [prf.gridsearch_params.shape[0]])
    params_output[mask] = prf.gridsearch_params.T
    params_out_nii = nb.Nifti1Image(
        params_output, affine=mask_nii.affine, header=mask_nii.header)
    params_out_nii.to_filename(
        os.path.join(base_dir, 'derivatives', 'out', 'pp', 'sub-{subject}'.format(subject=subject), 'ses-{ses}'.format(ses=settings['session']), 'grid_params.nii.gz'))


###################################################################
###
# iterative fitting stage
###
###################################################################

if args.iter == 1:
    rsq_grid_nii = nb.load(
        os.path.join(base_dir, 'derivatives', 'out', 'pp', 'sub-{subject}'.format(subject=subject), 'ses-{ses}'.format(ses=settings['session']), 'grid_rsq.nii.gz'))
    # mask this result with the rsq_threshold for the iterative fit.
    rsq_mask = rsq_grid_nii.get_data() > settings['rsq_threshold']
    prf.gridsearch_r2 = rsq_grid_nii.get_data()[rsq_mask]

    params_grid_nii = nb.load(
        os.path.join(base_dir, 'derivatives', 'out', 'pp', 'sub-{subject}'.format(subject=subject), 'ses-{ses}'.format(ses=settings['session']), 'grid_params.nii.gz'))
    prf.gridsearch_params = params_grid_nii.get_data()[rsq_mask]

    prf.data = input_data[rsq_mask]

    print('starting iterative prf fitting on {nvox} voxels'.format(
        nvox=rsq_mask.sum()))

    # actual fit
    prf.n_units = np.sum(rsq_mask)
    fit_results = prf.iterative_fit()

    output_data = np.zeros(
        (list(rsq_grid_nii.shape) + [fit_results.shape[-1]]))
    output_data[rsq_mask] = fit_results
    params_out_nii = nb.Nifti1Image(
        output_data, affine=mask_nii.affine, header=mask_nii.header)
    params_out_nii.to_filename(
        os.path.join(base_dir, 'derivatives', 'out', 'pp', 'sub-{subject}'.format(subject=subject), 'ses-{ses}'.format(ses=settings['session']), 'iter_params.nii.gz'))

    rsq_out_nii = nb.Nifti1Image(
        params_out_nii[..., -1], affine=mask_nii.affine, header=mask_nii.header)
    rsq_out_nii.to_filename(
        os.path.join(base_dir, 'derivatives', 'out', 'pp', 'sub-{subject}'.format(subject=subject), 'ses-{ses}'.format(ses=settings['session']), 'iter_rsq.nii.gz'))
