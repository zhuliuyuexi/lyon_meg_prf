import glob
import os
import shutil
from spynoza.filtering.nodes import savgol_filter
from spynoza.conversion.nodes import percent_signal_change
import nibabel as nb
import numpy as np

###################################################################
###
# start by parsing arguments
###
###################################################################

parser = argparse.ArgumentParser(
    description='Run the pRF fitting on given subject')
parser.add_argument('--subject', type=int, default=1,
                    help='BIDS integer for this subject')
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
for s in settings['system']:
    if settings['system'][s]['uname'] in uname:
        base_dir = settings['system'][s]['base_dir']
        system = settings['system'][s]['uname']

if base_dir == None:
    sys.exit()

###################################################################
###
# load stuff
###
###################################################################

all_nii = glob.glob(os.path.join(
    base_dir, settings['result_dir'], settings['fmriprep_output_wildcard'].format(subject=subject)))

# create parallel folder structure for post-preprocessing
os.makedirs(os.path.split(os.path.join(base_dir, settings['result_dir'], settings['fmriprep_output_wildcard'].format(subject=subject)))
            [0].replace('fmriprep', 'pp'))

for niif in all_nii:
    print('filtering the data')
    sg_nii = savgol_filter(in_file=niif, polyorder=settings['sg_filt_polyorder'],
                           deriv=settings['sg_filt_deriv'], window_length=settings['sg_filt_window_length'], tr=settings['TR'])
    print('percent-signal change converting the data')
    psc_nii = percent_signal_change(
        in_file=sg_nii, func=settings['psc_method'])

    print('moving the data')
    shutil.move(sg_nii, sg_nii.replace('fmriprep', 'pp'))
    shutil.move(psc_nii, psc_nii.replace('fmriprep', 'pp'))

# just going to average the T1w space data for now.
to_be_averaged = glob.glob(os.path.join(
    base_dir, settings['average_wildcard'].format(subject=subject)))

# set up the data
to_be_medianed_data = np.zeros(
    [len(to_be_averaged)] + list(nb.load(to_be_averaged[0]).shape))

# load the data
for i, tba in enumerate(to_be_averaged):
    to_be_medianed_data[i] = nb.load(tba).get_data()

median_data = np.median(to_be_medianed_data, axis=0)

# create output image
nii_img_median_data = nb.Nifti1Image(median_data,
    affine=nb.load(to_be_averaged[0]).affine,
    header=nb.load(to_be_averaged[0]).header)

# save
nii_img_median_data.to_filename(os.path.join(os.path.split(
    to_be_averaged[0])[0], settings['median_psc_filename'].format(subject=subject)))


# now, set up all the required files in that folder
# using the last item of the nii_files, and the to_be_averaged path.
shutil.copyfile(to_be_averaged.replace('_desc-preproc_bold', '_brain-mask'),
        os.path.join(os.path.split(to_be_averaged[0])[0], settings['brainmask_filename'].format(subject=subject))))
shutil.copyfile(to_be_averaged.replace('_desc-preproc_bold', '_boldref'),
        os.path.join(os.path.split(to_be_averaged[0])[0], settings['example_epi_filename'].format(subject=subject))))
