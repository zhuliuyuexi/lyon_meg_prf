import re
import os
import glob
import json
import sys
import argparse
import yaml

from utils import *

###################################################################
###
# start by parsing arguments
###
###################################################################

parser = argparse.ArgumentParser(
    description='Run the singularity images on a gives set of subjects. NOTE that paths in this file are hard-coded due to timeconstraints.')
parser.add_argument('--subjects', metavar='N', type=int, nargs='+',
                    help='BIDS integer for this subject')

parser.add_argument('--fmriprep', type=str2bool, nargs='?',
                    const=True, default=False,
                    help='whether to run fmriprep')
parser.add_argument('--mriqc_subject', type=str2bool, nargs='?',
                    const=True, default=False,
                    help='whether to run mriqc on the subject')
parser.add_argument('--mriqc_group', type=str2bool, nargs='?',
                    const=True, default=False,
                    help='whether to run mriqc on the group')


# parser.set_defaults(fmriprep=False, mriqc_subject=False, mriqc_group=False)

args = parser.parse_args()
subjects = [str(s).zfill(2) for s in args.subjects]

with open('settings.yml', 'r') as f:
    settings = yaml.load(f)

batchdir = '/home/rasa.gulbinaite/batch/'
os.chdir(batchdir)

#############################################################################################
#############################################################################################
# FMRIPREP
#############################################################################################
#############################################################################################

if args.fmriprep:
    for subject in subjects:

        batch_string = """#!/bin/bash
#SBATCH -p normal -t 20:00:00 -N 1
# job requires at most 100 hours, 0 minutes
#     and 0 seconds wallclock time

PYTHONPATH="" singularity run -B /mnt/data/rasa,/scratch \
/mnt/data/rasa/software/poldracklab_fmriprep_1.4.1-2019-07-09-412a69224405.simg \
/mnt/data/rasa/prf_lyon/bids/ /mnt/data/rasa/prf_lyon/derivatives/out/ participant \
--participant-label sub-$SJ_NR --output-spaces T1w MNI152NLin2009cAsym fsaverage fsnative \
--use-syn-sdc --write-graph --nthreads 15 --omp-nthreads 15 --low-mem --fs-license-file /home/rasa.gulbinaite/software/freesurfer/license.txt \
--ignore slicetiming -w /scratch --skip_bids_validation 

wait          # wait until programs are finished
        """

        working_string = batch_string.replace('$SJ_NR', subject)

        js_name = os.path.join(batchdir, str(
            subject).zfill(2) + '_lyon-prf_slurm_fmriprep.sh')
        of = open(js_name, 'w')
        of.write(working_string)
        of.close()

        print('submitting ' + js_name + ' to queue')
        print(working_string)
        os.system('sbatch ' + js_name)

#############################################################################################
#############################################################################################
# MRIQC PARTICIPANT
#############################################################################################
#############################################################################################

if args.mriqc_subject:
    for subject in subjects:

        batch_string = """#!/bin/bash
#SBATCH -p normal -t 1:00:00 -N 1

PYTHONPATH="" singularity run -B /mnt/data/rasa,/scratch \
/mnt/data/rasa/software/poldracklab_mriqc_latest-2019-04-05-f2009956414a.simg \
/mnt/data/rasa/prf_lyon/bids/ /mnt/data/rasa/prf_lyon/derivatives/out/ participant \
--participant-label sub-$SJ_NR --n_procs 15 -m bold --verbose-reports --mem_gb 32 \
--ants-nthreads 15 -w /scratch

wait          # wait until programs are finished

        """

        working_string = batch_string.replace('$SJ_NR', subject)

        working_string = batch_string

        js_name = os.path.join(batchdir, str(
            subject).zfill(2) + '_lyon-prf_slurm_mriqc.sh')
        of = open(js_name, 'w')
        of.write(working_string)
        of.close()

        print('submitting ' + js_name + ' to queue')
        print(working_string)
        os.system('sbatch ' + js_name)


#############################################################################################
#############################################################################################
# MRIQC GROUP
#############################################################################################
#############################################################################################

if args.mriqc_group:

    batch_string = """#!/bin/bash
#SBATCH -p normal -t 1:00:00 -N 1

PYTHONPATH="" singularity run -B /mnt/data/rasa,/scratch \
/mnt/data/rasa/software/poldracklab_mriqc_latest-2019-04-05-f2009956414a.simg \
/mnt/data/rasa/prf_lyon/bids/ /mnt/data/rasa/prf_lyon/derivatives/out/ participant \
--n_procs 15 -m bold --verbose-reports --mem_gb 32 --ants-nthreads 15 \
-w /scratch

wait          # wait until programs are finished
    """

    js_name = os.path.join(batchdir, 'prf_slurm_mriqc_group.sh')
    of = open(js_name, 'w')
    of.write(batch_string)
    of.close()

    print('submitting ' + js_name + ' to queue')
    print(batch_string)
    os.system('sbatch ' + js_name)
