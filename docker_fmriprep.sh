bids_folder=/scratch/2019/visual/prf_lyon/bids
derivatives_folder=/scratch/2019/visual/prf_lyon/derivatives
license_file=/tank/tkn219/software/license.txt
subject_nr=01

docker run -ti --rm \
-v ${bids_folder}:/data:ro \
-v ${derivatives_folder}:/out \
-v ${license_file}:/license.txt \
poldracklab/fmriprep:latest \
/data /out/out participant --fs-license-file /license.txt \
--participant_label ${subject_nr} \
--omp-nthreads 30 --nthreads 30 --skip_bids_validation \
--ignore slicetiming \
--output-spaces T1w MNI152NLin2009cAsym fsaverage fsnative --use-syn-sdc --write-graph