#!/usr/bin/env bash
VERSION=${1?Error: give version of fmriprep you want to build!}
docker run --privileged -t --rm \
    -v /var/run/docker.sock:/var/run/docker.sock \
    -v /home/shared/software/bids_apps:/output \
    singularityware/docker2singularity \
    poldracklab/fmriprep:$VERSION
