#!/bin/bash
## Author: Bansho Masutani<banmasutani@gmail.com>
set -ue

## Install h5py and ont_fast5_api
cd ../
python3 -m pip install h5py
git clone https://github.com/nanoporetech/ont_fast5_api.git
cd ont_fast5_api
make develop
cd ../

## Cloning dependent libraries
for repository in dtw utility knn_predictor histogram_minimizer squiggler fast5wrapper
do
    git clone git@github.com:ban-m/${repository}.git
done
cd dyss
git clone https://github.com/nanoporetech/kmer_models.git
python3 -m pip install -r requirements.txt
cargo build --release
