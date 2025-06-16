#!/bin/bash

#this pipeline refer to https://github.com/dpeerlab/SEACells
#Lin Lyu 2025.3.25

set -e

#create environment, mamba can be replaced by conda 
mamba create -n seacells python=3.8
mamba activate seacells

#set mirror for pip to accelerate downloadig, conda/mamba do not include seacells package
pip config set global.index-url https://mirrors.aliyun.com/pypi/simple/
pip install seacells

#install jupyterlab for learning of 'seacells'
# notice!!! jupyterlab must be installed as seacells need it to run
mamba install jupyterlab

#download notebook file for learning of 'seacells'
#https://github.com/dpeerlab/SEACells/blob/main/notebooks/SEACell_computation.ipynb
