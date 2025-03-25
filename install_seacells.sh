#!/bin/bash

set -e

#create environment
mamba create -n seacells python=3.8
mamba activate seacells

#set mirror for pip to accelerate downloadig, conda/mamba do not include seacells package
pip config set global.index-url https://mirrors.aliyun.com/pypi/simple/
pip install seacells
