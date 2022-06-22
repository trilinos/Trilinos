#!/bin/bash
#SBATCH -A cfd116
#SBATCH -J push-results
#SBATCH -t 00:05:00
#SBATCH -N 1
#SBATCH --mail-type=BEGIN,END,INVALID_DEPEND

#proxies to allow git operations from compute nodes
export all_proxy=socks://proxy.ccs.ornl.gov:3128/
export ftp_proxy=ftp://proxy.ccs.ornl.gov:3128/
export http_proxy=http://proxy.ccs.ornl.gov:3128/
export https_proxy=http://proxy.ccs.ornl.gov:3128/
export no_proxy='localhost,127.0.0.0/8,*.ccs.ornl.gov,*.ncrc.gov',*.olcf.ornl.gov

cd $HOME/crusher/sources/trilinos/Trilinos-muelu
SHA=`git rev-parse HEAD | head -c 8`

cd $HOME/spock/trilinos-kernel-performance/CrusherData
TODAY=`date +%F`
git pull
git add .
git commit -m "${TODAY} ${SHA}: crusher results"
git push
