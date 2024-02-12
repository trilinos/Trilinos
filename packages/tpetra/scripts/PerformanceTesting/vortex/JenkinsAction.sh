set -e

export WATCHR_PERF_DIR="$WORKSPACE/VortexData"
export TRILINOS_SRC="$WORKSPACE/Trilinos"

#Remove performance dirs, if they exist
rm -rf PerfScripts
rm -rf VortexData

#make all the new directories
mkdir PerfScripts
mkdir VortexData

cp tpetra_performance/scripts/vortex/* PerfScripts/
cp tpetra_performance/scripts/timingPercentages/* PerfScripts/

cd $WORKSPACE
#Create build dir, if not already there.
#Don't need to clear it out, since script will reconfigure everything.
#Even with reconfigure, this is an incremental build, saving a lot of time.
mkdir -p build
#Configure and build on a compute node, then run run tests on 4 compute nodes.
$WORKSPACE/PerfScripts/launch.sh

#calculate percentages then stick XMLs into the repo and push
python $WORKSPACE/PerfScripts/timing_extractor.py -p VortexData

if [ "$DRYRUN" != "ON" ]
then
  cp VortexData/*.xml tpetra_performance/VortexData/
  cd tpetra_performance
  git add VortexData/*
  git commit -m "updating VortexData results" 
  git pull --rebase
  git push
fi
