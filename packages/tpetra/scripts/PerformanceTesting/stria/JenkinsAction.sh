#Set up other env variables
#$DATASTORE is the path to SEMS data store on sherlock
export DATASTORE="/storage/elements/sems-data-store/trilinos-performance"
export DATASTORE_MIRROR="/storage/elements/sems-data-store/trilinos-performance-son"
export WATCHR_PERF_DIR="$WORKSPACE/StriaData"
export TRILINOS_SRC="$WORKSPACE/Trilinos"
#Remove script dir, if it exists
rm -rf StriaPerfScripts
#Get a fresh copy
scp -rp jenkins@sherlock.sandia.gov:$DATASTORE/Scripts/stria StriaPerfScripts || true
scp -rp jenkins@sherlock.sandia.gov:$DATASTORE/Scripts/timingPercentages StriaPerfScripts || true
#Copy the newest performance report tarball from data store
scp -p jenkins@sherlock.sandia.gov:$DATASTORE/StriaData.tar.gz StriaData.tar.gz || true
#Unpack it, is directory "VortexData"
rm -rf StriaData
tar -xf StriaData.tar.gz
#Create build dir, if not already there.
#Don't need to clear it out, since script will reconfigure everything.
#Even with reconfigure, this is an incremental build, saving a lot of time.
mkdir -p build
#Configure and build on a compute node, then run run tests on 4 compute nodes.
$WORKSPACE/StriaPerfScripts/launch.sh
#put the time percentages in each xml file
python $WORKSPACE/StriaPerfScripts/timing_extractor.py -p StriaData
#Repack the up to date tarball of data, and put back in datastore
tar -czf StriaData.tar.gz StriaData
scp -p StriaData.tar.gz jenkins@sherlock.sandia.gov:$DATASTORE/StriaData.tar.gz || true
scp -p StriaData.tar.gz jenkins@watson.sandia.gov:$DATASTORE_MIRROR/StriaData.tar.gz || true

