#Set up other env variables
#$DATASTORE is the path to SEMS data store on sherlock
export DATASTORE="/storage/elements/sems-data-store/trilinos-performance"
export DATASTORE_MIRROR="/storage/elements/sems-data-store/trilinos-performance-son"
export WATCHR_PERF_DIR="$WORKSPACE/EclipseData"
export TRILINOS_SRC="$WORKSPACE/Trilinos"
#Remove script dir, if it exists
rm -rf EclipsePerfScripts
#Get a fresh copy
scp -rp jenkins@sherlock.sandia.gov:$DATASTORE/Scripts/eclipse EclipsePerfScripts || true
scp -rp jenkins@sherlock.sandia.gov:$DATASTORE/Scripts/timingPercentages EclipsePerfScripts || true
#Copy the newest performance report tarball from data store
scp -p jenkins@sherlock.sandia.gov:$DATASTORE/EclipseData.tar.gz EclipseData.tar.gz || true
#Unpack it
rm -rf EclipseData
tar -xf EclipseData.tar.gz
#Create build dir, if not already there.
#Don't need to clear it out, since script will reconfigure everything.
#Even with reconfigure, this is an incremental build, saving a lot of time.
mkdir -p build
#Configure and build on a compute node, then run run tests on 4 compute nodes.
$WORKSPACE/EclipsePerfScripts/launch.sh
#put the time percentages in each xml file
python $WORKSPACE/EclipsePerfScripts/timing_extractor.py -p EclipseData
#Repack the up to date tarball of data, and put back in both datastores
tar -czf EclipseData.tar.gz EclipseData
scp -p EclipseData.tar.gz jenkins@sherlock.sandia.gov:$DATASTORE/EclipseData.tar.gz || true
scp -p EclipseData.tar.gz jenkins@watson.sandia.gov:$DATASTORE_MIRROR/EclipseData.tar.gz || true

