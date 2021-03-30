#Set up other env variables
#$DATASTORE is the path to SEMS data store on sherlock
export DATASTORE="/storage/elements/sems-data-store/trilinos-performance"
export DATASTORE_MIRROR="/storage/elements/sems-data-store/trilinos-performance-son"
export WATCHR_PERF_DIR="$WORKSPACE/VortexData"
export TRILINOS_SRC="$WORKSPACE/Trilinos"
#Remove script dir, if it exists
rm -rf VortexPerfScripts
#Get a fresh copy
scp -rp jenkins@sherlock.sandia.gov:$DATASTORE/Scripts/vortex VortexPerfScripts || true
#Copy the newest performance report tarball from data store
scp -p jenkins@sherlock.sandia.gov:$DATASTORE/VortexData.tar.gz VortexData.tar.gz || true
#Unpack it, is directory "VortexData"
rm -rf VortexData
tar -xf VortexData.tar.gz
#Create build dir, if not already there.
#Don't need to clear it out, since script will reconfigure everything.
#Even with reconfigure, this is an incremental build, saving a lot of time.
mkdir -p build
#Configure and build on a compute node, then run run tests on 4 compute nodes.
$WORKSPACE/VortexPerfScripts/launch.sh
#Repack the up to date tarball of data, and put back in datastore
tar -czf VortexData.tar.gz VortexData
scp -p VortexData.tar.gz jenkins@sherlock.sandia.gov:$DATASTORE/VortexData.tar.gz || true
scp -p VortexData.tar.gz jenkins@watson.sandia.gov:$DATASTORE_MIRROR/VortexData.tar.gz || true
