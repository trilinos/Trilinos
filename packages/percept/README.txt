
****************************
* Building with ATDM scripts
****************************

mkdir BUILD
cd BUILD
export TRILINOS_DIR=/gpfs1/bcarnes/percept_trilinos
. $TRILINOS_DIR/packages/percept/scripts/set_env_gcc.sh 
$TRILINOS_DIR/packages/percept/scripts/build_gcc.sh 
ninja

****************************
* Snapshot in a new Percept
****************************

Temporarily, we currently just snapshot the src directory from percept
on github into this package's src directory. There is a snapshot.py
script in the scripts directory.

cd Trilinos/packages/percept
./scripts/snapshot.py ~/percept/src .

## for experimenting in non-clean repo
## ./scripts/snapshot.py -n --no-validate-repo  ~/junk/percept/src .

