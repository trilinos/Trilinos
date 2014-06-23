#
# Set TRIBITS_DIR in the env then call this!
#

TRIBITS_DIR_ABS=$(readlink -f $TRIBITS_DIR)

EXTRA_ARGS=$@
cmake \
-DTribitsExProj_TRIBITS_DIR=$TRIBITS_DIR_ABS \
-DTribitsExProj_ENABLE_TESTS=ON \
$EXTRA_ARGS \
$TRIBITS_DIR_ABS/doc/examples/TribitsExampleProject

# -DCTEST_DROP_SITE=testing.sandia.gov \
