#
# Purge the modules and unload the current Trilinos SEMS development
# environment
#
# USAGE:
#
#   $ source unload_sems_dev_env.sh
#
# NOTE: It would be nice to just unload the loaded SEMS modules but the way
# the SEMS modules are set up makes that impossible.  This is because loading
# a given SEMS module can trigger the loading of another module which will
# linger.  Then that lingering module can conflict with the loading of other
# modules.  One example was loading 'sems-clang/3.6.1' would load
# 'sems-gcc/4.8.4' but unloading 'sems-clang/3.6.1' would not unload
# 'sems-gcc/4.8.4'.  That was fixed but there are other examples like that.
# Therefore, the solution was to just purge the modules.
#

module purge 
