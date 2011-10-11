#!/bin/bash

# Load the appropriate bash environment
. /opt/modules/default/init/bash

# Unload GNU modules for catamount
module load PrgEnv-pgi
module load xt-papi

# These two specific modules are required (defaults don't work)
module load fftw/3.2.2.1
module swap pgi/10.9.0

