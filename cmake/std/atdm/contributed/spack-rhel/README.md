
* <a href="#spack-rhel-environment">Spack RHEL Environment</a>


### Spack RHEL Environment

The env 'spack-rhel' should work on any Red Hat Enterprise Linux (RHEL) (and
perhaps many other Linux systems) that have the SNL ATDM Spack modules
installed on them.  See the [installation
documentation](https://gitlab.sandia.gov/atdm/atdm-spack-scripts/blob/master/README.md).
**WARNING:** This Spack env is still under development and may change in the
future.

Once logged onto a Linux machine with the SNL ATDM Spack modules installed,
one can directly configure, build, and run tests using the `spack-rhel` env.
For example, to configure, build and run the tests for `MueLu` one would clone
Trilinos on the `develop` branch and then do the following:


```
$ cd <some_build_dir>/

$ source <spack-install-base-dir>/setup-env.sh

$ source $TRILINOS_DIR/cmake/std/atdm/load-env.sh spack-rhel-gnu-openmp-opt

$ cmake \
  -GNinja \
  -DTrilinos_CONFIGURE_OPTIONS_FILE:STRING=cmake/std/atdm/ATDMDevEnv.cmake \
  -DTrilinos_ENABLE_TESTS=ON -DTrilinos_ENABLE_MueLu=ON \
  $TRILINOS_DIR

$ make NP=16

$ ctest -j8
```

One can also run the same build and tests using the <a
href="#checkin-test-atdmsh">checkin-test-atdm.sh</a> script as:

```
$ cd <some_build_dir>/
$ ln -s $TRILINOS_DIR/cmake/std/atdm/checkin-test-atdm.sh .
$ env ATDM_CHT_DEFAULT_ENV=spack-rhel-default \
  ./checkin-test-atdm.sh spack-rhel-gnu-openmp-opt \
  --enable-packages=MueLu \
  --local-do-all
```

NOTE: Above one must set `ATDM_CHT_DEFAULT_ENV=spack-rhel-default` in the env
when passing in `all` in order for it to select the correct set of supported
builds for the `spack-rhel` env and also to load the correct env to find
Python, etc.


* `spack-rhel/`: RHEL (and likely other Linux) systems with the SNL ATDM Spack modules installed.
