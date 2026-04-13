# Simple Trilinos build configuration using Trilinos clang-19 container

This build configuration enables a lot of Trilinos without using GenConfig.  The
purpose is to be a very simple configuration that is easy to understand and work
with.

To do a configuration of Trilinos, create a build directory and sym-link in these files:

```bash
cd Trilinos/
mkdir -p BUILDS/clang-19-simple
cd BUILDS/clang-19-simple/
ln -s ../../sampleScripts/simple-genconfig-scripts/clang-simple/* .
```

Then to configure Trilinos like:

```bash
. load-env.sh
./do-configure -DTrilinos_ENABLE_Tpetra=ON  # Or whatever packages you want
```

Then you are ready to build, run tests, etc. with:

```bash
ninja -j $(($(nproc)/2))
ctest -j $(($(nproc)/2))
```
