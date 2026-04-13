# Official Trilinos clang GenConfig configuration

This configuration matches the exact configuration that Trilinos uses for its
clang PR and nightly builds that run out of this container.

To do a configuration of Trilinos, on the host, get the GenCongif files with:

```bash
cd Trilinos/
./packages/framework/get_dependencies.sh
```

Then, inside of the container, create the build dir with:

```bash
cd /mounted_from_host/Trilinos/
mkdir -p BUILDS/clang-19-genconfig
cd BUILDS/clang-19-genconfig/
ln -s ../../sampleScripts/simple-genconfig-scripts/clang-genconfig/* .
```

Then, inside of the container, you can configure and build any Trilinos package with:

```bash
cd /mounted_from_host/Trilinos/BUILDS/clang-19-genconfig/
. load-env.sh [--ci-mode]  # Will drop you into a new shell unless you pass in --ci-mode
./do-configure -DTrilinos_ENABLE_Tpetra=ON  # Or any other Trilinos packages you want
ninja -j $(($(nproc)/2))
ctest -j $(($(nproc)/2))
```

You can also build and submit to the default internal SNL CDash site with:

```bash
make dashboard
```

Or, if you are working on an external machine using an externally built Trilinos container, you can use:

```bash
./do-configure-mycdash -DTrilinos_ENABLE_Tpetra=ON
make dashboard
```

(see the STDOUT for where that is going and the link to the build).
