# TribitsOldSimpleExampleApp

Simple example project `TribitsOldSimpleExampleApp` is a raw CMake project
that pulls in libraries from a few packages from `TribitsExampleProject` using
just `find_package(TribitsExProj REQUIRED COMPONENTS ...)` but uses the old
specification for using an installed TriBITS project (i.e. through CMake
variables `<Project>_LIBRARIES`, `<Project>_INCLUDE_DIRS`, and
`<Project>_TPL_INCLUDE_DIRS`).

After building and installing TribitsExampleProject under
`<upstreamInstallDir>`, then configure, build, and test
`TribitsSimpleExampleApp` with:

```
  cmake \
    -DCMAKE_PREFIX_PATH=<upstreamInstallDir> \
    <base-dir>/TribitsOldSimpleExampleApp

  make

  ctest
```

This project can be instructed to use the old deprecated targets by adding the
CMake cache var:

```
  -D TribitsOldSimpleExApp_USE_DEPRECATED_TARGETS=ON
```

This will generated deprecated warnings with newer versions of TriBITS.
