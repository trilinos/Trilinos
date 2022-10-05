# TribitsOldSimpleExampleApp

This simple example project `TribitsOldSimpleExampleApp` is a raw CMake
project that pulls in libraries from a few packages from
[`TribitsExampleProject`](../TribitsExampleProject/README.md) using just
`find_package(TribitsExProj REQUIRED COMPONENTS ...)` but uses the old
specification for using an installed TriBITS project (i.e. through CMake
variables `<Project>_LIBRARIES`, `<Project>_INCLUDE_DIRS`, and
`<Project>_TPL_INCLUDE_DIRS`).

After building and installing TribitsExampleProject under
`<upstreamInstallDir>`, then configure, build, and test
`TribitsOldSimpleExampleApp` with:

```
  cmake \
    -DCMAKE_PREFIX_PATH=<upstreamInstallDir> \
    <base-dir>/TribitsOldSimpleExampleApp

  make

  ctest
```

This project can be instructed to use the non-namespaced old deprecated
targets by adding the CMake cache var:

```
  -D TribitsOldSimpleExApp_USE_DEPRECATED_TARGETS=ON
```

This will generated deprecated warnings with newer versions of TriBITS.  (This
is to demonstrate and test these deprecated non-namespaced targets and to
demonstrate how to remove the deprecated warnings in a downstream customer
CMake project.)

NOTE: The modern version of this project using the modern CMake targets
produced by `TribitsExampleProject` is
[`TribitsSimpleExampleApp`](../TribitsSimpleExampleApp/README.md).  Comparing
the `CMakeLists.txt` files in these two projects (e.g. [here](CMakeLists.txt)
and [here](../TribitsSimpleExampleApp/CMakeLists.txt)) demonstrates the
simplifications afforded by the namespaced modern CMake targets.
