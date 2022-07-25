# TribitsSimpleExampleApp

This simple example project `TribitsSimpleExampleApp` is a raw CMake project
that pulls in libraries from a few packages from
[`TribitsExampleProject`](../TribitsExampleProject/README.md) using just
`find_package(TribitsExProj REQUIRED COMPONENTS ...)` and also links a smaller
program using a subset of the libraries from one of the packages from
`TribitsExampleProject`.

After building and installing `TribitsExampleProject` under
`<upstreamInstallDir>`, then configure, build, and test
`TribitsSimpleExampleApp` with:

```
  cmake \
    -DCMAKE_PREFIX_PATH=<upstreamInstallDir> \
    <base-dir>/TribitsSimpleExampleApp

  make

  ctest
```

NOTE: The version of this project that demonstrates how to use the old interface, 
with variables that work with much older versions of TriBITS, is given in
[`TribitsOldSimpleExampleApp`](../TribitsOldSimpleExampleApp/README.md).
