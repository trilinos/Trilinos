# TribitsSimpleExampleApp

Simple example project `TribitsSimpleExampleApp` is a raw CMake project that
pulls in libraries from a few packages from `TribitsExampleProject` using just
`find_package(TribitsExProj REQUIRED COMPONENTS ...)`.

After building and installing TribitsExampleProject under
`<upstreamInstallDir>`, then configure, build, and test
`TribitsSimpleExampleApp` with:

```
  cmake \
    -DCMAKE_PREFIX_PATH=<upstreamInstallDir> \
    <base-dir>/TribitsSimpleExampleApp

  make

  ctest
```
