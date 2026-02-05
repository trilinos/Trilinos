# Speeding up Trilinos builds

Trilinos contains a lot of functionality. Build times can be
significant. There are however different approaches to achieve
speedups.


## CMake configuration options

The simplest way to reduce overall build time is to disable unneeded
targets. The output of the CMake invocation shows which packages were
enabled and why they were enabled.

- Unneeded packages should be disabled.

  Packages have required and optional dependencies on other packages.
  By enabling a package, e.g.
  
  ```
  -D Trilinos_ENABLE_MueLu:BOOL=ON
  ```
  
  all its required dependencies also get enabled. Optional
  dependencies get automatically enabled when
  
  ```
  -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES:BOOL=ON
  ```
  
  which is the default behavior.

- Depending on the package, building tests and examples can take
  significantly longer than the libraries themselves. Instead of
  enabling all tests or examples with
  
  ```
  -D Trilinos_ENABLE_TESTS:BOOL=ON \
  -D Trilinos_ENABLE_EXAMPLES:BOOL=ON \
  ```
  
  one should prefer to only enable tests and examples for certain
  packages. For example, to enable tests and examples for the MueLu
  package one can set
  
  ```
  -D MueLu_ENABLE_TESTS:BOOL=ON \
  -D MueLu_ENABLE_EXAMPLES:BOOL=ON \
  ```

- Unneeded scalar types should be disabled. By default only `double`
  scalars are enabled.
  
- Explicit template instantiation (ETI) is crucial for fast builds. It
  is enabled by default:

  ```
  -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON
  ```


## Ninja build system

The [Ninja](https://ninja-build.org/) build system generator can be
used as a replacement for the default UNIX Makefiles. It offers much
faster dependency checking.

Install ninja:
- Debian/Ubuntu: `sudo apt-get install ninja-build`
- RedHat: `sudo yum install ninja-build`
- MacOS: `brew install ninja`

and configure Trilinos with the CMake flag `-G Ninja`.

Trilinos will generate Makefiles in all sub-directories that call
Ninja internally. This implies that `make` can still be used to build
Trilinos. One difference is that instead of calling

```sh
make -j 10
```

we have to switch to

```sh
NP=10 make
```

to specify the level of build parallelism.


## Ccache

Repeated build (e.g. in a CI pipeline or during development) can be
sped up significantly by caching the results of compilation using
[ccache](https://ccache.dev/).

Install Ccache:
- Debian/Ubuntu: `sudo apt-get install ccache`
- RedHat: `sudo yum install ccache`
- MacOS: `brew install ccache`

and configure Trilinos with the CMake flags

```
-D CMAKE_CXX_COMPILER_LAUNCHER:STRING="ccache" \
-D CMAKE_C_COMPILER_LAUNCHER:STRING="ccache"
```

Statistics about cache hits and misses can be displayed by calling

```sh
ccache --show-stats
```

This allows to verify that the build cache is being used.

When ccache is used in CI pipelines one has to make sure that the
cache directory is preserved between subsequent runs.


## Linker

If linker performance is an issue, it might make sense to explore
alternative linkers such as LLVM's [LLD](https://lld.llvm.org/) or
[mold](https://github.com/rui314/mold).
