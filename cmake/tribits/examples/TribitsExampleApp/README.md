# TribitsExampleApp

The example project `TribitsExampleApp` is a raw CMake project that pulls in
libraries from packages from `TribitsExampleProject`. If `WithSubpackageC` 
subpackage or the `MixedLang` packages are built, then a Fortran compiler 
is required to build `TribitsExampleApp`. To build against all of the 
installed packages from an upstream `TribitsExampleProject`, configure, 
build, and run the tests with:

```
  cmake \
    -DCMAKE_PREFIX_PATH=<upstreamInstallDir> \
    <base-dir>/TribitsExampleApp

  make

  ctest
```

That will call `find_package(TribitsExProj)` and will set the compilers
pulled in from the found `TribitsExProjConfig.cmake` file.


## Pulling in only a subset of packages

To configure and build `TribitsExampleApp` against only a subset of the
installed `TribitsExampleProject` packages, configure `TribitsExampleApp`
using, for example:

```
  cmake \
    -DCMAKE_PREFIX_PATH=<upstreamInstallDir> \
    -DTribitsExApp_USE_COMPONENTS=SimpleCxx,WithSubpackages \
    <base-dir>/TribitsExampleApp
```

Internally, that that will call:

```
  find_package(TribitsExApp REQUIRED COMPONENTS ${TribitsExApp_USE_COMPONENTS})
```

(where `,` in `TribitsExApp_USE_COMPONENTS` is replaced with `;` internally to
create a proper CMake list internally before calling `find_package()`).


## Pulling in only a subset of packages by finding each package individually

`TribitsExampleApp` is also set up to demonstrate finding the individual
packages using separate calls to `find_package(<Package>)` by configuring
with, for example:

```
  cmake \
    -DCMAKE_PREFIX_PATH=<upstreamInstallDir> \
    -DTribitsExApp_USE_COMPONENTS=SimpleCxx,WithSubpackages \
    -DTribitsExApp_FIND_INDIVIDUAL_PACKAGES=ON \
    <base-dir>/TribitsExampleApp
```  

That essentially results in calling:

```
  find_package(SimpleCxx REQUIRED)
  find_package(WithSubpackages REQUIRED)
  ...
  target_link_libraries(app
    PRIVATE SimpleCxx::all_libs
    PRIVATE WithSubpackages::all_libs
    )
```

(but does so with a loop and a list internally).


## Pulling in packages from the build tree

`TribitsExampleApp` is additionally set up to demonstrate finding the
individual packages from the build directory (where `<upstreamBuildDir>` is
the build directory for `TribitsExampleProject`) instead of the install tree
using separate calls to `find_package(<Package>)` by configuring with, for
example:

```
  cmake \
    -DCMAKE_PREFIX_PATH=<upstreamBuildDir>/cmake_packages \
    -DTribitsExApp_USE_COMPONENTS=SimpleCxx,WithSubpackages \
    -DTribitsExApp_FIND_INDIVIDUAL_PACKAGES=ON \
    <base-dir>/TribitsExampleApp
```  

This is identical to the case to finding the individual packages under the
install tree except here the packages and libraries are found under the build
tree and the include directories point into the source tree.  (There is no
need for an install directory in this case.)  This is to simulate more
advanced use cases such as where `TribitsExampleApp` and
`TribitsExampleProject` may be part of a larger super-build that works out of
the build tree without needing to install.

NOTE: There is no `TribitsExProjConfig.cmake` file generated in the build tree
to find so it is not possible to call `find_package(TribitsExProjConfig)`
pointing into the build tree.


## Tests demonstrating usage of TribitsExampleApp

The TriBITS project contains automated tests that run the various use cases
described above and check their result.  These tests are shown defined
[here](https://github.com/TriBITSPub/TriBITS/blob/master/test/core/ExamplesUnitTests/TribitsExampleApp_Tests.cmake)
and are run locally with CTest and are submitted to
[CDash](https://github.com/TriBITSPub/TriBITS/wiki/TriBITS-CDash-Dashboard) as
part of regular testing.
