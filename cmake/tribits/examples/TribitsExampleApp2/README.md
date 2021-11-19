# TribitsExampleApp2

The example project `TribitsExampleApp2` is a raw CMake project that pulls in
libraries from packages from `TribitsExampleProject2`.  To build against all of
the installed packages from an upstream `TribitsExampleProject2`, configure,
build, and run the tests with:

```
  cmake \
    -DCMAKE_PREFIX_PATH=<upstreamInstallDir> \
    <base-dir>/TribitsExampleApp2

  make

  ctest
```

That will call `find_package(TribitsExProj2)` and will set the compilers
pulled in from the found `TribitsExProj2Config.cmake` file.


## Pulling in only a subset of packages

To configure and build `TribitsExampleApp2` against only a subset of the
installed `TribitsExampleProject2` packages, configure `TribitsExampleApp2`
using, for example:

```
  cmake \
    -DCMAKE_PREFIX_PATH=<upstreamInstallDir> \
    -DTribitsExApp_USE_COMPONENTS=Package1,Package2 \
    <base-dir>/TribitsExampleApp2
```

Internally, that that will call:

```
  find_package(TribitsExApp REQUIRED COMPONENTS ${TribitsExApp_USE_COMPONENTS})
```

(where `,` in `TribitsExApp_USE_COMPONENTS` is replaced with `;` internally to
create a proper CMake list internally before calling `find_package()`).


## Pulling in only a subset of packages by finding each package individually

`TribitsExampleApp2` is also set up to demonstrate finding the individual
packages using separate calls to `find_package(<Package>)` by configuring
with, for example:

```
  cmake \
    -DCMAKE_PREFIX_PATH=<upstreamInstallDir> \
    -DTribitsExApp_USE_COMPONENTS=Package1,Package2 \
    -DTribitsExApp_FIND_INDIVIDUAL_PACKAGES=ON \
    <base-dir>/TribitsExampleApp2
```  

That essentially results in calling:

```
  find_package(Package1 REQUIRED)
  find_package(Package2 REQUIRED)
  ...
  target_link_libraries(app
    PRIVATE Package1::all_libs
    PRIVATE Package2::all_libs
    )
```

(but does so with a loop and a list internally).


## Pulling in packages from the build tree

`TribitsExampleApp2` is additionally set up to demonstrate finding the
individual packages from the build directory (where `<upstreamBuildDir>` is
the build directory for `TribitsExampleProject2`) instead of the install tree
using separate calls to `find_package(<Package>)` by configuring with, for
example:

```
  cmake \
    -DCMAKE_PREFIX_PATH=<upstreamBuildDir>/cmake_packages \
    -DTribitsExApp_USE_COMPONENTS=Package1,Package2 \
    -DTribitsExApp_FIND_INDIVIDUAL_PACKAGES=ON \
    <base-dir>/TribitsExampleApp2
```  

This is identical to the case to finding the individual packages under the
install tree except here the packages and libraries are found under the build
tree and the include directories point into the source tree.  (There is no
need for an install directory in this case.)  This is to simulate more
advanced use cases such as where `TribitsExampleApp2` and
`TribitsExampleProject2` may be part of a larger super-build that works out of
the build tree without needing to install.

NOTE: There is no `TribitsExProj2Config.cmake` file generated in the build tree
to find so it is not possible to call `find_package(TribitsExProj2Config)`
pointing into the build tree.


## Tests demonstrating usage of TribitsExampleApp2

The TriBITS project contains automated tests that run the various use cases
described above and check their result.  These tests are shown defined
[here](https://github.com/TriBITSPub/TriBITS/blob/master/test/core/ExamplesUnitTests/TribitsExampleApp2_Tests.cmake)
and are run locally with CTest and are submitted to
[CDash](https://github.com/TriBITSPub/TriBITS/wiki/TriBITS-CDash-Dashboard) as
part of regular testing.
