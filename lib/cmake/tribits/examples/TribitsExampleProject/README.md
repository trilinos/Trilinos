# TribitsExampleProject Documentation

The project `TribitsExampleProject` defines a TriBITS CMake project designed to
provide a simple example to demonstrate how to use the TriBITS system to
create a CMake build, test, and deployment system using a package-based
architecture. To build all of the packages from `TribitsExampleProject`,
a Fortran compiler is needed.

To build and test the project, one must first create a build directory and
configure by pointing to the TriBITS source dir `<tribits-dir>`
(i.e. `TriBITS/tribits`) and the `TribitsExampleProject` project source dir
with:

```
  mkdir <build-dir>
  
  cd <build-dir>
  
  cmake \
    -DTribitsExProj_TRIBITS_DIR=<tribits-dir> \
    -DCMAKE_INSTALL_PREFIX=<TribitsExampleProject-install-dir>
    -DTribitsExProj_ENABLE_TESTS=ON \
    -DTribitsExProj_ENABLE_ALL_PACKAGES=ON
    -DTribitsExProj_ENABLE_SECONDARY_TESTED_CODE=ON \
    -DCMAKE_CXX_COMPILER=g++ \
    <path-to-TribitsExampleProject>
```

then build and test with:

```
  make -j4 && ctest -j4
```

and then install:

```
  make install
```

`TribitsExampleProject` will be installed to 
`<TribitsExampleProject-install-dir>` and this path should be then used 
to point `-DCMAKE_PREFIX_PATH=<TribitsExampleProject-install-dir>` when
building e.g. `TribitsExampleApp`

The layout of a TriBITS project is described in:

* https://tribits.org/doc/TribitsUsersGuide.html#tribits-project-structure

Otherwise, this example TriBITS project is simple enough that it should be
enough to get started as a template.
