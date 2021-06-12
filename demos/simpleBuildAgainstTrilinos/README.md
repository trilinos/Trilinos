# Simple example of building against an install Trilinos with CMake

```
******************************************************************
*** WARNING: THIS CODE IS NOT AUTOMATICALLY TESTED IN TRILINOS ***
***          BECAUSE IT REQUIRES AN INSTALLED TRILINOS         ***
******************************************************************
```

This is a small demonstration of how to build an application code using CMake
against an installation of Trilinos.  This example also demonstrates how to
use CTest.

To run this example, do the following steps.

## 1. Configure Trilinos to be installed:

```
  $ cd <any-base-dir>/
  $ mkdir <tril-build-dir>
  $ cd <tril-build-dir>/
  $ cmake \
    -D CMAKE_INSTALL_PREFIX=<tril-install-prefix> \
    -D Trilinos_ENABLE_Teuchos=ON \
    -D Trilinos_ENABLE_Tpetra=ON \
    [other options] \
    <trilinos-src>
```

NOTE: Must enable at least the Trilinos packages Teuchos and Tpetra in order
to build this example.  But in a real case, enable and install any packages
needed by the application.

## 2. Build and install Trilinos:

```
  $ make -j4 install
```

This will put a file called `TrilinosConfig.cmake` under a subdirectory of
`<tril-install-prefix>/`.

## 3. Configure the example project:

```
  $ cd <any-base-dir>/
  $ mkdir <app-build-dir>
  $ cd <app-build-dir>/
  $ cmake \
    -D CMAKE_PREFIX_PATH=<tril-install-prefix> \
    <trilinos-src>/demos/simpleBuildAgainstTrilinos
```

## 4. Build the application

```
  $ make -j4
```

## 5. Run the application test(s):

```
  $ ctest
```

NOTE: That test will not pass unless Trilinos is propertly built with MPI and
MPI is working.

## More info and questions

Look into comments in the `CMakeLists.txt` file for more info.

Send questions to trilinos-framework@software.sandia.gov.
