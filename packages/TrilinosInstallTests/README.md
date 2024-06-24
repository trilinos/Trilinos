# TrilinosInstallTests package

The purpose of the TrilinosInstallTests package is to test the installation of
Trilinos and related functionalities such as reduced tarball creation and
generated and installed `<Package>Config.cmake` files.  It also builds and
tests the demo project `simpleBuildAgainstTrilinos` which depends on the
package Tpetra.

## The role of the TrilinosInstallTests package in CI Testing

The role of the package `TrilinosInstallTests` in Trilinos CI testing is a bit
tricky and the details are described below.

Firs, the script `get-changed-trilinos-packages.sh` is set up to enable the
`TrilinosInstallTests` package if **any** Trilinos package is changed.  The
tests `TrilinosInstallTests_doinstall` and
`TrilinosInstallTests_find_package_Trilinos` are always run no matter what
Trilinos packages are enabled (and will pass even if no other Trilinos
packages are enabled).  This ensures that any changed package will have some
aspect of its configuration tested; that is, its install and its generated and
installed `<Package>Config.cmake` file are tested.

If the Tpetra package is enabled in the CI build (e.g. as a consequence of
setting `Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES=ON`), then the
`TrilinosInstallTests` package and the tests based on
`simpleBuildAgainstTrilinos` will be enabled.  This ensures that any changes
to Tpetra or its upstream dependencies will trigger the testing with
`simpleBuildAgainstTrilinos` (and thereby provide some installation testing
for Tpetra and its upstream dependencies such as Kokkos.)

However, if the directory `demos/simpleBuildAgainstTrilinos` is changed and
Tpetra is not enabled, then this will **not** trigger the enable of the tests
based on `simpleBuildAgainstTrilinos`.  This is a gap in the CI testing
process and therefore testing for changes to `simpleBuildAgainstTrilinos` must
be done locally.

Also, based on the above-described logic, not that if a package downstream
from Tpetra is changed, then this will trigger the enable of the tests based
on `simpleBuildAgainstTrilinos`, even though changing such a package can not
break `simpleBuildAgainstTrilinos` or the tests based on it.  But the CMake
project `simpleBuildAgainstTrilinos` configures, builds, and runs its tests
very quickly, so this should not negatively impact the CI testing process.

See more detailed comments in the file `./CMakeLists.txt` where the tests are
defined.

## Copyright and License
See TrilinosInstallTests/COPYRIGHT, TrilinosInstallTests/LICENSE, https://trilinos.github.io/license.html and individual file headers for additional information.

## Questions? 
Contact lead developers:

* Framework team   (GitHub handle: @trilinos/framework)
* Sam Browne       (GitHub handle: [sebrowne](https://github.com/sebrowne) or sebrown@sandia.gov)
