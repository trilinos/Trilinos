
INCLUDE(PackageListHelpers)

#
# Define the Trilinos package names, directories, and classification.
#
# Package classifications are:
#
#   PS: Primary Stable Package
#
#     Primary Stable Packages have at least some Primary Stable Code which is
#     expected to be fully tested before every push to the global repo.  The
#     default enable for PS packages is empty "" which allows the PS package
#     to be enabled implicitly based on other criteria.  The option
#     Trilinos_ENABLE_ALL_PACKAGES=ON will cause all PS packages to be enabled
#     unless they are explicitly disabled.
#
#   SS: Secondary Stable Package
#
#     Secondary Stable Packages have no PS code or they would be classified as
#     PS packages.  A package must be classified as SS if it has a required
#     dependency on another SS package or SS TPL.  A package may also be
#     declared SS to avoid requiring it to be tested before every push to the
#     global repo.  For example, a package that does not provide any
#     significant functionally like Didasko is classified as a SS package even
#     through it could be classified as PS just based on its required package
#     and TPL dependencies.  SS packages will have their default enables set
#     to empty "".  This allows them to be enabled implicilty.  When
#     Trilinos_ENABLE_ALL_PACKAGES=ON but
#     Trilinos_ENABLE_SECONDARY_STABLE_CODE=OFF, the SS packages will not be
#     enabled.  However, when Trilinos_ENABLE_ALL_PACKAGES=ON and
#     Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON, then SS packages will be
#     enabled if they are not explicitly disabled.  Packages that are SS but
#     not PS must be disabled in pre-push testing.  However, SS packages are
#     tested by the post-push CI and nightly testing processes.
#
#   EX: Experimental Package
#
#     Experimental packages are those packages that contain no PS or SS
#     code. The default enable for EX packages is always OFF which requires
#     that they be explicitly enabled in order to be turned on. EX packages
#     must be disabled in pre-push testring and are not tested as part of the
#     post-push CI or nightly testing processes.  However, package developers
#     of EX pacakges are encouraged to set up their own nightly testing for
#     thier EX packages.
#
# NOTE: These packages must be listed in strictly assending order in terms of
# package dependencies.  If you get the order wrong, then an error message
# will be printed during configuration with CMake.
#

SET( Trilinos_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
  TrilinosFramework     ../cmake                       PS # Only tests, no libraries/capabilities!
  Teuchos               teuchos                        PS
  ThreadPool            ThreadPool                     PS # Depends on ptheads system library
  Sacado                sacado                         PS
  RTOp                  rtop                           PS
  Kokkos                kokkos                         PS
  Epetra                epetra                         PS
  Zoltan                zoltan                         PS
  Shards                shards                         PS
  GlobiPack             globipack                      PS
  Triutils              triutils                       PS
  Tpetra                tpetra                         PS
  EpetraExt             epetraext                      PS
  Thyra                 thyra                          PS
  OptiPack              optipack                       PS
  Isorropia             isorropia                      PS
  Pliris                pliris                         PS
  Claps                 claps                          SS
  AztecOO               aztecoo                        PS
  Galeri                galeri                         PS
  Amesos                amesos                         PS
  Pamgen                pamgen                         PS
  SEACAS                seacas                         SS # Depends on netcdf, optionally hdf5, xdmf, pamgen
  Trios                 Trios                          SS 
  Ifpack                ifpack                         PS
  Komplex               komplex                        PS
  ML                    ml                             PS
  Belos                 belos                          PS
  Ifpack2               ifpack2                        PS
  Stratimikos           stratimikos                    PS
  Meros                 meros                          EX # no longer released
  FEI                   fei                            PS
  Anasazi               anasazi                        PS
  Teko                  teko                           SS
  RBGen                 rbgen                          PS
  TriKota               TriKota                        SS
  Intrepid              intrepid                       PS
  STK                   stk                            SS # Depends on boost
  Phalanx               phalanx                        SS
  Phdmesh               phdmesh                        EX # to be replaced by STK
  NOX                   nox                            PS
  Moertel               moertel                        PS
  TrilinosCouplings     trilinoscouplings              SS
  Rythmos               rythmos                        PS
  MOOCHO                moocho                         PS
  Aristos               aristos                        EX
  Stokhos               stokhos                        SS
  Piro                  piro                           SS
  Sundance              Sundance                       SS # Could be PS based on deps (BUG: 4669)
  CTrilinos             CTrilinos                      SS # Switched to SS to speed up checkin testing
  ForTrilinos           ForTrilinos                    EX
  PyTrilinos            PyTrilinos                     SS
  WebTrilinos           WebTrilinos                    EX # Should be SS
  Didasko               didasko                        SS
  NewPackage            new_package                    EX # Should be SS
  Optika		optika			       SS
  Mesquite              mesquite                       PS
  FEApp                 ../demos/FEApp                 SS # Capability demonstration package
  )


#
# Disable certain packages on certain platforms.
#
# NOTE: This just makes the packages experimental 'EX' and therefore still
# allows the user to enable the package explicitly but the package will not
# get enabled implicitly.
#

PACKAGE_DISABLE_ON_PLATFORMS(MOOCHO Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Phalanx Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Phdmesh Windows)
PACKAGE_DISABLE_ON_PLATFORMS(PyTrilinos Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Sundance Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Tpetra Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Ifpack2 Windows)
PACKAGE_DISABLE_ON_PLATFORMS(TriKota Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Pamgen Windows)
PACKAGE_DISABLE_ON_PLATFORMS(STK Windows)
PACKAGE_DISABLE_ON_PLATFORMS(SEACAS Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Trios Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Anasazi Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Zoltan Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Isorropia Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Teko Windows)
PACKAGE_DISABLE_ON_PLATFORMS(Mesquite AIX)
