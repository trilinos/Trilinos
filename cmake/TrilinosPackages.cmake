
#
# Define the Trilinos package names, directories, and classification.
#
# Package classifications are:
#
#   PS: Primary Stable Package
#
#     Primary Stable Packages have at least some Primary Stable Code which is
#     expected to be fully tested before every checkin.  The default enable
#     for PS packages is empty "" which allows the PS package to be enabled
#     implicitly based on other criteria.  The option
#     Trilinos_ENABLE_ALL_PACKAGES=ON will cause all PS packages to be enabled
#     unless they are explicitly disabled.
#
#   SS: Secondary Stable Package
#
#     Secondary Stable Packages have no PS code or they would be classified as
#     PS packages.  A package must be classified as SS if it has a required
#     dependency on another SS package or SS TPL.  A package may also be
#     declared SS to avoid requiring it to be tested before every checkin.
#     For example, a package that does not provide any significant
#     functionally like Didasko is classified as a SS package even through it
#     could be classified as PS just based on its required package and TPL
#     dependencies.  SS packages will have their default enables set to empty
#     "".  This allows them to be enabled implicilty.  When
#     Trilinos_ENABLE_ALL_PACKAGES=ON but
#     Trilinos_ENABLE_SECONDARY_STABLE_CODE=OFF, the SS packages will not be
#     enabled.  However, when Trilinos_ENABLE_ALL_PACKAGES=ON and
#     Trilinos_ENABLE_SECONDARY_STABLE_CODE=ON, then SS packages will be
#     enabled if they are not explicitly disabled.  Packages that are SS but
#     not PS must be disabled in precheckin testing.  However, SS packages are
#     tested by the nightly testing process.
#
#   EX: Experimental Package
#
#     Experimental packages are those packages that contain no PS or SS
#     code. The default enable for EX packages is always OFF which requires
#     that they be explicitly enabled in order to be turned on. EX packages
#     must be disabled in precheckin testring and are not tested as part of
#     the nightly testing process.  However, package developers of EX pacakges
#     are encouraged to set up their own nightly testing for ther EX packages.

#
# NOTE: These packages must be listed in strictly assending order in
# terms of package dependencies.  If you get the order wrong, then an
# error message will be printed.
#

SET( Trilinos_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
  Teuchos               teuchos                        PS
  ThreadPool            ThreadPool                     SS # Depends on Ptheads TPL
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
  Ifpack                ifpack                         PS
  Komplex               komplex                        PS
  ML                    ml                             PS
  Belos                 belos                          PS
  Stratimikos           stratimikos                    PS
  Meros                 meros                          PS
  FEI                   fei                            PS
  RBGen                 rbgen                          PS
  Anasazi               anasazi                        PS
  TriKota               TriKota                        SS
  Stokhos               stokhos                        SS
  Sacado                sacado                         PS
  Intrepid              intrepid                       PS
  Phalanx               phalanx                        SS
  Phdmesh               phdmesh                        PS
  NOX                   nox                            PS
  Moertel               moertel                        PS
  TrilinosCouplings     trilinoscouplings              SS
  Rythmos               rythmos                        PS
  MOOCHO                moocho                         PS
  Aristos               aristos                        EX
  Sundance              Sundance                       PS
  CTrilinos             CTrilinos                      EX
  ForTrilinos           ForTrilinos                    EX
  PyTrilinos            PyTrilinos                     SS
  WebTrilinos           WebTrilinos                    EX # Should be SS
  Didasko               didasko                        SS
  NewPackage            new_package                    EX # Should be SS
  Optika		optika			       SS
  Mesquite              mesquite                       EX
  )
