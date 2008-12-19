#
# Define the Trilinos package names, directories, and default
# enables/disables
#
# NOTE: These packages must be listed in strictly assending order in
# terms of package dependencies.  If you get the order wrong, then an
# error message will be printed.
#
# NOTE: Packages that are not CMakeified yet or are experimental and
# should not be build by default should be marked with OFF instead of
# empty "" in the third 'ENABLES' column.  All packages that are
# marked with the empty enable "" are left to be enabled or disabled
# based on whatever logic the user or the scripts want to use.  Note
# that the user can explicitly override any enable value by setting
# the value in the cache.


SET(Trilinos_PACKAGES_AND_DIRS_AND_ENABLES
  Teuchos               teuchos                        ""
  RTOp                  rtop                           ""
  Kokkos                kokkos                         ""
  Epetra                epetra                         ""
  Stokhos               stokhos                        OFF
  Sacado                sacado                         ""
  Zoltan                zoltan                         ""
  Shards                shards                         ""
  Intrepid              intrepid                       OFF
  Triutils              triutils                       ""
  Tpetra                tpetra                         ""
  EpetraExt             epetraext                      ""
  Thyra                 thyra                          ""
  Isorropia             isorropia                      ""
  Pliris                pliris                         ""
  Claps                 claps                          ""
  AztecOO               aztecoo                        ""
  Galeri                galeri                         ""
  Amesos                amesos                         ""
  Ifpack                ifpack                         ""
  Komplex               komplex                        ""
  ML                    ml                             ""
  Belos                 belos                          ""
  Stratimikos           stratimikos                    ""
  Meros                 meros                          OFF
  FEI                   fei                            OFF
  RBGen                 rbgen                          ""
  Anasazi               anasazi                        ""
  ThreadPool            ThreadPool                     OFF
  Phalanx               phalanx                        OFF
  Pamgen                pamgen                         OFF
  Phdmesh               phdmesh                        OFF
  NOX                   nox                            OFF
  Moertel               moertel                        OFF
  TrilinosCouplings     trilinoscouplings              OFF
  Rythmos               rythmos                        OFF
  MOOCHO                moocho                         OFF
  Aristos               aristos                        OFF
  Sundance              Sundance                       OFF
  TriKota               TriKota                        OFF
  CTrilinos             CTrilinos                      OFF
  ForTrilinos           ForTrilinos                    OFF
  PyTrilinos            PyTrilinos                     OFF
  WebTrilinos           WebTrilinos                    OFF
  Didasko               didasko                        OFF
  NewPackage            new_package                    OFF
  )
