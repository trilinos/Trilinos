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
  GlobiPack             globipack                      ""
  Triutils              triutils                       ""
  Tpetra                tpetra                         ""
  EpetraExt             epetraext                      ""
  Thyra                 thyra                          ""
  OptiPack              optipack                       ""
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
  Meros                 meros                          ""
  FEI                   fei                            ""
  RBGen                 rbgen                          ""
  Anasazi               anasazi                        ""
  ThreadPool            ThreadPool                     OFF
  Phalanx               phalanx                        ""
  Pamgen                pamgen                         ""
  Phdmesh               phdmesh                        OFF
  NOX                   nox                            ""
  Moertel               moertel                        ""
  TrilinosCouplings     trilinoscouplings              ""
  Rythmos               rythmos                        ""
  MOOCHO                moocho                         OFF
  Aristos               aristos                        OFF
  Sundance              Sundance                       OFF
  TriKota               TriKota                        OFF
  CTrilinos             CTrilinos                      OFF
  ForTrilinos           ForTrilinos                    OFF
  PyTrilinos            PyTrilinos                     ""
  WebTrilinos           WebTrilinos                    OFF
  Didasko               didasko                        OFF
  NewPackage            new_package                    OFF
  )
