
# This list is just used for unit testing the dependency handling
# CMake code.  The reason that we have a separate list is so that we
# can keep very stable unit tests.

SET( Trilinos_PACKAGES_AND_DIRS_AND_ENABLES
  Teuchos          teuchos                ""
  Epetra           epetra                 "" 
  Zoltan           zoltan                 ""
  RTOp             rtop                   ""
  Triutils         triutils               ""
  EpetraExt        epetraext              ""
  Isorropia        isorropia              ""
  Thyra            thyra                  ""
  AztecOO          aztecoo                ""
  Galeri           galeri                 ""
  Amesos           amesos                 ""
  Ifpack           ifpack                 ""
  ML               ml                     ""
  Belos            belos                  ""
  RBGen            rbgen                  ""
  Stratimikos      stratimikos            ""
  Stokhos          stokhos                OFF
  )
