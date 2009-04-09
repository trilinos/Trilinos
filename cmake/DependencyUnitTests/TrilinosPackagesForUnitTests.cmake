
# This list is just used for unit testing the dependency handling
# CMake code.  The reason that we have a separate list is so that we
# can keep very stable unit tests.


SET( Trilinos_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
  Teuchos          teuchos                PS
  Epetra           epetra                 PS 
  Zoltan           zoltan                 PS
  RTOp             rtop                   PS
  Triutils         triutils               PS
  Tpetra           tpetra                 PS
  EpetraExt        epetraext              PS
  Isorropia        isorropia              PS
  Thyra            thyra                  PS
  AztecOO          aztecoo                PS
  Galeri           galeri                 PS
  Amesos           amesos                 PS
  Ifpack           ifpack                 PS
  ML               ml                     PS
  Belos            belos                  PS
  RBGen            rbgen                  PS
  Stratimikos      stratimikos            PS
  Stokhos          stokhos                EX
  Sacado           sacado                 PS
  Phalanx          phalanx                SS
  )

# NOTE: Sacado is really PS but for testing purpose it is make SS
