INCLUDE(PackageListHelpers)


# This list is just used for unit testing the dependency handling
# CMake code.  The reason that we have a separate list is so that we
# can keep very stable unit tests.


SET( Trilinos_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
  TrilinosFramework ../cmake              PS
  Teuchos          teuchos                PS
  RTOp             rtop                   PS
  Epetra           epetra                 PS 
  Zoltan           zoltan                 PS
  Shards           shards                 PS
  Triutils         triutils               PS
  Tpetra           tpetra                 PS
  EpetraExt        epetraext              PS
  Stokhos          stokhos                EX
  Sacado           sacado                 SS
  Thyra            thyra                  PS
  Isorropia        isorropia              PS
  AztecOO          aztecoo                PS
  Galeri           galeri                 PS
  Amesos           amesos                 PS
  Intrepid         intrepid               PS
  Ifpack           ifpack                 PS
  ML               ml                     PS
  Belos            belos                  SS
  Stratimikos      stratimikos            PS
  RBGen            rbgen                  PS
  Phalanx          phalanx                SS
  )

# NOTE: Sacado is really PS but for testing purpose it is made SS
# NOTE: Belos is really PS but for testing purpose it is made SS

PACKAGE_DISABLE_ON_PLATFORMS(ML BadSystem1)
PACKAGE_DISABLE_ON_PLATFORMS(Ifpack BadSystem1 BadSystem2)
