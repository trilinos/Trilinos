TRIBITS_REPOSITORY_DEFINE_PACKAGES(
  TrilinosFramework   cmake                           PT
  Teuchos             packages/teuchos                PT
  RTOp                packages/rtop                   PT
  Epetra              packages/epetra                 PT
  Zoltan              packages/zoltan                 PT
  Shards              packages/shards                 PT
  Triutils            packages/triutils               PT
  Tpetra              packages/tpetra                 PT
  EpetraExt           packages/epetraext              PT
  Stokhos             packages/stokhos                EX
  Sacado              packages/sacado                 ST
  Thyra               packages/thyra                  PT
  Isorropia           packages/isorropia              PT
  AztecOO             packages/aztecoo                PT
  Galeri              packages/galeri                 PT
  Amesos              packages/amesos                 PT
  Intrepid            packages/intrepid               PT
  Ifpack              packages/ifpack                 PT
  ML                  packages/ml                     PT
  Belos               packages/belos                  ST
  Stratimikos         packages/stratimikos            PT
  RBGen               packages/rbgen                  PT
  Phalanx             packages/phalanx                ST
  Panzer              packages/panzer                 ST
  )

# NOTE: Sacado was really PT but for testing purpose it is made ST
# NOTE: Belos was really PT but for testing purpose it is made ST

TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(ML BadSystem1)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Ifpack BadSystem1 BadSystem2)
