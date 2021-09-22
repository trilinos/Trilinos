tribits_repository_define_packages(
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
  Phalanx             packages/phalanx                ST
  Panzer              packages/panzer                 ST
  AlwaysMissing       AlwaysMissing                   PT
  )

# NOTE: Sacado was really PT but for testing purpose it is made ST
# NOTE: Belos was really PT but for testing purpose it is made ST

tribits_allow_missing_external_packages(AlwaysMissing)
tribits_disable_package_on_platforms(ML BadSystem1)
tribits_disable_package_on_platforms(Ifpack BadSystem1 BadSystem2)
