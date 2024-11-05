tribits_package_define_dependencies(
  LIB_REQUIRED_PACKAGES Teuchos Sacado Phalanx Intrepid Thyra
    Tpetra Epetra EpetraExt
  LIB_OPTIONAL_PACKAGES Stokhos
  TEST_OPTIONAL_PACKAGES Stratimikos
  LIB_REQUIRED_TPLS MPI Boost
  )
