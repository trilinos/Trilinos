TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    Core                core                PT  REQUIRED
    EpetraAdapters      adapters/epetra     ST  OPTIONAL
    EpetraExtAdapters   adapters/epetraext  ST  OPTIONAL
    TpetraAdapters      adapters/tpetra     PT  OPTIONAL
  )
