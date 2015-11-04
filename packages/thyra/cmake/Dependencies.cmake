TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    Core                core                PT  REQUIRED
    EpetraAdapters      adapters/epetra     PT  OPTIONAL
    EpetraExtAdapters   adapters/epetraext  PT  OPTIONAL
    TpetraAdapters      adapters/tpetra     PT  OPTIONAL
  )
