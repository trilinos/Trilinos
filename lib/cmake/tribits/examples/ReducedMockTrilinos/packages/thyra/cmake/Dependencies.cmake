tribits_package_define_dependencies(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    CoreLibs      src                  PT  REQUIRED
    GoodStuff     good_stuff           ST  OPTIONAL
    CrazyStuff    crazy_stuff          EX  OPTIONAL
    Epetra        adapters/epetra      PT  OPTIONAL
    EpetraExt     adapters/epetraext   ST  OPTIONAL
  )
