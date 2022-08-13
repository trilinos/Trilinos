tribits_package_define_dependencies(
  LIB_OPTIONAL_PACKAGES  Tpetra)
# NOTE: By declaring a dependence on Tpetra, if that package or any package
# upstream from it changes, then it will trigger the enable of this package
# and the running of installation tests and the demo example.  But since
# TrilinosInstallTests does not have any libs or execs of its own (all of its
# work is done in tests), there is no link-time dependency here in the build
# dir.
