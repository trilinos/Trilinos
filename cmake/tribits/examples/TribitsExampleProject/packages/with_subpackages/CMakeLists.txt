tribits_package_decl(WithSubpackages)

tribits_add_debug_option()

tribits_process_subpackages()

tribits_exclude_files(
  b/ExcludeFromRelease.txt
  b/src/AlsoExcludeFromTarball.txt
  )

tribits_package_def()

tribits_package_postprocess()
