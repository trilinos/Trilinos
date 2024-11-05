tribits_repository_define_tpls(
  MPI  "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/FindTPLMPI.cmake"   PT
  HeaderOnlyTpl  "cmake/tpls/"   PT
  SimpleTpl      "cmake/tpls/"   ST
  )
