
macro(seacas_disable_subpackages_depending_on_fortran)
  if (NOT ${PROJECT_NAME}_ENABLE_Fortran)
    seacas_disable_subpackage_since_no_fortran(Mapvarlib)
    seacas_disable_subpackage_since_no_fortran(Exodus_for)
    seacas_disable_subpackage_since_no_fortran(ExoIIv2for32)
    seacas_disable_subpackage_since_no_fortran(Supes)
    seacas_disable_subpackage_since_no_fortran(Suplib)
    seacas_disable_subpackage_since_no_fortran(PLT)
    seacas_disable_subpackage_since_no_fortran(Blot)
    seacas_disable_subpackage_since_no_fortran(Fastq)
    seacas_disable_subpackage_since_no_fortran(SVDI)
    seacas_disable_subpackage_since_no_fortran(Algebra)
    seacas_disable_subpackage_since_no_fortran(Exotxt)
    seacas_disable_subpackage_since_no_fortran(Gjoin)
    seacas_disable_subpackage_since_no_fortran(Gen3D)
    seacas_disable_subpackage_since_no_fortran(Genshell)
    seacas_disable_subpackage_since_no_fortran(Grepos)
    seacas_disable_subpackage_since_no_fortran(Explore)
    seacas_disable_subpackage_since_no_fortran(Mapvar)
    seacas_disable_subpackage_since_no_fortran(Mapvar-kd)
    seacas_disable_subpackage_since_no_fortran(Numbers)
    seacas_disable_subpackage_since_no_fortran(Txtexo)
    seacas_disable_subpackage_since_no_fortran(Ex2ex1v2)
    seacas_disable_subpackage_since_no_fortran(Ex1ex2v2)
  endif()
endmacro()


macro(seacas_disable_subpackage_since_no_fortran subpackage)
  if (${PROJECT_NAME}_ENABLE_SEACAS${subpackage})
    message("-- "
      "WARNING: Setting ${PROJECT_NAME}_ENABLE_SEACAS${subpackage}=OFF"
      " even though it was set to ${${PROJECT_NAME}_ENABLE_SEACAS${subpackage}}"
      " because ${PROJECT_NAME}_ENABLE_Fortran=${${PROJECT_NAME}_ENABLE_Fortran}!"
      )
  elseif("${${PROJECT_NAME}_ENABLE_SEACAS${subpackage}}" STREQUAL "")
     message("-- "
      "NOTE: Setting ${PROJECT_NAME}_ENABLE_SEACAS${subpackage}=OFF"
      " because ${PROJECT_NAME}_ENABLE_Fortran=${${PROJECT_NAME}_ENABLE_Fortran}!"
       )
  endif()
  set(${PROJECT_NAME}_ENABLE_SEACAS${subpackage} OFF)
endmacro()
