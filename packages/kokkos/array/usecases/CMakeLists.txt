INCLUDE(TribitsAddOptionAndDefine)
TRIBITS_ADD_OPTION_AND_DEFINE(
  Kokkos_ENABLE_Implicit_Thermal_GMRES
  Kokkos_ENABLE_IMPLICIT_THERMAL_GMRES
  "Enable the implicit thermal gmres example"
  OFF
)


IF(Kokkos_ENABLE_Implicit_Thermal_GMRES)
  ADD_SUBDIRECTORY(implicit_thermal_gmres)
ENDIF(Kokkos_ENABLE_Implicit_Thermal_GMRES)

