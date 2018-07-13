TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  Exodus      libraries/exodus        PT  REQUIRED
  Exodus_for  libraries/exodus_for    PT  OPTIONAL
  ExoIIv2for32 libraries/exoIIv2for32 PT  OPTIONAL
  Nemesis     libraries/nemesis       PT  OPTIONAL
  Ioss        libraries/ioss          PT  REQUIRED
  Chaco       libraries/chaco         PT  OPTIONAL
  Aprepro_lib libraries/aprepro_lib   PT  OPTIONAL
  Supes       libraries/supes         PT  OPTIONAL
  Suplib      libraries/suplib        PT  OPTIONAL
  SuplibC     libraries/suplib_c      PT  OPTIONAL
  SuplibCpp   libraries/suplib_cpp    PT  OPTIONAL
  SVDI        libraries/svdi          ST  OPTIONAL
  PLT         libraries/plt           ST  OPTIONAL
  Algebra     applications/algebra    PT  OPTIONAL
  Aprepro     applications/aprepro    PT  OPTIONAL
  Blot        applications/blot       ST  OPTIONAL
  Conjoin     applications/conjoin    PT  OPTIONAL
  Ejoin       applications/ejoin      PT  OPTIONAL
  Epu         applications/epu        PT  OPTIONAL
  Exo2mat     applications/exo2mat    ST  OPTIONAL
  Exodiff     applications/exodiff    PT  OPTIONAL
  Exomatlab   applications/exomatlab  ST  OPTIONAL
  Exotxt      applications/exotxt     PT  OPTIONAL
  Exo_format  applications/exo_format PT  OPTIONAL
  Ex1ex2v2    applications/ex1ex2v2   PT  OPTIONAL
  Fastq       applications/fastq      ST  OPTIONAL
  Gjoin       applications/gjoin      PT  OPTIONAL
  Gen3D       applications/gen3d      PT  OPTIONAL
  Genshell    applications/genshell   PT  OPTIONAL
  Grepos      applications/grepos     PT  OPTIONAL
  Grope       applications/grope      PT  OPTIONAL
  Mapvarlib   libraries/mapvarlib     PT  OPTIONAL
  Mapvar      applications/mapvar     PT  OPTIONAL
  Mapvar-kd   applications/mapvar-kd  PT  OPTIONAL
  Mat2exo     applications/mat2exo    ST  OPTIONAL
  Nemslice    applications/nem_slice  PT  OPTIONAL
  Nemspread   applications/nem_spread PT  OPTIONAL
  Numbers     applications/numbers    PT  OPTIONAL
#  Slice	      applications/slice      ST  OPTIONAL
  Txtexo      applications/txtexo     PT  OPTIONAL
  Ex2ex1v2    applications/ex2ex1v2   PT  OPTIONAL
)

SET(LIB_OPTIONAL_DEP_TPLS MPI)

