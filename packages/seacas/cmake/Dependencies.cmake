TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  Exodus      libraries/exodus        PT  REQUIRED
  Exodus_for  libraries/exodus_for    PT  REQUIRED
  ExoIIv2for32 libraries/exoIIv2for32 PT  REQUIRED
  Nemesis     libraries/nemesis       PT  REQUIRED
  Ioss        libraries/ioss          PT  REQUIRED
  Chaco       libraries/chaco         PT  REQUIRED
  Aprepro_lib libraries/aprepro_lib   PT  REQUIRED
  Supes       libraries/supes         PT  REQUIRED
  Suplib      libraries/suplib        PT  REQUIRED
  SuplibC     libraries/suplib_c      PT  REQUIRED
  SuplibCpp   libraries/suplib_cpp    PT  REQUIRED
  SVDI        libraries/svdi          ST  OPTIONAL
  PLT         libraries/plt           ST  OPTIONAL
  Algebra     applications/algebra    PT  REQUIRED
  Aprepro     applications/aprepro    PT  REQUIRED
  Blot        applications/blot       ST  OPTIONAL
  Conjoin     applications/conjoin    PT  REQUIRED
  Ejoin       applications/ejoin      PT  REQUIRED
  Epu         applications/epu        PT  REQUIRED
  Exo2mat     applications/exo2mat    ST  OPTIONAL
  Exodiff     applications/exodiff    PT  REQUIRED
  Exomatlab   applications/exomatlab  ST  REQUIRED
  Exotxt      applications/exotxt     PT  REQUIRED
  Exo_format  applications/exo_format PT  REQUIRED
  Ex1ex2v2    applications/ex1ex2v2   PT  OPTIONAL
  Fastq       applications/fastq      ST  OPTIONAL
  Gjoin       applications/gjoin      PT  REQUIRED
  Gen3D       applications/gen3d      PT  REQUIRED
  Genshell    applications/genshell   PT  REQUIRED
  Grepos      applications/grepos     PT  REQUIRED
  Grope       applications/grope      PT  REQUIRED
  Mapvarlib   libraries/mapvarlib     PT  REQUIRED
  Mapvar      applications/mapvar     PT  REQUIRED
  Mapvar-kd   applications/mapvar-kd  PT  REQUIRED
  Mat2exo     applications/mat2exo    ST  OPTIONAL
  Nemslice    applications/nem_slice  PT  REQUIRED
  Nemspread   applications/nem_spread PT  REQUIRED
  Numbers     applications/numbers    PT  REQUIRED
  Txtexo      applications/txtexo     PT  REQUIRED
  Ex2ex1v2    applications/ex2ex1v2   PT  OPTIONAL
)

SET(LIB_OPTIONAL_DEP_TPLS MPI)

