TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  Exodus      libraries/exodus        SS  REQUIRED
  Exodus_for  libraries/exodus_for    SS  REQUIRED
  ExoIIv2for32 libraries/exoIIv2for32 SS  REQUIRED
  Nemesis     libraries/nemesis       SS  REQUIRED
  Ioss        libraries/ioss          SS  REQUIRED
  Chaco       libraries/chaco         SS  REQUIRED
  Aprepro_lib libraries/aprepro_lib   SS  REQUIRED
  Supes       libraries/supes         SS  REQUIRED
  Suplib      libraries/suplib        SS  REQUIRED
  SuplibC     libraries/suplib_c      SS  REQUIRED
  SuplibCpp   libraries/suplib_cpp    SS  REQUIRED
  SVDI        libraries/svdi          SS  OPTIONAL
  PLT         libraries/plt           SS  OPTIONAL
  Algebra     applications/algebra    SS  REQUIRED
  Aprepro-orig applications/aprepro-orig    SS  OPTIONAL
  Aprepro     applications/aprepro    SS  REQUIRED
  Blot        applications/blot       SS  OPTIONAL
  Conjoin     applications/conjoin    SS  REQUIRED
  Ejoin       applications/ejoin      SS  REQUIRED
  Epu         applications/epu        SS  REQUIRED
  Exo2mat     applications/exo2mat    SS  OPTIONAL
  Exodiff     applications/exodiff    SS  REQUIRED
  Exomatlab   applications/exomatlab  SS  REQUIRED
  Exotxt      applications/exotxt     SS  REQUIRED
  Exo_format  applications/exo_format SS  REQUIRED
  Ex1ex2v2    applications/ex1ex2v2   SS  OPTIONAL
  Fastq       applications/fastq      SS  OPTIONAL
  Gjoin       applications/gjoin      SS  REQUIRED
  Gen3D       applications/gen3d      SS  REQUIRED
  Genshell    applications/genshell   SS  REQUIRED
  Grepos      applications/grepos     SS  REQUIRED
  Grope       applications/grope      SS  REQUIRED
  Mapvarlib   libraries/mapvarlib     SS  REQUIRED
  Mapvar      applications/mapvar     SS  REQUIRED
  Mapvar-kd   applications/mapvar-kd  SS  REQUIRED
  Mat2exo     applications/mat2exo    SS  OPTIONAL
  Nemslice    applications/nem_slice  SS  REQUIRED
  Nemspread   applications/nem_spread SS  REQUIRED
  Numbers     applications/numbers    SS  REQUIRED
  Txtexo      applications/txtexo     SS  REQUIRED
  Ex2ex1v2    applications/ex2ex1v2   SS  OPTIONAL
)

SET(LIB_OPTIONAL_DEP_TPLS MPI)

