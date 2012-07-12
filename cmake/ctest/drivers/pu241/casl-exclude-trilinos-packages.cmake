# These are all packages that are not needed by CASL development.  This list
# is the only list that need to be maintained to exclude and disable Trilinos
# packages.  This list is read in and used in a variety of places.

SET(${PROJECT_NAME}_EXCLUDE_PACKAGES
  GlobiPack
  OptiPack
  Pliris
  Claps
  Galeri
  Pamgen
  Komplex
  RBGen
  TriKota
  Phdmesh
  Moertel
  TrilinosCouplings
  MOOCHO
  Stokhos
  Sundance
  CTrilinos
  ForTrilinos
  PyTrilinos
  Didasko
  Optika
  Mesquite
  FEApp
  Zoltan2
  )
