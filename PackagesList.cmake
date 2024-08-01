#
# Define the Trilinos packages
#
TRIBITS_REPOSITORY_DEFINE_PACKAGES(
  TrilinosFrameworkTests  commonTools/framework           PT
  TrilinosATDMConfigTests cmake/std/atdm                  PT
  Gtest                 commonTools/gtest                 PT
  Kokkos                packages/kokkos                   PT
  Teuchos               packages/teuchos                  PT
  KokkosKernels         packages/kokkos-kernels           PT
  RTOp                  packages/rtop                     PT
  Sacado                packages/sacado                   PT
  MiniTensor            packages/minitensor               PT
  Epetra                packages/epetra                   ST
  SCOREClion            SCOREC/lion                       ST
  SCORECpcu             SCOREC/pcu                        ST
  SCORECgmi             SCOREC/gmi                        ST
  SCORECgmi_sim         SCOREC/gmi_sim                    ST
  SCORECapf             SCOREC/apf                        ST
  SCORECapf_sim         SCOREC/apf_sim                    ST
  SCORECmds             SCOREC/mds                        ST
  SCORECparma           SCOREC/parma                      ST
  SCORECspr             SCOREC/spr                        ST
  AvatarT               packages/avatart                  EX
  Zoltan                packages/zoltan                   PT
  Shards                packages/shards                   PT
  Triutils              packages/triutils                 ST
  EpetraExt             packages/epetraext                ST
  Tpetra                packages/tpetra                   PT
  TrilinosSS            packages/common/auxiliarySoftware/SuiteSparse PT # Auxiliary software.
  Domi                  packages/domi                     PT
  Thyra                 packages/thyra                    PT
  Xpetra                packages/xpetra                   PT
  Isorropia             packages/isorropia                ST
  Pliris                packages/pliris                   ST
  AztecOO               packages/aztecoo                  ST
  Galeri                packages/galeri                   PT
  Amesos                packages/amesos                   ST
  Pamgen                packages/pamgen                   PT
  Zoltan2Core           packages/zoltan2/core             PT
  Ifpack                packages/ifpack                   ST
  ML                    packages/ml                       ST
  Belos                 packages/belos                    PT
  ShyLU_Node            packages/shylu/shylu_node         PT
  Amesos2               packages/amesos2                  PT
  SEACAS                packages/seacas                   PT # Depends on netcdf, optionally hdf5, xdmf, pamgen
  Komplex               packages/komplex                  ST
  Anasazi               packages/anasazi                  PT
  Ifpack2               packages/ifpack2                  PT
  Stratimikos           packages/stratimikos              PT
  FEI                   packages/fei                      PT
  Teko                  packages/teko                     PT
  TriKota               packages/TriKota                  ST
  Intrepid              packages/intrepid                 ST
  Intrepid2             packages/intrepid2                PT
  Compadre              packages/compadre                 ST
  STK                   packages/stk                      PT # Depends on boost
  Percept               packages/percept                  PT # Depends on boost
  Krino                 packages/krino                    PT # Depends on boost
  SCORECapf_zoltan      SCOREC/zoltan                     ST
  SCORECapf_stk         SCOREC/stk                        ST
  SCORECma              SCOREC/ma                         ST
  SCORECpumi            SCOREC/pumi                       ST
  SCOREC                SCOREC                            ST
  Phalanx               packages/phalanx                  PT
  NOX                   packages/nox                      PT
  Moertel               packages/moertel                  ST
  MueLu                 packages/muelu                    PT
  TrilinosLinearSolvers packages/trilinos_linear_solvers  PT
  Zoltan2Sphynx         packages/zoltan2/sphynx           PT
  Zoltan2               packages/zoltan2                  PT
  ShyLU_DD              packages/shylu/shylu_dd           PT
  ShyLU                 packages/shylu                    PT
  Rythmos               packages/rythmos                  PT
  Tempus                packages/tempus                   PT
  MOOCHO                packages/moocho                   ST
  Stokhos               packages/stokhos                  PT
  ROL                   packages/rol                      PT
  Piro                  packages/piro                     PT
  SGM                   packages/sgm                      ST
  UMR                   packages/umr                      ST
  Panzer                packages/panzer                   PT
  CTrilinos             packages/CTrilinos                ST # Switched to ST to speed up checkin testing
  PyTrilinos            packages/PyTrilinos               ST
  PyTrilinos2           packages/PyTrilinos2              EX
  WebTrilinos           packages/WebTrilinos              EX # Should be ST
  NewPackage            packages/new_package              EX # Should be ST
  Optika		packages/optika		          EX
  Adelus                packages/adelus                   PT
  TrilinosCouplings     packages/trilinoscouplings        PT
  Pike                  packages/pike                     PT
  xSDKTrilinos          packages/xSDKTrilinos             ST
  TrilinosBuildStats    commonTools/build_stats           PT
  TrilinosInstallTests  packages/TrilinosInstallTests     PT
  )

# Allow builds even if some packages are missing
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(AvatarT)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCOREC)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCOREClion)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECgmi)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECgmi_sim)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECpcu)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECapf)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECapf_sim)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECmds)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECparma)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECspr)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECapf_stk)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECapf_zoltan)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECma)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SCORECpumi)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Avatar)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(MOOCHO)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Sundance)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(CTrilinos)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Optika)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Mesquite)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(WebTrilinos)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(xSDKTrilinos)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(SGM)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(UMR)
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(TrilinosLinearSolvers)

# TRILFRAME-500
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Rythmos)    # 27115 targets
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Pike)       # 27048 targets
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Komplex)    # 27030 targets
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Moertel)    # 26995 targets
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(TriKota)    # 26995 targets
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(Domi)       # 26946 targets
TRIBITS_ALLOW_MISSING_EXTERNAL_PACKAGES(FEI)

#
# Disable certain packages on certain platforms.
#
# NOTE: This just makes the packages experimental 'EX' and therefore still
# allows the user to enable the package explicitly but the package will not
# get enabled implicitly.
#

TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(MOOCHO Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Phalanx Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(PyTrilinos Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Sundance Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Tpetra Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Ifpack2 Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(TriKota Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Pamgen Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(STK Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Anasazi Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Isorropia Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Zoltan Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Teko Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Panzer Windows)
TRIBITS_DISABLE_PACKAGE_ON_PLATFORMS(Compadre Windows)
