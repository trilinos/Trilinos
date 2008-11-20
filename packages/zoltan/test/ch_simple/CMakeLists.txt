
INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Copy_Files_To_Binary_Dir)


TRILINOS_ADD_TEST(
  ../../src/driver/zdrive
  NOEXEPREFIX
  NAME zdrive_ch_simple
  ARGS zdrive.inp.rcb
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "fatal: "
  )
# 2008/11/18: rabartl: ToDo: Above, we are *not* running the MPI test because
# it fails.  We need to talk with the Zoltan developers about this!

SET(INPUTFILES
  zdrive.inp.adaptive-partless-v3        zdrive.inp.partkway-order-metis-v3
  zdrive.inp.adaptive-partlocal1-v3      zdrive.inp.partkway-partless-v3
  zdrive.inp.adaptive-partlocal2-v3      zdrive.inp.partkway-partlocal1-v3
  zdrive.inp.adaptive-partlocal4-v3      zdrive.inp.partkway-partlocal2-v3
  zdrive.inp.adaptive-partlocal5-v3      zdrive.inp.partkway-partlocal4-v3
  zdrive.inp.adaptive-partmore-v3        zdrive.inp.partkway-partmore-v3
  zdrive.inp.adaptive-v3                 zdrive.inp.partkway-v3
  zdrive.inp.block                       zdrive.inp.partkway-v3.addobjwgt
  zdrive.inp.block-oneproc               zdrive.inp.phg
  zdrive.inp.block-partless              zdrive.inp.phg.addobjwgt
  zdrive.inp.block-partmore              zdrive.inp.random
  zdrive.inp.color                       zdrive.inp.random-oneproc
  zdrive.inp.color-after-phg             zdrive.inp.random-partless
  zdrive.inp.color-dist2                 zdrive.inp.random-partmore
  zdrive.inp.color-orderB                zdrive.inp.rcb
  zdrive.inp.color-orderU                zdrive.inp.rcb-avg
  zdrive.inp.diffusion-v3                zdrive.inp.rcb-avg-ts
  zdrive.inp.graph                       zdrive.inp.rcb-oneproc
  zdrive.inp.graph-parmetis              zdrive.inp.rcb-partless
  zdrive.inp.graph-parmetis-partition    zdrive.inp.rcb-partless-ts
  zdrive.inp.graph-parmetis-refine       zdrive.inp.rcb-partlocal1
  zdrive.inp.graph-parmetis-repartition  zdrive.inp.rcb-partlocal1-ts
  zdrive.inp.graph-partition             zdrive.inp.rcb-partlocal2
  zdrive.inp.graph-refine                zdrive.inp.rcb-partlocal2-ts
  zdrive.inp.graph-repartition           zdrive.inp.rcb-partlocal3
  zdrive.inp.hier-parmetis2-hsfc         zdrive.inp.rcb-partlocal3-ts
  zdrive.inp.hier-parmetis2-rcb          zdrive.inp.rcb-partlocal4
  zdrive.inp.hier-rcb2-parmetis          zdrive.inp.rcb-partlocal4-ts
  zdrive.inp.hsfc                        zdrive.inp.rcb-partmore
  zdrive.inp.hsfc-drum                   zdrive.inp.rcb-partmore-ts
  zdrive.inp.hsfc-oneproc                zdrive.inp.rcb-ts
  zdrive.inp.hsfc-partless               zdrive.inp.rib
  zdrive.inp.hsfc-partlocal1             zdrive.inp.rib-avg
  zdrive.inp.hsfc-partlocal2             zdrive.inp.rib-avg-ts
  zdrive.inp.hsfc-partlocal3             zdrive.inp.rib-oneproc
  zdrive.inp.hsfc-partlocal4             zdrive.inp.rib-partless
  zdrive.inp.hsfc-partmore               zdrive.inp.rib-partless-ts
  zdrive.inp.oct0                        zdrive.inp.rib-partlocal1
  zdrive.inp.oct0-partlocal3             zdrive.inp.rib-partlocal1-ts
  zdrive.inp.oct1                        zdrive.inp.rib-partlocal2
  zdrive.inp.oct1-partlocal3             zdrive.inp.rib-partlocal2-ts
  zdrive.inp.oct2                        zdrive.inp.rib-partlocal3
  zdrive.inp.oct2-partlocal3             zdrive.inp.rib-partlocal3-ts
  zdrive.inp.order-metis-v3              zdrive.inp.rib-partlocal4
  zdrive.inp.order-parmetis-v3           zdrive.inp.rib-partlocal4-ts
  zdrive.inp.partgeomkway-v3             zdrive.inp.rib-partmore
  zdrive.inp.partgeom-v3                 zdrive.inp.rib-partmore-ts
  zdrive.inp.partkway-cyclic-v3          zdrive.inp.rib-ts
  zdrive.inp.partkway-oneproc-v3
  simple.graph
  )
# 2008/11/18: rabartl: Above, we most likely do not need all of these files
# :-)  We should get Zoltan developers to work this out.


COPY_FILES_TO_BINARY_DIR(ZDriveChSimpleCopyFiles
  SOURCE_FILES ${INPUTFILES}
  DEST_FILES ${INPUTFILES}
  )
