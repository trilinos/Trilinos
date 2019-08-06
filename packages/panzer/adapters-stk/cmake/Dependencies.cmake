# NOTE: Pamgen is an indirect dependency brought in through
# SEACASIOSS. We do NOT need to declare it as a direct library
# dependency. The pamgen stk_interface factory unit test does need
# Pamgen enabled to run, so Pamgen is only declared an optional test
# dependency.

SET(LIB_REQUIRED_DEP_PACKAGES STKUtil STKTopology STKMesh STKIO Zoltan Stratimikos Piro NOX Rythmos PanzerCore PanzerDiscFE)
SET(LIB_OPTIONAL_DEP_PACKAGES SEACASIoss SEACASExodus Percept Teko MueLu Ifpack2)
SET(TEST_REQUIRED_DEP_PACKAGES SEACASIoss SEACASExodus Teko MueLu Ifpack2)
SET(TEST_OPTIONAL_DEP_PACKAGES Pamgen)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)
