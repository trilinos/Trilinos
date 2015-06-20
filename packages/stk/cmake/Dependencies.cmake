SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  Classic         stk_classic          SS  OPTIONAL
  Util            stk_util             SS  OPTIONAL
  Topology        stk_topology         SS  OPTIONAL
  Mesh            stk_mesh             SS  OPTIONAL
  IO              stk_io               SS  OPTIONAL
  Search          stk_search           SS  OPTIONAL
  SearchUtil      stk_search_util      SS  OPTIONAL
  Transfer        stk_transfer         SS  OPTIONAL
  Unit_test_utils stk_unit_test_utils  SS  OPTIONAL
  Unit_tests      stk_unit_tests       SS  OPTIONAL
  Doc_tests       stk_doc_tests        SS  OPTIONAL
  Exp             stk_exp              EX  OPTIONAL
)

SET(LIB_REQUIRED_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_PACKAGES)
SET(TEST_REQUIRED_DEP_PACKAGES Gtest)
SET(TEST_OPTIONAL_DEP_PACKAGES KokkosCore)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS MPI)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)
