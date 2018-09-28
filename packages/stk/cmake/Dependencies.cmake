SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  Util            stk_util             PT  OPTIONAL
  Simd            stk_simd             PT  OPTIONAL
  Topology        stk_topology         PT  OPTIONAL
  Mesh            stk_mesh             PT  OPTIONAL
  NGP             stk_ngp              PT  OPTIONAL
  IO              stk_io               PT  OPTIONAL
  Unit_test_utils stk_unit_test_utils  PT  OPTIONAL
  Math            stk_math             PT  OPTIONAL
  Search          stk_search           PT  OPTIONAL
  SearchUtil      stk_search_util      PT  OPTIONAL
  Transfer        stk_transfer         PT  OPTIONAL
  Tools           stk_tools            PT  OPTIONAL
  Unit_tests      stk_unit_tests       PT  OPTIONAL
  Doc_tests       stk_doc_tests        PT  OPTIONAL
  Exp             stk_exp              EX  OPTIONAL
  ExprEval        stk_expreval         PT  OPTIONAL
)

SET(LIB_REQUIRED_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_PACKAGES)
SET(TEST_REQUIRED_DEP_PACKAGES Gtest)
SET(TEST_OPTIONAL_DEP_PACKAGES KokkosCore)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS MPI)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)
