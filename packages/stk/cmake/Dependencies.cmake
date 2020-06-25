SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  Math            stk_math             PT  OPTIONAL
  Util            stk_util             PT  OPTIONAL
  Simd            stk_simd             PT  OPTIONAL
  Topology        stk_topology         PT  OPTIONAL
  Mesh            stk_mesh             PT  OPTIONAL
  NGP             stk_ngp              EX  OPTIONAL
  IO              stk_io               PT  OPTIONAL
  NGP_TEST        stk_ngp_test         PT  OPTIONAL
  Unit_test_utils stk_unit_test_utils  PT  OPTIONAL
  Search          stk_search           PT  OPTIONAL
  SearchUtil      stk_search_util      PT  OPTIONAL
  Transfer        stk_transfer         PT  OPTIONAL
  Tools           stk_tools            PT  OPTIONAL
  Balance         stk_balance          PT  OPTIONAL
  Unit_tests      stk_unit_tests       PT  OPTIONAL
  Doc_tests       stk_doc_tests        PT  OPTIONAL
  ExprEval        stk_expreval         PT  OPTIONAL
)

SET(LIB_REQUIRED_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_PACKAGES)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS MPI)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)
