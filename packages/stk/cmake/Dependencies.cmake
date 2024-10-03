# the order of the packages below is important.  Stk packages can
# depend on packages listed above them, and not below them.
SET(SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
  Util             stk_util             PT  OPTIONAL
  Coupling         stk_coupling         PT  OPTIONAL
  Math             stk_math             PT  OPTIONAL
  Simd             stk_simd             PT  OPTIONAL
  NGP_TEST         stk_ngp_test         PT  OPTIONAL
  ExprEval         stk_expreval         PT  OPTIONAL
  Topology         stk_topology         PT  OPTIONAL
  Search           stk_search           PT  OPTIONAL
  Middle_mesh      stk_middle_mesh      EX  OPTIONAL
  Transfer         stk_transfer         PT  OPTIONAL
  Mesh             stk_mesh             PT  OPTIONAL
  SearchUtil       stk_search_util      PT  OPTIONAL
  TransferUtil     stk_transfer_util    EX  OPTIONAL
  IO               stk_io               PT  OPTIONAL
  Middle_mesh_util stk_middle_mesh_util EX  OPTIONAL
  Tools            stk_tools            PT  OPTIONAL
  Balance          stk_balance          PT  OPTIONAL
  Unit_test_utils  stk_unit_test_utils  PT  OPTIONAL
  Unit_tests       stk_unit_tests       PT  OPTIONAL
  Doc_tests        stk_doc_tests        PT  OPTIONAL
  Emend            stk_emend            PT  OPTIONAL
)

SET(LIB_REQUIRED_DEP_PACKAGES)
SET(LIB_OPTIONAL_DEP_PACKAGES)
SET(TEST_REQUIRED_DEP_PACKAGES)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS MPI)
SET(TEST_REQUIRED_DEP_TPLS)
SET(TEST_OPTIONAL_DEP_TPLS)
