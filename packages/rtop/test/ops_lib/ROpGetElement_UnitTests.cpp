
#include "RTOpPack_ROpGetElement.hpp"
#include "RTOpPack_TOpSetAssendingValues.hpp"
#include "Teuchos_ExplicitInstantiationHelpers.hpp"
#include "Teuchos_implicit_cast.hpp"

// Must come last!
#include "supportUnitTestsHelpers.hpp"


namespace RTOpPack {


const double g_tol = 100.0 * ScalarTraits<double>::eps();


template<class Scalar>
void basicTest(const int size, const int stride, FancyOStream &out, bool &success)
{
  using Teuchos::as;
  using Teuchos::null;
  using Teuchos::tuple;
  typedef ScalarTraits<Scalar> ST;

  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(size, stride, ST::zero());

  TOpSetAssendingValues<Scalar> setAssendingOp;
  setAssendingOp.apply_op(null, tuple<SubVectorView<Scalar> >(sv)(), null);

  ROpGetElement<Scalar> getEleOp(0);
  const RCP<ReductTarget> reduct_obj = getEleOp.reduct_obj_create();

  getEleOp.globalIndex(0);
  getEleOp.reduct_obj_reinit(reduct_obj.ptr());
  getEleOp.apply_op(tuple<ConstSubVectorView<Scalar> >(sv)(), null, reduct_obj.ptr());
  TEST_EQUALITY_CONST(getEleOp(*reduct_obj), as<Scalar>(1));

  getEleOp.globalIndex(1);
  getEleOp.reduct_obj_reinit(reduct_obj.ptr());
  getEleOp.apply_op(tuple<ConstSubVectorView<Scalar> >(sv)(), null, reduct_obj.ptr());
  TEST_EQUALITY_CONST(getEleOp(*reduct_obj), as<Scalar>(2));

  getEleOp.globalIndex(size-2);
  getEleOp.reduct_obj_reinit(reduct_obj.ptr());
  getEleOp.apply_op(tuple<ConstSubVectorView<Scalar> >(sv)(), null, reduct_obj.ptr());
  TEST_EQUALITY_CONST(getEleOp(*reduct_obj), as<Scalar>(size-1));

  getEleOp.globalIndex(size-1);
  getEleOp.reduct_obj_reinit(reduct_obj.ptr());
  getEleOp.apply_op(tuple<ConstSubVectorView<Scalar> >(sv)(), null, reduct_obj.ptr());
  TEST_EQUALITY_CONST(getEleOp(*reduct_obj), as<Scalar>(size));

}


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpGetElement, basicTest_4_1, Scalar )
{
  basicTest<Scalar>(4, 1, out, success);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpGetElement, basicTest_4_1)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpGetElement, basicTest_4_3, Scalar )
{
  basicTest<Scalar>(4, 3, out, success);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpGetElement, basicTest_4_3)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpGetElement, basicTest_6_1, Scalar )
{
  basicTest<Scalar>(6, 1, out, success);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpGetElement, basicTest_6_1)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpGetElement, basicTest_6_3, Scalar )
{
  basicTest<Scalar>(6, 3, out, success);
}
TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpGetElement, basicTest_6_3)


} // namespace RTOpPack
