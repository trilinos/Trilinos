
#include "RTOpPack_ROpNorm2.hpp"
#include "opsUnitTestsHelpers.hpp"


namespace {


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNorm2, unitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;
  const Scalar v = ST::random();
  out << "v="<<v<<"\n";
  SubVectorView<Scalar> sv = newSubVectorView<Scalar>(n, v);
  RTOpPack::ROpNorm2<Scalar> norm2Op;
  RCP<RTOpPack::ReductTarget> norm2 = norm2Op.reduct_obj_create();
  norm2Op.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    norm2.ptr()
    );
  TEST_FLOATING_EQUALITY( norm2Op(*norm2),
    SMT::squareroot(n)*ST::magnitude(v),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNorm2, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNorm2, nonunitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;
  const Scalar v = ST::random();
  out << "v="<<v<<"\n";
  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(n, 3, v);
  RTOpPack::ROpNorm2<Scalar> norm2Op;
  RCP<RTOpPack::ReductTarget> norm2 = norm2Op.reduct_obj_create();
  norm2Op.apply_op(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    norm2.ptr()
    );
  TEST_FLOATING_EQUALITY( norm2Op(*norm2),
    SMT::squareroot(n)*ST::magnitude(v),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNorm2, nonunitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNorm2, reduct, Scalar )
{
  using Teuchos::as;
  using Teuchos::dyn_cast;
  using RTOpPack::ReductTargetScalar;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef ScalarTraits<ScalarMag> SMT;

  const ScalarMag three = as<ScalarMag>(3.0);
  const ScalarMag four = as<ScalarMag>(4.0);
  const ScalarMag two = as<ScalarMag>(2.0);

  RTOpPack::ROpNorm2<Scalar> norm2Op;

  RCP<RTOpPack::ReductTarget> reduct1 = norm2Op.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = norm2Op.reduct_obj_create();

  ReductTargetScalar<Scalar> &scalarReduct1 =
    dyn_cast<ReductTargetScalar<Scalar> >(*reduct1); 
  ReductTargetScalar<Scalar> &scalarReduct2 =
    dyn_cast<ReductTargetScalar<Scalar> >(*reduct2); 

  scalarReduct1.set(three);
  scalarReduct2.set(four);
  norm2Op.reduct_reduct_objs( *reduct1, reduct2.ptr() );

  scalarReduct1.set(two);
  norm2Op.reduct_reduct_objs( *reduct1, reduct2.ptr() );

  TEST_FLOATING_EQUALITY( norm2Op(*reduct2), SMT::squareroot(three+four+two),
    as<ScalarMag>(ST::eps() * errorTolSlack) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNorm2, reduct )


} // namespace
