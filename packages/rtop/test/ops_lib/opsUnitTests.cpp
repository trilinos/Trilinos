
#include "RTOpPack_ROpDotProd.hpp"
#include "RTOpPack_ROpNormInf.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_Tuple.hpp"
#include "Teuchos_implicit_cast.hpp"
#include "Teuchos_as.hpp"


namespace {


// Size of the vectors
int n = 4;
// Cushion off of machine eps
double errorTolSlack = 1e-2;


class UnitTestSetup {
public:
  UnitTestSetup()
    {
      Teuchos::CommandLineProcessor &clp =
        Teuchos::UnitTestRepository::getCLP();
      clp.setOption(
        "n", &n, "Number of elements in the local vectors" );
      clp.setOption(
        "error-tol-slack", &errorTolSlack,
        "Slack off of machine epsilon used to check test results" );
    }
} unitTestSetup;


//
// Helpers
//

using Teuchos::RCP;
using Teuchos::as;
using Teuchos::tuple;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::ScalarTraits;
using RTOpPack::SubVectorView;
using RTOpPack::ConstSubVectorView;


template<class Scalar>
SubVectorView<Scalar>
newSubVectorView(const int n, const Scalar &val)
{
  ArrayRCP<Scalar> vals = Teuchos::arcp<Scalar>(n);
  std::fill(vals.begin(), vals.end(), val);
  return SubVectorView<Scalar>(
    0, n, vals, 1);
}


template<class Scalar>
SubVectorView<Scalar>
newStridedSubVectorView(const int n, const int stride, const Scalar &val)
{
  ArrayRCP<Scalar> vals = Teuchos::arcp<Scalar>(n*stride);
  std::fill(vals.begin(), vals.end(), Teuchos::ScalarTraits<Scalar>::nan());
  for (
    typename ArrayRCP<Scalar>::iterator itr = vals.begin();
    itr != vals.end();
    itr += stride
    )
  {
    *itr = val;
  }
  return SubVectorView<Scalar>(
    0, n, vals, stride);
}


//
// ROpDotProd
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpDotProd, unitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v1 = ST::random();
  const Scalar v2 = ST::random();
  out << "v1="<<v1<<", v2="<<v2<<"\n";
  SubVectorView<Scalar> sv1 = newSubVectorView<Scalar>(n, v1);
  SubVectorView<Scalar> sv2 = newSubVectorView<Scalar>(n, v2);
  RTOpPack::ROpDotProd<Scalar> dotProdOp;
  RCP<RTOpPack::ReductTarget> dotProd = dotProdOp.reduct_obj_create();
  dotProdOp.apply_op_new(
    tuple<ConstSubVectorView<Scalar> >(sv1, sv2)(),
    Teuchos::null,
    dotProd.ptr()
    );
  TEST_FLOATING_EQUALITY( dotProdOp(*dotProd), ST::conjugate(v1)*v2*as<Scalar>(n),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpDotProd, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpDotProd, nonunitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v1 = ST::random();
  const Scalar v2 = ST::random();
  out << "v1="<<v1<<", v2="<<v2<<"\n";
  SubVectorView<Scalar> sv1 = newStridedSubVectorView<Scalar>(n, 2, v1);
  SubVectorView<Scalar> sv2 = newStridedSubVectorView<Scalar>(n, 3, v2);
  RTOpPack::ROpDotProd<Scalar> dotProdOp;
  RCP<RTOpPack::ReductTarget> dotProd = dotProdOp.reduct_obj_create();
  dotProdOp.apply_op_new(
    tuple<ConstSubVectorView<Scalar> >(sv1, sv2)(),
    Teuchos::null,
    dotProd.ptr()
    );
  TEST_FLOATING_EQUALITY( dotProdOp(*dotProd), ST::conjugate(v1)*v2*as<Scalar>(n),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpDotProd, nonunitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpDotProd, reduct, Scalar )
{
  using Teuchos::as;
  using Teuchos::dyn_cast;
  using RTOpPack::ReductTargetScalar;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;

  const Scalar v1 = ST::random();
  const Scalar v2 = ST::random();
  out << "v1="<<v1<<", v2="<<v2<<"\n";

  RTOpPack::ROpDotProd<Scalar> dotProdOp;

  RCP<RTOpPack::ReductTarget> reduct1 = dotProdOp.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = dotProdOp.reduct_obj_create();

  ReductTargetScalar<Scalar> &scalarReduct1 =
    dyn_cast<ReductTargetScalar<Scalar> >(*reduct1); 
  ReductTargetScalar<Scalar> &scalarReduct2 =
    dyn_cast<ReductTargetScalar<Scalar> >(*reduct2); 

  scalarReduct1.set(v1);
  scalarReduct2.set(v2);

  dotProdOp.reduce_reduct_objs_new( *reduct1, reduct2.ptr() );

  TEST_FLOATING_EQUALITY( scalarReduct2.get(), v1+v2, as<ScalarMag>(ST::eps()*errorTolSlack) );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpDotProd, reduct )


//
// ROpNormInf
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNormInf, unitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v = ST::random();
  out << "v="<<v<<"\n";
  SubVectorView<Scalar> sv = newSubVectorView<Scalar>(n, v);
  RTOpPack::ROpNormInf<Scalar> normInfOp;
  RCP<RTOpPack::ReductTarget> normInf = normInfOp.reduct_obj_create();
  normInfOp.apply_op_new(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    normInf.ptr()
    );
  TEST_FLOATING_EQUALITY( normInfOp(*normInf), ST::magnitude(v),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNormInf, unitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNormInf, nonunitStride, Scalar )
{
  using Teuchos::as;
  typedef ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  const Scalar v = ST::random();
  out << "v="<<v<<"\n";
  SubVectorView<Scalar> sv = newStridedSubVectorView<Scalar>(n, 3, v);
  RTOpPack::ROpNormInf<Scalar> normInfOp;
  RCP<RTOpPack::ReductTarget> normInf = normInfOp.reduct_obj_create();
  normInfOp.apply_op_new(
    tuple<ConstSubVectorView<Scalar> >(sv)(),
    Teuchos::null,
    normInf.ptr()
    );
  TEST_FLOATING_EQUALITY( normInfOp(*normInf), ST::magnitude(v),
    as<ScalarMag>(ST::eps() * errorTolSlack) );
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNormInf, nonunitStride )


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( ROpNormInf, reduct, Scalar )
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

  RTOpPack::ROpNormInf<Scalar> normInfOp;

  RCP<RTOpPack::ReductTarget> reduct1 = normInfOp.reduct_obj_create();
  RCP<RTOpPack::ReductTarget> reduct2 = normInfOp.reduct_obj_create();

  ReductTargetScalar<ScalarMag> &scalarReduct1 =
    dyn_cast<ReductTargetScalar<ScalarMag> >(*reduct1); 
  ReductTargetScalar<ScalarMag> &scalarReduct2 =
    dyn_cast<ReductTargetScalar<ScalarMag> >(*reduct2); 

  scalarReduct1.set(three);
  scalarReduct2.set(four);
  normInfOp.reduce_reduct_objs_new( *reduct1, reduct2.ptr() );

  scalarReduct1.set(two);
  normInfOp.reduce_reduct_objs_new( *reduct1, reduct2.ptr() );

  TEST_EQUALITY( scalarReduct2.get(), four );

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( ROpNormInf, reduct )


} // namespace
