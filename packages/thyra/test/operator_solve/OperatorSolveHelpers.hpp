#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"


namespace Thyra {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::inOutArg;


/** \brief Generate a very well scaled unsymmetic dense matrix for
 * use in testing a dense solve.
*/
template<class Scalar>
RCP<MultiVectorBase<Scalar> >
createNonsingularMultiVector(const RCP<const VectorSpaceBase<Scalar> > &vs)
{
  using Teuchos::as;
  const RCP<MultiVectorBase<Scalar> > M = createMembers(vs, vs);
  assign(M.ptr(), as<Scalar>(0.0));
  {
    DetachedMultiVectorView<Scalar> M_detached_view(*M);
    const RTOpPack::SubMultiVectorView<Scalar> M_smvv = M_detached_view.smv(); 
    const Scalar
      one = as<Scalar>(1.0),
      two = as<Scalar>(2.0),
      four = as<Scalar>(4.0);
    for (int i = 0; i < M_smvv.subDim(); ++i) {
      if (i > 0) M_smvv(i, i-1) = one;
      M_smvv(i, i) = four;
      if (i < M_smvv.subDim() - 1) M_smvv(i, i+1) = two;
    }
  }
  return M;
}


} // namespace Thyra
