
#ifndef THYRA_DIAGONAL_SCALAR_PROD_DECL_HPP
#define THYRA_DIAGONAL_SCALAR_PROD_DECL_HPP


#include "Thyra_ScalarProdBase.hpp"


namespace Thyra {


/** \brief Concrete implementation of a scalar product using a diagonal
 * vector.
 *
 * This test class really shows how to create an application-defined scalar
 * product.
 */
template<class Scalar>
class DiagonalScalarProd : public ScalarProdBase<Scalar> {
public:
  
  /** @name Consturctors/Initializers/Accessors */
  //@{

  /** \brief . */
  DiagonalScalarProd();

  /** \brief . */
  void initialize( const RCP<const VectorBase<Scalar> > &s_diag );

  //@}

protected:
  
  /** @name Overridden protected virtual functions from ScalarProdBase */
  //@{

  /** \brief Returns <tt>false</tt>. */
  virtual bool isEuclideanImpl() const;
  
  /** \brief . */
  virtual void scalarProdsImpl(
    const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
    const ArrayView<Scalar> &scalarProds_out ) const;

  /** \brief . */
  RCP<const LinearOpBase<Scalar> > getLinearOpImpl() const;

  //@}

private:

  RCP<const VectorBase<Scalar> > s_diag_;

};


/** \brief Nonmember constructor.
 *
 * \relates DiagonalScalarProd
 */
template<class Scalar>
RCP<DiagonalScalarProd<Scalar> >
diagonalScalarProd(const RCP<const VectorBase<Scalar> > &s_diag)
{
  const RCP<DiagonalScalarProd<Scalar> > scalarProd =
    Teuchos::rcp(new DiagonalScalarProd<Scalar>());
  scalarProd->initialize(s_diag);
  return scalarProd;
}



} // end namespace Thyra


#endif  // THYRA_DIAGONAL_SCALAR_PROD_DECL_HPP
