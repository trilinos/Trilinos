// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_EXAMPLE_TRIDIAG_SERIAL_LINEAR_OP_HPP
#define THYRA_EXAMPLE_TRIDIAG_SERIAL_LINEAR_OP_HPP

#include "Thyra_LinearOpDefaultBase.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Teuchos_Assert.hpp"


/** \brief Simple example subclass for serial tridiagonal matrices.
 *
 * This subclass form a linear operator for tridiagonal matrices
 * of the form:
 *
 \f[

 M=
 \left[\begin{array}{ccccc}
 d_{(1)} & u_{(1)}  \\
 l_{(1)} & d_{(2)} & u_{(2)} \\
         & \ddots  & \ddots    & \ddots \\
         &         & l_{(n-2)} & d_{(n-1)}  & u_{(n-1)} \\
         &         &           & l_{(n-1)}  & d_{(n)}
 \end{array}\right]
 \f]
 *
 * given the arrays <tt>lower[]</tt>, <tt>diag[]</tt>, and
 * <tt>upper[]</tt> of dimension <tt>dim-1</tt>, <tt>dim</tt> and <tt>dim-1</tt>
 * respectively (see <tt>initialize()</tt>).
 *
 * See the source code for this simple example by clicking on the
 * link to the definition below.
 *
 * \ingroup Thyra_Op_Vec_examples_power_method_serial_grp
 * \ingroup Thyra_Op_Vec_examples_cg_serial_grp
 */
template<class Scalar>
class ExampleTridiagSerialLinearOp : public Thyra::LinearOpDefaultBase<Scalar>
{
public:

  /** \brief Construct to uninitialized. */
  ExampleTridiagSerialLinearOp() {}

  /** \brief <tt>initialize()</tt>. */
  ExampleTridiagSerialLinearOp(
    const Thyra::Ordinal dim,
    const Teuchos::ArrayView<const Scalar> &lower,
    const Teuchos::ArrayView<const Scalar> &diag,
    const Teuchos::ArrayView<const Scalar> &upper
    )
    { this->initialize(dim, lower, diag, upper);	}
  
  /** Initialize given lower, diagonal and upper arrays of data.
   *
   * \param dim [in] Dimension of this matrix (must be >= 2).
   *
   * \param lower [in] Array (length <tt>dim-1</tt>) of the lower diagonal
   * elements
   *
   * \param diag [in] Array (length <tt>dim</tt>) of the central diagonal
   * elements
   *
   * \param upper [in] Array (length <tt>dim-1</tt>) of the upper diagonal
   * elements
   *
   * Preconditions:<ul>
   * <li><tt>dim >= 2</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li>Should be obvious!
   * </ul>
   */
  void initialize(
    const Thyra::Ordinal dim,
    const Teuchos::ArrayView<const Scalar> &lower,
    const Teuchos::ArrayView<const Scalar> &diag,
    const Teuchos::ArrayView<const Scalar> &upper
    )
    {
      TEST_FOR_EXCEPT( dim < 2 );
      space_ = Thyra::defaultSpmdVectorSpace<Scalar>(dim);
      lower_ = lower;
      diag_ = diag;
      upper_ = upper;
    }

protected:

  // Overridden from LinearOpBase

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > range() const
    { return space_; }

  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > domain() const
    { return space_; }

  /** \brief . */
  bool opSupportedImpl(Thyra::EOpTransp M_trans) const
    { return true; }  // This class supports everything!

  /** \brief . */
  void applyImpl(
    const Thyra::EOpTransp M_trans,
    const Thyra::MultiVectorBase<Scalar> &X_in,
    const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> > &Y_inout,
    const Scalar alpha,
    const Scalar beta
    ) const;

private:

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > space_;
  Teuchos::Array<Scalar> lower_;   // size = dim - 1    
  Teuchos::Array<Scalar> diag_;    // size = dim
  Teuchos::Array<Scalar> upper_;   // size = dim - 1

};	// end class ExampleTridiagSerialLinearOp


template<class Scalar>
void ExampleTridiagSerialLinearOp<Scalar>::applyImpl(
  const Thyra::EOpTransp M_trans,
  const Thyra::MultiVectorBase<Scalar> &X_in,
  const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> > &Y_inout,
  const Scalar alpha,
  const Scalar beta
  ) const
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::Ordinal Ordinal;

  const Ordinal dim = space_->dim();
      
  // Loop over the input columns

  const Ordinal m = X_in.domain()->dim();

  for (Ordinal col_j = 0; col_j < m; ++col_j) {

    // Get access the the elements of column j
    Thyra::ConstDetachedVectorView<Scalar> x_vec(X_in.col(col_j));
    Thyra::DetachedVectorView<Scalar> y_vec(Y_inout->col(col_j));
    const Teuchos::ArrayRCP<const Scalar> x = x_vec.sv().values();
    const Teuchos::ArrayRCP<Scalar> y = y_vec.sv().values();
        
    // Perform y = beta*y (being careful to set y=0 if beta=0 in case y is
    // uninitialized on input!)
    if( beta == ST::zero() ) {
      for( Ordinal k = 0; k < dim; ++k ) y[k] = ST::zero();
    }
    else if( beta != ST::one() ) {
      for( Ordinal k = 0; k < dim; ++k ) y[k] *= beta;
    }

    // Perform y = alpha*op(M)*x 
    Ordinal k = 0;
    if( M_trans == Thyra::NOTRANS ) {
      y[k] += alpha * ( diag_[k]*x[k] + upper_[k]*x[k+1] );  // First row
      for( k = 1; k < dim - 1; ++k )   // Middle rows
        y[k] += alpha * ( lower_[k-1]*x[k-1] + diag_[k]*x[k] + upper_[k]*x[k+1] );
      y[k] += alpha * ( lower_[k-1]*x[k-1] + diag_[k]*x[k] ); // Last row
    }
    else if( M_trans == Thyra::CONJ ) {
      y[k] += alpha * ( ST::conjugate(diag_[k])*x[k] + ST::conjugate(upper_[k])*x[k+1] );
      for( k = 1; k < dim - 1; ++k )
        y[k] += alpha * ( ST::conjugate(lower_[k-1])*x[k-1]
          + ST::conjugate(diag_[k])*x[k] + ST::conjugate(upper_[k])*x[k+1] );
      y[k] += alpha * ( ST::conjugate(lower_[k-1])*x[k-1] + ST::conjugate(diag_[k])*x[k] );
    }
    else if( M_trans == Thyra::TRANS ) {
      y[k] += alpha * ( diag_[k]*x[k] + lower_[k]*x[k+1] );
      for( k = 1; k < dim - 1; ++k )
        y[k] += alpha * ( upper_[k-1]*x[k-1] + diag_[k]*x[k] + lower_[k]*x[k+1] );
      y[k] += alpha * ( upper_[k-1]*x[k-1] + diag_[k]*x[k] );
    }
    else if( M_trans == Thyra::CONJTRANS ) {
      y[k] += alpha * ( ST::conjugate(diag_[k])*x[k] + ST::conjugate(lower_[k])*x[k+1] );
      for( k = 1; k < dim - 1; ++k )
        y[k] += alpha * ( ST::conjugate(upper_[k-1])*x[k-1]
          + ST::conjugate(diag_[k])*x[k] + ST::conjugate(lower_[k])*x[k+1] );
      y[k] += alpha * ( ST::conjugate(upper_[k-1])*x[k-1] + ST::conjugate(diag_[k])*x[k] );
    }
    else {
      TEST_FOR_EXCEPT(true); // Throw exception if we get here!
    }
  }

}


#endif	// THYRA_EXAMPLE_TRIDIAG_SERIAL_LINEAR_OP_HPP
