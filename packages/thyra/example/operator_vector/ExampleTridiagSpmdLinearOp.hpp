// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_EXAMPLE_TRIDIAG_SPMD_LINEAR_OP_HPP
#define THYRA_EXAMPLE_TRIDIAG_SPMD_LINEAR_OP_HPP

#include "Thyra_LinearOpDefaultBase.hpp"
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DetachedSpmdVectorView.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"


/** \brief Simple example subclass for Spmd tridiagonal matrices.
 *
 * This subclass represents a linear operator for tridiagonal matrices
 * of the global form:
 *
 \f[

 M=
 \left[\begin{array}{ccccc}
 d_{(1)} & u_{(1)}  \\
 l_{(1)} & d_{(2)} & u_{(2)} \\
         & \ddots  & \ddots    & \ddots \\
         &         & l_{(n-2)} & d_{(n-1)}  & u_{(n-1)} \\
         &         &           & l_{(n-1)}  & d_{(n)}
 \end{array}\right].
 \f]
 *
 * If there is only \f$N = 1\f$ processes then the arrays <tt>lower[]</tt>,
 * <tt>diag[]</tt>, and <tt>upper[]</tt> of dimension <tt>localDim-1</tt>,
 * <tt>localDim</tt> and <tt>localDim-1</tt> respectively are stored (see
 * <tt>initialize()</tt>).
 *
 * If there \f$N > 1\f$ processes then locally this class stores
 * slightly different arrays of data depending on which process this
 * is and process-to-process communication is required.
 *
 * <ul>
 *
 * <li> On process 0 with \f$n_0\f$ local vector elements the
 * following sub-matrix is stored:
 *
 \f[

 \left[\begin{array}{cccccc}
 d_{(1)} & u_{(1)}  \\
 l_{(1)} & d_{(2)} & u_{(2)} \\
         & \ddots  & \ddots    & \ddots \\
         &         & l_{(n_0-2)} & d_{(n_0-1)}  & u_{(n_0-1)} \\
         &         &           & l_{(n_0-1)}  & d_{(n_0)}  & u_{(n_0)}
 \end{array}\right].
 \f]
 *
 * In this case, arrays <tt>lower[]</tt>, <tt>diag[]</tt>, and
 * <tt>upper[]</tt> of dimension <tt>localDim-1</tt>,
 * <tt>localDim</tt> and <tt>localDim</tt> respectively are stored
 * (see <tt>initialize()</tt>).
 *
 * <li> On process \f$i\f$, where \f$0 < i < N-1\f$, with local
 * offset \f$g_i\f$ and \f$n_i\f$ local vector elements the following
 * sub-matrix is stored:
 *
 \f[

 \left[\begin{array}{cccccc}
 l_{(g_i)} & d_{(g_i+1)} & u_{(g_i+1)} \\
         & \ddots  & \ddots    & \ddots \\
         &         &           &  l_{(g_i+n_i-1)} & d_{(g_i+n_i)} & u_{(g_i+n_i)}
 \end{array}\right].
 \f]
 *
 * In this case, arrays <tt>lower[]</tt>, <tt>diag[]</tt>, and
 * <tt>upper[]</tt> of dimension <tt>localDim</tt>,
 * <tt>localDim</tt> and <tt>localDim</tt> respectively are stored
 * (see <tt>initialize()</tt>).
 *
 * <li> On process \f$N-1\f$ with local offset \f$g_{N-1}\f$ and
 * \f$n_{N-1}\f$ local vector elements the following sub-matrix is
 * stored:
 *
 \f[

 \left[\begin{array}{cccccc}
 l_{(g_{N-1})} & d_{(g_{N-1}+1)} & u_{(g_{N-1}+1)} \\
         & \ddots  & \ddots    & \ddots \\
         &         &           &  l_{(g_{N-1}+n_{N-1}-1)} & d_{(g_{N-1}+n_{N-1})}
 \end{array}\right].
 \f]
 *
 * In this case, arrays <tt>lower[]</tt>, <tt>diag[]</tt>, and
 * <tt>upper[]</tt> of dimension <tt>localDim</tt>,
 * <tt>localDim</tt> and <tt>localDim-1</tt> respectively are stored
 * (see <tt>initialize()</tt>).
 *
 * </ul>
 *
 * See the source code for this simple example by clicking on the
 * link to the definition below.
 *
 * \ingroup Thyra_Op_Vec_examples_cg_Spmd_grp
 */
template<class Scalar>
class ExampleTridiagSpmdLinearOp : public Thyra::LinearOpDefaultBase<Scalar> {
public:

  /** Construct to uninitialized. */
  ExampleTridiagSpmdLinearOp() {}

  /** Calls <tt>initialize()</tt>. */
  ExampleTridiagSpmdLinearOp(
    const Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > &comm,
    const Thyra::Ordinal localDim,
    const Teuchos::ArrayView<const Scalar> &lower,
    const Teuchos::ArrayView<const Scalar> &diag,
    const Teuchos::ArrayView<const Scalar> &upper
    )
    { this->initialize(comm, localDim, lower, diag, upper); }
  
  /** Initialize given lower, diagonal and upper arrays of data.
   *
   * \param comm [in] Communicator (allowed to be Teuchos::null)
   *
   * \param localDim [in] Dimension of this matrix (must be >= 2).
   *
   * \param lower [in] Array (length <tt>( procRank == 0 ? localDim - 1 :
   * localDim )</tt>) of the lower diagonal elements \param diag [in] Array
   * (length <tt>localDim</tt>) of the central diagonal elements
   *
   * \param upper [in] Array (length <tt>( procRank == numProc-1 ? localDim -
   * 1 : localDim )</tt>) of the upper diagonal elements
   *
   * Preconditions:<ul>
   * <li><tt>localDim >= 2</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li>Should be obvious!
   * </ul>
   */
  void initialize(
    const Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > &comm,
    const Thyra::Ordinal localDim,
    const Teuchos::ArrayView<const Scalar> &lower,
    const Teuchos::ArrayView<const Scalar> &diag,
    const Teuchos::ArrayView<const Scalar> &upper
    );

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
    {
      return (M_trans == Thyra::NOTRANS);
    }

  /** \brief . */
  void applyImpl(
    const Thyra::EOpTransp M_trans,
    const Thyra::MultiVectorBase<Scalar> &X_in,
    const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> > &Y_inout,
    const Scalar alpha,
    const Scalar beta
    ) const;

private:

  Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > comm_;
  Teuchos::RCP<const Thyra::SpmdVectorSpaceBase<Scalar> > space_;
  Teuchos::Array<Scalar> lower_; // size = ( procRank == 0 ? localDim - 1 : localDim ) 
  Teuchos::Array<Scalar> diag_; // size = localDim
  Teuchos::Array<Scalar> upper_; // size = ( procRank == numProc-1 ? localDim - 1 : localDim )

  void communicate( const bool first, const bool last,
    const Teuchos::ArrayView<const Scalar> &x,
    const Teuchos::Ptr<Scalar> &x_km1,
    const Teuchos::Ptr<Scalar> &x_kp1 ) const;

};	// end class ExampleTridiagSpmdLinearOp


template<class Scalar>
void ExampleTridiagSpmdLinearOp<Scalar>::initialize(
  const Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> > &comm,
  const Thyra::Ordinal localDim,
  const Teuchos::ArrayView<const Scalar> &lower,
  const Teuchos::ArrayView<const Scalar> &diag,
  const Teuchos::ArrayView<const Scalar> &upper
  )
{
  TEST_FOR_EXCEPT( localDim < 2 );
  comm_ = comm;
  space_ = Thyra::defaultSpmdVectorSpace<Scalar>(comm, localDim, -1);
  lower_ = lower;
  diag_ = diag;
  upper_ = upper;
}


template<class Scalar>
void ExampleTridiagSpmdLinearOp<Scalar>::applyImpl(
  const Thyra::EOpTransp M_trans,
  const Thyra::MultiVectorBase<Scalar> &X_in,
  const Teuchos::Ptr<Thyra::MultiVectorBase<Scalar> > &Y_inout,
  const Scalar alpha,
  const Scalar beta
  ) const
{

  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef Thyra::Ordinal Ordinal;
  using Teuchos::outArg;

  TEUCHOS_ASSERT(this->opSupported(M_trans));

  // Get constants
  const Scalar zero = ST::zero();
  const Ordinal localDim = space_->localSubDim();
  const Ordinal procRank = comm_->getRank();
  const Ordinal numProcs = comm_->getSize();
  const Ordinal m = X_in.domain()->dim();
      
  // Loop over the input columns

  for (Ordinal col_j = 0; col_j < m; ++col_j) {
    
    // Get access the the elements of column j
    Thyra::ConstDetachedSpmdVectorView<Scalar> x_vec(X_in.col(col_j));
    Thyra::DetachedSpmdVectorView<Scalar> y_vec(Y_inout->col(col_j));
    const Teuchos::ArrayRCP<const Scalar> x = x_vec.sv().values();
    const Teuchos::ArrayRCP<Scalar> y = y_vec.sv().values();

    // Determine what process we are
    const bool first = (procRank == 0), last = ( procRank == numProcs-1 );
    
    // Communicate ghost elements
    Scalar x_km1, x_kp1;
    communicate( first, last, x(), outArg(x_km1), outArg(x_kp1) );

    // Perform operation (if beta==0 then we must be careful since y could
    // be uninitialized on input!)
    Thyra::Ordinal k = 0, lk = 0;
    if( beta == zero ) {
      // First local row
      y[k] = alpha * ( (first?zero:lower_[lk]*x_km1) + diag_[k]*x[k]
        + upper_[k]*x[k+1] );
      if(!first) ++lk;
      // Middle local rows
      for( k = 1; k < localDim - 1; ++lk, ++k )
        y[k] = alpha * ( lower_[lk]*x[k-1] + diag_[k]*x[k] + upper_[k]*x[k+1] );
      // Last local row
      y[k] = alpha * ( lower_[lk]*x[k-1] + diag_[k]*x[k]
        + (last?zero:upper_[k]*x_kp1) );
    }
    else {
      // First local row
      y[k] = alpha * ( (first?zero:lower_[lk]*x_km1) + diag_[k]*x[k]
        + upper_[k]*x[k+1] ) + beta*y[k];
      if(!first) ++lk;
      // Middle local rows
      for( k = 1; k < localDim - 1; ++lk, ++k )
        y[k] = alpha * ( lower_[lk]*x[k-1] + diag_[k]*x[k] + upper_[k]*x[k+1] )
          + beta*y[k];
      // Last local row
      y[k] = alpha * ( lower_[lk]*x[k-1] + diag_[k]*x[k]
        + (last?zero:upper_[k]*x_kp1) ) + beta*y[k];
    }
    
  }

}


template<class Scalar>
void ExampleTridiagSpmdLinearOp<Scalar>::communicate(
  const bool first, const bool last,
  const Teuchos::ArrayView<const Scalar> &x,
    const Teuchos::Ptr<Scalar> &x_km1,
    const Teuchos::Ptr<Scalar> &x_kp1
  ) const
{
  typedef Thyra::Ordinal Ordinal;
  const Ordinal localDim = space_->localSubDim();
  const Ordinal procRank = comm_->getRank();
  const Ordinal numProcs = comm_->getSize();
  const bool isEven = (procRank % 2 == 0);
  const bool isOdd = !isEven;
  // 1) Send x[localDim-1] from each even process forward to the next odd
  // process where it is received in x_km1.
  if(isEven && procRank+1 < numProcs) send(*comm_, x[localDim-1], procRank+1);
  if(isOdd && procRank-1 >= 0) receive(*comm_, procRank-1, &*x_km1);
  // 2) Send x[0] from each odd process backward to the previous even
  // process where it is received in x_kp1.
  if( isOdd && procRank-1 >= 0 )         send(*comm_, x[0], procRank-1);
  if( isEven && procRank+1 < numProcs ) receive(*comm_, procRank+1, &*x_kp1);
  // 3) Send x[localDim-1] from each odd process forward to the next even
  // process where it is received in x_km1.
  if (isOdd && procRank+1 < numProcs) send(*comm_, x[localDim-1], procRank+1);
  if (isEven && procRank-1 >= 0) receive(*comm_, procRank-1, &*x_km1);
  // 4) Send x[0] from each even process backward to the previous odd
  // process where it is received in x_kp1.
  if (isEven && procRank-1 >= 0) send(*comm_, x[0], procRank-1);
  if (isOdd && procRank+1 < numProcs) receive(*comm_, procRank+1, &*x_kp1);
}


#endif	// THYRA_EXAMPLE_TRIDIAG_SPMD_LINEAR_OP_HPP
