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

#include "Thyra_SpmdLinearOpBase.hpp"
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
 * Note, this is just an example class and for the sake of simplified
 * presentation the private members are listed first and in-class
 * declarations are used which are not a good idea in production code.
 * However, in this case, they make the example code easier to read
 * and maintaining encapsulation and a well defined interface are
 * unnecessary here.
 *
 * See the source code for this simple example by clicking on the
 * link to the definition below.
 *
 * \ingroup Thyra_Op_Vec_examples_cg_Spmd_grp
 */
template<class Scalar>
class ExampleTridiagSpmdLinearOp : public Thyra::SpmdLinearOpBase<Scalar> {
private:

  Teuchos::RCP<const Teuchos::Comm<Thyra::Index> >  comm_;
  int                  procRank_;
  int                  numProcs_;
  Thyra::Index         localDim_;
  std::vector<Scalar>  lower_;   // size = ( procRank == 0         ? localDim - 1 : localDim )    
  std::vector<Scalar>  diag_;    // size = localDim
  std::vector<Scalar>  upper_;   // size = ( procRank == numProc-1 ? localDim - 1 : localDim )

  void communicate( const bool first, const bool last, const Scalar x[], Scalar *x_km1, Scalar *x_kp1 ) const;

public:

  /** \brief . */
  using Thyra::SpmdLinearOpBase<Scalar>::euclideanApply;

  /// Construct to uninitialized
  ExampleTridiagSpmdLinearOp() : procRank_(0), numProcs_(0) {}

  /// Calls <tt>initialize()</tt>.
  ExampleTridiagSpmdLinearOp(
    const Teuchos::RCP<const Teuchos::Comm<Thyra::Index> > &comm
    ,const Thyra::Index localDim, const Scalar lower[], const Scalar diag[], const Scalar upper[] )
    { this->initialize(comm,localDim,lower,diag,upper);	}
  
  /** Initialize given lower, diagonal and upper arrays of data.
   *
   * @param  comm      [in] Communicator (allowed to be Teuchos::null)
   * @param  localDim  [in] Dimension of this matrix (must be >= 2).
   * @param  lower     [in] Array (length <tt>( procRank == 0 ? localDim - 1 : localDim )</tt>)
   *                   of the lower diagonal elements
   * @param  diag      [in] Array (length <tt>localDim</tt>) of the central diagonal elements
   * @param  upper     [in] Array (length <tt>( procRank == numProc-1 ? localDim - 1 : localDim )</tt>)
   *                   of the upper diagonal elements
   *
   * Preconditions:<ul>
   * <li><tt>localDim >= 2</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li>Should be obvious!
   * <li>See <tt>setLocalDimensions()</tt>
   * </ul>
   */
  void initialize(
    const Teuchos::RCP<const Teuchos::Comm<Thyra::Index> > &comm
    ,const Thyra::Index             localDim // >= 2
    ,const Scalar                   lower[]  // size == ( procRank == 0         ? localDim - 1 : localDim )
    ,const Scalar                   diag[]   // size == localDim
    ,const Scalar                   upper[]  // size == ( procRank == numProc-1 ? localDim - 1 : localDim )
    )
    {
      TEST_FOR_EXCEPT( localDim < 2 );
      this->setLocalDimensions(comm,localDim,localDim); // Needed to set up range() and domain()
      comm_  = Teuchos::DefaultComm<Thyra::Index>::getDefaultSerialComm(comm);
      localDim_ = localDim;
      numProcs_ = comm->getSize();
      procRank_ = comm->getRank();
      const Thyra::Index
        lowerDim = ( procRank_ == 0           ? localDim - 1 : localDim ),
        upperDim = ( procRank_ == numProcs_-1 ? localDim - 1 : localDim );
      lower_.resize(lowerDim);  for( int k = 0; k < lowerDim; ++k ) lower_[k] = lower[k];
      diag_.resize(localDim);   for( int k = 0; k < localDim; ++k ) diag_[k]  = diag[k];
      upper_.resize(upperDim);  for( int k = 0; k < upperDim; ++k ) upper_[k] = upper[k];
    }

  // Overridden form Teuchos::Describable */

  std::string description() const
    {
      return (std::string("ExampleTridiagSpmdLinearOp<") + Teuchos::ScalarTraits<Scalar>::name() + std::string(">"));
    }

protected:


  // Overridden from SingleScalarEuclideanLinearOpBase

  bool opSupported( Thyra::EOpTransp M_trans ) const
    {
      typedef Teuchos::ScalarTraits<Scalar> ST;
      return (M_trans == Thyra::NOTRANS || (!ST::isComplex && M_trans == Thyra::CONJ) );
    }

  // Overridden from SerialLinearOpBase

  void euclideanApply(
    const Thyra::EOpTransp                         M_trans
    ,const RTOpPack::ConstSubVectorView<Scalar>  &local_x_in
    ,const RTOpPack::SubVectorView<Scalar>       *local_y_out
    ,const Scalar                                alpha
    ,const Scalar                                beta
    ) const
    {
      typedef Teuchos::ScalarTraits<Scalar> ST;
      TEST_FOR_EXCEPTION( M_trans != Thyra::NOTRANS, std::logic_error, "Error, can not handle transpose!" );
      // Get constants
      const Scalar zero = ST::zero();
      // Get raw pointers to vector data to make me feel better!
      const Scalar *x = local_x_in.values().get();
      Scalar       *y = local_y_out->values().get();
      // Determine what process we are
      const bool first = ( procRank_ == 0 ), last = ( procRank_ == numProcs_-1 );
      // Communicate ghost elements
      Scalar x_km1, x_kp1;
      communicate( first, last, x, &x_km1, &x_kp1 );
      // Perform operation (if beta==0 then we must be careful since y could be uninitialized on input!)
      Thyra::Index k = 0, lk = 0;
      if( beta == zero ) {
        y[k] = alpha * ( (first?zero:lower_[lk]*x_km1) + diag_[k]*x[k] + upper_[k]*x[k+1] ); if(!first) ++lk;             // First local row
        for( k = 1; k < localDim_ - 1; ++lk, ++k )
          y[k] = alpha * ( lower_[lk]*x[k-1] + diag_[k]*x[k] + upper_[k]*x[k+1] );                                        // Middle local rows
        y[k] = alpha * ( lower_[lk]*x[k-1] + diag_[k]*x[k] + (last?zero:upper_[k]*x_kp1) );                               // Last local row
      }
      else {
        y[k] = alpha * ( (first?zero:lower_[lk]*x_km1) + diag_[k]*x[k] + upper_[k]*x[k+1] ) + beta*y[k]; if(!first) ++lk; // First local row
        for( k = 1; k < localDim_ - 1; ++lk, ++k )
          y[k] = alpha * ( lower_[lk]*x[k-1] + diag_[k]*x[k] + upper_[k]*x[k+1] ) + beta*y[k];                            // Middle local rows
        y[k] = alpha * ( lower_[lk]*x[k-1] + diag_[k]*x[k] + (last?zero:upper_[k]*x_kp1) ) + beta*y[k];                   // Last local row
      }
      //std::cout << "\ny = ["; for(k=0;k<localDim_;++k) { std::cout << y[k]; if(k<localDim_-1) std::cout << ","; } std::cout << "]\n";
    }

};	// end class ExampleTridiagSpmdLinearOp

// private

template<class Scalar>
void ExampleTridiagSpmdLinearOp<Scalar>::communicate(
  const bool first, const bool last, const Scalar x[], Scalar *x_km1, Scalar *x_kp1
  ) const
{
  const bool isEven = (procRank_ % 2 == 0);
  const bool isOdd = !isEven;
  // 1) Send x[localDim_-1] from each even process forward to the next odd
  // process where it is received in x_km1.
  if( isEven && procRank_+1 < numProcs_ ) send(*comm_,x[localDim_-1],procRank_+1);
  if( isOdd && procRank_-1 >= 0 )         receive(*comm_,procRank_-1,x_km1);
  // 2) Send x[0] from each odd process backward to the previous even
  // process where it is received in x_kp1.
  if( isOdd && procRank_-1 >= 0 )         send(*comm_,x[0],procRank_-1);
  if( isEven && procRank_+1 < numProcs_ ) receive(*comm_,procRank_+1,x_kp1);
  // 3) Send x[localDim_-1] from each odd process forward to the next even
  // process where it is received in x_km1.
  if( isOdd && procRank_+1 < numProcs_ )  send(*comm_,x[localDim_-1],procRank_+1);
  if( isEven && procRank_-1 >= 0 )        receive(*comm_,procRank_-1,x_km1);
  // 4) Send x[0] from each even process backward to the previous odd
  // process where it is received in x_kp1.
  if(isEven && procRank_-1 >= 0 )         send(*comm_,x[0],procRank_-1);
  if(isOdd && procRank_+1 < numProcs_ )   receive(*comm_,procRank_+1,x_kp1);
}

#endif	// THYRA_EXAMPLE_TRIDIAG_SPMD_LINEAR_OP_HPP
