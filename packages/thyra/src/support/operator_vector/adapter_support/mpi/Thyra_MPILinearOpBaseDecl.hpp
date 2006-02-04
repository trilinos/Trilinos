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

#ifndef THYRA_MPI_LINEAR_OP_BASE_DECL_HPP
#define THYRA_MPI_LINEAR_OP_BASE_DECL_HPP

#include "Thyra_SingleScalarEuclideanLinearOpBaseDecl.hpp"
#include "RTOp_MPI_config.h"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace Thyra {

template<class Scalar> class MPIVectorSpaceBase;

/** \brief Base subclass for simplistic MPI SPMD linear operators.
 *
 * This subclass defines machinery for developing concrete
 * <tt>LinearOpBase</tt> subclasses for MPI SPMD vectors where it is
 * assumed that all of the local elements in associated vectors and
 * multi-vectors are immediately and cheaply available on each
 * processor.
 *
 * This base subclass derives from <tt>EuclideanLinearOpBase</tt> and
 * therefore any application-specific scalar products can easily be
 * incorporated.
 *
 * <b>Notes to subclass developers</b>
 *
 * The only function that a subclass must override in order to provide
 * a concrete implementation is the explicit single-vector version
 * \ref apply_expl_vec "euclideanApply()".
 *
 * This function is called on the subclass implementation passing in
 * views of explicit data.  The raw pointers to the local input and
 * input/output arrays are passed in simple templated classes
 * <tt>RTOpPack::SubVectorT</tt> and
 * <tt>RTOpPack::MutableSubVectorT</tt>.  Getting raw pointers out of
 * these objects is easy.
 *
 * It is very easy to create concrete subclasses of
 * <tt>MPILinearOpBase</tt> (see <tt>MPITridiagLinearOp</tt> for
 * a concrete example).  All one has to do is to create a derived base class
 * (<tt>MyMPILinearOp</tt> for example) of the form:
 *
 \code
template<class Scalar>
class MyMPILinearOp : public MPILinearOpBase<Scalar> {
private:
  // Declare your classes private data
  ...
public:
  // Declare you classes constructors, destructor and other initialization functions
  ...
protected:
  // Override the version of euclideanApply() that takes explicit data
  void euclideanApply(
    const ETransp                                M_trans
    ,const RTOpPack::SubVectorT<Scalar>          &local_x_in
    ,const RTOpPack::MutableSubVectorT<Scalar>   *local_y_out
    ,const Scalar                                alpha
    ,const Scalar                                beta
    ) const
    {
      // Get raw pointers to local vector data to make me feel better!
      const Scalar *local_x     = local_x_in.values();
      const Index  local_x_dim  = local_x_in.subDim();
      Scalar       *local_y     = local_y_out->values();
      const Index  local_y_dim  = local_y_out->subDim();
      // Perform operation the operation
      if( real_trans(M_trans) == ::Thyra::NOTRANS ) {
        // Perform the non-transposed operator: y = alpha*M*x + beta*y
        ...
      }
      else {
        // Perform the transposed operation: y = alpha*M'*x + beta*y
        ...
      }
  };
 \endcode
 *
 * Of course the above function will have to perform some type of
 * processor-to-processor communication in order to apply any
 * non-trivial distributed-memory linear operator but that is always
 * the case.
 *
 * If you do not need to handle arbitrary scalar data types then you
 * do not have to support them.  For example, to define a subclass
 * that only supports <b><tt>double</tt></b> you would declare a
 * non-templated version of the form:
 *
 \code
class MyMPILinearOp : public MPILinearOpBase<double> {
private:
  // Declare your classes private data
  ...
public:
  // Declare you classes constructors, destructor and other initialization functions
  ...
protected:
  // Override the version of euclideanApply() that takes explicit data
  void euclideanApply(
    const ETransp                                M_trans
    ,const RTOpPack::SubVectorT<double>          &local_x_in
    ,const RTOpPack::MutableSubVectorT<double>   *local_y_out
    ,const double                                alpha
    ,const double                                beta
    ) const
    {
      // Get raw pointers to vector data to make me feel better!
      const double *local_x     = local_x_in.values();
      const Index  local_x_dim  = local_x_in.subDim();
      double       *local_y     = local_y_out->values();
      const Index  local_y_dim  = local_y_out->subDim();
      // Perform operation the operation
      if( real_trans(M_trans) == ::Thyra::NOTRANS ) {
        // Perform the non-transposed operator: y = alpha*M*x + beta*y
        ...
      }
      else {
        // Perform the transposed operation: y = alpha*M'*x + beta*y
        ...
      }
  };
 \endcode
 *
 * By default, pointers to explicit data returned from
 * <tt>local_x.values()</tt> and <tt>local_y->values()</tt> above are
 * forced to have unit stride to simplify things.  However, if your
 * subclass can efficiently handle non-unit stride vector data (as the
 * BLAS can for example) then you can allow this by calling the
 * function <tt>this->forceUnitStride()</tt> and passing in
 * <tt>false</tt>.  The function <tt>this->forceUnitStride()</tt> can
 * only be called by your subclasses as it is declared
 * <tt>protected</tt> so do not worry about silly users messing with
 * this, it is none of their business.
 *
 * The explicit multi-vector version of \ref apply_expl_multi_vec "euclideanApply()"
 * has a default implementation that calls the explicit
 * single-vector version (that a subclass must supply) one column at a
 * time.  A subclass should only override this default multi-vector version
 * if it can do something more efficient.
 *
 * \ingroup Thyra_Op_Vec_adapters_MPI_support_grp
 */
template<class Scalar>
class MPILinearOpBase : virtual public SingleScalarEuclideanLinearOpBase<Scalar> {
public:

  /** \brief . */
  using SingleScalarEuclideanLinearOpBase<Scalar>::euclideanApply;

  /** @name Overridden from EuclideanLinearOpBase */
  //@{
  /** \brief . */
  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> > rangeScalarProdVecSpc() const;
  /** \brief . */
  Teuchos::RefCountPtr<const ScalarProdVectorSpaceBase<Scalar> > domainScalarProdVecSpc() const;
  /** \brief Calls protected <tt>euclideanApply()</tt> function.
   *
   * \anchor apply_multi_vec
   */
  void euclideanApply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;
  //@}

protected:

  /** @name Protected constructors/initializers/accessors */
  //@{

  /** \brief Set if unit stride is forced for vector data views or not
   *
   * @param forceUnitStride  [in]
   *
   * Postconditions:<ul>
   * <li><tt>this->forceUnitStride() == forceUnitStride</tt>
   * </ul>
   */
  STANDARD_MEMBER_COMPOSITION_MEMBERS( bool, forceUnitStride )

  /** Construct to uninitialized
   *
   * Postconditions:<ul>
   * <li><tt>this->domain().get() == NULL</tt>
   * <li><tt>this->range().get() == NULL</tt>
   * </ul>
   */
  MPILinearOpBase();

  /** \brief Initialize vector spaces using pre-formed <tt>MPIVectorSpaceBase</tt> objects.
   *
   * @param  domain   [in] Smart pointer to domain space
   * @param  range    [in] Smart pointer to range space
   *
   * Precondition:<ul>
   * <li><tt>domain.get() != NULL</tt>
   * <li><tt>range.get()  != NULL</tt>
   * </ul>
   *
   * Postcondition:<ul>
   * <li><tt>this->domain().get() == domain.get()</tt>
   * <li><tt>this->range().get()  == range.get()</tt>
   * </ul>
   */
  virtual void setSpaces(
    const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >      &range
    ,const Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >     &domain
    );

  /** \brief Initialize vector spaces given local dimensions (uses <tt>MPIVectorSpaceStd</tt>).
   *
   * @param  mpiComm         [in] MPI Communicator
   * @param  localDimRange   [in] The local number of vector elements in the domain space.
   * @param  localDimDomain  [in] The local number of vector elements in the domain space.
   *
   * Precondition:<ul>
   * <li><tt>localDimRange  > 0</tt>
   * <li><tt>localDimDomain > 0</tt>
   * </ul>
   *
   * Postcondition:<ul>
   * <li><tt>dynamic_cast<MPIVectorSpace<Scalar>&>(*this->range()).localSubDim() == localDimRange</tt>
   * <li><tt>dynamic_cast<MPIVectorSpace<Scalar>&>(*this->domain()).localSubDim() == localDimDomain</tt>
   * <li><tt>dynamic_cast<const MPIVectorSpaceStd<Scalar>*>(this->range().get())  != NULL</tt>
   * <li><tt>dynamic_cast<const MPIVectorSpaceStd<Scalar>*>(this->domain().get()) != NULL</tt>
   * </ul>
   */
  virtual void setLocalDimensions(
    MPI_Comm                                                    mpiComm 
    ,const Index                                                localDimRange
    ,const Index                                                localDimDomain
    );

  //@}

  /** @name Protected virtual functions to be overridden by subclasses */
  //@{

  /** \brief Apply the operator to explicit vector data.
   *
   * \anchor apply_expl_vec
   *
   * See <tt>LinearOpBase::euclideanApply()</tt> for a discussion of the
   * arguments to this function.  What differentiates this function is
   * that <tt>local_x</tt> and <tt>local_y</tt> are passed as objects
   * with explicit pointers to local vector data.
   *
   * Since this function is protected and does not get directly called by a client.
   * Instead, this function is called by the vector version of \ref apply_vec "euclideanApply()".
   */
  virtual void euclideanApply(
    const ETransp                                M_trans
    ,const RTOpPack::SubVectorT<Scalar>          &local_x
    ,const RTOpPack::MutableSubVectorT<Scalar>   *local_y
    ,const Scalar                                alpha
    ,const Scalar                                beta
    ) const = 0;

  /** \brief Apply the operator to explicit multi-vector data.
   *
   * \anchor apply_expl_multi_vec
   *
   * See <tt>LinearOpBase::euclideanApply()</tt> for a discussion of the
   * arguments to this function.  What differentiates this function is
   * that <tt>local_X</tt> and <tt>local_Y</tt> are passed as objects
   * with explicit pointers to local multi-vector data.
   *
   * Since this function is protected and does not get directly called by a client.
   * Instead, this function is called by the multi-vector version of \ref apply_multi_vec "euclideanApply()".
   *
   * The default implementation just calls the above vector version
   * one column at a time.  A subclass should only override this
   * function if it can provide a cache-smart version.  At any rate,
   * one can get up and going very quickly by just providing an
   * override for the simpler single-vector version.  Then latter, if
   * profiling data justifies it, one can provide a specialized
   * override for this function in an attempt to improve performance.
   */
  virtual void euclideanApply(
    const ETransp                                     M_trans
    ,const RTOpPack::SubMultiVectorT<Scalar>          &local_X
    ,const RTOpPack::MutableSubMultiVectorT<Scalar>   *local_Y
    ,const Scalar                                     alpha
    ,const Scalar                                     beta
    ) const;

  //@}

private:

  Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    range_;
  Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> >    domain_;

};	// end class LinearOpBase

}	// end namespace Thyra

#endif	// THYRA_MPI_LINEAR_OP_BASE_DECL_HPP
