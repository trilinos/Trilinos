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


#ifndef THYRA_DEFAULT_SPMD_VECTOR_SPACE_DECL_HPP
#define THYRA_DEFAULT_SPMD_VECTOR_SPACE_DECL_HPP


#include "Thyra_SpmdVectorSpaceDefaultBase_decl.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Teuchos_Handleable.hpp"


namespace Thyra {


/** \brief Concrete implementation of an SPMD vector space subclass which
 * creates <tt>DefaultSpmdVector</tt> and <tt>DefaultSpmdMultiVector</tt>
 * objects.
 *
 * This is a simple but yet fully general and functional concrete subclass of
 * <tt>SpmdVectorSpaceBase</tt> that returns <tt>DefaultSpmdMultiVector</tt>
 * objects from <tt>createMembers()</tt> and <tt>DefaultSpmdVector</tt>
 * objects from <tt>createMember()</tt>.
 *
 * See the function <tt>initialize()</tt> that describes the different kinds
 * of distributions this class can handle.
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_concrete_std_grp
 */
template<class Scalar>
class DefaultSpmdVectorSpace
  : public SpmdVectorSpaceDefaultBase<Scalar>,
    public Teuchos::Handleable<VectorSpaceBase<Scalar> >
{
public:
  /* handleable interface */
  TEUCHOS_GET_RCP(VectorSpaceBase<Scalar>);

  /** @name Constructors and initializers */
  //@{

  /** \brief Create with weak ownership to self. */
  static RCP<DefaultSpmdVectorSpace<Scalar> > create();

  /** \brief Initialize a serial space.
   *
   * \param dim
   * [in] Gives the dimension of the vector space.
   *
   * Equivalent to calling <tt>this->initialize(Teuchos::null,dim,dim)</tt>
   */
  void initialize(
    const Ordinal dim
    );

  /** \brief Initialize an SPMD space.
   *
   * \param comm
   * [in] The communicator. This object must be maintained
   * by the client the entire time that <tt>this</tt> is in use.
   *
   * \param localSubDim
   * [in] The number of elements in the local process. This number
   * can be different in every process.
   *
   * \param globalDim
   * [in] Gives the number of global elements in the vector
   * if <tt>globalDim > 0</tt>. If <tt>globalDim < 0</tt>
   * then the global dimension is determined by the above argument
   * <tt>localSubDim</tt> but requires a global communication to
   * do so (i.e. <tt>Spmd_Allreduce()</tt>).
   *
   * Preconditions:<ul>
   * <li><tt>localSubDim > 0</tt>
   * <li><tt>globalDim != 0</tt>
   * <li>[<tt>comm.get() != NULL && globalDim > 0<tt>] <tt>globalDim >= localSubDim</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>this->getComm().get() == comm.get()</tt>
   * <li><tt>this->localSubDim() == localSubDim</tt>
   * <li>[<tt>comm.get() == NULL</tt>] <tt>this->dim() == localSubDim</tt>
   * <li>[<tt>comm.get() != NULL && globalDim > 0</tt>] <tt>this->dim() == globalDim</tt>
   * <li>[<tt>comm.get() != NULL && globalDim < 0</tt>] <tt>this->dim() == sum(localSubDim[i],i=0...numProc-1)</tt>
   * </ul>
   *
   * This function supports three different types of use-cases:
   *
   * <ul>
   *
   * <li><tt>comm.get()==NULL</tt> : Serial (i.e. single process) vectors
   * where <tt>this->dim() == localSubDim</tt>.
   *
   * <li><tt>comm.get()!=NULL && globalDim < 0</tt> : Distributed-memory vectors
   * where <tt>this->dim()</tt> is equal to the sum of the <tt>localSubDim</tt>
   * arguments in each process. This will result in a call to <tt>Spmd_Allreduce()</tt>
   * inside of this function.
   *
   * <li><tt>comm.get()!=NULL && globalDim > 0</tt> : Distributed-memory vectors
   * where <tt>this->dim()</tt> returns <tt>globalDim</tt>. This will not result
   * in a call to <tt>Teuchos::reduceAll()</tt> inside this function and therefore the client had better
   * be sure that <tt>globalDim</tt> is consistent with <tt>localSubDim</tt>
   * in each process.
   *
   * <li><tt>comm.get()!=NULL && globalDim == localSubDim</tt> : Locally-replicated
   * distributed-memory vectors where <tt>this->dim() == globalDim == localSubDim</tt>.
   *
   * </ul>
   */
  void initialize(
    const RCP<const Teuchos::Comm<Ordinal> > &comm,
    const Ordinal localSubDim, const Ordinal globalDim
    );

  /** \brief Set to an uninitialized state.
   *
   * Postconditions:<ul>
   * <li>this->getComm().get() == NULL</tt>
   * <li>this->localSubDim() == 0</tt>
   * <li>this->dim() == 0</tt>
   * </ul>
   */
  void uninitialize();

  //@}

  /** @name Public overridden from VectorSpaceBase */
  //@{
  /** \brief Returns true if all the elements in <tt>rng</tt> are in this
   * process.
   */
  bool hasInCoreView(
    const Range1D& rng, const EViewType viewType, const EStrideType strideType
    ) const;
  /** \brief . */
  RCP< const VectorSpaceBase<Scalar> > clone() const;
  //@}

protected:

  /** @name Protected overridden from VectorSpaceBase */
  //@{

  /** \brief . */
  RCP<VectorBase<Scalar> >
  createMember() const;
  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  createMembers(int numMembers) const;
  /** \brief . */
  RCP<VectorBase<Scalar> >
  createMemberView( const RTOpPack::SubVectorView<Scalar> &raw_v ) const;
  /** \brief . */
  RCP<const VectorBase<Scalar> >
  createMemberView( const RTOpPack::ConstSubVectorView<Scalar> &raw_v ) const;
  /** \brief . */
  RCP<MultiVectorBase<Scalar> >
  createMembersView( const RTOpPack::SubMultiVectorView<Scalar> &raw_mv ) const;
  /** \brief . */
  RCP<const MultiVectorBase<Scalar> >
  createMembersView( const RTOpPack::ConstSubMultiVectorView<Scalar> &raw_mv ) const;

  //@}

public:

  /** @name Public overridden from SpmdVectorSpaceDefaultBase */
  //@{

  /** \brief . */
  RCP<const Teuchos::Comm<Ordinal> > getComm() const;
  /** \brief . */
  Ordinal localSubDim() const;

  //@}

private:

  // //////////////////////////////////////
  // Private data members

  RCP<const Teuchos::Comm<Ordinal> > comm_;
  Ordinal localSubDim_;
  int numProc_;
  int procRank_;
  RCP<DefaultSpmdVectorSpace<Scalar> > weakSelfPtr_;

  // /////////////////////////////////////
  // Private member functions
 
  DefaultSpmdVectorSpace();

public:

  /** \name Deprecated */
  //@{

  /** \brief Deprecated. */
  DefaultSpmdVectorSpace(
    const Ordinal dim
    );

  /** \brief Deprecated. */
  DefaultSpmdVectorSpace(
    const RCP<const Teuchos::Comm<Ordinal> > &comm,
    const Ordinal localSubDim, const Ordinal globalDim
    );

  //@}
 
}; // end class DefaultSpmdVectorSpace


/** \brief Nonmember consturctor that creats an uninitialized vector space.
 *
 * \relates DefaultSpmdVectorSpace
 */
template<class Scalar>
RCP<DefaultSpmdVectorSpace<Scalar> >
defaultSpmdVectorSpace()
{
  return DefaultSpmdVectorSpace<Scalar>::create();
}


/** \brief Nonmember consturctor that creats a serial vector space.
 *
 * \relates DefaultSpmdVectorSpace
 */
template<class Scalar>
RCP<DefaultSpmdVectorSpace<Scalar> >
defaultSpmdVectorSpace( const Ordinal dim )
{
  RCP<DefaultSpmdVectorSpace<Scalar> > vs =
    DefaultSpmdVectorSpace<Scalar>::create();
  vs->initialize(dim);
  return vs;
}


/** \brief Nonmember consturctor function that creates a distributed or
 * locally-replicated parallel vector space.
 *
 * \relates DefaultSpmdVectorSpace
 */
template<class Scalar>
RCP<DefaultSpmdVectorSpace<Scalar> >
defaultSpmdVectorSpace(
  const RCP<const Teuchos::Comm<Ordinal> > &comm,
  const Ordinal localSubDim, const Ordinal globalDim
  )
{
  RCP<DefaultSpmdVectorSpace<Scalar> > vs =
    DefaultSpmdVectorSpace<Scalar>::create();
  vs->initialize(comm, localSubDim, globalDim);
  return vs;
}


} // end namespace Thyra


#endif // THYRA_DEFAULT_SPMD_VECTOR_SPACE_DECL_HPP
