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

#ifndef THYRA_MPI_VECTOR_SPACE_STD_DECL_HPP
#define THYRA_MPI_VECTOR_SPACE_STD_DECL_HPP

#include "Thyra_MPIVectorSpaceBaseDecl.hpp"

namespace Thyra {

/** \brief Concrete implementation of a MPI-based SPMD vector space subclass which creates
 * <tt>MPIMultiVectorStd</tt> objects.
 *
 * This is a simple but yet fully general and functional concrete
 * subclass of <tt>MPIVectorSpace</tt> that creates <tt>MPIMultiVectorStd</tt>
 * objects from <tt>createMembers()</tt> and <tt>VectorMultiVector</tt>
 * wrapped <tt>MPIMulitVectorStd</tt> objects from <tt>createMember()</tt>.
 *
 * See the function <tt>initialize()</tt> that describes the different
 * kinds of distributions this class can handle.
 *
 * \ingroup Thyra_Op_Vec_adapters_MPI_concrete_std_grp
 */
template<class Scalar>
class MPIVectorSpaceStd : public MPIVectorSpaceBase<Scalar> {
public:

  /** @name Constructors and initializers */
  //@{
 
  /// Construct to uninitialized (see postconditions from <tt>uninitialize()</tt>)
  MPIVectorSpaceStd();

  /// Calls <tt>initialize()</tt>.
  MPIVectorSpaceStd( MPI_Comm mpiComm, const Index localSubDim, const Index globalDim );

  /** \brief Initialize.
   *
   * @param  mpiComm      [in] The MPI communicator.  This object must be maintained
   *                      by the client the entire time that <tt>this</tt> is in use.
   * @param  localSubDim  [in] The number of elements on local processor.  This number
   *                      can be different on every processor.
   * @param  globalDim    [in] Gives the number of global elements in the vector
   *                      if <tt>globalDim > 0</tt>.  If <tt>globalDim < 0</tt>
   *                      then the global dimension is determined by the above argument
   *                      <tt>localSubDim</tt> but requires a global communication to
   *                      do so (i.e. <tt>MPI_Allreduce()</tt>).
   *
   * Preconditions:<ul>
   * <li><tt>localSubDim > 0</tt>
   * <li><tt>globalDim != 0</tt>
   * <li>[<tt>mpiComm != MPI_COMM_NULL && globalDim > 0<tt>]  <tt>globalDim >= localSubDim</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li><tt>this->mpiComm() == mpiComm</tt>
   * <li><tt>this->localSubDim() == localSubDim</tt>
   * <li>[<tt>mpiComm == MPI_COMM_NULL</tt>] <tt>this->dim() == localSubDim</tt>
   * <li>[<tt>mpiComm != MPI_COMM_NULL && globalDim > 0</tt>] <tt>this->dim() == globalDim</tt>
   * <li>[<tt>mpiComm != MPI_COMM_NULL && globalDim < 0</tt>] <tt>this->dim() == sum(localSubDim[i],i=0...numProc-1)</tt>
   * </ul>
   *
   * This function supports three different types of use-cases:
   * <ul>
   * <li><tt>mpiComm==MPI_COMM_NULL</tt> : Serial (i.e. single process) vectors
   *     where <tt>this->dim() == localSubDim</tt>.
   * <li><tt>mpiComm!=MPI_COMM_NULL && globalDim < 0</tt> : Distributed-memory vectors
   *     where <tt>this->dim()</tt> is equal to the sum of the <tt>localSubDim</tt>
   *     arguments on each processor.  This will result in a call to <tt>MPI_Allreduce()</tt>
   *     inside of this function.
   * <li><tt>mpiComm!=MPI_COMM_NULL && globalDim > 0</tt> : Distributed-memory vectors
   *     where <tt>this->dim()</tt> returns <tt>globalDim</tt>.  This will not result
   *     in a call to <tt>MPI_Allreduce()</tt> inside this function and therefore the client had better
   *     be sure that <tt>globalDim</tt> is consistent with <tt>localSubDim</tt>
   *     on each processor.
   * <li><tt>mpiComm!=MPI_COMM_NULL && globalDim == localSubDim</tt> : Locally-replicated
   *     distributed-memory vectors where <tt>this->dim() == globalDim == localSubDim</tt>.
   * </ul>
   */
  void initialize( MPI_Comm mpiComm, const Index localSubDim, const Index globalDim );

  /** \brief Set to an uninitialized state.
   *
   * Postconditions:<ul>
   * <li>this->mpiComm() == MPI_COMM_NULL</tt>
   * <li>this->localSubDim() == 0</tt>
   * <li>this->dim() == 0</tt>
   * </ul>
   */
  void uninitialize( MPI_Comm *mpiComm = NULL, Index *localSubDim = NULL, Index *globalDim=NULL );

  //@}

  /** @name Overridden form Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  //@}

  /** @name Public overridden from VectorSpace */
  //@{
  /** \brief . */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > clone() const;
  //@}

protected:

  /** @name Protected overridden from VectorSpace */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<VectorBase<Scalar> > createMember() const;
  /** \brief . */
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > createMembers(int numMembers) const;
  /** \brief . */
  Teuchos::RefCountPtr<VectorBase<Scalar> > createMemberView( const RTOpPack::MutableSubVectorT<Scalar> &raw_v ) const;
  /** \brief . */
  Teuchos::RefCountPtr<const VectorBase<Scalar> > createMemberView( const RTOpPack::SubVectorT<Scalar> &raw_v ) const;
  /** \brief . */
  Teuchos::RefCountPtr<MultiVectorBase<Scalar> > createMembersView( const RTOpPack::MutableSubMultiVectorT<Scalar> &raw_mv ) const;
  /** \brief . */
  Teuchos::RefCountPtr<const MultiVectorBase<Scalar> > createMembersView( const RTOpPack::SubMultiVectorT<Scalar> &raw_mv ) const;

  //@}

public:

  /** @name Public overridden from MPIVectorSpaceBase */
  //@{

  /** \brief . */
  MPI_Comm mpiComm() const;
  /** \brief . */
   Index localSubDim() const;

  //@}

private:

  // //////////////////////////////////////
  // Private data members

  MPI_Comm           mpiComm_;
  Index              localSubDim_;
  int                numProc_;
  int                procRank_;
  
}; // end class MPIVectorSpaceStd

} // end namespace Thyra

#endif // THYRA_MPI_VECTOR_SPACE_STD_DECL_HPP
