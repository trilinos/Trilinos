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

#ifndef THYRA_SPMD_VECTOR_SPACE_BASE_DECL_HPP
#define THYRA_SPMD_VECTOR_SPACE_BASE_DECL_HPP

#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_ScalarProdVectorSpaceBase_decl.hpp"

namespace Thyra {

/** \brief Base <tt>%VectorSpaceBase</tt> class for all SPMD vector spaces
 * with contiguous local-to-global indexing.
 *
 * See <tt>SpmdVectorSpaceBase</tt> for details on what this class represents
 * in an abstract way.
 *
 * <b>Notes to subclass developers:</b>
 *
 * The pure virtual methods <tt>mpiComm()</tt>, <tt>localSubDim()</tt> and
 * <tt>createMember()</tt> are the only functions that must be overridden.
 *
 * If <tt>this</tt> this is in an uninitialized state then
 * <tt>localSubDim()</tt> should return <tt>0</tt>.
 *
 * It should never be necessary to override the virtual functions
 * <tt>mapCode()</tt> and <tt>isCompatible()</tt> as these functions
 * have very good and very general implementations.
 *
 * If optimized implementations of multi-vectors can be supported,
 * then the <tt>createMembers()</tt> method should also be overridden.
 *
 * This class defines a very general default implementation for
 * <tt>smallVecSpcFcty()</tt> that returns a
 * <tt>DefaultSpmdVectorSpaceFactory</tt> object which in turn can be used to
 * create <tt>DefaultSpmdVectorSpace</tt> objects.  These are the vector space
 * objects that are used by the domain space of multi-vectors created by
 * <tt>createMembers()</tt>.  The class <tt>DefaultSpmdVectorSpace</tt>
 * creates <tt>DefaultSpmdVector</tt> and <tt>DefaultSpmdMultiVector</tt>
 * objects.  This implementation of <tt>smallVecSpcFcty()</tt> is very general
 * should be very appropriate for many different concrete implementations.
 *
 * <b>Note:</b> It is very important that subclasses call the
 * <tt>updateState()</tt> function whenever the state of <tt>*this</tt>
 * changes in a way that might affect the behavior of any of the public member
 * functions.  For example, if a different value of <tt>localSubDim()</tt>
 * will be returned the next time it is called by a client, then
 * <tt>%updateState()</tt> needs to be called by the subclass.  External
 * clients should never need to worry about this function and that is why
 * <tt>%updateState()</tt> is declared protected.
 *
 * \ingroup Thyra_Op_Vec_adapters_Spmd_support_grp
 */
template<class Scalar>
class SpmdVectorSpaceDefaultBase
  : virtual public SpmdVectorSpaceBase<Scalar>
  , virtual public ScalarProdVectorSpaceBase<Scalar>
{
public:

  /** \brief . */
  SpmdVectorSpaceDefaultBase();

  /** @name Overridden from SpmdVectorSpaceBase */
  //@{

  /** \brief .
   *
   * This method has a default implementation which just assigns this offset
   * based on counting up <tt>localSubDim()</tt> on each process and then
   * setting <tt>localOffset()</tt> by the rank of the process.  For
   * example, if there are 5 elements in process 0 and 4 elements in process
   * rank, then <tt>localOffset</tt> on each of these processes will be set
   * as: <tt>localOffset=0</tt> on process 0, <tt>localOffset=5</tt> on
   * process 1, <tt>localOffset=9</tt> on process 2 and so on.
   */
  Ordinal localOffset() const;
  
  /** \brief .
   *
   * This method takes the data <tt>getComm()</tt>, <tt>numProc</tt> (where
   * <tt>numProc=this->getComm()->size()</tt>), <tt>localOffset()</tt> or
   * <tt>localSubDim()</tt> on each process and then uses it to compute a
   * value for <tt>mapCode</tt> (using a single global reduction if
   * <tt>numProc > 1</tt>) which is returned from this function.
   *
   * The value returned from this default implementation of this
   * method must not be changed or this approach breaks down.  The
   * only reason for overriding this method is for the subclass to
   * be alerted of <em>when</em> this method is called but not
   * <em>what</em> is returned from this method.  If a subclass
   * developer does not understand what this means then <b>don't</b>
   * override this method!
   *
   * The default implementation will always return <tt>return >
   * 0</tt> (unless <tt>this</tt> is uninitialized) so that if this
   * method is overridden to return <tt>return <= </tt> then this is
   * a flag that the underlying vector map does not satisfy the
   * assumptions of this vector space interface and vectors that are
   * in <tt>*this</tt> vector space can not collaborate with other
   * Spmd-based vector implementations.
   */
  Ordinal mapCode() const;

  //@}

  /** @name Overridden form Teuchos::Describable */
  //@{
  /** \brief . */
  std::string description() const;
  //@}

  /** @name Overridden from VectorSpaceBase */
  //@{

  /** \brief Returns the sum of the local number of elements on every process.
   */
  Ordinal dim() const;

  /** \brief Returns a <tt>DefaultSpmdVectorSpaceFactory</tt> object that has been given <tt>getComm()</tt>.
   */
  Teuchos::RCP< const VectorSpaceFactoryBase<Scalar> > smallVecSpcFcty() const;

  /** \brief Checks the general compatibility of parallel (or serial on one
   * process) Spmd-based vector spaces.
   *
   * @return Returns true if <tt>*this</tt> and <tt>vecSpace</tt>
   * are both serial in-core vectors or if <tt>vecSpc</tt> is of
   * type <tt>SpmdVectorSpaceDefaultBase<Scalar></tt> and both <tt>*this</tt>
   * and <tt>vecSpc</tt> have the same Spmd communicators and the same
   * mapping of elements to processes.
   *
   * Postconditions:<ul>
   * <li> <tt>return = ( this->hasInCoreView() &&
   *      vecSpc->hasInCoreView() ) || ( (mpiVecSpc = dynamic_cast<const
   *      SpmdVectorSpaceDefaultBase<Scalar>*>(&vecSpc)) && this->getComm() ==
   *      mpiVecSpc->getComm() && this->mapCode() ==
   *      mpiVecSpc->getComm())</tt>.
   *
   * </ul>
   *
   * If the mapping of vector elements to processes is not as
   * described
   * <A HREF="classThyra_1_1SpmdVectorSpaceBase.html#SpmdVectorSpaceBase_Vector_layout">above</A>
   * then this method should be overridden in a way that is specific
   * to the vector implementation.
   */
   bool isCompatible(const VectorSpaceBase<Scalar>& vecSpc) const;
  
  //@}

protected:

  /** \brief This function must be called whenever the state of
   * <tt>this</tt> changes and some internal state must be updated.
   *
   * @param  globalDim
   *             [in] If <tt>globalDim > 0</tt> then this determines
   *             the global dimension of the vector space.  If <tt>globalDim==this->localSubDim()</tt>
   *             then this is a locally replicated vector space.  If <tt>globalDim < 0</tt> then
   *             the global dimension is computed using a global reduction.
   *             If <tt>this->getComm()->size()==1</tt> then this argument is ignored.
   *
   * Note that calling this function may involve one or more global reductions
   * being called if this is parallel vector space so it should only be called
   * when needed by subclasses.
   *
   * Usually, this operation only needs to be called once for every
   * *new* parallel vector space constructed and very few parallel
   * vector spaces will be created per application usually.
   */
  virtual void updateState( const Ordinal globalDim );

private:

  // //////////////////////////////////////
  // Private data members

  Ordinal     mapCode_;    // < 0 is a flag that everything needs initialized
  Ordinal     defaultLocalOffset_;
  Ordinal     defaultGlobalDim_;
  Ordinal     localSubDim_;

  Teuchos::RCP< const VectorSpaceFactoryBase<Scalar> >  smallVecSpcFcty_;
  
}; // end class SpmdVectorSpaceDefaultBase

} // end namespace Thyra

#endif // THYRA_SPMD_VECTOR_SPACE_BASE_DECL_HPP
