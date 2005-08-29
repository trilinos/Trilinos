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

#ifndef THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DECL_HPP
#define THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DECL_HPP

#include "Thyra_LinearOpWithSolveBaseDecl.hpp"

namespace Thyra {

/** \brief Factory interface for creating <tt>LinearOpWithSolveBase</tt>
 * objects from <tt>LinearOpBase</tt> objects.
 *
 * This strategy interface allows a client to take one or more "compatible"
 * <tt>LinearOpBase</tt> objects and then create one or more
 * <tt>LinearOpWithSolveBase</tt> objects that can then be used to solve for
 * linear systems.  This interface carefully separates the construction from
 * the initialization of a <tt>LinearOpWithSolveBase</tt> object.
 *
 * Note that the non-member functions defined
 * \ref Thyra_LinearOpWithSolveFactoryBase_helper_grp "here" provide for
 * simpler use cases and are recommended as a way to access the capabilities
 * of this interface.
 *
 * This interface can be implemented by both direct and iterative linear
 * solvers.
 *
 * \ingroup Thyra_Op_Vec_Interoperability_Extended_Interfaces_grp
 */
template <class RangeScalar, class DomainScalar = RangeScalar>
class LinearOpWithSolveFactoryBase : virtual public Teuchos::Describable {
public:

  /** \brief . */
  virtual ~LinearOpWithSolveFactoryBase() {}

  /** @name Pure virtual functions that must be overridden in subclasses */
  //@{

  /** \brief Check that a <tt>LinearOpBase</tt> objects is compatible with
   * <tt>*this</tt> factory object.
   */
  virtual bool isCompatible( const LinearOpBase<RangeScalar,DomainScalar> &fwdOp ) const = 0;

  /** \brief Create an (uninitialized) <tt>LinearOpWithSolveBase</tt> object
   * to be initialized later in <tt>this->initializeOp()</tt>.
   *
   * Note that on output <tt>return->domain().get()==NULL</tt> may be true
   * which means that the operator is not fully initialized.  In fact, the
   * output operator object is not guaranteed to be fully initialized until
   * after it is passed through <tt>this->initializeOp()</tt>.
   */
  virtual Teuchos::RefCountPtr<LinearOpWithSolveBase<RangeScalar,DomainScalar> > createOp() const = 0;

  /** \brief Initialize a pre-created <tt>LinearOpWithSolveBase</tt> object
   * given a "compatible" <tt>LinearOpBase</tt> object.
   *
   * \param  fwdOp  [in] The forward linear operator that will be used to create
   *                the output <tt>LinearOpWithSolveBase</tt> object.
   * \param  Op     [out] The output <tt>LinearOpWithSolveBase</tt> object.  This object must have
   *                be created first by <tt>this->createOp()</tt>.  The object may have also
   *                already been passed through this function several times.  Note that subclasses
   *                should always first strip off the transpose and scaling by calling <tt>unwrap()</tt>
   *                before attempting to dynamic cast the object.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->isCompatible(*fwdOp)==true</tt>
   * <li><tt>Op!=NULL</tt>
   * <li><tt>*Op</tt> must have been created by <tt>this->createOp()</tt> prior to calling
   *     this function.
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>Throws <tt>CatastrophicSolveFailure</tt> if the underlying linear solver could
   *     not be created sucessfully (do to a factorization failure or some other cause).
   * <li><tt>Op->range()->isCompatible(*fwdOp->range())==true</tt>
   * <li><tt>Op->domain()->isCompatible(*fwdOp->domain())==true</tt>
   * <li><tt>Op->apply()</tt> and <tt>Op->applyTranspose()</tt> must behave
   *     exactly the same as <tt>fwdOp->apply()</tt> and <tt>fwdOp->applyTranspose()</tt>
   * <li><tt>fwdOp.count()</tt> after output is greater than <tt>fwdOp.count()</tt>
   *     just before this call and therefore the client can assume that the <tt>*fwdOp</tt> object will 
   *     be remembered by the <tt>*Op</tt> object.  The client must be careful
   *     not to modify the <tt>*fwdOp</tt> object or else the <tt>*Op</tt> object may also
   *     be modified.
   * </ul>
   */
  virtual void initializeOp(
    const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &fwdOp
    ,LinearOpWithSolveBase<RangeScalar,DomainScalar>                             *Op
    ) const = 0;

  /** \brief Uninitialize a <tt>LinearOpWithSolveBase</tt> object and return its
   * remembered forward linear operator.
   *
   * \param  Op     [in/out] On input, <tt>*Op</tt> is an initialized or uninitialized
   *                object and on output is uninitialized.  Note that "uninitialized"
   *                does not mean that <tt>Op</tt> is completely stateless.  It may still
   *                remember some aspect of the matrix <tt>fwdOp</tt> that will allow
   *                for a more efficient initialization next time through
   *                <tt>this->initializeOp()</tt>.
   * \param  fwdOp  [in/out] If <tt>fwdOp!=NULL</tt> on input, the on output, the
   *                same forward operator passed into <tt>this->initailzeOp()</tt> will be
   *                returned.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>*Op</tt> must have been created by <tt>this->createOp()</tt> prior to calling
   *     this function.
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>If <tt>*Op</tt> on input was initialized through a call to <tt>this->initializeOp()</tt>
   *     then <tt>return.get()!=NULL</tt>.
   * <li>If <tt>*Op</tt> was uninitialized on input and <tt>fwdOp!=NULL</ttt> then <tt>fwdOp->get()==NULL</tt>.
   * <li>On output, <tt>*Op</tt> can be considered to be uninitialized and
   *     it is safe to modify the forward operator object <tt>*(*fwdOp)</tt> returned in <tt>fwdOp</tt>.
   *     The default is <tt>fwdOp==NULL</tt> in which case the forward operator will not be returned in <tt>*fwdOp</tt>.
   * </ul>
   */
  virtual void uninitializeOp(
    LinearOpWithSolveBase<RangeScalar,DomainScalar>                       *Op
    ,Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >  *fwdOp = NULL
    ) const = 0;
  
  //@}


  /** @name Virtual functions with default implementations */
  //@{


  /** \brief Initialize a pre-created <tt>LinearOpWithSolveBase</tt> object
   * given a "compatible" <tt>LinearOpBase</tt> object but allow for reuse of
   * any preprocessing that is in <tt>*Op</tt>..
   *
   * \param  fwdOp  [in] The forward linear operator that will be used to create
   *                the output <tt>LinearOpWithSolveBase</tt> object.
   * \param  Op     [in/out] The output <tt>LinearOpWithSolveBase</tt> object.  This object must have
   *                be created first by <tt>this->createOp()</tt> and may have already been through
   *                at least one previous set of calls to <tt>this->initializeOp()</tt> and
   *                <tt>this->uninitializeOp()</tt>.  Note that subclasses
   *                should always first strip off the transpose and scaling by calling <tt>unwrap()</tt>
   *                before attempting to dynamic cast the object.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>this->isCompatible(*fwdOp)==true</tt>
   * <li><tt>Op!=NULL</tt>
   * <li><tt>*Op</tt> must have been created by <tt>this->createOp()</tt> prior to calling
   *     this function.
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>Throws <tt>CatastrophicSolveFailure</tt> if the underlying linear solver could
   *     not be created sucessfully (do to a factorization failure or some other cause).
   * <li><tt>Op->range()->isCompatible(*fwdOp->range())==true</tt>
   * <li><tt>Op->domain()->isCompatible(*fwdOp->domain())==true</tt>
   * <li><tt>Op->apply()</tt> and <tt>Op->applyTranspose()</tt> must behave
   *     exactly the same as <tt>fwdOp->apply()</tt> and <tt>fwdOp->applyTranspose()</tt>
   * <li><tt>fwdOp.count()</tt> after output is greater than <tt>fwdOp.count()</tt>
   *     just before this call and therefore the client can assume that the <tt>*fwdOp</tt> object will 
   *     be remembered by the <tt>*Op</tt> object.  The client must be careful
   *     not to modify the <tt>*fwdOp</tt> object or else the <tt>*Op</tt> object may also
   *     be modified.
   * </ul>
   *
   * The purpose of this function is to allow the reuse of old factorizations
   * and/or preconditioners that may go into the initialization of the
   * <tt>*Op</tt> objects.  Note that by calling this function, the peformance
   * <tt>Op->solve(...)</tt> may not be as good as when calling the function
   * <tt>this->initializeOp(...,Op)</tt>.
   *
   * The default implemenation of this function just calls
   * <tt>this->initializeOp(fwdOp,Op)</tt>.
   */
  virtual void initializeAndReuseOp(
    const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &fwdOp
    ,LinearOpWithSolveBase<RangeScalar,DomainScalar>                             *Op
    ) const;

  //@}

};

//@}

/** \defgroup Thyra_LinearOpWithSolveFactoryBase_helper_grp Non-member LinearOpWithSolveFactoryBase helper functions.
 *
 * These functions provide for simpler use cases for the use of
 * <tt>LinearOpWithSolveFactoryBase</tt> objects. and provide some
 * documentation for the various use cases.
 *
 * \ingroup Thyra_Op_Vec_Interoperability_Extended_Interfaces_grp
 */

/** \brief Create and initialize a <tt>LinearOpWithSolveBase</tt> object from
 * a <tt>LinearOpBase</tt> object using a
 * <tt>LinearOpWithSolveFactoryBase</tt> strategy object.
 *
 * \ingroup Thyra_LinearOpWithSolveFactoryBase_helper_grp
 */
template<class RangeScalar, class DomainScalar>
Teuchos::RefCountPtr<LinearOpWithSolveBase<RangeScalar,DomainScalar> >
createAndInitializeLinearOpWithSolve(
  const LinearOpWithSolveFactoryBase<RangeScalar,DomainScalar>                  &factory
  ,const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &fwdOp
  )
{
  Teuchos::RefCountPtr<LinearOpWithSolveBase<RangeScalar,DomainScalar> >
    Op = factory.createOp();
  factory.initializeOp(fwdOp,&*Op);
  return Op;
}

// /////////////////////////
// Implementations

template<class RangeScalar, class DomainScalar>
void LinearOpWithSolveFactoryBase<RangeScalar,DomainScalar>::initializeAndReuseOp(
  const Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >    &fwdOp
  ,LinearOpWithSolveBase<RangeScalar,DomainScalar>                             *Op
  ) const
{
  this->initializeOp(fwdOp,Op);
}

} // namespace Thyra

#endif // THYRA_LINEAR_OP_WITH_SOLVE_FACTORY_BASE_DECL_HPP
