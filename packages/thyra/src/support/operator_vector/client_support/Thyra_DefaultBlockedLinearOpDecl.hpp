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

#ifndef THYRA_DEFAULT_BLOCKED_LINEAR_OP_DECL_HPP
#define THYRA_DEFAULT_BLOCKED_LINEAR_OP_DECL_HPP

#include "Thyra_MultipliedLinearOpBase.hpp"
#include "Thyra_SingleScalarLinearOpBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Teuchos_arrayArg.hpp"

namespace Thyra {

/** \brief Concrete composite <tt>LinearOpBase</tt> subclass that creates
 * single linear operator object out of a set of constituent <tt>LinearOpBase</tt>
 * blocks.
 *
 * This class represents a blocked linear operator <tt>M</tt> of the form:

 \verbatim
 
  M =  [ Op[0,0], Op[0,1], ... , Op[0,N];
         Op[1,0], Op[1,1], ... , Op[1,N];
         .        .              .
         Op[M,0], Op[M,1], ... , Op[M,N]; ]

 \endverbatim

 * where <tt>Op[]</tt> is a logical 2D array of <tt>LinearOpBase</tt> objects
 * and <tt>M=this->productRange()->getNumBlocks()</tt> and
 * <tt>N=this->productDomain()->getNumBlocks().  Of course the operator
 * <tt>M</tt> is not constructed explicitly but instead just applies the
 * constituent linear operators with each set of blocks.
 *
 * ToDo: Finish Documentation!
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class DefaultBlockedLinearOp
  : virtual public PhysicallyBlockedLinearOpBase<Scalar>  // Public interface
  , virtual protected SingleScalarLinearOpBase<Scalar>    // Implementation detail
{
public:

  /** @name Overridden from PhysicallyBlockedLinearOpBase */
  //@{

  /** \brief . */
  void beginBlockFill();
  /** \brief . */
  void beginBlockFill(
    const Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> >  &productRange
    ,const Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> > &productDomain
    );
  /** \brief . */
  bool blockFillIsActive() const;
  /** \brief . */
  bool acceptsBlock(const int i, const int j) const;
  /** \brief . */
  void setNonconstBlock(
    const int i, const int j
    ,const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &block
    );
  /** \brief . */
  void setBlock(
    const int i, const int j
    ,const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &block
    );
  /** \brief . */
  void endBlockFill();

  //@}

  /** @name Overridden from BlockedLinearOpBase */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const ProductVectorSpaceBase<RangeScalar> >
  /** \brief . */
  productRange() const;
  /** \brief . */
  Teuchos::RefCountPtr<const ProductVectorSpaceBase<DomainScalar> >
  /** \brief . */
  productDomain() const;
  /** \brief . */
  bool blockExists(const int i, const int j) const; 
  /** \brief . */
  bool blockIsConst(const int i, const int j) const; 
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<RangeScalar,DomainScalar> >
  getNonconstBlock(const int i, const int j); 
  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpBase<RangeScalar,DomainScalar> >
  getBlock(const int i, const int j) const; 

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > range() const;
  /** \brief . */
  Teuchos::RefCountPtr< const VectorSpaceBase<Scalar> > domain() const;
  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpBase<Scalar> > clone() const;

  //@}


  /** @name Overridden from Teuchos::Describable */
  //@{
                                                
  /** \brief Prints just the name <tt>DefaultBlockedLinearOp</tt> along with
   * the overall dimensions and the number of constituent operators.
   */
  std::string description() const;

  /** \brief Prints the details about the constituent linear operators.
   *
   * This function outputs different levels of detail based on the value passed in
   * for <tt>verbLevel</tt>:
   *
   * ToDo: Finish documentation!
   */
  void describe(
    Teuchos::FancyOStream                &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ) const;

  //@}

protected:

  /** @name Overridden from SingleScalarLinearOpBase */
  //@{

  /** \brief Returns <tt>true</tt> only if all constituent operators support
   * <tt>M_trans</tt>.
   */
  bool opSupported(ETransp M_trans) const;

  /** \brief . */
  void apply(
    const ETransp                     M_trans
    ,const MultiVectorBase<Scalar>    &X
    ,MultiVectorBase<Scalar>          *Y
    ,const Scalar                     alpha
    ,const Scalar                     beta
    ) const;
  
  //@}

private:

  // ///////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<LinearOpBase<Scalar> > CNCLO;

  // /////////////////////////
  // Private data members
 
  std::vector<CNCLO> Ops_; // M x N

  // ///////////////////////////
  // Private member functions
  
  void assertInitialized() const;

  // Not defined and not to be called
  DefaultBlockedLinearOp(const DefaultBlockedLinearOp&);
  DefaultBlockedLinearOp& operator=(const DefaultBlockedLinearOp&);

};

} // namespace Thyra

#endif	// THYRA_DEFAULT_BLOCKED_LINEAR_OP_DECL_HPP
