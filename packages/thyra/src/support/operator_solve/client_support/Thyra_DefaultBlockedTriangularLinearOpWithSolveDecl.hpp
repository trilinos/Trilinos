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

#ifndef THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_DECL_HPP
#define THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_DECL_HPP


#include "Thyra_PhysicallyBlockedLinearOpWithSolveBase.hpp"
#include "Thyra_SingleScalarLinearOpWithSolveBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Teuchos_Array.hpp"


namespace Thyra {


/** \brief Concrete composite <tt>LinearOpWithSolveBase</tt> subclass that
 * creates single upper or lower block triangular LOWSB object out of a set of
 * LOWSB objects along the diagonal with LOB object on the off diagonal..
 *
 * ToDo: Finish Documentatioin!
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class DefaultBlockedTriangularLinearOpWithSolve
  : virtual public PhysicallyBlockedLinearOpWithSolveBase<Scalar>, // Public interface
    virtual protected SingleScalarLinearOpWithSolveBase<Scalar> // Implementation detail
{
public:

  /** @name Constructors */
  //@{

  /** \brief . */
  DefaultBlockedTriangularLinearOpWithSolve();

  //@}

  /** @name Overridden from PhysicallyBlockedLinearOpWithSolveBase */
  //@{
  
  /** \brief . */
  bool acceptsLOWSBlock(const int i, const int j) const;
  /** \brief . */
  void setNonconstLOWSBlock(
    const int i, const int j,
    const Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> > &block
    );
  /** \brief . */
  void setLOWSBlock(
    const int i, const int j,
    const Teuchos::RefCountPtr<const LinearOpWithSolveBase<Scalar> > &block
    );

  //@}

  /** @name Overridden from PhysicallyBlockedLinearOpBase */
  //@{

  /** \brief . */
  void beginBlockFill();
  /** \brief . */
  void beginBlockFill(
    const int numRowBlocks, const int numColBlocks
    );
  /** \brief . */
  void beginBlockFill(
    const Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> >  &productRange,
    const Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> > &productDomain
    );
  /** \brief . */
  bool blockFillIsActive() const;
  /** \brief . */
  bool acceptsBlock(const int i, const int j) const;
  /** \brief . */
  void setNonconstBlock(
    const int i, const int j,
    const Teuchos::RefCountPtr<LinearOpBase<Scalar> > &block
    );
  /** \brief . */
  void setBlock(
    const int i, const int j,
    const Teuchos::RefCountPtr<const LinearOpBase<Scalar> > &block
    );
  /** \brief . */
  void endBlockFill();
  /** \brief . */
  void uninitialize();

  //@}

  /** @name Overridden from BlockedLinearOpWithSolveBase */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<LinearOpWithSolveBase<Scalar> >
  getNonconstLOWSBlock(const int i, const int j); 
  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpWithSolveBase<Scalar> >
  getLOWSBlock(const int i, const int j) const; 

  //@}

  /** @name Overridden from BlockedLinearOpBase */
  //@{

  /** \brief . */
  Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> >
  productRange() const;
  /** \brief . */
  Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> >
  productDomain() const;
  /** \brief . */
  bool blockExists(const int i, const int j) const; 
  /** \brief . */
  bool blockIsConst(const int i, const int j) const; 
  /** \brief . */
  Teuchos::RefCountPtr<LinearOpBase<Scalar> >
  getNonconstBlock(const int i, const int j); 
  /** \brief . */
  Teuchos::RefCountPtr<const LinearOpBase<Scalar> >
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
                                                
  /** \brief Prints just the name <tt>DefaultBlockedTriangularLinearOpWithSolve</tt> along with
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
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

protected:

  /** @name Overridden from SingleScalarLinearOpWithSolveBase */
  //@{

  /** \brief . */
  bool solveSupportsTrans(ETransp M_trans) const;
  /** \brief . */
  bool solveSupportsSolveMeasureType(
    ETransp M_trans, const SolveMeasureType& solveMeasureType ) const;
  /** \brief . */
  void solve(
    const ETransp M_trans,
    const MultiVectorBase<Scalar> &B,
    MultiVectorBase<Scalar> *X,
    const int numBlocks,
    const BlockSolveCriteria<Scalar> blockSolveCriteria[],
    SolveStatus<Scalar> blockSolveStatus[]
    ) const;
  
  //@}

  /** @name Overridden from SingleScalarLinearOpBase */
  //@{

  /** \brief Returns <tt>true</tt> only if all constituent operators support
   * <tt>M_trans</tt>.
   */
  bool opSupported(ETransp M_trans) const;

  /** \brief . */
  void apply(
    const ETransp M_trans,
    const MultiVectorBase<Scalar> &X,
    MultiVectorBase<Scalar> *Y,
    const Scalar alpha,
    const Scalar beta
    ) const;
  
  //@}

private:

  // //////////////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<LinearOpWithSolveBase<Scalar> >
  CNCLOWS;

  // /////////////////////////
  // Private data members
  
  bool blockFillIsActive_;

  Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> > productRange_;
  Teuchos::RefCountPtr<const ProductVectorSpaceBase<Scalar> > productDomain_;
  int numDiagBlocks_;
  
  Teuchos::Array<CNCLOWS> diagonalBlocks_;

  // ToDo: Add members to support off-diagonal blocks!

  // /////////////////////////
  // Private member functions

  void assertBlockFillIsActive(bool) const;

  void assertBlockRowCol(const int i, const int j) const;

  template<class LinearOpWithSolveType>
  void setLOWSBlockImpl(
    const int i, const int j
    ,const Teuchos::RefCountPtr<LinearOpWithSolveType> &block
    );
  
  // Not defined and not to be called
  DefaultBlockedTriangularLinearOpWithSolve(
    const DefaultBlockedTriangularLinearOpWithSolve& );
  DefaultBlockedTriangularLinearOpWithSolve&
  operator=(const DefaultBlockedTriangularLinearOpWithSolve&);
  
};


} // namespace Thyra


#endif	// THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_DECL_HPP
