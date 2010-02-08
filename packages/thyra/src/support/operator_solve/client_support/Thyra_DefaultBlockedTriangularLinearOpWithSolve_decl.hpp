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
#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Teuchos_Array.hpp"


namespace Thyra {


/** \brief Concrete composite <tt>LinearOpWithSolveBase</tt> subclass that
 * creates single upper or lower block triangular LOWSB object out of a set of
 * LOWSB objects along the diagonal with LOB objects off diagonal.
 *
 * This subclass implements a strictly block upper or lower triangular LOWSB
 * object.  With LOWSB objects on the diagonal, the block system can be solved
 * by solving for each of the diagonal blocks and moving lower (or upper)
 * blocks to the RHS on each successive solve.
 *
 * For example, consider the lower block triangular linear operator:

 \verbatim

       [ M(0,0)                   ]
   M = [ M(1,0)   M(1,1)          ]
       [ M(2,0)   M(2,1)   M(2,2) ]  

 \endverbatim
 
 * A linear system of the form:

 \verbatim

   M * x = b

   =>

   [ M(0,0)                   ]   [ x(0) ]   [ b(0) ]
   [ M(1,0)   M(1,1)          ] * [ x(1) ] = [ b(1) ]
   [ M(2,0)   M(2,1)   M(2,2) ]   [ x(2) ]   [ b(2) ]

 \endverbatim

 * is solved as:

 \verbatim

  x(0) = inv(M(0,0)) * b(0)
  x(1) = inv(M(1,1)) * ( b(1) - M(1,0)*x(0) )
  x(2) = inv(M(2,2)) * ( b(2) - M(2,0)*x(0) - M(2,1)*x(1) )

 \endverbatim

 * The same approach can be used for block upper triangular linear operators
 * as well of course.
 *
 * See the class <tt>DefaultBlockedTriangularLinearOpWithSolveFactory</tt> for
 * an example of how one of these objects can be created from any
 * <tt>PhysicallyBlockedLinearOpBase</tt> object and compatible
 * <tt>LinearWithSolveBase</tt> objects.
 *
 * ToDo: Finish Documentation!
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class DefaultBlockedTriangularLinearOpWithSolve
  : virtual public PhysicallyBlockedLinearOpWithSolveBase<Scalar>
{
public:

  /** @name Constructors/Initializers/Accessors */
  //@{

  /** \brief . */
  DefaultBlockedTriangularLinearOpWithSolve();

  /** \brief . */
  void setNonconstBlocks( const RCP<PhysicallyBlockedLinearOpBase<Scalar> > &blocks );

  /** \brief . */
  void setBlocks( const RCP<const PhysicallyBlockedLinearOpBase<Scalar> > &blocks );

  /** \brief . */
  RCP<PhysicallyBlockedLinearOpBase<Scalar> > getNonconstBlocks();

  /** \brief . */
  RCP<const PhysicallyBlockedLinearOpBase<Scalar> > getBlocks();

  //@}

  /** @name Overridden from PhysicallyBlockedLinearOpWithSolveBase */
  //@{
  
  /** \brief . */
  bool acceptsLOWSBlock(const int i, const int j) const;
  /** \brief . */
  void setNonconstLOWSBlock(
    const int i, const int j,
    const RCP<LinearOpWithSolveBase<Scalar> > &block
    );
  /** \brief . */
  void setLOWSBlock(
    const int i, const int j,
    const RCP<const LinearOpWithSolveBase<Scalar> > &block
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
    const RCP<const ProductVectorSpaceBase<Scalar> > &productRange,
    const RCP<const ProductVectorSpaceBase<Scalar> > &productDomain
    );
  /** \brief . */
  bool blockFillIsActive() const;
  /** \brief . */
  bool acceptsBlock(const int i, const int j) const;
  /** \brief . */
  void setNonconstBlock(
    const int i, const int j,
    const RCP<LinearOpBase<Scalar> > &block
    );
  /** \brief . */
  void setBlock(
    const int i, const int j,
    const RCP<const LinearOpBase<Scalar> > &block
    );
  /** \brief . */
  void endBlockFill();
  /** \brief . */
  void uninitialize();

  //@}

  /** @name Overridden from BlockedLinearOpWithSolveBase */
  //@{

  /** \brief . */
  RCP<LinearOpWithSolveBase<Scalar> >
  getNonconstLOWSBlock(const int i, const int j); 
  /** \brief . */
  RCP<const LinearOpWithSolveBase<Scalar> >
  getLOWSBlock(const int i, const int j) const; 

  //@}

  /** @name Overridden from BlockedLinearOpBase */
  //@{

  /** \brief . */
  RCP<const ProductVectorSpaceBase<Scalar> >
  productRange() const;
  /** \brief . */
  RCP<const ProductVectorSpaceBase<Scalar> >
  productDomain() const;
  /** \brief . */
  bool blockExists(const int i, const int j) const; 
  /** \brief . */
  bool blockIsConst(const int i, const int j) const; 
  /** \brief . */
  RCP<LinearOpBase<Scalar> >
  getNonconstBlock(const int i, const int j); 
  /** \brief . */
  RCP<const LinearOpBase<Scalar> >
  getBlock(const int i, const int j) const; 

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > range() const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > domain() const;
  /** \brief . */
  RCP<const LinearOpBase<Scalar> > clone() const;

  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{
                                                
  /** \brief Prints just the name
   * <tt>DefaultBlockedTriangularLinearOpWithSolve</tt> along with the overall
   * dimensions and the number of constituent operators.
   */
  std::string description() const;

  /** \brief Prints the details about the constituent linear operators.
   *
   * This function outputs different levels of detail based on the value
   * passed in for <tt>verbLevel</tt>:
   *
   * ToDo: Finish documentation!
   */
  void describe(
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

protected:
  
  /** @name Overridden from LinearOpBase */
  //@{
  /** \brief . */
  bool opSupportedImpl(EOpTransp M_trans) const;
  /** \brief . */
 void applyImpl(
   const EOpTransp M_trans,
   const MultiVectorBase<Scalar> &X,
   const Ptr<MultiVectorBase<Scalar> > &Y,
   const Scalar alpha,
   const Scalar beta
   ) const;
  //@}

  /** @name Overridden from LinearOpWithSolveBase */
  //@{
  /** \brief . */
  bool solveSupportsImpl(EOpTransp M_trans) const;
  /** \brief . */
  bool solveSupportsSolveMeasureTypeImpl(
    EOpTransp M_trans, const SolveMeasureType& solveMeasureType) const;
  /** \brief . */
  SolveStatus<Scalar> solveImpl(
    const EOpTransp transp,
    const MultiVectorBase<Scalar> &B,
    const Ptr<MultiVectorBase<Scalar> > &X,
    const Ptr<const SolveCriteria<Scalar> > solveCriteria
    ) const;
  //@}

private:

  // //////////////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<LinearOpWithSolveBase<Scalar> >
  CNCLOWS;

  typedef Teuchos::ConstNonconstObjectContainer<PhysicallyBlockedLinearOpBase<Scalar> >
  CNCPBLOB;

  // /////////////////////////
  // Private data members
  
  bool blockFillIsActive_;

  RCP<const ProductVectorSpaceBase<Scalar> > productRange_;
  RCP<const ProductVectorSpaceBase<Scalar> > productDomain_;
  int numDiagBlocks_;
  
  Array<CNCLOWS> diagonalBlocks_;

  // All blocks (including LOB form of diagonals)
  CNCPBLOB blocks_;

  // /////////////////////////
  // Private member functions

  void assertBlockFillIsActive(bool) const;

  void assertBlockRowCol(const int i, const int j) const;

  template<class LinearOpWithSolveType>
  void setLOWSBlockImpl( const int i, const int j,
    const RCP<LinearOpWithSolveType> &block );

  void assertAndSetBlockStructure(const PhysicallyBlockedLinearOpBase<Scalar>& blocks);
  
  // Not defined and not to be called

  DefaultBlockedTriangularLinearOpWithSolve(
    const DefaultBlockedTriangularLinearOpWithSolve&);

  DefaultBlockedTriangularLinearOpWithSolve&
  operator=(const DefaultBlockedTriangularLinearOpWithSolve&);
  
};


/** \brief Nonmember constructor.
 *
 * \relates DefaultBlockedTriangularLinearOpWithSolve
 */
template<class Scalar>
RCP<DefaultBlockedTriangularLinearOpWithSolve<Scalar> >
defaultBlockedTriangularLinearOpWithSolve()
{
  return Teuchos::rcp(new DefaultBlockedTriangularLinearOpWithSolve<Scalar>);
}


} // namespace Thyra


#endif	// THYRA_DEFAULT_BLOCKED_TRIANGULAR_LINEAR_OP_WITH_SOLVE_DECL_HPP
