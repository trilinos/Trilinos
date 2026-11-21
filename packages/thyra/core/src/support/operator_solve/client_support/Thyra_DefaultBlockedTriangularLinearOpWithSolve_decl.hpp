// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
 * The diagonal blocks are the LOWSB objects set in the block fill.  The
 * off-diagonal blocks are those of the blocked linear operator set with
 * <tt>setBlocks()</tt> or <tt>setNonconstBlocks()</tt> after the block fill
 * (whose diagonal blocks are not used); they can not be set in the block fill.
 * This blocked operator must have as many block rows and columns as there are
 * diagonal blocks, and be block lower or block upper triangular, not both.
 * The apply and the solves (block forward or backward substitution, also for
 * the transposed operator) include the off-diagonal blocks.
 * <tt>uninitialize()</tt> and a new block fill remove the off-diagonal blocks.
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

  /** \brief Throws <tt>std::logic_error</tt> if <tt>blocks</tt> does not
   * have as many block rows as block columns, or is not block triangular
   * (has nonzero blocks both below and above the diagonal).
   *
   * <tt>setBlocks()</tt> and <tt>setNonconstBlocks()</tt> also check this
   * (and the number of blocks); a factory can check it before reinitializing
   * the diagonal blocks, so that a rejected operator leaves the object
   * unchanged.
   */
  static void assertBlockTriangular(
    const PhysicallyBlockedLinearOpBase<Scalar> &blocks );

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
  /** \brief Returns false: the diagonal blocks are set with
   * <tt>setLOWSBlock()</tt> or <tt>setNonconstLOWSBlock()</tt>, and the
   * off-diagonal blocks with <tt>setBlocks()</tt> or
   * <tt>setNonconstBlocks()</tt> after the block fill. */
  bool acceptsBlock(const int i, const int j) const;
  /** \brief Throws, see <tt>acceptsBlock()</tt>. */
  void setNonconstBlock(
    const int i, const int j,
    const RCP<LinearOpBase<Scalar> > &block
    );
  /** \brief Throws, see <tt>acceptsBlock()</tt>. */
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
  /** \brief The diagonal blocks are the LOWSB blocks, and the off-diagonal
   * blocks those of the blocked operator of <tt>setBlocks()</tt>. */
  bool blockExists(const int i, const int j) const; 
  /** \brief . */
  bool blockIsConst(const int i, const int j) const; 
  /** \brief See <tt>blockExists()</tt>.  The off-diagonal blocks are
   * nonconst only if the blocked operator was set with
   * <tt>setNonconstBlocks()</tt> (throws otherwise). */
  RCP<LinearOpBase<Scalar> >
  getNonconstBlock(const int i, const int j); 
  /** \brief See <tt>blockExists()</tt>. */
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
  // Whether blocks_ has nonzero blocks below (above) the diagonal
  bool hasLowerBlocks_;
  bool hasUpperBlocks_;

  // /////////////////////////
  // Private member functions

  void assertBlockFillIsActive(bool) const;

  void assertBlockRowCol(const int i, const int j) const;

  template<class LinearOpWithSolveType>
  void setLOWSBlockImpl( const int i, const int j,
    const RCP<LinearOpWithSolveType> &block );

  void assertAndSetBlockStructure(const PhysicallyBlockedLinearOpBase<Scalar>& blocks);

  // Check the numbers of blocks and the structure of blocks (see
  // assertBlockTriangular()), and return whether it has nonzero blocks below
  // (hasLower) and above (hasUpper) the diagonal.  hasLower and hasUpper are
  // set only if the checks pass (they are unchanged if this throws).
  static void getBlockStructure(
    const PhysicallyBlockedLinearOpBase<Scalar> &blocks,
    bool &hasLower, bool &hasUpper );

  // Remove the off-diagonal blocks
  void removeOffDiagonalBlocks();

  // The off-diagonal block M[i,j], i != j, of blocks_ (null if zero or if
  // blocks_ is not set)
  RCP<const LinearOpBase<Scalar> >
  getOffDiagonalBlock(const int i, const int j) const;

  // The off-diagonal block op(M)[i,j] of op(M), to be applied with M_trans:
  // M[i,j] for M_trans = NOTRANS or CONJ, and M[j,i] for TRANS or CONJTRANS
  RCP<const LinearOpBase<Scalar> >
  getOpOffDiagonalBlock(const EOpTransp M_trans, const int i, const int j) const;

  // Whether all the off-diagonal blocks support op(), as needed by the apply
  // and the solves with op()
  bool offDiagonalBlocksSupport(const EOpTransp M_trans) const;
  
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
