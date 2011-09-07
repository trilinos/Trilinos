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


#ifndef THYRA_DEFAULT_BLOCKED_LINEAR_OP_DECL_HPP
#define THYRA_DEFAULT_BLOCKED_LINEAR_OP_DECL_HPP


#include "Thyra_PhysicallyBlockedLinearOpBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"



namespace Thyra {


template<class Scalar> class DefaultProductVectorSpace;


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
 * <tt>N=this->productDomain()->getNumBlocks()</tt>.  Of course the operator
 * <tt>M</tt> is not constructed explicitly but instead just applies the
 * constituent linear operators with each set of blocks.
 *
 * ToDo: Finish Documentation!
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class DefaultBlockedLinearOp
  : virtual public PhysicallyBlockedLinearOpBase<Scalar>

{
public:

  /** @name Constructors */
  //@{

  /** \brief . */
  DefaultBlockedLinearOp();

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
    const Teuchos::RCP<const ProductVectorSpaceBase<Scalar> > &productRange,
    const Teuchos::RCP<const ProductVectorSpaceBase<Scalar> > &productDomain
    );
  /** \brief . */
  bool blockFillIsActive() const;
  /** \brief . */
  bool acceptsBlock(const int i, const int j) const;
  /** \brief . */
  void setNonconstBlock(
    const int i, const int j,
    const Teuchos::RCP<LinearOpBase<Scalar> > &block
    );
  /** \brief . */
  void setBlock(
    const int i, const int j
    ,const Teuchos::RCP<const LinearOpBase<Scalar> > &block
    );
  /** \brief . */
  void endBlockFill();
  /** \brief . */
  void uninitialize();

  //@}

  /** @name Overridden from BlockedLinearOpBase */
  //@{

  /** \brief . */
  Teuchos::RCP<const ProductVectorSpaceBase<Scalar> >
  productRange() const;
  /** \brief . */
  Teuchos::RCP<const ProductVectorSpaceBase<Scalar> >
  productDomain() const;
  /** \brief . */
  bool blockExists(const int i, const int j) const; 
  /** \brief . */
  bool blockIsConst(const int i, const int j) const; 
  /** \brief . */
  Teuchos::RCP<LinearOpBase<Scalar> >
  getNonconstBlock(const int i, const int j); 
  /** \brief . */
  Teuchos::RCP<const LinearOpBase<Scalar> >
  getBlock(const int i, const int j) const; 

  //@}

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief . */
  Teuchos::RCP< const VectorSpaceBase<Scalar> > range() const;
  /** \brief . */
  Teuchos::RCP< const VectorSpaceBase<Scalar> > domain() const;
  /** \brief . */
  Teuchos::RCP<const LinearOpBase<Scalar> > clone() const;

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
    Teuchos::FancyOStream &out,
    const Teuchos::EVerbosityLevel verbLevel
    ) const;

  //@}

protected:

  /** @name Overridden from LinearOpBase */
  //@{

  /** \brief Returns <tt>true</tt> only if all constituent operators support
   * <tt>M_trans</tt>.
   */
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

private:

  // ///////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<LinearOpBase<Scalar> > CNCLO;
  typedef Teuchos::Array<Teuchos::RCP<const VectorSpaceBase<Scalar> > > vec_array_t;

  template<class Scalar2>
  struct BlockEntry {
    BlockEntry() : i(-1), j(-1) {}
    BlockEntry( const int i_in, const int j_in, const CNCLO &block_in )
      :i(i_in),j(j_in),block(block_in)
      {}
    int i;
    int j;
    CNCLO block;
  };

  // /////////////////////////
  // Private data members

  Teuchos::RCP<const ProductVectorSpaceBase<Scalar> > productRange_;
  Teuchos::RCP<const ProductVectorSpaceBase<Scalar> > productDomain_;
  Teuchos::RCP<const DefaultProductVectorSpace<Scalar> > defaultProductRange_;
  Teuchos::RCP<const DefaultProductVectorSpace<Scalar> > defaultProductDomain_;
  int numRowBlocks_; // M
  int numColBlocks_; // N
 
  std::vector<CNCLO> Ops_; // Final M x N storage
 
  vec_array_t rangeBlocks_;
  vec_array_t domainBlocks_;
  std::vector<BlockEntry<Scalar> > Ops_stack_; // Temp stack of ops begin filled (if Ops_.size()==0).
  bool blockFillIsActive_;

  // ///////////////////////////
  // Private member functions
 
  void resetStorage( const int numRowBlocks, const int numColBlocks );
  void assertBlockFillIsActive(bool) const;
  void assertBlockRowCol(const int i, const int j) const;
  void setBlockSpaces(
    const int i, const int j, const LinearOpBase<Scalar> &block
    );
  template<class LinearOpType>
  void setBlockImpl(
    const int i, const int j,
    const Teuchos::RCP<LinearOpType> &block
    );
  void adjustBlockSpaces();

  // Not defined and not to be called
  DefaultBlockedLinearOp(const DefaultBlockedLinearOp&);
  DefaultBlockedLinearOp& operator=(const DefaultBlockedLinearOp&);

};


/** \brief Nonmember default constructor.
 *
 * \relates DefaultBlockedLinearOp
 */
template<class Scalar>
RCP<DefaultBlockedLinearOp<Scalar> > defaultBlockedLinearOp();


/** \brief Form an implicit block 1x1 linear operator <tt>[ A00 ]</tt>.
 *
 * \relates DefaultBlockedLinearOp
 */
template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
block1x1(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &A00,
  const std::string &label = ""
  );


/** \brief Form an implicit block 1x2 linear operator <tt>[ A00, A01 ]</tt>.
 *
 * \relates DefaultBlockedLinearOp
 */
template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
block1x2(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &A00,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &A01,
  const std::string &label = ""
  );


/** \brief Form an implicit block 2x1 linear operator <tt>[ A00; A10 ]</tt>.
 *
 * \relates DefaultBlockedLinearOp
 */
template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
block2x1(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &A00,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &A10,
  const std::string &label = ""
  );


/** \brief Form an implicit block 2x2 linear operator <tt>[ A00, A01; A10, A11 ]</tt>.
 *
 * \relates DefaultBlockedLinearOp
 */
template<class Scalar>
Teuchos::RCP<const LinearOpBase<Scalar> >
block2x2(
  const Teuchos::RCP<const LinearOpBase<Scalar> > &A00,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &A01,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &A10,
  const Teuchos::RCP<const LinearOpBase<Scalar> > &A11,
  const std::string &label = ""
  );


/** \brief Form an implicit block 1x1 linear operator <tt>[ A00 ]</tt>.
 *
 * \relates DefaultBlockedLinearOp
 */
template<class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
nonconstBlock1x1(
  const Teuchos::RCP<LinearOpBase<Scalar> > &A00,
  const std::string &label = ""
  );


/** \brief Form an implicit block 1x2 linear operator <tt>[ A00, A01 ]</tt>.
 *
 * \relates DefaultBlockedLinearOp
 */
template<class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
nonconstBlock1x2(
  const Teuchos::RCP<LinearOpBase<Scalar> > &A00,
  const Teuchos::RCP<LinearOpBase<Scalar> > &A01,
  const std::string &label = ""
  );


/** \brief Form an implicit block 2x1 linear operator <tt>[ A00; A10 ]</tt>.
 *
 * \relates DefaultBlockedLinearOp
 */
template<class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
nonconstBlock2x1(
  const Teuchos::RCP<LinearOpBase<Scalar> > &A00,
  const Teuchos::RCP<LinearOpBase<Scalar> > &A10,
  const std::string &label = ""
  );


/** \brief Form an implicit block 2x2 linear operator <tt>[ A00, A01; A10, A11 ]</tt>.
 *
 * \relates DefaultBlockedLinearOp
 */
template<class Scalar>
Teuchos::RCP<LinearOpBase<Scalar> >
nonconstBlock2x2(
  const Teuchos::RCP<LinearOpBase<Scalar> > &A00,
  const Teuchos::RCP<LinearOpBase<Scalar> > &A01,
  const Teuchos::RCP<LinearOpBase<Scalar> > &A10,
  const Teuchos::RCP<LinearOpBase<Scalar> > &A11,
  const std::string &label = ""
  );


} // namespace Thyra


#endif	// THYRA_DEFAULT_BLOCKED_LINEAR_OP_DECL_HPP
