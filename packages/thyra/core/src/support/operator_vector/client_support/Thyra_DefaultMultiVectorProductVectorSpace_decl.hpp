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

#ifndef THYRA_DEFAULT_MULTI_VECTOR_PRODUCT_VECTOR_SPACE_DECL_HPP
#define THYRA_DEFAULT_MULTI_VECTOR_PRODUCT_VECTOR_SPACE_DECL_HPP

#include "Thyra_ProductVectorSpaceBase.hpp" // Interface
#include "Thyra_DefaultProductVectorSpace.hpp" // Implementation
#include "Thyra_VectorSpaceDefaultBase.hpp"


namespace Thyra {


/** \brief Standard concrete implementation of a product vector space that
 * creates product vectors fromed implicitly from the columns of a
 * multi-vector.
 *
 * The default copy constructor is allowed since it has just the right
 * behavior (i.e. shallow copy).
 *
 * ToDo: Finish Documentation!
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
class DefaultMultiVectorProductVectorSpace
  : virtual public ProductVectorSpaceBase<Scalar>
  , virtual protected VectorSpaceDefaultBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to an uninitialized state. */
  DefaultMultiVectorProductVectorSpace();
  
  /** \brief Initialize with a list of constituent vector spaces.
   *
   * \param space [in,persisting] The vector space used to create the
   * multi-vectors.
   *
   * \param numColunns [in] The number of columns to create in the
   * multi-vector represented as a product vector.
   *
   * Preconditions:<ul>
   * <li> <tt>!is_null(space)</tt>
   * <li> <tt>numColumns > 0</tt>
   * </ul>
   *
   * Postconditions:<ul>
   * <li> <tt>this->dim() == space->dim() * numColumns</tt>
   * <li> <tt>this->numBlocks() == numColumns</tt>
   * <li> <tt>getBlock(i).get() == space.get(), i=0,...,numColumns-1</tt>
   * </ul>
   */
  void initialize(
    const RCP<const VectorSpaceBase<Scalar> > &space,
    const int numColumns
    );

  /** \brief . */
  RCP<const DefaultProductVectorSpace<Scalar> >
  getDefaultProductVectorSpace() const;

  /** \brief Uninitialize.
   *
   * \param numBlocks [out] If <tt>numBlocks!=NULL</tt> then on output
   * <tt>*numBlocks</tt> will be set to <tt>this->numBlocks()</tt>.
   *
   * \param vecSpaces [out] If <tt>vecSpaces!=NULL</tt> then
   * <tt>vecSpaces</tt> must point to an array of length
   * <tt>this->numBlocks</tt> and on output <tt>vecSpace[i]</tt> will be set
   * to <tt>this->vecSpaces()[i]</tt> for <tt>i=0,..,this->numBlocks()-1</tt>.
   *
   * Postconditions:<ul>
   * <li> <tt>this->numBlocks()==0</tt>
   * <li> <tt>vecSpaces()==NULL</tt>
   * </ul>
   *
   * <b>Warning!</b> If <tt>this->hasBeenCloned()==true</tt> then the client
   * had better not mess with the constituent vector spaces returned in
   * <tt>vecSpaces[]</tt> since another
   * <tt>DefaultMultiVectorProductVectorSpace</tt> object is still using them.
   */
  void uninitialize(
    RCP<const VectorSpaceBase<Scalar> > *space = 0,
    int *numColumns = 0
    );
    
  //@}
    
  /** @name Overridden from DefaultMultiVectorProductVectorSpace */
  //@{
    
  /** \brief . */
  int numBlocks() const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > getBlock(const int k) const; 

  //@}

  /** @name Overridden from VectorSpaceBase */
  //@{

  /** \brief . */
  Ordinal dim() const;
  /** \brief . */
  bool isCompatible( const VectorSpaceBase<Scalar>& vecSpc ) const;
  /** \brief Returns a <tt>DefaultMultiVectorProductVector</tt> object. */
  RCP< VectorBase<Scalar> > createMember() const;
  /** \brief Returns the sum of the scalar products of the constituent vectors. */
  Scalar scalarProd( const VectorBase<Scalar>& x, const VectorBase<Scalar>& y ) const;
  /** \brief Returns the sum of the scalar products of each of the columns of
   * the constituent multi-vectors.
   */
  void scalarProdsImpl(
    const MultiVectorBase<Scalar>& X, const MultiVectorBase<Scalar>& Y,
    const ArrayView<Scalar> &scalarProds ) const;
  /** \brief Returns true if all of the constituent vector spaces return
   * true.
   */
  bool hasInCoreView( const Range1D& rng, const EViewType viewType,
    const EStrideType strideType ) const;
  /** \brief Returns <tt>getBlock(0)->smallVecSpcFcty()</tt>. */
  RCP< const VectorSpaceFactoryBase<Scalar> > smallVecSpcFcty() const;
  /** \brief Returns a <tt>DefaultColumnwiseMultiVector</tt> object.
   *
   * ToDo: It is possible an general and well optimized multi-vector
   * implementation called something like MultiVectorProducMultiVector.  This
   * would require that you create an underlying multi-vector with
   * numBlocks*numMembers total columns.  However, this class is not needed at this
   * time so it is not provided.
   */
  RCP< MultiVectorBase<Scalar> > createMembers(int numMembers) const;
  /** \brief Clones the object as promised. */
  RCP< const VectorSpaceBase<Scalar> > clone() const;

  //@}

  /** @name Overridden from Teuchos::Describable */
  //@{
                                                
  /** \brief Prints just the name <tt>DefaultMultiVectorProductVectorSpace</tt> along
   * with the overall dimension and the number of blocks.
   */
  std::string description() const;

  /** \brief Prints the details about the constituent vector space.
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

private:
 
  // ///////////////////////////////////
  // Private data members

  RCP<const VectorSpaceBase<Scalar> > space_;
  int numColumns_;
  RCP<const DefaultProductVectorSpace<Scalar> > defaultProdVecSpc_;

  // ///////////////////////////////////
  // Private member functions

  void assertInitialized() const;

};


/** \brief Nonmember constructor function.
 *
 * \relates DefaultMultiVectorProductVectorSpace
 */
template<class Scalar>
inline
RCP<DefaultMultiVectorProductVectorSpace<Scalar> >
multiVectorProductVectorSpace()
{
  return Teuchos::rcp(new DefaultMultiVectorProductVectorSpace<Scalar>());
}


/** \brief Nonmember constructor function.
 *
 * \relates DefaultMultiVectorProductVectorSpace
 */
template<class Scalar>
inline
RCP<DefaultMultiVectorProductVectorSpace<Scalar> >
multiVectorProductVectorSpace(
  const RCP<const VectorSpaceBase<Scalar> > &space,
  const int numColumns
  )
{
  RCP<DefaultMultiVectorProductVectorSpace<Scalar> > multiVecProdVecSpace =
    multiVectorProductVectorSpace<Scalar>();
  multiVecProdVecSpace->initialize(space,numColumns);
  return multiVecProdVecSpace;
}


// /////////////////////////////////
// Inline members


template<class Scalar>
inline
RCP<const DefaultProductVectorSpace<Scalar> >
DefaultMultiVectorProductVectorSpace<Scalar>::getDefaultProductVectorSpace() const
{
  return defaultProdVecSpc_;
}


template<class Scalar>
inline
void DefaultMultiVectorProductVectorSpace<Scalar>::assertInitialized() const
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT( is_null(space_) );
#endif
}


} // namespace Thyra


#endif // THYRA_DEFAULT_MULTI_VECTOR_PRODUCT_VECTOR_SPACE_DECL_HPP
