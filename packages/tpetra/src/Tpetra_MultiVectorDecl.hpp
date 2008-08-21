// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
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

// FINISH: some of these arrayview objects should be something else, like Ptr

#ifndef TPETRA_MULTIVECTOR_DECL_HPP
#define TPETRA_MULTIVECTOR_DECL_HPP

#include <Teuchos_CompObject.hpp>
#include <Teuchos_Object.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_BLAS_types.hpp>
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of MultiVectorData, needed to prevent circular inclusions
  template<typename Ordinal, typename Scalar> class MultiVectorData;
#endif

  /*! multivector */
  template<class Ordinal, class Scalar>
  class MultiVector : public Teuchos::CompObject, public DistObject<Ordinal,Scalar>{

    public:

    /* FINISH: what with these?
       int replaceMap (const Map<Ordinal> &map);
       int reduce ();
     */

    //@{ \name Constructor/Destructor Methods

    //! Basic MultiVector constuctor.
    MultiVector(const Map<Ordinal> &map, Ordinal numVectors, bool zeroOut=true);

    //! MultiVector copy constructor.
    MultiVector(const MultiVector<Ordinal,Scalar> &source);

    //! Set multi-vector values from two-dimensional array. (copy)
    MultiVector(const Map<Ordinal> &map, const Teuchos::ArrayView<const Scalar> &A, Ordinal LDA, Ordinal numVectors);

    //! Set multi-vector values from array of pointers. (copy)
    MultiVector(const Map<Ordinal> &map, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &arrayOfArrays, Ordinal numVectors);

    //! Set multi-vector values from a list of vectors in an existing MultiVector. (copy)
    MultiVector(const MultiVector<Ordinal,Scalar> &source, const Teuchos::ArrayView<const Ordinal> &indices);

    //! Set multi-vector values from a range of vectors in an existing MultiVector. (copy)
    MultiVector(const MultiVector<Ordinal,Scalar> &source, Ordinal startIndex, Ordinal numVectors);

    //! MultiVector destructor.
    virtual ~MultiVector();

    //@}

    //@{ \name Post-construction modification routines

    //! Replace current value at the specified (GlobalRow, VectorIndex) location with ScalarValue.
    void replaceGlobalValue(Ordinal globalRow, Ordinal vectorIndex, const Scalar &value);

    //! Adds ScalarValue to existing value at the specified (GlobalRow, VectorIndex) location.
    int sumIntoGlobalValue(Ordinal globalRow, Ordinal vectorIndex, const Scalar &value);

    //! Replace current value at the specified (MyRow, VectorIndex) location with ScalarValue.
    int replaceMyValue(int MyRow, int VectorIndex, const Scalar &ScalarValue);

    //! Adds ScalarValue to existing value at the specified (MyRow, VectorIndex) location.
    int sumIntoMyValue(int MyRow, int VectorIndex, const Scalar &ScalarValue);

    //! Initialize all values in a multi-vector with constant value.
    int putScalar(const Scalar &ScalarConstant);

    //! Set multi-vector values to random numbers.
    int random();

    //@} 

    //@{ \name Extraction methods

    // FINISH: should these be const or not?

    //! Return multi-vector values in user-provided two-dimensional array.
    void extractCopy(const Teuchos::ArrayView<Scalar> &A, Ordinal &MyLDA) const;

    //! Return multi-vector values in user-provided array of pointers.
    void extractCopy(const Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > &arrayOfArrays) const;

    //! Return non-const pointers to multi-vector values in user-provided two-dimensional array.
    void extractView(Teuchos::ArrayRCP<Scalar> &A, Ordinal &MyLDA);

    //! Return const pointers to multi-vector values in user-provided two-dimensional array.
    void extractConstView(Teuchos::ArrayRCP<const Scalar> &A, Ordinal &MyLDA) const;

    //! Return non-const pointers to multi-vector values in user-provided array of pointers.
    void extractView(Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > &arrayOfArrays);

    //! Return const pointers to multi-vector values in user-provided array of pointers.
    void extractConstView(Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > &arrayOfArrays) const;

    //@} 

    //@{ \name Mathematical methods

    // FINISH: expand documentation of these functions

    //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
    void dot(const MultiVector<Ordinal,Scalar> &A, Teuchos::Array<Scalar> &dots) const;

    //! Puts element-wise absolute values of input Multi-vector in target, this = abs(A).
    void abs(const MultiVector<Ordinal,Scalar> &A);

    //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    void reciprocal(const MultiVector<Ordinal,Scalar> &A);

    //! Scale the current values of a multi-vector, this = alpha*this.
    void scale(const Scalar &alpha);

    //! Replace multi-vector values with scaled values of A, this = alpha*A.
    void scale(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A);

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    void update(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const Scalar &beta);

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    void update(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const Scalar &beta, const MultiVector<Ordinal,Scalar> &B, const Scalar &gamma);

    //! Compute 1-norm of each vector in multi-vector.
    void norm1(Teuchos::Array<Scalar> &norms) const;

    //! Compute 2-norm of each vector in multi-vector.
    void norm2(Teuchos::Array<Scalar> &norms) const;

    //! Compute Inf-norm of each vector in multi-vector.
    void normInf(Teuchos::Array<Scalar> &norms) const;

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    void normWeighted(const MultiVector<Ordinal,Scalar> &weights, Teuchos::Array<Scalar> &norms) const;

    //! Compute minimum value of each vector in multi-vector.
    void minValue(Teuchos::Array<Scalar> &mins) const;

    //! Compute maximum value of each vector in multi-vector.
    void maxValue(Teuchos::Array<Scalar> &maxs) const;

    //! Compute mean (average) value of each vector in multi-vector.
    void meanValue(Teuchos::Array<Scalar> &means) const;

    //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
    void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const MultiVector<Ordinal,Scalar> &B, const Scalar &beta);

    //! Multiply a MultiVector with another, element-by-element: this(i,j) = beta*this(i,j) + alpha*A(i,j)*B(i,j)
    void multiply(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const MultiVector<Ordinal,Scalar> &B, const Scalar &beta);

    //! Multiply a MultiVector by the reciprocal of another, element-by-element. this(i,j) = beta*this(i,j) + alpha*B(i,j)/A(i,j)
    void reciprocalMultiply(const Scalar &alpha, const MultiVector<Ordinal,Scalar> &A, const MultiVector<Ordinal,Scalar> &B, const Scalar &beta);

    //@} 

    //@{ \name Random number utilities

    //! Set seed for Random function.
    int setSeed(unsigned int Seed_in);

    //! Get seed from Random function.
    unsigned getSeed();

    //@} 

    //@{ \name Overloaded operators

    //! = Operator.
    /*! \param In 
      A - Multivector to copy
     */
    MultiVector<Ordinal,Scalar> operator=(const MultiVector<Ordinal,Scalar> &source);

    //! Vector access function.
    /*! ArrayRCP to the local values in the ith vector of this multi-vector.
     */
    const Teuchos::ArrayRCP<Scalar> & operator[](Ordinal i);

    //! Vector access function.
    /** ArrayRCP to the local values in the ith vector of this multi-vector.
     */
    const Teuchos::ArrayRCP<const Scalar> & operator[](Ordinal i) const;

    //! Vector access function.
    Vector<Ordinal,Scalar> & operator()(Ordinal i);

    //! Vector access function.
    const Vector<Ordinal,Scalar> & operator() (Ordinal i) const;

    //@} 

    //@{ \name Attribute access functions

    //! Returns the number of vectors in the multi-vector.
    int numVectors() const;

    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
    int myLength() const;

    //! Returns the global vector length of vectors in the multi-vector.
    int globalLength() const;

    //! Returns the stride between vectors in the multi-vector (only meaningful if ConstantStride() is true).
    int stride() const;

    //! Returns true if this multi-vector has constant stride between vectors.
    bool constantStride() const;

    //@} 

    //@{ \name I/O methods

    //! Print method.
    void print(std::ostream &os) const;
    void printValues(std::ostream &os) const;

    //@} 

    //@{ \name Expert-only unsupported methods

    //! Reset the view of an existing multivector to point to new user data.
    int resetView(const Teuchos::ArrayRCP<const Teuchos::ArrayRCP<Scalar> > &arrayOfArrays);

    //! Get pointer to MultiVector values.
    const Teuchos::ArrayRCP<Scalar> & values();

    //! Get pointer to MultiVector values.
    const Teuchos::ArrayRCP<const Scalar> & valuesConst() const;

    //! Get pointer to individual vector pointers.
    const Teuchos::ArrayRCP<const Teuchos::ArrayRCP<Scalar> > & pointers();

    //! Get pointer to individual vector pointers.
    const Teuchos::ArrayRCP<const Teuchos::ArrayRCP<const Scalar> > & pointersConst() const;

    //@}

    protected:

    Teuchos::RCP<MultiVectorData<Ordinal,Scalar> > MVData_;

    // four functions needed for DistObject derivation
    bool checkSizes(const DistObject<Ordinal,Scalar> & sourceObj);

    int copyAndPermute(const DistObject<Ordinal,Scalar> & sourceObj,
               Ordinal numSameIDs,
               Ordinal numPermuteIDs,
               const std::vector<Ordinal> & permuteToLIDs,
               const std::vector<Ordinal> & permuteFromLIDs);

    int packAndPrepare(const DistObject<Ordinal,Scalar> & sourceObj,
               Ordinal numExportIDs,
               const std::vector<Ordinal> & exportLIDs,
               std::vector<Scalar>& exports,
               Ordinal &packetSize,
               Distributor<Ordinal> &distor);
  
    int unpackAndCombine(Ordinal numImportIDs,
               const std::vector<Ordinal> & importLIDs,
               const std::vector<Scalar> & imports,
               Distributor<Ordinal> &distor,
               CombineMode CM);


  }; // class MultiVector

} // namespace Tpetra


#endif // TPETRA_MULTIVECTOR_DECL_HPP
