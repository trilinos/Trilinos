// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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

#ifndef TPETRA_MULTIVECTOR_DECL_HPP
#define TPETRA_MULTIVECTOR_DECL_HPP

#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_MultiVector.hpp>
#include <Teuchos_LabeledObject.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayView.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_Range1D.hpp>
#include "Tpetra_DistObject.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CombineMode.hpp"

// TODO: add principal use case instructions for memory management interfaces (view/copy extraction)
// TODO: expand user-visible documentation 

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of Vector, needed to prevent circular inclusions
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node> class Vector;
#endif

  //! \brief A class for constructing and using dense, distributors multivectors.
  /*!
     This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
     The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
     type, if omitted, defaults to the \c LocalOrdinal type.
   */
  template<class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal, class Node=Kokkos::DefaultNode::DefaultNodeType>
  class MultiVector : public DistObject<Scalar,LocalOrdinal,GlobalOrdinal> {

    public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Basic MultiVector constuctor.
    MultiVector(Node &node, const Map<LocalOrdinal,GlobalOrdinal> &map, Teuchos_Ordinal NumVectors, bool zeroOut=true);

    //! MultiVector copy constructor.
    MultiVector(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source);

    //! Set multi-vector values from two-dimensional array using Teuchos memory management classes. (copy)
    /*! Post-condition: constantStride() == true */
    MultiVector(Node &node, const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayView<const Scalar> &A, Teuchos_Ordinal LDA, Teuchos_Ordinal NumVectors);

    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy)
    /*! Post-condition: constantStride() == true */
    MultiVector(Node &node, const Map<LocalOrdinal,GlobalOrdinal> &map, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &ArrayOfPtrs, Teuchos_Ordinal NumVectors);

    //! MultiVector destructor.
    virtual ~MultiVector();

    //@}

    //! @name Post-construction modification routines
    //@{ 

    //! Replace current value at the specified (globalRow, vectorIndex) location with specified value.
    void replaceGlobalValue(GlobalOrdinal globalRow, Teuchos_Ordinal vectorIndex, const Scalar &value);

    //! Adds specified value to existing value at the specified (globalRow, vectorIndex) location.
    void sumIntoGlobalValue(GlobalOrdinal globalRow, Teuchos_Ordinal vectorIndex, const Scalar &value);

    //! Replace current value at the specified (myRow, vectorIndex) location with specified value.
    void replaceMyValue(LocalOrdinal myRow, Teuchos_Ordinal vectorIndex, const Scalar &value);

    //! Adds specified value to existing value at the specified (myRow, vectorIndex) location.
    void sumIntoMyValue(LocalOrdinal myRow, Teuchos_Ordinal vectorIndex, const Scalar &value);

    //! Initialize all values in a multi-vector with specified value.
    void putScalar(const Scalar &value);

    //! Set multi-vector values to random numbers.
    void randomize();

    //! Replace the underlying Map with a compatible one.
    void replaceMap(const Map<LocalOrdinal,GlobalOrdinal> &map);

    //! Instruct a local (non-distributed) MultiVector to sum values across all nodes.
    void reduce();

    //! = Operator.
    /*! \param In A - Multivector to copy
     */
    MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& operator=(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source);

    //@}

    //! @name Data Copy and View get methods
    /** These methods are used to get the data underlying the MultiVector. They return data in one of three forms: 
      - a MultiVector with a subset of the columns of the target MultiVector
      - a raw C pointer or array of raw C pointers
      - one of the Teuchos memory management classes
      Not all of these methods are valid for a particular MultiVector. For instance, calling a method that accesses a 
      view of the data in a 1-D format (i.e., get1dView) requires that the target MultiVector has constant stride.
     */
    //@{

    //! Returns a MultiVector with copies of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subCopy(const Teuchos::Range1D &colRng) const;

    //! Returns a MultiVector with copies of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subCopy(const Teuchos::ArrayView<const Teuchos_Ordinal> &cols) const;

    //! Returns a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subView(const Teuchos::Range1D &colRng) const;

    //! Returns a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subView(Teuchos::ArrayView<const Teuchos_Ordinal> cols) const;

    //! Returns a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subViewNonConst(const Teuchos::Range1D &colRng);

    //! Returns a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subViewNonConst(Teuchos::ArrayView<const Teuchos_Ordinal> cols);

    //! Const Vector access function.
    Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getVector(Teuchos_Ordinal j) const;

    //! Vector access function.
    Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getVectorNonConst(Teuchos_Ordinal j);

    //! Const Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    Teuchos::ArrayRCP<const Scalar> getData(Teuchos_Ordinal j) const;

    //! Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    Teuchos::ArrayRCP<Scalar> getDataNonConst(Teuchos_Ordinal j);

    //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
    void get1dCopy(Teuchos::ArrayView<Scalar> A, Teuchos_Ordinal LDA) const;

    //! Return multi-vector values in user-provided two-dimensional array (using a C pointer).
    void get1dCopy(Scalar *A, Teuchos_Ordinal LDA) const;

    //! Return multi-vector values in user-provided array of pointers (using Teuchos memory management classes).
    void get2dCopy(Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > ArrayOfPtrs) const;

    //! Return multi-vector values in user-provided array of pointers (using C pointers).
    void get2dCopy(Scalar * const * ArrayOfPtrs) const;

    //! Return const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.
    Teuchos::ArrayRCP<const Scalar> get1dView() const;

    //! Return const persisting pointers to values.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > get2dView() const;

    //! Return non-const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.  Teuchos::ArrayRCP<Scalar> get1dViewNonConst();
    Teuchos::ArrayRCP<Scalar> get1dViewNonConst();

    //! Return non-const persisting pointers to values.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > get2dViewNonConst();

    //@}

    //! @name Mathematical methods
    //@{ 

    //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
    void dot(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Teuchos::ArrayView<Scalar> &dots) const;

    //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
    void abs(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A);

    //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    void reciprocal(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A);

    //! Scale the current values of a multi-vector, this = alpha*this.
    void scale(const Scalar &alpha);

    //! Scale the current values of a multi-vector, this = alpha*this.
    void scale(Teuchos::ArrayView<const Scalar> alpha);

    //! Replace multi-vector values with scaled values of A, this = alpha*A.
    void scale(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A);

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    void update(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Scalar &beta);

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    void update(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Scalar &beta, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &gamma);

    //! Compute 1-norm of each vector in multi-vector.
    void norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute 2-norm of each vector in multi-vector.
    void norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute Inf-norm of each vector in multi-vector.
    void normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    void normWeighted(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const;

    //! Compute mean (average) value of each vector in multi-vector.
    void meanValue(const Teuchos::ArrayView<Scalar> &means) const;

    //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
    void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &beta);

    //! Multiply a MultiVector with another, element-by-element: this(i,j) = beta*this(i,j) + alpha*A(i,j)*B(i,j)
    void multiply(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &beta);

    //! Multiply a MultiVector by the reciprocal of another, element-by-element. this(i,j) = beta*this(i,j) + alpha*B(i,j)/A(i,j)
    void reciprocalMultiply(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &beta);

    //@} 

    //! @name Attribute access functions
    //@{ 

    //! Returns the number of vectors in the multi-vector.
    Teuchos_Ordinal getNumVectors() const;

    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
    LocalOrdinal getMyLength() const;

    //! Returns the global vector length of vectors in the multi-vector.
    GlobalOrdinal getGlobalLength() const;

    //! Returns the stride between vectors in the multi-vector (only meaningful if ConstantStride() is true). WARNING: this may vary from node to node.
    Teuchos_Ordinal getStride() const;

    //! Returns true if this multi-vector has constant stride between vectors. WARNING: This may vary from node to node.
    bool isConstantStride() const;

    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    std::string description() const;

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

    //@}

    protected:

    typedef Kokkos::MultiVector<Scalar,Node>  KMV;
    typedef Kokkos::DefaultArithmetic<KMV>   DMVA;
    KMV lclMV_;
    Teuchos::Array<Teuchos_Ordinal> whichVectors_;

    //! Advanced constructor for non-contiguous views.
    MultiVector(Node &node, const Map<LocalOrdinal,GlobalOrdinal> &map,
                Teuchos::ArrayRCP<Scalar> data, Teuchos_Ordinal LDA, Teuchos::ArrayView<const Teuchos_Ordinal> whichVectors);

    //! Advanced constructor for contiguous views.
    MultiVector(Node &node, const Map<LocalOrdinal,GlobalOrdinal> &map,
                Teuchos::ArrayRCP<Scalar> data, Teuchos_Ordinal LDA, Teuchos_Ordinal NumVectors);

    // four functions needed for DistObject derivation
    bool checkSizes(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> &sourceObj, Teuchos_Ordinal &packetSize);

    void copyAndPermute(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> &sourceObj,
                              Teuchos_Ordinal numSameIDs,
                        const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                        const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs);

    void packAndPrepare(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal> &sourceObj,
                        const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
                        const Teuchos::ArrayView<Scalar> &exports,
                        Distributor &distor);

    void unpackAndCombine(const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                          const Teuchos::ArrayView<const Scalar> &imports,
                          Distributor &distor,
                          CombineMode CM);

  }; // class MultiVector

} // namespace Tpetra


#endif // TPETRA_MULTIVECTOR_DECL_HPP
