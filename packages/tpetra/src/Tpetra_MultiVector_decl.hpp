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

#include <Teuchos_LabeledObject.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_Range1D.hpp>

#include <Kokkos_MultiVector.hpp>
#include <Kokkos_DefaultArithmetic.hpp>

#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_DistObject.hpp"
#include "Tpetra_Map.hpp"

// TODO: add principal use case instructions for memory management interfaces (view/copy extraction)
// TODO: expand user-visible documentation 

namespace Tpetra {

#ifndef DOXYGEN_SHOULD_SKIP_THIS
  // forward declaration of Vector, needed to prevent circular inclusions
  template<class S, class LO, class GO, class N> class Vector;
#endif

  //! \brief A class for constructing and using dense, distributors multivectors.
  /*!
     This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
     The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
     type, if omitted, defaults to the \c LocalOrdinal type.
   */
  template <class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal, class Node=Kokkos::DefaultNode::DefaultNodeType>
  class MultiVector : public DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

    public:
      typedef Scalar        scalar_type;
      typedef LocalOrdinal  local_ordinal_type;
      typedef GlobalOrdinal global_ordinal_type;
      typedef Node          node_type;

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Basic MultiVector constuctor.
    MultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, size_t NumVectors, bool zeroOut=true);

    //! MultiVector copy constructor.
    MultiVector(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source);

    //! Set multi-vector values from two-dimensional array using Teuchos memory management classes. (copy)
    /*! Post-condition: constantStride() == true */
    MultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, const Teuchos::ArrayView<const Scalar> &A, size_t LDA, size_t NumVectors);

    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy)
    /*! Post-condition: constantStride() == true */
    MultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &ArrayOfPtrs, size_t NumVectors);

    //! MultiVector destructor.
    virtual ~MultiVector();

    //@}

    //! @name Post-construction modification routines
    //@{ 

    //! Replace current value at the specified (globalRow, vectorIndex) location with specified value.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
      */
    void replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value);

    //! Adds specified value to existing value at the specified (globalRow, vectorIndex) location.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
      */
    void sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value);

    //! Replace current value at the specified (myRow, vectorIndex) location with specified value.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
      */
    void replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value);

    //! Adds specified value to existing value at the specified (myRow, vectorIndex) location.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
      */
    void sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value);

    //! Initialize all values in a multi-vector with specified value.
    void putScalar(const Scalar &value);

    //! Set multi-vector values to random numbers.
    void randomize();

    //! Replace the underlying Map with a compatible one.
    void replaceMap(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map);

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
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subCopy(const Teuchos::ArrayView<const size_t> &cols) const;

    //! Returns a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subView(const Teuchos::Range1D &colRng) const;

    //! Returns a const MultiVector with const views of selected columns.
    Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subView(const Teuchos::ArrayView<const size_t> &cols) const;

    //! Returns a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subViewNonConst(const Teuchos::Range1D &colRng);

    //! Returns a MultiVector with views of selected columns.
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subViewNonConst(const Teuchos::ArrayView<const size_t> &cols);

    //! \brief Returns a const MultiVector view of a subset of rows.
    /** 
        Returns a const view of this MultiVector consisting of a subset of the rows, as specified by an offset and a sub-Map.

        \param In subMap - The row map for the new MultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre  <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
     */
    Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > offsetView(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &subMap, size_t offset) const;

    //! \brief Returns a non-const MultiVector view of a subset of rows.
    /** 
        Returns a non-const view of this MultiVector consisting of a subset of the rows, as specified by an offset and a sub-Map.

        \param In subMap - The row map for the new MultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre  <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
     */
    Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > offsetViewNonConst(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &subMap, size_t offset);

    //! Const Vector access function.
    Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getVector(size_t j) const;

    //! Vector access function.
    Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getVectorNonConst(size_t j);

    //! Const Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    Teuchos::ArrayRCP<const Scalar> getData(size_t j) const;

    //! Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    Teuchos::ArrayRCP<Scalar> getDataNonConst(size_t j);

    //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
    void get1dCopy(Teuchos::ArrayView<Scalar> A, size_t LDA) const;

    //! Return multi-vector values in user-provided array of pointers (using Teuchos memory management classes).
    void get2dCopy(Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > ArrayOfPtrs) const;

    //! Return const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.
    Teuchos::ArrayRCP<const Scalar> get1dView() const;

    //! Return const persisting pointers to values.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > get2dView() const;

    //! Return non-const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.  Teuchos::ArrayRCP<Scalar> get1dViewNonConst();
    Teuchos::ArrayRCP<Scalar> get1dViewNonConst();

    //! Return non-const persisting pointers to values.
    Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > get2dViewNonConst();

    //! Return a const reference to the underlying Kokkos::MultiVector object (advanced use only)
    const Kokkos::MultiVector<Scalar,Node> & getLocalMV() const;

    //! Return a non-const reference to the underlying Kokkos::MultiVector object (advanced use only)
    Kokkos::MultiVector<Scalar,Node> & getLocalMVNonConst();

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

    //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
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

    //! Element-wise multiply of a Vector A with a MultiVector B.
    /** Forms this = scalarThis * this + scalarAB * B @ A
     *  where @ denotes element-wise multiplication.
     *  B must be the same shape (size and num-vectors) as this, while
     *  A is the same size but a single vector (column).
     */
    void elementWiseMultiply(Scalar scalarAB, const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, Scalar scalarThis);
    //@} 

    //! @name Attribute access functions
    //@{ 

    //! Returns the number of vectors in the multi-vector.
    size_t getNumVectors() const;

    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
    size_t getLocalLength() const;

    //! Returns the global vector length of vectors in the multi-vector.
    global_size_t getGlobalLength() const;

    //! Returns the stride between vectors in the multi-vector (only meaningful if ConstantStride() is true). WARNING: this may vary from node to node.
    size_t getStride() const;

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

    typedef Kokkos::MultiVector<Scalar,Node> KMV;
    typedef Kokkos::DefaultArithmetic<KMV>   MVT;
  
    inline bool vectorIndexOutOfRange(size_t VectorIndex) const {
      return (VectorIndex < 1 && VectorIndex != 0) || VectorIndex >= getNumVectors();
    }

    KMV lclMV_;
    Teuchos::Array<size_t> whichVectors_;

    template <class T>
    //! Get persisting view of j-th column in given ArrayRCP, considering isConstantStride(). ArrayRCP may correspond to a compute buffer or host view.
    Teuchos::ArrayRCP<T> getSubArrayRCP(Teuchos::ArrayRCP<T> arr, size_t j) const;

    //! Advanced constructor for non-contiguous views.
    MultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map,
                Teuchos::ArrayRCP<Scalar> data, size_t LDA, Teuchos::ArrayView<const size_t> whichVectors);

    //! Advanced constructor for contiguous views.
    MultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map,
                Teuchos::ArrayRCP<Scalar> data, size_t LDA, size_t NumVectors);

    // four functions needed for DistObject derivation
    bool checkSizes(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> &sourceObj);

    void copyAndPermute(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> &sourceObj,
                        size_t numSameIDs,
                        const Teuchos::ArrayView<const LocalOrdinal> &permuteToLIDs,
                        const Teuchos::ArrayView<const LocalOrdinal> &permuteFromLIDs);

    void packAndPrepare(const DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> &sourceObj,
                        const Teuchos::ArrayView<const LocalOrdinal> &exportLIDs,
                        Teuchos::Array<Scalar> &exports,
                        const Teuchos::ArrayView<size_t> &numExportPacketsPerLID,
                        size_t& constantNumPackets,
                        Distributor &distor);

    void unpackAndCombine(const Teuchos::ArrayView<const LocalOrdinal> &importLIDs,
                          const Teuchos::ArrayView<const Scalar> &imports,
                          const Teuchos::ArrayView<size_t> &numPacketsPerLID,
                          size_t constantNumPackets,
                          Distributor &distor,
                          CombineMode CM);

  mutable Teuchos::ArrayRCP<Scalar> ncview_;
  mutable Teuchos::ArrayRCP<const Scalar> cview_;

  void createViews() const;
  void createViewsNonConst(Kokkos::ReadWriteOption rwo);
  void releaseViews() const;

  }; // class MultiVector

  /** \brief Non-member function to create a MultiVector from a specified Map.
      \relates MultiVector
   */
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
  createMultiVector(const Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, size_t numVectors);


} // namespace Tpetra


#endif // TPETRA_MULTIVECTOR_DECL_HPP
