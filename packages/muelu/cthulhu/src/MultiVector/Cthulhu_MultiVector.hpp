#ifndef CTHULHU_MULTIVECTOR_DECL_HPP
#define CTHULHU_MULTIVECTOR_DECL_HPP

#include <Teuchos_LabeledObject.hpp>
#include <Teuchos_DataAccess.hpp>
#include <Teuchos_BLAS_types.hpp>
#include <Teuchos_Range1D.hpp>

#include <Kokkos_MultiVector.hpp>
#include <Kokkos_DefaultArithmetic.hpp>

#include "Cthulhu_ConfigDefs.hpp"
//#include "Cthulhu_DistObject.hpp" TODO
//#include "Cthulhu_MultiVector.hpp"
#include "Cthulhu_Map.hpp"

#include "Cthulhu_Debug.hpp"

namespace Cthulhu {

  // TODO
  // #ifndef DOXYGEN_SHOULD_SKIP_THIS
  //   // forward declaration of Vector, needed to prevent circular inclusions
  //   template<class S, class LO, class GO, class N> class Vector;
  // #endif

  //! \brief A class for constructing and using dense, distributors multivectors.
  /*!
    This class is templated on \c Scalar, \c LocalOrdinal and \c GlobalOrdinal. 
    The \c LocalOrdinal type, if omitted, defaults to \c int. The \c GlobalOrdinal 
    type, if omitted, defaults to the \c LocalOrdinal type.
  */
  template <class Scalar, class LocalOrdinal=int, class GlobalOrdinal=LocalOrdinal, class Node=Kokkos::DefaultNode::DefaultNodeType>
  class MultiVector { //: public DistObject<Scalar,LocalOrdinal,GlobalOrdinal,Node> { // TODO

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! MultiVector destructor.
    virtual ~MultiVector() { CTHULHU_DEBUG_ME; }

    //@}

    //! @name Post-construction modification routines
    //@{ 

    //! Replace current value at the specified (globalRow, vectorIndex) location with specified value.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
#ifdef CTHULHU_TPETRA_ONLY
    virtual void replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Adds specified value to existing value at the specified (globalRow, vectorIndex) location.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
#ifdef CTHULHU_TPETRA_ONLY
    virtual void sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Replace current value at the specified (myRow, vectorIndex) location with specified value.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
#ifdef CTHULHU_TPETRA_ONLY
    virtual void replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Adds specified value to existing value at the specified (myRow, vectorIndex) location.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
#ifdef CTHULHU_TPETRA_ONLY
    virtual void sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Initialize all values in a multi-vector with specified value.
#ifdef CTHULHU_TPETRA_ONLY
    virtual void putScalar(const Scalar &value) =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Set multi-vector values to random numbers.
    virtual void randomize() =0;

    //! Replace the underlying Map with a compatible one.
#ifdef CTHULHU_NOT_IMPLEMENTED
    virtual void replaceMap(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map) =0;
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Instruct a local (non-distributed) MultiVector to sum values across all nodes.
#ifdef CTHULHU_TPETRA_ONLY
    virtual void reduce() =0;
#endif // CTHULHU_TPETRA_ONLY

    //! = Operator.
    /*! \param In A - Multivector to copy
     */
#ifdef CTHULHU_NOT_IMPLEMENTED
    virtual MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& operator=(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source) =0;
#endif // CTHULHU_NOT_IMPLEMENTED

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

#ifdef CTHULHU_NOT_IMPLEMENTED

    //! Returns a MultiVector with copies of selected columns.
    virtual Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subCopy(const Teuchos::Range1D &colRng) const =0;

    //! Returns a MultiVector with copies of selected columns.
    virtual Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subCopy(const Teuchos::ArrayView<const size_t> &cols) const =0;

    //! Returns a const MultiVector with const views of selected columns.
    virtual Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subView(const Teuchos::Range1D &colRng) const =0;

    //! Returns a const MultiVector with const views of selected columns.
    virtual Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subView(const Teuchos::ArrayView<const size_t> &cols) const =0;

    //! Returns a MultiVector with views of selected columns.
    virtual Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subViewNonConst(const Teuchos::Range1D &colRng) =0;

    //! Returns a MultiVector with views of selected columns.
    virtual Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subViewNonConst(const Teuchos::ArrayView<const size_t> &cols) =0;

    //! \brief Returns a const MultiVector view of a subset of rows.
    /** 
        Returns a const view of this MultiVector consisting of a subset of the rows, as specified by an offset and a sub-Map.

        \param In subMap - The row map for the new MultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre  <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
    */
    virtual Teuchos::RCP<const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > offsetView(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &subMap, size_t offset) const =0;

    //! \brief Returns a non-const MultiVector view of a subset of rows.
    /** 
        Returns a non-const view of this MultiVector consisting of a subset of the rows, as specified by an offset and a sub-Map.

        \param In subMap - The row map for the new MultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre  <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
    */
    virtual Teuchos::RCP<MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > offsetViewNonConst(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &subMap, size_t offset) =0;

    //! Const Vector access function.
    virtual Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getVector(size_t j) const =0;

    //! Vector access function.
    virtual Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getVectorNonConst(size_t j) =0;
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Const Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
#ifdef CTHULHU_TPETRA_ONLY
    virtual Teuchos::ArrayRCP<const Scalar> getData(size_t j) const =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    virtual Teuchos::ArrayRCP<Scalar> getDataNonConst(size_t j) =0;

    //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
#ifdef CTHULHU_TPETRA_ONLY
    virtual void get1dCopy(Teuchos::ArrayView<Scalar> A, size_t LDA) const =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Return multi-vector values in user-provided array of pointers (using Teuchos memory management classes).
#ifdef CTHULHU_TPETRA_ONLY
    virtual void get2dCopy(Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > ArrayOfPtrs) const =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Return const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.
#ifdef CTHULHU_TPETRA_ONLY
    virtual Teuchos::ArrayRCP<const Scalar> get1dView() const =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Return const persisting pointers to values.
#ifdef CTHULHU_TPETRA_ONLY
    virtual Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > get2dView() const =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Return non-const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.  Teuchos::ArrayRCP<Scalar> get1dViewNonConst() =0;
#ifdef CTHULHU_TPETRA_ONLY
    virtual Teuchos::ArrayRCP<Scalar> get1dViewNonConst() =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Return non-const persisting pointers to values.
#ifdef CTHULHU_TPETRA_ONLY
    virtual Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > get2dViewNonConst() =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Return a const reference to the underlying Kokkos::MultiVector object (advanced use only)
#ifdef CTHULHU_TPETRA_ONLY
    virtual const Kokkos::MultiVector<Scalar,Node> & getLocalMV() const =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Return a non-const reference to the underlying Kokkos::MultiVector object (advanced use only)
#ifdef CTHULHU_TPETRA_ONLY
    virtual Kokkos::MultiVector<Scalar,Node> & getLocalMVNonConst() =0;
#endif // CTHULHU_TPETRA_ONLY

    //@}

    //! @name Mathematical methods
    //@{ 

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
    virtual void dot(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Teuchos::ArrayView<Scalar> &dots) const =0;

    //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
    virtual void abs(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) =0;

    //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    virtual void reciprocal(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) =0;
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Scale the current values of a multi-vector, this = alpha*this.
#ifdef CTHULHU_TPETRA_ONLY
    virtual void scale(const Scalar &alpha) =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
#ifdef CTHULHU_TPETRA_ONLY
    virtual void scale(Teuchos::ArrayView<const Scalar> alpha) =0;
#endif // CTHULHU_TPETRA_ONLY

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Replace multi-vector values with scaled values of A, this = alpha*A.
    virtual void scale(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) =0;

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    virtual void update(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Scalar &beta) =0;

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    virtual void update(const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Scalar &beta, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &gamma) =0;
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Compute 1-norm of each vector in multi-vector.
#ifdef CTHULHU_TPETRA_ONLY
    virtual void norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Compute 2-norm of each vector in multi-vector.
#ifdef CTHULHU_TPETRA_ONLY
    virtual void norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Compute Inf-norm of each vector in multi-vector.
#ifdef CTHULHU_TPETRA_ONLY
    virtual void normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
#ifdef CTHULHU_NOT_IMPLEMENTED
    virtual void normWeighted(const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const =0;
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Compute mean (average) value of each vector in multi-vector.
#ifdef CTHULHU_TPETRA_ONLY
    virtual void meanValue(const Teuchos::ArrayView<Scalar> &means) const =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
#ifdef CTHULHU_NOT_IMPLEMENTED
    virtual void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &beta) =0;
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Element-wise multiply of a Vector A with a MultiVector B.
    /** Forms this = scalarThis * this + scalarAB * B @ A
     *  where @ denotes element-wise multiplication.
     *  B must be the same shape (size and num-vectors) as this, while
     *  A is the same size but a single vector (column).
     */
#ifdef CTHULHU_NOT_IMPLEMENTED
    void elementWiseMultiply(Scalar scalarAB, const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, Scalar scalarThis) =0;
#endif // CTHULHU_NOT_IMPLEMENTED
    //@} 

    //! @name Attribute access functions
    //@{ 

    //! Returns the number of vectors in the multi-vector.
#ifdef CTHULHU_TPETRA_ONLY
    virtual size_t getNumVectors() const =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
#ifdef CTHULHU_TPETRA_ONLY
    virtual size_t getLocalLength() const =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Returns the global vector length of vectors in the multi-vector.
#ifdef CTHULHU_TPETRA_ONLY
    virtual global_size_t getGlobalLength() const =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Returns the stride between vectors in the multi-vector (only meaningful if ConstantStride() is true). WARNING: this may vary from node to node.
#ifdef CTHULHU_TPETRA_ONLY
    virtual size_t getStride() const =0;
#endif // CTHULHU_TPETRA_ONLY

    //! Returns true if this multi-vector has constant stride between vectors. WARNING: This may vary from node to node.
#ifdef CTHULHU_TPETRA_ONLY
    virtual bool isConstantStride() const =0;
#endif // CTHULHU_TPETRA_ONLY

    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
#ifdef CTHULHU_TPETRA_ONLY
    virtual std::string description() const =0;
#endif // CTHULHU_TPETRA_ONLY

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
#ifdef CTHULHU_TPETRA_ONLY
    virtual void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const =0;
#endif // CTHULHU_TPETRA_ONLY

    //@}

  }; // class MultiVector

} // namespace Cthulhu


#endif // CTHULHU_MULTIVECTOR_DECL_HPP
