#ifndef CTHULHU_TPETRAMULTIVECTOR_DECL_HPP
#define CTHULHU_TPETRAMULTIVECTOR_DECL_HPP

#include "Cthulhu_MultiVector.hpp"

#include "Tpetra_MultiVector.hpp"

namespace Cthulhu {

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
  class TpetraMultiVector : public Cthulhu::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> {

    public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Basic TpetraMultiVector constuctor.
    TpetraMultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, size_t NumVectors, bool zeroOut=true) : vec_(rcp(Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, NumVectors, zeroOut))) {}

    //! TpetraMultiVector copy constructor.
    TpetraMultiVector(const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source) : vec_(rcp(Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(source))) {}

    //! Set multi-vector values from two-dimensional array using Teuchos memory management classes. (copy)
    /*! Post-condition: constantStride() == true */
    TpetraMultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, const Teuchos::ArrayView<const Scalar> &A, size_t LDA, size_t NumVectors) : vec_(rcp(Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, A, LDA, NumVectors))) {}

    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy)
    /*! Post-condition: constantStride() == true */
    TpetraMultiVector(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> > &ArrayOfPtrs, size_t NumVectors) : vec_(rcp(Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>(map, ArrayOfPtrs, NumVectors))) {}

    TpetraMultiVector(const Teuchos::RCP<const Tpetra::MultiVector<LocalOrdinal, GlobalOrdinal, Node> > &vec) : vec_(vec) {}

    //! TpetraMultiVector destructor.
    virtual ~TpetraMultiVector();

    //@}

    //! @name Post-construction modification routines
    //@{ 

    //! Replace current value at the specified (globalRow, vectorIndex) location with specified value.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
      */
    inline void replaceGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) { vec_->replaceGlobalValue(globalRow, vectorIndex, value); }

    //! Adds specified value to existing value at the specified (globalRow, vectorIndex) location.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
      */
    inline void sumIntoGlobalValue(GlobalOrdinal globalRow, size_t vectorIndex, const Scalar &value) { vec_->sumIntoGlobalValue(globalRow, vectorIndex, value); }

    //! Replace current value at the specified (myRow, vectorIndex) location with specified value.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
      */
    inline void replaceLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) { vec_->replaceLocalValue(myRow, vectorIndex, value); }

    //! Adds specified value to existing value at the specified (myRow, vectorIndex) location.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
      */
    inline void sumIntoLocalValue(LocalOrdinal myRow, size_t vectorIndex, const Scalar &value) { vec_->sumIntoLocalValue(myRow, vectorIndex, value); }

    //! Initialize all values in a multi-vector with specified value.
    inline void putScalar(const Scalar &value) { vec_->putScalar(value); }

    //! Set multi-vector values to random numbers.
    inline void randomize() { vec_->randomize(); }

    //! Replace the underlying Map with a compatible one.
    inline void replaceMap(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &map) { vec_->replaceMap(map); }

    //! Instruct a local (non-distributed) TpetraMultiVector to sum values across all nodes.
    inline void reduce() { vec_->reduce(); }

    //! = Operator.
    /*! \param In A - Multivector to copy
     */
    //TODO inline TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& operator=(const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &source) { return vec_->(source); }

    //@}

    //! @name Data Copy and View get methods
    /** These methods are used to get the data underlying the TpetraMultiVector. They return data in one of three forms: 
      - a TpetraMultiVector with a subset of the columns of the target TpetraMultiVector
      - a raw C pointer or array of raw C pointers
      - one of the Teuchos memory management classes
      Not all of these methods are valid for a particular TpetraMultiVector. For instance, calling a method that accesses a 
      view of the data in a 1-D format (i.e., get1dView) requires that the target TpetraMultiVector has constant stride.
     */
    //@{

    //! Returns a TpetraMultiVector with copies of selected columns.
    inline Teuchos::RCP<TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subCopy(const Teuchos::Range1D &colRng) const { return vec_->subCopy(colRng); }

    //! Returns a TpetraMultiVector with copies of selected columns.
    inline Teuchos::RCP<TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subCopy(const Teuchos::ArrayView<const size_t> &cols) const { return vec_->subCopy(cols); }

    //! Returns a const TpetraMultiVector with const views of selected columns.
    inline Teuchos::RCP<const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subView(const Teuchos::Range1D &colRng) const { return vec_->subView(colRng); }

    //! Returns a const TpetraMultiVector with const views of selected columns.
    inline Teuchos::RCP<const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subView(const Teuchos::ArrayView<const size_t> &cols) const { return vec_->subView(cols); }

    //! Returns a TpetraMultiVector with views of selected columns.
    inline Teuchos::RCP<TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subViewNonConst(const Teuchos::Range1D &colRng) { return vec_->subViewNonConst(colRng); }

    //! Returns a TpetraMultiVector with views of selected columns.
    inline Teuchos::RCP<TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > subViewNonConst(const Teuchos::ArrayView<const size_t> &cols) { return vec_->subViewNonConst(cols); }

    //! \brief Returns a const TpetraMultiVector view of a subset of rows.
    /** 
        Returns a const view of this TpetraMultiVector consisting of a subset of the rows, as specified by an offset and a sub-Map.

        \param In subMap - The row map for the new TpetraMultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre  <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
     */
    inline Teuchos::RCP<const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > offsetView(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &subMap, size_t offset) const { return vec_->offsetView(subMap, offset); }

    //! \brief Returns a non-const TpetraMultiVector view of a subset of rows.
    /** 
        Returns a non-const view of this TpetraMultiVector consisting of a subset of the rows, as specified by an offset and a sub-Map.

        \param In subMap - The row map for the new TpetraMultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre  <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
     */
    inline Teuchos::RCP<TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > offsetViewNonConst(const Teuchos::RCP<const Map<LocalOrdinal,GlobalOrdinal,Node> > &subMap, size_t offset) { return vec_->offsetViewNonConst(subMap, offset); }

    //! Const Vector access function.
    // inline Teuchos::RCP<const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getVector(size_t j) const { return vec_->getVector(j); }

    //! Vector access function.
    // inline Teuchos::RCP<Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > getVectorNonConst(size_t j) { return vec_->getVectorNonConst(j); }

    //! Const Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    inline Teuchos::ArrayRCP<const Scalar> getData(size_t j) const { return vec_->getData(j); }

    //! Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    inline Teuchos::ArrayRCP<Scalar> getDataNonConst(size_t j) { return vec_->getDataNonConst(j); }

    //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
    inline void get1dCopy(Teuchos::ArrayView<Scalar> A, size_t LDA) const { vec_->get1dCopy(A, LDA); }

    //! Return multi-vector values in user-provided array of pointers (using Teuchos memory management classes).
    inline void get2dCopy(Teuchos::ArrayView<const Teuchos::ArrayView<Scalar> > ArrayOfPtrs) const { vec_->get2dCopy(ArrayOfPtrs); }

    //! Return const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.
    inline Teuchos::ArrayRCP<const Scalar> get1dView() const { return vec_->get1dView(); }

    //! Return const persisting pointers to values.
    inline Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > get2dView() const { return vec_->get2dView(); }

    //! Return non-const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.  Teuchos::ArrayRCP<Scalar> get1dViewNonConst() { return vec_->(); }
    inline Teuchos::ArrayRCP<Scalar> get1dViewNonConst() { return vec_->get1dViewNonConst(); }

    //! Return non-const persisting pointers to values.
    inline Teuchos::ArrayRCP<Teuchos::ArrayRCP<Scalar> > get2dViewNonConst() { return vec_->get2dViewNonConst(); }

    //! Return a const reference to the underlying Kokkos::TpetraMultiVector object (advanced use only)
    inline const Kokkos::MultiVector<Scalar,Node> & getLocalMV() const { return vec_->getLocalMV(); }

    //! Return a non-const reference to the underlying Kokkos::TpetraMultiVector object (advanced use only)
    inline Kokkos::MultiVector<Scalar,Node> & getLocalMVNonConst() { return vec_->getLocalMVNonConst(); }

    //@}

    //! @name Mathematical methods
    //@{ 

    //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
    inline void dot(const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Teuchos::ArrayView<Scalar> &dots) const { vec_->dot(A, dots); }

    //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
    inline void abs(const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) { vec_->abs(A); }

    //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    inline void reciprocal(const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) { vec_->reciprocal(A); }

    //! Scale the current values of a multi-vector, this = alpha*this.
    inline void scale(const Scalar &alpha) { vec_->scale(alpha); }

    //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
    inline void scale(Teuchos::ArrayView<const Scalar> alpha) { vec_->scale(alpha); }

    //! Replace multi-vector values with scaled values of A, this = alpha*A.
    inline void scale(const Scalar &alpha, const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A) { vec_->scale(A); }

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    inline void update(const Scalar &alpha, const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Scalar &beta) { vec_->update(A, beta); }

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    inline void update(const Scalar &alpha, const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const Scalar &beta, const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &gamma) { vec_->update(A, beta); }

    //! Compute 1-norm of each vector in multi-vector.
    inline void norm1(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const { vec_->norm1(norms); }

    //! Compute 2-norm of each vector in multi-vector.
    inline void norm2(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const { vec_->norm2(norms); }

    //! Compute Inf-norm of each vector in multi-vector.
    inline void normInf(const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const { vec_->normInf(norms); }

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    inline void normWeighted(const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &weights, const Teuchos::ArrayView<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &norms) const { vec_->normWeighted(weights, norms); }

    //! Compute mean (average) value of each vector in multi-vector.
    inline void meanValue(const Teuchos::ArrayView<Scalar> &means) const { vec_->meanValue(means); }

    //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
    inline void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const Scalar &alpha, const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, const Scalar &beta) { vec_->multiply(transA, transB, alpha, A, B, beta); }

    //! Element-wise multiply of a Vector A with a TpetraMultiVector B.
    /** Forms this = scalarThis * this + scalarAB * B @ A
     *  where @ denotes element-wise multiplication.
     *  B must be the same shape (size and num-vectors) as this, while
     *  A is the same size but a single vector (column).
     */
    // inline void elementWiseMultiply(Scalar scalarAB, const Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &A, const TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &B, Scalar scalarThis) { vec_->elementWiseMultiply(scalarAB, A , B , scalarThis); }
    //@} 

    //! @name Attribute access functions
    //@{ 

    //! Returns the number of vectors in the multi-vector.
    inline size_t getNumVectors() const { return vec_->getNumVectors(); }

    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
    inline size_t getLocalLength() const { return vec_->getLocalLength(); }

    //! Returns the global vector length of vectors in the multi-vector.
    inline global_size_t getGlobalLength() const { return vec_->getGlobalLength(); }

    //! Returns the stride between vectors in the multi-vector (only meaningful if ConstantStride() is true). WARNING: this may vary from node to node.
    inline size_t getStride() const { return vec_->getStride(); }

    //! Returns true if this multi-vector has constant stride between vectors. WARNING: This may vary from node to node.
    inline bool isConstantStride() const { return vec_->isConstantStride(); }

    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    inline std::string description() const { return vec_->description(); }

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { vec_->describe(out, verbLevel); }

    //@}

    RCP< const Tpetra::MultiVector<LocalOrdinal, GlobalOrdinal, Node> > getTpetra_MultiVector() const { return vec_; }
    
  private:
    const RCP< const Tpetra::MultiVector<LocalOrdinal, GlobalOrdinal, Node> > vec_;


  }; // class TpetraMultiVector

  /** \brief Non-member function to create a TpetraMultiVector from a specified Map.
      \relates TpetraMultiVector
   */
//   template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
//   Teuchos::RCP< TpetraMultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
//   createTpetraMultiVector(const Teuchos::RCP< const Map<LocalOrdinal,GlobalOrdinal,Node> > &map, size_t numVectors) 


} // namespace Cthulhu


#endif // CTHULHU_MULTIVECTOR_DECL_HPP
