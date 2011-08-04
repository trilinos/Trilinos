#ifndef CTHULHU_EPETRAINTVECTOR_HPP
#define CTHULHU_EPETRAINTVECTOR_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_EPETRA
#error This file should be included only if HAVE_CTHULHU_EPETRA is defined.
#endif

#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_MultiVector.hpp"
#include "Cthulhu_Vector.hpp"
#include "Cthulhu_Exceptions.hpp"

#include "Cthulhu_EpetraMap.hpp"
#include "Cthulhu_EpetraMultiVector.hpp"
#include "Epetra_IntVector.h"

namespace Cthulhu {

  //! \brief A class for constructing and using dense, distributors vectors.
  class EpetraIntVector : public Vector<int,int,int> {
    
  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Sets all vector entries to zero.
    explicit EpetraIntVector(const Teuchos::RCP<const Map<int,int> > &map, bool zeroOut=true) 
      
    {
      CTHULHU_DEBUG_ME;
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, map, eMap, "Cthulhu::EpetraCrsMatrix constructors only accept Cthulhu::EpetraMap as input arguments.");
      vec_ = rcp(new Epetra_IntVector(eMap->getEpetra_BlockMap(), zeroOut));
    }
    
#ifdef CTHULHU_NOT_IMPLEMENTED
    //! EpetraIntVector copy constructor.
    EpetraIntVector(const EpetraIntVector &source);
#endif

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! \brief Set multi-vector values from an array using Teuchos memory management classes. (copy)
    EpetraIntVector(const Teuchos::RCP<const Map<int,int> > &map, const Teuchos::ArrayView<const int> &A);
#endif

    //! Destructor.  
    inline ~EpetraIntVector() { CTHULHU_DEBUG_ME; };

    //@}

    //! @name Post-construction modification routines
    //@{ 
#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Replace current value at the specified location with specified value.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void replaceGlobalValue(int globalRow, const int &value) { CTHULHU_DEBUG_ME; this->EpetraMultiVector::replaceGlobalValue(globalRow,0,value); };

    //! Adds specified value to existing value at the specified location.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void sumIntoGlobalValue(int globalRow, const int &value) { CTHULHU_DEBUG_ME; this->EpetraMultiVector::sumIntoGlobalValue(globalRow, 0, value); };

    //! Replace current value at the specified location with specified values.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void replaceLocalValue(int myRow, const int &value) { CTHULHU_DEBUG_ME; this->EpetraMultiVector::replaceLocalValue(myRow, 0, value); };

    //! Adds specified value to existing value at the specified location.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void sumIntoLocalValue(int myRow, const int &value) { CTHULHU_DEBUG_ME; this->EpetraMultiVector::sumIntoLocalValue(myRow, 0, value); };
#endif
    //@}

#ifdef CTHULHU_TODO
    //! @name Extraction methods
    //@{
    //using MultiVector<int,int,int>::get1dCopy; // overloading, not hiding
    //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
    void get1dCopy(Teuchos::ArrayView<int> A) const;
    //@}
#endif

    //! @name Mathematical methods
    //@{ 
    //using EpetraMultiVector::dot; // overloading, not hiding
    //! Computes dot product of this Vector against input Vector x.
    int dot(const Vector<int,int,int> &a) const { 
      CTHULHU_DEBUG_ME; 
      TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO");
      return -1;
    }

    //using EpetraMultiVector::norm1; // overloading, not hiding
    //! Return 1-norm of this Vector.
    Teuchos::ScalarTraits<int>::magnitudeType norm1() const { CTHULHU_DEBUG_ME; TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO"); return -1; }

    //using EpetraMultiVector::norm2; // overloading, not hiding
    //! Compute 2-norm of this Vector.
    Teuchos::ScalarTraits<int>::magnitudeType norm2() const { CTHULHU_DEBUG_ME; TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO"); return -1; }

    //using EpetraMultiVector::normInf; // overloading, not hiding
    //! Compute Inf-norm of this Vector.
    Teuchos::ScalarTraits<int>::magnitudeType normInf() const { CTHULHU_DEBUG_ME; TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO"); return -1; }

    //using EpetraMultiVector::normWeighted; // overloading, not hiding
    //! Compute Weighted 2-norm (RMS Norm) of this Vector.
    Teuchos::ScalarTraits<int>::magnitudeType normWeighted(const Vector<int,int,int> &weights) const { CTHULHU_DEBUG_ME;
      TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO"); return -1;
    }

    //using EpetraMultiVector::meanValue; // overloading, not hiding
    //! Compute mean (average) value of this Vector.
    int meanValue() const { CTHULHU_DEBUG_ME; TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO"); return -1; }

    //using EpetraMultiVector::maxValue; // overloading, not hiding
    //! Compute max value of this Vector.
    int maxValue() const { CTHULHU_DEBUG_ME; TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO"); return -1; }
    //@} 

    //@}

    // protected:

    //     //! Advanced constructor accepting parallel buffer view.
    //     Vector(const Teuchos::RCP<const Map<int,int> > &map, Teuchos::ArrayRCP<double> data) { CTHULHU_DEBUG_ME; vec_->(); };

    //TODO wrap RCP etc.   RCP<Epetra_Vector> getEpetra_Vector() const { CTHULHU_DEBUG_ME; this->EpetraMultiVector::getEpetra_MultiVector()->getVector(0); }

    // **********************************
    // **********************************
    // **********************************
    // **********************************
    // IMPLEMENTING MULTIVECTOR INTERFACE
    // **********************************
    // **********************************
    // **********************************
    // **********************************

    //! @name Post-construction modification routines
    //@{ 

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Replace current value at the specified (globalRow, vectorIndex) location with specified value.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void replaceGlobalValue(int globalRow, size_t vectorIndex, const int &value) { CTHULHU_DEBUG_ME; vec_->replaceGlobalValue(globalRow, vectorIndex, value); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Adds specified value to existing value at the specified (globalRow, vectorIndex) location.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void sumIntoGlobalValue(int globalRow, size_t vectorIndex, const int &value) { CTHULHU_DEBUG_ME; vec_->sumIntoGlobalValue(globalRow, vectorIndex, value); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Replace current value at the specified (myRow, vectorIndex) location with specified value.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void replaceLocalValue(int myRow, size_t vectorIndex, const int &value) { CTHULHU_DEBUG_ME; vec_->replaceLocalValue(myRow, vectorIndex, value); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Adds specified value to existing value at the specified (myRow, vectorIndex) location.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void sumIntoLocalValue(int myRow, size_t vectorIndex, const int &value) { CTHULHU_DEBUG_ME; vec_->sumIntoLocalValue(myRow, vectorIndex, value); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! Initialize all values in a multi-vector with specified value.
    inline void putScalar(const int &value) { CTHULHU_DEBUG_ME; vec_->PutValue(value); }

    //! Set multi-vector values to random numbers.
    inline void randomize() { CTHULHU_DEBUG_ME;  TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "Cthulhu::EpetraIntVector::randomize(): Functionnality not available in Epetra"); }

    //! Set seed for Random function.
    /** Note: this method does not exist in Tpetra interface. Added for MueLu. */
    inline void setSeed(unsigned int seed) {
      TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "Cthulhu::EpetraIntVector::setSeed(): Functionnality not available in Epetra");
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Replace the underlying Map with a compatible one.
    inline void replaceMap(const Teuchos::RCP<const Map<int,int> > &map) { CTHULHU_DEBUG_ME; vec_->replaceMap(map); }
#endif // CTHULHU_NOT_IMPLEMENTED

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Instruct a local (non-distributed) EpetraMultiVector to sum values across all nodes.
    inline void reduce() { CTHULHU_DEBUG_ME; vec_->reduce(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! = Operator.
    /*! \param In A - Multivector to copy
     */
    inline EpetraMultiVector<int,int,int>& operator=(const EpetraMultiVector<int,int,int> &source) { CTHULHU_DEBUG_ME; return vec_->(source); }
#endif // CTHULHU_NOT_IMPLEMENTED

    //@}

    //! @name Data Copy and View get methods
    /** These methods are used to get the data underlying the EpetraMultiVector. They return data in one of three forms: 
        - a EpetraMultiVector with a subset of the columns of the target EpetraMultiVector
        - a raw C pointer or array of raw C pointers
        - one of the Teuchos memory management classes
        Not all of these methods are valid for a particular EpetraMultiVector. For instance, calling a method that accesses a 
        view of the data in a 1-D format (i.e., get1dView) requires that the target EpetraMultiVector has constant stride.
    */
    //@{
#ifdef CTHULHU_NOT_IMPLEMENTED

    //! Returns a MultiVector with copies of selected columns.
    inline Teuchos::RCP<MultiVector<int,int,int> > subCopy(const Teuchos::Range1D &colRng) const { CTHULHU_DEBUG_ME; return vec_->subCopy(colRng); }

    //! Returns a EpetraMultiVector with copies of selected columns.
    inline Teuchos::RCP<MultiVector<int,int,int> > subCopy(const Teuchos::ArrayView<const size_t> &cols) const { CTHULHU_DEBUG_ME; return vec_->subCopy(cols); }

    //! Returns a const MultiVector with const views of selected columns.
    inline Teuchos::RCP<const MultiVector<int,int,int> > subView(const Teuchos::Range1D &colRng) const { CTHULHU_DEBUG_ME; return vec_->subView(colRng); }

    //! Returns a const MultiVector with const views of selected columns.
    inline Teuchos::RCP<const MultiVector<int,int,int> > subView(const Teuchos::ArrayView<const size_t> &cols) const { CTHULHU_DEBUG_ME; return vec_->subView(cols); }

    //! Returns a MultiVector with views of selected columns.
    inline Teuchos::RCP<MultiVector<int,int,int> > subViewNonConst(const Teuchos::Range1D &colRng) { CTHULHU_DEBUG_ME; return vec_->subViewNonConst(colRng); }

    //! Returns a MultiVector with views of selected columns.
    inline Teuchos::RCP<MultiVector<int,int,int> > subViewNonConst(const Teuchos::ArrayView<const size_t> &cols) { CTHULHU_DEBUG_ME; return vec_->subViewNonConst(cols); }

    //! \brief Returns a const MultiVector view of a subset of rows.
    /** 
        Returns a const view of this MultiVector consisting of a subset of the rows, as specified by an offset and a sub-Map.

        \param In subMap - The row map for the new MultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre  <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
    */
    inline Teuchos::RCP<const MultiVector<int,int,int> > offsetView(const Teuchos::RCP<const Map<int,int> > &subMap, size_t offset) const { CTHULHU_DEBUG_ME; return vec_->offsetView(subMap, offset); }

    //! \brief Returns a non-const MultiVector view of a subset of rows.
    /** 
        Returns a non-const view of this MultiVector consisting of a subset of the rows, as specified by an offset and a sub-Map.

        \param In subMap - The row map for the new MultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre  <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
    */
    inline Teuchos::RCP<MultiVector<int,int,int> > offsetViewNonConst(const Teuchos::RCP<const Map<int,int> > &subMap, size_t offset) { CTHULHU_DEBUG_ME; return vec_->offsetViewNonConst(subMap, offset); }

    //! Const Vector access function.
    inline Teuchos::RCP<const Vector<int,int,int> > getVector(size_t j) const { CTHULHU_DEBUG_ME; return vec_->getVector(j); }

    //! Vector access function.
    inline Teuchos::RCP<Vector<int,int,int> > getVectorNonConst(size_t j) { CTHULHU_DEBUG_ME; return vec_->getVectorNonConst(j); }
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Const Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    inline Teuchos::ArrayRCP<const int> getData(size_t j) const { CTHULHU_DEBUG_ME; 
      CTHULHU_DEBUG_ME; 

      int * data = vec_->Values();
      int localLength = vec_->MyLength();
      
      return ArrayRCP<int>(data, 0, localLength, false); // not ownership
    }

    //! Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    inline Teuchos::ArrayRCP<int> getDataNonConst(size_t j) { 
      CTHULHU_DEBUG_ME; 

      int * data = vec_->Values();
      int localLength = vec_->MyLength();
      
      return ArrayRCP<int>(data, 0, localLength, false); // not ownership
    }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
    inline void get1dCopy(Teuchos::ArrayView<int> A, size_t LDA) const { CTHULHU_DEBUG_ME; vec_->get1dCopy(A, LDA); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return multi-vector values in user-provided array of pointers (using Teuchos memory management classes).
    inline void get2dCopy(Teuchos::ArrayView<const Teuchos::ArrayView<int> > ArrayOfPtrs) const { CTHULHU_DEBUG_ME; vec_->get2dCopy(ArrayOfPtrs); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.
    inline Teuchos::ArrayRCP<const int> get1dView() const { CTHULHU_DEBUG_ME; return vec_->get1dView(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return const persisting pointers to values.
    inline Teuchos::ArrayRCP<Teuchos::ArrayRCP<const int> > get2dView() const { CTHULHU_DEBUG_ME; return vec_->get2dView(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return non-const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.  Teuchos::ArrayRCP<int> get1dViewNonConst() { CTHULHU_DEBUG_ME; return vec_->(); }
    inline Teuchos::ArrayRCP<int> get1dViewNonConst() { CTHULHU_DEBUG_ME; return vec_->get1dViewNonConst(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return non-const persisting pointers to values.
    inline Teuchos::ArrayRCP<Teuchos::ArrayRCP<int> > get2dViewNonConst() { CTHULHU_DEBUG_ME; return vec_->get2dViewNonConst(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return a const reference to the underlying Kokkos::MultiVector object (advanced use only)
    inline const Kokkos::MultiVector<int> & getLocalMV() const { CTHULHU_DEBUG_ME; return vec_->getLocalMV(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return a non-const reference to the underlying Kokkos::MultiVector object (advanced use only)
    inline Kokkos::MultiVector<int> & getLocalMVNonConst() { CTHULHU_DEBUG_ME; return vec_->getLocalMVNonConst(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //@}

    //! @name Mathematical methods
    //@{ 
    //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
    inline void dot(const MultiVector<int,int,int> &A, const Teuchos::ArrayView<int> &dots) const { 
      CTHULHU_DEBUG_ME; 
      //CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO");
    }

    //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
    inline void abs(const MultiVector<int,int,int> &A) { 
      CTHULHU_DEBUG_ME; 
      //CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");

      TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO");
    }

    //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    inline void reciprocal(const MultiVector<int,int,int> &A) { 
      CTHULHU_DEBUG_ME; 
      //CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");

      TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO");
    }

    //! Scale the current values of a multi-vector, this = alpha*this.
    inline void scale(const int &alpha) { CTHULHU_DEBUG_ME; TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO"); }

//     //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
//     inline void scale(Teuchos::ArrayView<const int> alpha) { CTHULHU_DEBUG_ME; vec_->Scale(alpha.getRawPtr()); }

//     //! Replace multi-vector values with scaled values of A, this = alpha*A.
//     inline void scale(const int &alpha, const MultiVector<int,int,int> &A) { 
//       CTHULHU_DEBUG_ME; 
//       CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
//       vec_->Scale(*eA.getEpetra_MultiVector()); 
//     }

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    inline void update(const int &alpha, const MultiVector<int,int,int> &A, const int &beta) { 
      CTHULHU_DEBUG_ME; 
      // CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO");
    }

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    inline void update(const int &alpha, const MultiVector<int,int,int> &A, const int &beta, const MultiVector<int,int,int> &B, const int &gamma) {
      CTHULHU_DEBUG_ME; 
      //CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      //CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, B, eB, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO");
    }

    //! Compute 1-norm of each vector in multi-vector.
    inline void norm1(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const { CTHULHU_DEBUG_ME; TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO"); }

    //! Compute 2-norm of each vector in multi-vector.
    inline void norm2(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const { CTHULHU_DEBUG_ME; TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO"); }

    //! Compute Inf-norm of each vector in multi-vector.
    inline void normInf(const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const { CTHULHU_DEBUG_ME; TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO"); }

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    inline void normWeighted(const MultiVector<int,int,int> &weights, const Teuchos::ArrayView<Teuchos::ScalarTraits<int>::magnitudeType> &norms) const { 
      CTHULHU_DEBUG_ME; 
      TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO"); 
    }

    //! Compute mean (average) value of each vector in multi-vector.
    inline void meanValue(const Teuchos::ArrayView<int> &means) const { CTHULHU_DEBUG_ME;       TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO"); }

    //! Compute max value of each vector in multi-vector.
    inline void maxValue(const Teuchos::ArrayView<int> &maxs) const { CTHULHU_DEBUG_ME;       TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO"); }

    //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
    inline void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const int &alpha, const MultiVector<int,int,int> &A, const MultiVector<int,int,int> &B, const int &beta) { 
      CTHULHU_DEBUG_ME; 
      TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "Not available in Epetra"); 
    }


    //! Element-wise multiply of a Vector A with a EpetraMultiVector B.
    /** Forms this = scalarThis * this + scalarAB * B @ A
     *  where @ denotes element-wise multiplication.
     *  B must be the same shape (size and num-vectors) as this, while
     *  A is the same size but a single vector (column).
     */
    inline void elementWiseMultiply(int scalarAB, const Vector<int,int,int> &A, const MultiVector<int,int,int> &B, int scalarThis) {
      CTHULHU_DEBUG_ME;
      TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "Cthulhu_EpetraIntVector: elementWiseMultiply not implemented because Epetra_IntVector does not support this operation");

    }
    //@} 

    //! @name Attribute access functions
    //@{ 

    //! Returns the number of vectors in the multi-vector.
    inline size_t getNumVectors() const { CTHULHU_DEBUG_ME; return 1; }

    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
    inline size_t getLocalLength() const { CTHULHU_DEBUG_ME; return vec_->MyLength(); }

    //! Returns the global vector length of vectors in the multi-vector.
    inline global_size_t getGlobalLength() const { CTHULHU_DEBUG_ME; return vec_->GlobalLength(); }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns the stride between vectors in the multi-vector (only meaningful if ConstantStride() is true). WARNING: this may vary from node to node.
    inline size_t getStride() const { CTHULHU_DEBUG_ME; return vec_->getStride(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns true if this multi-vector has constant stride between vectors. WARNING: This may vary from node to node.
    inline bool isConstantStride() const { CTHULHU_DEBUG_ME; return vec_->isConstantStride(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    inline std::string description() const { CTHULHU_DEBUG_ME; 
      CTHULHU_DEBUG_ME; 

      // This implementation come from Epetra_Vector_def.hpp (without modification)
      std::ostringstream oss;
      oss << Teuchos::Describable::description();
      oss << "{length="<<this->getGlobalLength()
      << "}";
      return oss.str();
    }

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { 
      CTHULHU_DEBUG_ME; 
      
      typedef Kokkos::MultiVector<double> KMV;
      typedef Kokkos::DefaultArithmetic<KMV>   MVT;

      // This implementation come from Tpetra_Vector_def.hpp (without modification) // JG: true?
      using std::endl;
      using std::setw;
      using Teuchos::VERB_DEFAULT;
      using Teuchos::VERB_NONE;
      using Teuchos::VERB_LOW;
      using Teuchos::VERB_MEDIUM;
      using Teuchos::VERB_HIGH;
      using Teuchos::VERB_EXTREME;

      if (verbLevel > Teuchos::VERB_NONE)
        vec_->Print(out);
    }

    //@}

    RCP< Epetra_IntVector > getEpetra_IntVector() const { CTHULHU_DEBUG_ME; return vec_; }

    // From DistObject
    const Teuchos::RCP<const Map<int,int> > getMap() const { 
      CTHULHU_DEBUG_ME; 
      
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(vec_->Map()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    const RCP<const Comm<int> > getComm() const {
      CTHULHU_DEBUG_ME; 
      
      TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO getComm Epetra MultiVector not implemented");
    }

    // end of "Implementing Epetra interface"

    // Implementing DistObject
    
    inline void doImport(const Vector<int, int, int> &source, 
                         const Import<int, int> &importer, CombineMode CM) {
      CTHULHU_DEBUG_ME;

      CTHULHU_DYNAMIC_CAST(const EpetraIntVector, source, tSource, "Cthulhu::EpetraIntVector::doImport only accept Cthulhu::EpetraIntVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Cthulhu::EpetraIntVector::doImport only accept Cthulhu::EpetraImport as input arguments.");

      const Epetra_IntVector & v = *tSource.getEpetra_IntVector();
      int err = vec_->Import(v, *tImporter.getEpetra_Import(), Cthulhu2Epetra_CombineMode(CM)); 
      TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
    }

    void doExport(const Vector<int, int, int> &dest,
                  const Import<int, int>& importer, CombineMode CM) {
      CTHULHU_DEBUG_ME;

      CTHULHU_DYNAMIC_CAST(const EpetraIntVector, dest, tDest, "Cthulhu::EpetraIntVector::doImport only accept Cthulhu::EpetraIntVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Cthulhu::EpetraIntVector::doImport only accept Cthulhu::EpetraImport as input arguments.");

      const Epetra_IntVector & v = *tDest.getEpetra_IntVector();
      int err = vec_->Import(v, *tImporter.getEpetra_Import(), Cthulhu2Epetra_CombineMode(CM)); 
      TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
    }

    void doExport(const Vector<int, int, int> &dest,
                  const Export<int, int>& exporter, CombineMode CM) {
      CTHULHU_DEBUG_ME;

      CTHULHU_DYNAMIC_CAST(const EpetraIntVector, dest, tDest, "Cthulhu::EpetraIntVector::doImport only accept Cthulhu::EpetraIntVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Cthulhu::EpetraIntVector::doImport only accept Cthulhu::EpetraImport as input arguments.");

      const Epetra_IntVector & v = *tDest.getEpetra_IntVector();
      int err = vec_->Export(v, *tExporter.getEpetra_Export(), Cthulhu2Epetra_CombineMode(CM)); 
      TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
    }


    inline void doImport(const MultiVector<int, int, int> &source, 
                         const Import<int, int> &importer, CombineMode CM) {

      TEST_FOR_EXCEPTION(1, std::runtime_error, "NOT IMPLEMENTED");
    }

    void doExport(const MultiVector<int, int, int> &dest,
                  const Import<int, int>& importer, CombineMode CM) {
 
      TEST_FOR_EXCEPTION(1, std::runtime_error, "NOT IMPLEMENTED");
    }


    void doExport(const MultiVector<int, int, int> &dest,
                  const Export<int, int>& exporter, CombineMode CM) {
 
      TEST_FOR_EXCEPTION(1, std::runtime_error, "NOT IMPLEMENTED");
    }


  private:
    RCP< Epetra_IntVector > vec_;
    
  }; // class EpetraIntVector

#ifdef CTHULHU_NOT_IMPLEMENTED
  /** \brief Non-member function to create a Vector from a specified Map.
  
  \relates Vector
  */
  template <class int, class int, class int, class Node>
  Teuchos::RCP< Vector<int,int,int> >
  createVector(const Teuchos::RCP< const Map<int,int> > &map);
#endif
} // namespace Cthulhu

#endif // CTHULHU_VECTOR_DECL_HPP
