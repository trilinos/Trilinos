#ifndef CTHULHU_EPETRAMULTIVECTOR_HPP
#define CTHULHU_EPETRAMULTIVECTOR_HPP

#include "Cthulhu_ConfigDefs.hpp"

#ifndef HAVE_CTHULHU_EPETRA
#error This file should be included only if HAVE_CTHULHU_EPETRA is defined.
#endif

#include "Cthulhu_MultiVector.hpp"

#include "Cthulhu_EpetraMap.hpp" 
#include "Cthulhu_EpetraExport.hpp"
#include "Epetra_MultiVector.h"


#include "Cthulhu_Trans.hpp"

namespace Cthulhu {

  // #ifndef DOXYGEN_SHOULD_SKIP_THIS
  //   // forward declaration of Vector, needed to prevent circular inclusions
  //   template<class S, class LO, class GO, class N> class Vector;
  // #endif

  //! \brief A class for constructing and using dense, distributors multivectors.
  /*!
    This class is templated on \c double, \c int and \c GlobalOrdinal. 
    The \c int type, if omitted, defaults to \c int. The \c GlobalOrdinal 
    type, if omitted, defaults to the \c int type.
  */
  class EpetraMultiVector : public virtual Cthulhu::MultiVector<double,int,int> {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 
    //! Basic EpetraMultiVector constuctor.
    EpetraMultiVector(const Teuchos::RCP<const Map<int,int> > &map, size_t NumVectors, bool zeroOut=true) {
      
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, map, eMap, "Cthulhu::TpetraMultiVector constructors only accept Cthulhu::TpetraMap as input arguments.");
      vec_ = rcp(new Epetra_MultiVector(eMap->getEpetra_BlockMap(), NumVectors, zeroOut));
    }

#ifdef CTHULHU_NOT_IMPLEMENTED    
    //! EpetraMultiVector copy constructor.
    EpetraMultiVector(const EpetraMultiVector<double,int,int> &source){  } 
#endif // CTHULHU_NOT_IMPLEMENTED

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Set multi-vector values from two-dimensional array using Teuchos memory management classes. (copy)
    /*! Post-condition: constantStride() == true */
    EpetraMultiVector(const Teuchos::RCP<const Map<int,int> > &map, const Teuchos::ArrayView<const double> &A, size_t LDA, size_t NumVectors) {  } 
#endif // CTHULHU_NOT_IMPLEMENTED

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Set multi-vector values from array of pointers using Teuchos memory management classes. (copy)
    /*! Post-condition: constantStride() == true */
    EpetraMultiVector(const Teuchos::RCP<const Map<int,int> > &map, const Teuchos::ArrayView<const Teuchos::ArrayView<const double> > &ArrayOfPtrs, size_t NumVectors) {  } 
#endif // CTHULHU_NOT_IMPLEMENTED
  
    EpetraMultiVector(const Teuchos::RCP<Epetra_MultiVector> &vec) : vec_(vec) {  }

    //! EpetraMultiVector destructor.
    virtual ~EpetraMultiVector() {  }

    //@}

    //! @name Post-construction modification routines
    //@{ 

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Replace current value at the specified (globalRow, vectorIndex) location with specified value.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void replaceGlobalValue(int globalRow, size_t vectorIndex, const double &value) {  vec_->replaceGlobalValue(globalRow, vectorIndex, value); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Adds specified value to existing value at the specified (globalRow, vectorIndex) location.
    /** \pre \c globalRow must be a valid global element on this node, according to the row map.
     */
    inline void sumIntoGlobalValue(int globalRow, size_t vectorIndex, const double &value) {  vec_->sumIntoGlobalValue(globalRow, vectorIndex, value); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Replace current value at the specified (myRow, vectorIndex) location with specified value.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void replaceLocalValue(int myRow, size_t vectorIndex, const double &value) {  vec_->replaceLocalValue(myRow, vectorIndex, value); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Adds specified value to existing value at the specified (myRow, vectorIndex) location.
    /** \pre \c localRow must be a valid local element on this node, according to the row map.
     */
    inline void sumIntoLocalValue(int myRow, size_t vectorIndex, const double &value) {  vec_->sumIntoLocalValue(myRow, vectorIndex, value); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //! Initialize all values in a multi-vector with specified value.
    inline void putScalar(const double &value) {  vec_->PutScalar(value); }

    //! Set multi-vector values to random numbers.
    inline void randomize() {  vec_->Random(); }

    //! Set seed for Random function.
    /** Note: this method does not exist in Tpetra interface. Added for MueLu. */
    inline void setSeed(unsigned int seed) {
      vec_->SetSeed(seed);
    }

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! Replace the underlying Map with a compatible one.
    inline void replaceMap(const Teuchos::RCP<const Map<int,int> > &map) {  vec_->replaceMap(map); }
#endif // CTHULHU_NOT_IMPLEMENTED

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Instruct a local (non-distributed) EpetraMultiVector to sum values across all nodes.
    inline void reduce() {  vec_->reduce(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED
    //! = Operator.
    /*! \param In A - Multivector to copy
     */
    inline EpetraMultiVector<double,int,int>& operator=(const EpetraMultiVector<double,int,int> &source) {  return vec_->(source); }
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
    inline Teuchos::RCP<MultiVector<double,int,int> > subCopy(const Teuchos::Range1D &colRng) const {  return vec_->subCopy(colRng); }

    //! Returns a EpetraMultiVector with copies of selected columns.
    inline Teuchos::RCP<MultiVector<double,int,int> > subCopy(const Teuchos::ArrayView<const size_t> &cols) const {  return vec_->subCopy(cols); }

    //! Returns a const MultiVector with const views of selected columns.
    inline Teuchos::RCP<const MultiVector<double,int,int> > subView(const Teuchos::Range1D &colRng) const {  return vec_->subView(colRng); }

    //! Returns a const MultiVector with const views of selected columns.
    inline Teuchos::RCP<const MultiVector<double,int,int> > subView(const Teuchos::ArrayView<const size_t> &cols) const {  return vec_->subView(cols); }

    //! Returns a MultiVector with views of selected columns.
    inline Teuchos::RCP<MultiVector<double,int,int> > subViewNonConst(const Teuchos::Range1D &colRng) {  return vec_->subViewNonConst(colRng); }

    //! Returns a MultiVector with views of selected columns.
    inline Teuchos::RCP<MultiVector<double,int,int> > subViewNonConst(const Teuchos::ArrayView<const size_t> &cols) {  return vec_->subViewNonConst(cols); }

    //! \brief Returns a const MultiVector view of a subset of rows.
    /** 
        Returns a const view of this MultiVector consisting of a subset of the rows, as specified by an offset and a sub-Map.

        \param In subMap - The row map for the new MultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre  <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
    */
    inline Teuchos::RCP<const MultiVector<double,int,int> > offsetView(const Teuchos::RCP<const Map<int,int> > &subMap, size_t offset) const {  return vec_->offsetView(subMap, offset); }

    //! \brief Returns a non-const MultiVector view of a subset of rows.
    /** 
        Returns a non-const view of this MultiVector consisting of a subset of the rows, as specified by an offset and a sub-Map.

        \param In subMap - The row map for the new MultiVector.
        \param In offset - The offset into the data of <tt>(*this)</tt>.

        \pre  <tt>subMap->getNodeNumElements() + offset < this->getLocalLength()</tt>
    */
    inline Teuchos::RCP<MultiVector<double,int,int> > offsetViewNonConst(const Teuchos::RCP<const Map<int,int> > &subMap, size_t offset) {  return vec_->offsetViewNonConst(subMap, offset); }

    //! Const Vector access function.
    inline Teuchos::RCP<const Vector<double,int,int> > getVector(size_t j) const {  return vec_->getVector(j); }

    //! Vector access function.
    inline Teuchos::RCP<Vector<double,int,int> > getVectorNonConst(size_t j) {  return vec_->getVectorNonConst(j); }
#endif // CTHULHU_NOT_IMPLEMENTED

    //! Const Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    inline Teuchos::ArrayRCP<const double> getData(size_t j) const { 
       

      double ** arrayOfPointers;
      
      vec_->ExtractView(&arrayOfPointers);
     
      double * data = arrayOfPointers[j];
      int localLength = vec_->MyLength();
      
      return ArrayRCP<double>(data, 0, localLength, false); // not ownership
    }

    //! Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    inline Teuchos::ArrayRCP<double> getDataNonConst(size_t j) { 
       

      double ** arrayOfPointers;
      
      vec_->ExtractView(&arrayOfPointers);
     
      double * data = arrayOfPointers[j];
      int localLength = vec_->MyLength();
      
      return ArrayRCP<double>(data, 0, localLength, false); // not ownership
    }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return multi-vector values in user-provided two-dimensional array (using Teuchos memory management classes).
    inline void get1dCopy(Teuchos::ArrayView<double> A, size_t LDA) const {  vec_->get1dCopy(A, LDA); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return multi-vector values in user-provided array of pointers (using Teuchos memory management classes).
    inline void get2dCopy(Teuchos::ArrayView<const Teuchos::ArrayView<double> > ArrayOfPtrs) const {  vec_->get2dCopy(ArrayOfPtrs); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.
    inline Teuchos::ArrayRCP<const double> get1dView() const {  return vec_->get1dView(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return const persisting pointers to values.
    inline Teuchos::ArrayRCP<Teuchos::ArrayRCP<const double> > get2dView() const {  return vec_->get2dView(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return non-const persisting view of values in a one-dimensional array. Throws std::runtime_error if the underlying data is non-contiguous.  Teuchos::ArrayRCP<double> get1dViewNonConst() {  return vec_->(); }
    inline Teuchos::ArrayRCP<double> get1dViewNonConst() {  return vec_->get1dViewNonConst(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return non-const persisting pointers to values.
    inline Teuchos::ArrayRCP<Teuchos::ArrayRCP<double> > get2dViewNonConst() {  return vec_->get2dViewNonConst(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return a const reference to the underlying Kokkos::MultiVector object (advanced use only)
    inline const Kokkos::MultiVector<double> & getLocalMV() const {  return vec_->getLocalMV(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Return a non-const reference to the underlying Kokkos::MultiVector object (advanced use only)
    inline Kokkos::MultiVector<double> & getLocalMVNonConst() {  return vec_->getLocalMVNonConst(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //@}

    //! @name Mathematical methods
    //@{ 
    //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
    inline void dot(const MultiVector<double,int,int> &A, const Teuchos::ArrayView<double> &dots) const { 
       
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      vec_->Dot(*eA.getEpetra_MultiVector(), dots.getRawPtr());
    }

    //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
    inline void abs(const MultiVector<double,int,int> &A) { 
       
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      vec_->Abs(*eA.getEpetra_MultiVector());
    }

    //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    inline void reciprocal(const MultiVector<double,int,int> &A) { 
       
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      vec_->Reciprocal(*eA.getEpetra_MultiVector());
    }

    //! Scale the current values of a multi-vector, this = alpha*this.
    inline void scale(const double &alpha) {  vec_->Scale(alpha); }

//     //! Scale the current values of a multi-vector, this[j] = alpha[j]*this[j].
//     inline void scale(Teuchos::ArrayView<const double> alpha) {  vec_->Scale(alpha.getRawPtr()); }

//     //! Replace multi-vector values with scaled values of A, this = alpha*A.
//     inline void scale(const double &alpha, const MultiVector<double,int,int> &A) { 
//        
//       CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
//       vec_->Scale(*eA.getEpetra_MultiVector()); 
//     }

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    inline void update(const double &alpha, const MultiVector<double,int,int> &A, const double &beta) { 
       
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      vec_->Update(alpha, *eA.getEpetra_MultiVector(), beta); 
    }

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    inline void update(const double &alpha, const MultiVector<double,int,int> &A, const double &beta, const MultiVector<double,int,int> &B, const double &gamma) {
       
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, B, eB, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      vec_->Update(alpha, *eA.getEpetra_MultiVector(), beta, *eB.getEpetra_MultiVector(), gamma); 
    }

    //! Compute 1-norm of each vector in multi-vector.
    inline void norm1(const Teuchos::ArrayView<Teuchos::ScalarTraits<double>::magnitudeType> &norms) const {  vec_->Norm1(norms.getRawPtr()); }

    //! Compute 2-norm of each vector in multi-vector.
    inline void norm2(const Teuchos::ArrayView<Teuchos::ScalarTraits<double>::magnitudeType> &norms) const {  vec_->Norm2(norms.getRawPtr()); }

    //! Compute Inf-norm of each vector in multi-vector.
    inline void normInf(const Teuchos::ArrayView<Teuchos::ScalarTraits<double>::magnitudeType> &norms) const {  vec_->NormInf(norms.getRawPtr()); }

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    inline void normWeighted(const MultiVector<double,int,int> &weights, const Teuchos::ArrayView<Teuchos::ScalarTraits<double>::magnitudeType> &norms) const { 
       
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, weights, eWeights, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      vec_->NormWeighted(*eWeights.getEpetra_MultiVector(), norms.getRawPtr()); 
    }

    //! Compute mean (average) value of each vector in multi-vector.
    inline void meanValue(const Teuchos::ArrayView<double> &means) const {  vec_->MeanValue(means.getRawPtr()); } //TODO: modify ArrayView size ??

    //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
    inline void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const double &alpha, const MultiVector<double,int,int> &A, const MultiVector<double,int,int> &B, const double &beta) { 
       

      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "Cthulhu::EpetraMultiVector->multiply() only accept Cthulhu::EpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, B, eB, "Cthulhu::EpetraMultiVector->multiply() only accept Cthulhu::EpetraMultiVector as input arguments.");

      // TODO: Check if following exception useful here.
      TEST_FOR_EXCEPTION((transA != Teuchos::NO_TRANS) && (transA == Teuchos::TRANS), Cthulhu::Exceptions::NotImplemented, "Cthulhu::EpetraMultiVector->multiply() only accept transA == NO_TRANS or transA == TRANS");
      TEST_FOR_EXCEPTION((transB != Teuchos::NO_TRANS) && (transB == Teuchos::TRANS), Cthulhu::Exceptions::NotImplemented, "Cthulhu::EpetraMultiVector->multiply() only accept transB == NO_TRANS or transB == TRANS");
      bool eTransA = Teuchos2Epetra_Trans(transA);
      bool eTransB = Teuchos2Epetra_Trans(transB);

      CTHULHU_ERR_CHECK(vec_->Multiply(eTransA, eTransB, alpha, *eA.getEpetra_MultiVector(), *eB.getEpetra_MultiVector(), beta));
    }

    //! Element-wise multiply of a Vector A with a EpetraMultiVector B.
    /** Forms this = scalarThis * this + scalarAB * B @ A
     *  where @ denotes element-wise multiplication.
     *  B must be the same shape (size and num-vectors) as this, while
     *  A is the same size but a single vector (column).
     */
    inline void elementWiseMultiply(double scalarAB, const Vector<double,int,int> &A, const MultiVector<double,int,int> &B, double
scalarThis) {
      
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "Cthulhu::EpetraMultiVector->elementWiseMultiply() only accepts Cthulhu::EpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, B, eB, "Cthulhu::EpetraMultiVector->elementWiseMultiply() only accepts Cthulhu::EpetraMultiVector as input arguments.");
      vec_->Multiply(scalarAB, *eA.getEpetra_MultiVector() , *eB.getEpetra_MultiVector() , scalarThis);
    }
    //@} 

    //! @name Attribute access functions
    //@{ 

    //! Returns the number of vectors in the multi-vector.
    inline size_t getNumVectors() const {  return vec_->NumVectors(); }

    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
    inline size_t getLocalLength() const {  return vec_->MyLength(); }

    //! Returns the global vector length of vectors in the multi-vector.
    inline global_size_t getGlobalLength() const {  return vec_->GlobalLength(); }

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns the stride between vectors in the multi-vector (only meaningful if ConstantStride() is true). WARNING: this may vary from node to node.
    inline size_t getStride() const {  return vec_->getStride(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

#ifdef CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA
    //! Returns true if this multi-vector has constant stride between vectors. WARNING: This may vary from node to node.
    inline bool isConstantStride() const {  return vec_->isConstantStride(); }
#endif // CTHULHU_NOT_IMPLEMENTED_FOR_EPETRA

    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{

    /** \brief Return a simple one-line description of this object. */
    inline std::string description() const {  
      TEST_FOR_EXCEPTION(1, Cthulhu::Exceptions::NotImplemented, "TODO");
      return "TODO"; 
    }

    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    inline void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const { 
      
      vec_->Print(out);
    }

    //@}

    RCP< Epetra_MultiVector > getEpetra_MultiVector() const {  return vec_; }

    //@{
    // Implements DistObject interface

    const Teuchos::RCP<const Map<int,int> > getMap() const { 
       
      
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(vec_->Map()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    inline void doImport(const DistObject<double, int, int> &source, 
                         const Import<int, int> &importer, CombineMode CM) {
      

      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, source, tSource, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraImport as input arguments.");

      RCP<Epetra_MultiVector> v = tSource.getEpetra_MultiVector();
      int err = this->getEpetra_MultiVector()->Import(*v, *tImporter.getEpetra_Import(), toEpetra(CM));
      TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
    }

    void doExport(const DistObject<double, int, int> &dest,
                  const Import<int, int>& importer, CombineMode CM) {
      
      
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, dest, tDest, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraImport, importer, tImporter, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraImport as input arguments.");

      RCP<Epetra_MultiVector> v = tDest.getEpetra_MultiVector();
      int err = this->getEpetra_MultiVector()->Export(*v, *tImporter.getEpetra_Import(), toEpetra(CM)); 
      TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
    }

    void doImport(const DistObject<double,int,int> &source,
                  const Export<int, int>& exporter, CombineMode CM) {
      

      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, source, tSource, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraImport as input arguments.");

      RCP<Epetra_MultiVector> v = tSource.getEpetra_MultiVector();
      int err = this->getEpetra_MultiVector()->Import(*v, *tExporter.getEpetra_Export(), toEpetra(CM));
      TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
    }

    void doExport(const DistObject<double, int, int> &dest,
                  const Export<int, int>& exporter, CombineMode CM) {
      
      
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, dest, tDest, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraExport, exporter, tExporter, "Cthulhu::EpetraMultiVector::doImport only accept Cthulhu::EpetraImport as input arguments.");

      RCP<Epetra_MultiVector> v = tDest.getEpetra_MultiVector();
      int err = this->getEpetra_MultiVector()->Export(*v, *tExporter.getEpetra_Export(), toEpetra(CM)); 
      TEST_FOR_EXCEPTION(err != 0, std::runtime_error, "Catch error code returned by Epetra.");
    }

    //@}

  private:
    RCP< Epetra_MultiVector > vec_;


  }; // class EpetraMultiVector

#ifdef CTHULHU_NOT_IMPLEMENTED
  /** \brief Non-member function to create a EpetraMultiVector from a specified Map.
      \relates EpetraMultiVector
  */
  template <class double, class int, class int, class Node>
  Teuchos::RCP< EpetraMultiVector<double,int,int> >
  createEpetraMultiVector(const Teuchos::RCP< const Map<int,int> > &map, size_t numVectors) {  }
#endif // CTHULHU_NOT_IMPLEMENTED
  
} // namespace Cthulhu

#endif // ifndef CTHULHU_EPETRAMULTIVECTOR_HPP
