#ifndef CTHULHU_EPETRAMULTIVECTOR_HPP
#define CTHULHU_EPETRAMULTIVECTOR_HPP

#include "Cthulhu_EpetraConfigDefs.hpp"

#include "Cthulhu_MultiVector.hpp"
#include "Cthulhu_Vector.hpp"

#include "Cthulhu_EpetraMap.hpp" 
#include "Cthulhu_EpetraExport.hpp"
#include "Epetra_MultiVector.h"
#include "Cthulhu_CombineMode.hpp"

#include "Cthulhu_Trans.hpp"

namespace Cthulhu {

  // #ifndef DOXYGEN_SHOULD_SKIP_THIS
  //   // forward declaration of Vector, needed to prevent circular inclusions
  //   template<class S, class LO, class GO, class N> class Vector;
  // #endif

  class EpetraMultiVector : public virtual Cthulhu::MultiVector<double,int,int> {

  public:

    //! @name Constructor/Destructor Methods
    //@{ 

    //! Basic EpetraMultiVector constuctor.
    EpetraMultiVector(const Teuchos::RCP<const Map<int,int> > &map, size_t NumVectors, bool zeroOut=true) {
      
      CTHULHU_RCP_DYNAMIC_CAST(const EpetraMap, map, eMap, "Cthulhu::TpetraMultiVector constructors only accept Cthulhu::TpetraMap as input arguments.");
      vec_ = rcp(new Epetra_MultiVector(eMap->getEpetra_BlockMap(), NumVectors, zeroOut));
    }

    EpetraMultiVector(const Teuchos::RCP<Epetra_MultiVector> &vec) : vec_(vec) {  }

    //! EpetraMultiVector destructor.
    virtual ~EpetraMultiVector() {  }

    //@}

    //! @name Post-construction modification routines
    //@{ 

    //! Initialize all values in a multi-vector with specified value.
    void putScalar(const double &value) {  vec_->PutScalar(value); }

    //! Set multi-vector values to random numbers.
    void randomize() {  vec_->Random(); }

    //! Set seed for Random function.
    /** Note: this method does not exist in Tpetra interface. Added for MueLu. */
    void setSeed(unsigned int seed) {
      vec_->SetSeed(seed);
    }

    //@}

    //! @name Data Copy and View get methods
    //@{

    //! Const Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    Teuchos::ArrayRCP<const double> getData(size_t j) const;

    //! Local vector access function.
    //! View of the local values in a particular vector of this multi-vector.
    Teuchos::ArrayRCP<double> getDataNonConst(size_t j);

    //@}

    //! @name Mathematical methods
    //@{ 
    //! Computes dot product of each corresponding pair of vectors, dots[i] = this[i].dot(A[i])
    void dot(const MultiVector<double,int,int> &A, const Teuchos::ArrayView<double> &dots) const;

    //! Puts element-wise absolute values of input Multi-vector in target: A = abs(this)
    void abs(const MultiVector<double,int,int> &A) { 
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      vec_->Abs(*eA.getEpetra_MultiVector());
    }

    //! Puts element-wise reciprocal values of input Multi-vector in target, this(i,j) = 1/A(i,j).
    void reciprocal(const MultiVector<double,int,int> &A) { 
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      vec_->Reciprocal(*eA.getEpetra_MultiVector());
    }

    //! Scale the current values of a multi-vector, this = alpha*this.
    void scale(const double &alpha) {  vec_->Scale(alpha); }

    //! Update multi-vector values with scaled values of A, this = beta*this + alpha*A.
    void update(const double &alpha, const MultiVector<double,int,int> &A, const double &beta) { 
       
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      vec_->Update(alpha, *eA.getEpetra_MultiVector(), beta); 
    }

    //! Update multi-vector with scaled values of A and B, this = gamma*this + alpha*A + beta*B.
    void update(const double &alpha, const MultiVector<double,int,int> &A, const double &beta, const MultiVector<double,int,int> &B, const double &gamma) {
       
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, B, eB, "This Cthulhu::EpetraMultiVector method only accept Cthulhu::EpetraMultiVector as input arguments.");
      vec_->Update(alpha, *eA.getEpetra_MultiVector(), beta, *eB.getEpetra_MultiVector(), gamma); 
    }

    //! Compute 1-norm of each vector in multi-vector.
    void norm1(const Teuchos::ArrayView<Teuchos::ScalarTraits<double>::magnitudeType> &norms) const {  vec_->Norm1(norms.getRawPtr()); }

    //! Compute 2-norm of each vector in multi-vector.
    void norm2(const Teuchos::ArrayView<Teuchos::ScalarTraits<double>::magnitudeType> &norms) const {  vec_->Norm2(norms.getRawPtr()); }

    //! Compute Inf-norm of each vector in multi-vector.
    void normInf(const Teuchos::ArrayView<Teuchos::ScalarTraits<double>::magnitudeType> &norms) const {  vec_->NormInf(norms.getRawPtr()); }

    //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
    void normWeighted(const MultiVector<double,int,int> &weights, const Teuchos::ArrayView<Teuchos::ScalarTraits<double>::magnitudeType> &norms) const;

    //! Compute mean (average) value of each vector in multi-vector.
    void meanValue(const Teuchos::ArrayView<double> &means) const;

    //! Matrix-Matrix multiplication, this = beta*this + alpha*op(A)*op(B).
    void multiply(Teuchos::ETransp transA, Teuchos::ETransp transB, const double &alpha, const MultiVector<double,int,int> &A, const MultiVector<double,int,int> &B, const double &beta) { 
       

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
    void elementWiseMultiply(double scalarAB, const Vector<double,int,int> &A, const MultiVector<double,int,int> &B, double scalarThis) {
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, A, eA, "Cthulhu::EpetraMultiVector->elementWiseMultiply() only accepts Cthulhu::EpetraMultiVector as input arguments.");
      CTHULHU_DYNAMIC_CAST(const EpetraMultiVector, B, eB, "Cthulhu::EpetraMultiVector->elementWiseMultiply() only accepts Cthulhu::EpetraMultiVector as input arguments.");
      vec_->Multiply(scalarAB, *eA.getEpetra_MultiVector() , *eB.getEpetra_MultiVector() , scalarThis);
    }
    //@} 

    //! @name Attribute access functions
    //@{ 

    //! Returns the number of vectors in the multi-vector.
    size_t getNumVectors() const {  return vec_->NumVectors(); }

    //! Returns the local vector length on the calling processor of vectors in the multi-vector.
    size_t getLocalLength() const {  return vec_->MyLength(); }

    //! Returns the global vector length of vectors in the multi-vector.
    global_size_t getGlobalLength() const {  return vec_->GlobalLength(); }

    //@} 

    //! @name Overridden from Teuchos::Describable 
    //@{

    //! Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const;

    //@}

    RCP< Epetra_MultiVector > getEpetra_MultiVector() const {  return vec_; }

    //@{
    // Implements DistObject interface

    const Teuchos::RCP<const Map<int,int> > getMap() const { 
       
      
      RCP<const Epetra_BlockMap> map = rcp(new Epetra_BlockMap(vec_->Map()));
      return rcp ( new Cthulhu::EpetraMap(map) );
    }

    void doImport(const DistObject<double, int, int> &source, const Import<int, int> &importer, CombineMode CM) ;
    void doExport(const DistObject<double, int, int> &dest, const Import<int, int>& importer, CombineMode CM) ;
    void doImport(const DistObject<double,int,int> &source, const Export<int, int>& exporter, CombineMode CM) ;
    void doExport(const DistObject<double, int, int> &dest, const Export<int, int>& exporter, CombineMode CM) ;

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
