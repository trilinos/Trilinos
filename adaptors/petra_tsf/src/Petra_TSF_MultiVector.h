#ifndef _PETRA_TSF_MULTIVECTOR_H_
#define _PETRA_TSF_MULTIVECTOR_H_


#include "Petra_RDP_MultiVector.h"
#include "TSF_MultiVector.h"
namespace Petra_TSF {

//! Petra_TSF::Multivector:  The Petra implementation of the Trilinos Virtual MultiVector Class.
/*! The Petra_TSF::Multivector class is a wrapper that uses the Petra_RDP_MultiVector class to
    provide a real double precision implement the TSF_MultiVector pure virtual class.  
*/
template<class scalarType>
class MultiVector: public virtual TSF::MultiVector, public Petra_RDP_MultiVector {
    
  public:

  //@{ \name Constructors/Destructor 
  //! Basic MultiVector constuctor.
  /*! Creates a MultiVector object and fills with zero values.  

    \param In 
           Map - A Petra_LocalMap, Petra_Map or Petra_BlockMap.

	   \warning Note that, because Petra_LocalMap
	   derives from Petra_Map and Petra_Map derives from Petra_BlockMap, this constructor works
	   for all three types of Petra map classes.
    \param In 
           NumVectors - Number of vectors in multi-vector.

    \return Pointer to a MultiVector.

  */
  MultiVector(const Petra_BlockMap& Map, int NumVectors);

  //! MultiVector copy constructor.
  
  MultiVector(const Petra_TSF::MultiVector<scalarType>& Source);
  
  //! Set multi-vector values from two-dimensional array.
  /*!
    \param In 
           Petra_DataAccess - Enumerated type set to Copy or View.
    \param In 
           Map - A Petra_LocalMap, Petra_Map or Petra_BlockMap.
    \param In
           A - Pointer to an array of double precision numbers.  The first vector starts at A.
	   The second vector starts at A+MyLDA, the third at A+2*MyLDA, and so on.
    \param In
           MyLDA - The "Leading Dimension", or stride between vectors in memory.
	   \warning This value refers to the stride on the calling processor.  Thus it is a
	   local quantity, not a global quantity.
    \param In 
           NumVectors - Number of vectors in multi-vector.

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  MultiVector(Petra_DataAccess CV, const Petra_BlockMap& Map, 
			double *A, int MyLDA, int NumVectors);

  //! Set multi-vector values from array of pointers.
  /*!
    \param In 
           Petra_DataAccess - Enumerated type set to Copy or View.
    \param In 
           Map - A Petra_LocalMap, Petra_Map or Petra_BlockMap.
    \param In
           ArrayOfPointers - An array of pointers such that ArrayOfPointers[i] points to the memory
	   location containing ith vector to be copied.
    \param In 
           NumVectors - Number of vectors in multi-vector.

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  MultiVector(Petra_DataAccess CV, const Petra_BlockMap& Map, 
			 double **ArrayOfPointers, int NumVectors);

  //! Set multi-vector values from list of vectors in an existing Petra_RDP_MultiVector.
  /*!
    \param In 
           Petra_DataAccess - Enumerated type set to Copy or View.
    \param In
           Source - An existing fully constructed Petra_RDP_MultiVector.
    \param In
           Indices - Integer list of the vectors to copy.  
    \param In 
           NumVectors - Number of vectors in multi-vector.

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  MultiVector(Petra_DataAccess CV,  
			const Petra_TSF::MultiVector<scalarType>& Source, int *Indices, int NumVectors);

  //! Set multi-vector values from range of vectors in an existing Petra_RDP_MultiVector.
  /*!
    \param In 
           Petra_DataAccess - Enumerated type set to Copy or View.
    \param In
           Source - An existing fully constructed Petra_RDP_MultiVector.
    \param In
           StartIndex - First of the vectors to copy.  
    \param In 
           NumVectors - Number of vectors in multi-vector.

    \return Integer error code, set to 0 if successful.

	   See Detailed Description section for further discussion.
  */
  MultiVector(Petra_DataAccess CV, 
			const Petra_TSF::MultiVector<scalarType>& Source, int StartIndex, 
			int NumVectors);

  //! MultiVector Destructor.
  virtual ~MultiVector(void){};
  //@}

  //@{ \name Clone methods
  //! Create a new multivector of the same class as the \this multivector, fill with zeros.
  /*! cloneZero is essentially a call to the non-trivial default constructor of the multivector class of
      which the \this object is a member.  clone returns with a new multivector in that class, 
      with numVector columns.
    \param In
           numVectors - The number of vectors that will be in clone.
    \param Out
           clone - The new multivector clone.
	   
    \return Integer error code, set to 0 if successful, -1 if not able to create.
  */

  virtual int cloneZero ( MultiVector<scalarType> * & clone, const int numVectors){

    Petra_TSF::MultiVector * clone1 = new Petra_TSF::MultiVector(Map(), NumVectors); 
    if (clone1==0) PETRA_CHK_ERR(-1);
    clone = dynamic_cast<TSF::MultiVector *>(clone1);
    return(0);
  };

  //! Create a new multivector that is an exact full(deep) copy the \this multivector.
  /*! cloneCopy is essentially a call to the standard copy constructor of the multivector class of
      which the \this object is a member.  clone returns with a complete replication of \this.
    \param Out
           clone - The new multivector clone.
	   
    \return Integer error code, set to 0 if successful, -1 if not able to create.
  */

  virtual int cloneCopy ( MultiVector<scalarType> * & clone){

    Petra_TSF::MultiVector * clone1 = new Petra_TSF::MultiVector(*this); 
    if (clone1==0) PETRA_CHK_ERR(-1);
    clone = dynamic_cast<TSF::MultiVector *>(clone1);
    return(0);  
  };

  //! Create a new multivector that is a view (shallow copy) of all or some vectors in the \this multivector.
  /*! cloneView creates a multivector that shares data with the \this multivector.  Thus, changes to the
      values of one of these vectors also changes the values to the other.  
    \param Out
           clone - The new multivector.
	   
    \return Integer error code, set to 0 if successful, -1 if not able to create.
  */

  virtual int cloneView ( MultiVector<scalarType> * & clone, int * indices, int numVectors){

    Petra_TSF::MultiVector * clone1 = new Petra_TSF::MultiVector(View, *this, Indices, NumVectors); 
    if (clone1==0) PETRA_CHK_ERR(-1);
    clone = dynamic_cast<TSF::MultiVector *>(clone1);
    return(0);
  };

  //@}

  //@{ \name Value initialization methods
  //! Set multi-vector values to random numbers.
  /*!
    \return Integer error code, set to 0 if successful.

  */
   int random(){PETRA_CHK_ERR(Petra_RDP_MultiVector::Random());};

  //! Initialize all values in a multi-vector with constant value.
  /*!
    \param In
           Scalar - Value to use.

    \return Integer error code, set to 0 if successful.
  */
  virtual int putScalar (scalarType scalar) {PETRA_CHK_ERR(Petra_RDP_MultiVector::PutScalar(Scalar));};

  //! Copies a multivector into an existing multivector using = overloading.
  /*! Operator= copies the contents of one multivector to another existing multivector.  
    \param In
           source - Multivector to copy.
	   
    \return Contents of source in \e this.
  */

  virtual MultiVector<scalarType> & operator =  ( const Petra_TSF::MultiVector<scalarType> & source) {
    const Petra_RDP_MultiVector & source1 = dynamic_cast<const Petra_RDP_MultiVector &>(source);
    if (source1==0) PETRA_CHK_ERR(-1);
    return(dynamic_cast<MultiVector<scalarType>&> Petra_RDP_MultiVector::operator=(source));
  };
  //@}

  //@{ \name Element-by-element methods
  //! Puts element-wise absolute values of input Multi-vector in target.
  /*!
    \param In
           A - Input Multi-vector.
    \param Out
           \e this will contain the absolute values of the entries of A.

    \return Integer error code, set to 0 if successful.
    
    Note:  It is possible to use the same argument for A and \e this.
  */
  virtual int absoluteValue(const Petra_TSF::MultiVector<scalarType>& A){
    const Petra_RDP_MultiVector *A1 = dynamic_cast<const Petra_RDP_MultiVector *>(&A);
    if (A1==0) PETRA_CHK_ERR(-1);
    PETRA_CHK_ERR(Petra_RDP_MultiVector::Abs(*A1));
  };

  //! Puts element-wise reciprocal values of input Multi-vector in target.
  /*!
    \param In
           A - Input Multi-vector.
    \param Out
           \e this will contain the absolute values of the entries of A.

    \return Integer error code, set to 0 if successful.  Returns 2 if some entry
            is too small, but not zero.  Returns 1 if some entry is zero.
    
    Note:  It is possible to use the same argument for A and \e this.  Also, 
    if a given value of A is smaller than Petra_DoubleMin (defined in Petra_Petra.h),
    but nonzero, then the return code is 2.  If an entry is zero, the return code
    is 1.  However, in all cases the reciprocal value is still used, even
    if a NaN is the result.
  */
  virtual int reciprocal(const Petra_TSF::MultiVector<scalarType>& A){
    const Petra_RDP_MultiVector *A1 = dynamic_cast<const Petra_RDP_MultiVector *>(&A);
    if (A1==0) PETRA_CHK_ERR(-1);
    PETRA_CHK_ERR(Petra_RDP_MultiVector::Reciprocal(*A1));
  };

  //! Scale the current values of a multi-vector, \e this = scalar*\e this.
  /*!
    \param In
           scalar - Scale value.
    \param Out
           \e this - Multi-vector with scaled values.

    \return Integer error code, set to 0 if successful.
  */
  virtual int scale(scalarType scalar){PETRA_CHK_ERR(Petra_RDP_MultiVector::Scale(scalar));};

  //! Replace multi-vector values with scaled values of A, \e this = scalar*A.
  /*!
    \param In
           scalarA - scale value.
    \param In
           A - Multi-vector to copy.
    \param Out
           \e This - Multi-vector with values overwritten by scaled values of A.

    \return Integer error code, set to 0 if successful.
  */
  virtual int scale(scalarType scalarA, const MultiVector<scalarType>& A) {
    const Petra_RDP_MultiVector *A1 = dynamic_cast<const Petra_RDP_MultiVector *>(&A);
    if (A1==0) PETRA_CHK_ERR(-1);
    PETRA_CHK_ERR(Petra_RDP_MultiVector::Scale(scalarA, *A1));
  };

  //! Multiply a Petra_MultiVector with another, element-by-element.
  /*! This function supports diagonal matrix multiply.  A is usually a single vector
      while B and \e this may have one or more columns.  Note that B and \e this must
      have the same shape.  A can be one vector or have the same shape as B.  The actual
      computation is \e this = scalar * \e this + scalarAB * B @ A where @ denotes element-wise
      multiplication.
  */
  virtual int hadamardProduct(scalarType scalarAB, const MultiVector<scalarType>& A, 
		       const MultiVector<scalarType>& B, scalarType scalar ) {
    const Petra_RDP_MultiVector *A1 = dynamic_cast<const Petra_RDP_MultiVector *>(&A);
    if (A1==0) PETRA_CHK_ERR(-1);
    const Petra_RDP_MultiVector *B1 = dynamic_cast<const Petra_RDP_MultiVector *>(&B);
    if (B1==0) PETRA_CHK_ERR(-2);
    PETRA_CHK_ERR(Multiply(scalarAB, *A1, *B1, scalar));
  };


  //! Multiply a Petra_MultiVector by the reciprocal of another, element-by-element.
  /*! This function supports diagonal matrix scaling.  A is usually a single vector
      while B and \e this may have one or more columns.  Note that B and \e this must
      have the same shape.  A can be one vector or have the same shape as B. The actual
      computation is \e this = scalar * \e this + scalarAB * B @ A where @ denotes element-wise
      division.
  */
  virtual int hadamardReciprocalProduct(scalarType scalarAB, const MultiVector<scalarType>& A, 
				 const MultiVector<scalarType>& B, scalarType scalar ) {
    const Petra_RDP_MultiVector *A1 = dynamic_cast<const Petra_RDP_MultiVector *>(&A);
    if (A1==0) PETRA_CHK_ERR(-1);
    const Petra_RDP_MultiVector *B1 = dynamic_cast<const Petra_RDP_MultiVector *>(&B);
    if (B1==0) PETRA_CHK_ERR(-2);
    PETRA_CHK_ERR(ReciprocalMultiply(scalarAB, *A1, *B1, scalar));
  };
  //@}

  //@{ \name Vector-by-vector methods
  //! Update multi-vector values with scaled values of A, \e this = scalar*\e this + scalarA*A.
  /*!
    \param In
           scalarA - scale value for A.
    \param In
           A - Multi-vector to add.
    \param In
           scalar - scale value for \e this.
    \param Out
           \e This - Multi-vector with updatede values.

    \return Integer error code, set to 0 if successful.
  */
  virtual int update(scalarType scalarA, const MultiVector<scalarType>& A, scalarType scalar) {
    const Petra_RDP_MultiVector *A1 = dynamic_cast<const Petra_RDP_MultiVector *>(&A);
    if (A1==0) PETRA_CHK_ERR(-1);
    PETRA_CHK_ERR(Update(scalarA, *A1, scalar));
  };

  //! Update multi-vector with scaled values of A and B, \e this = scalar*\e this + scalarA*A + scalarB*B.
  /*!
    \param In
           scalarA - scale value for A.
    \param In
           A - Multi-vector to add.
    \param In
           scalarB - scale value for B.
    \param In
           B - Multi-vector to add.
    \param In
           scalar - scale value for \e this.
    \param Out
           \e This - Multi-vector with updatede values.

    \return Integer error code, set to 0 if successful.
  */
  virtual int update(scalarType scalarA, const MultiVector<scalarType>& A, 
		     scalarType scalarB, const MultiVector<scalarType>& B, scalarType scalar) {
    const Petra_RDP_MultiVector *A1 = dynamic_cast<const Petra_RDP_MultiVector *>(&A);
    if (A1==0) PETRA_CHK_ERR(-1);
    const Petra_RDP_MultiVector *B1 = dynamic_cast<const Petra_RDP_MultiVector *>(&B);
    if (B1==0) PETRA_CHK_ERR(-2);
    PETRA_CHK_ERR(update(scalarA, *A1, scalarB, *B1, scalar));
  };

  //! Computes dot product of each corresponding pair of vectors.
  /*!
    \param In
           A - Multi-vector to be used with the "\e this" multivector.
    \param Out
           result - result[i] will contain the ith dot product result.

    \return Integer error code, set to 0 if successful.
  */

  virtual int dotProduct(const MultiVector<scalarType>& A, scalarType * & result){
    const Petra_RDP_MultiVector *A1 = dynamic_cast<const Petra_RDP_MultiVector *>(&A);
    if (A1==0) PETRA_CHK_ERR(-1);
    PETRA_CHK_ERR(dotProduct( A1, result));
  };

  //! Compute 1-norm of each vector in multi-vector.
  /*!
    \param Out
           result - result[i] contains 1-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int norm1   (scalarType * & result){return(Petra_RDP_MultiVector::Norm1(result));};

  //! Compute 2-norm of each vector in multi-vector.
  /*!
    \param Out
           result - result[i] contains 2-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int norm2   (scalarType * & result){return(Petra_RDP_MultiVector::Norm2(Result));};

  //! Compute Inf-norm of each vector in multi-vector.
  /*!
    \param Out
           result - result[i] contains Inf-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int normInf (scalarType * & result){return(Petra_RDP_MultiVector::NormInf(Result));};

  //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
  /*!
    \param In
           Weights - Multi-vector of weights.  If Weights contains a single vector,
           that vector will be used as the weights for all vectors of \e this.  Otherwise,
           Weights should have the same number of vectors as \e this.
    \param Out
           result - result[i] contains the weighted 2-norm of ith vector.  Specifically
           if we denote the ith vector in the multivector by \f$x\f$, and the ith weight
           vector by \f$w\f$ and let j represent the jth entry of each vector, on return
           result[i] will contain the following result:
           \f[\sqrt{(1/n)\sum_{j=1}^n(w_jx_j)^2}\f],
           where \f$n\f$ is the global length of the vectors.

    \return Integer error code, set to 0 if successful.
  */
  virtual int normWeighted   (const MultiVector<scalarType>& weights, scalarType * & result) {
    const Petra_RDP_MultiVector *A1 = dynamic_cast<const Petra_RDP_MultiVector *>(&weights);
    if (weights1==0) PETRA_CHK_ERR(-1);
    PETRA_CHK_ERR(NormWeighted( weights1, result));
  };
  //! Compute minimum value of each vector in multi-vector.
  /*!
    \param Out
           result - result[i] contains minimum value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int minValue  (scalarType * & result){return(Petra_RDP_MultiVector::MinValue(Result));};

  //! Compute maximum value of each vector in multi-vector.
  /*!
    \param Out
           result - result[i] contains maximum value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int maxValue  (scalarType * & result){return(Petra_RDP_MultiVector::MaxValue(Result));};

  //! Compute mean (average) value of each vector in multi-vector.
  /*!
    \param Out
           result - result[i] contains mean value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int meanValue (scalarType * & result){return(Petra_RDP_MultiVector::MeanValue(Result));};
  //@}  

  //@{ \name Multivector block operations
  //! Computes a multivector update of another multivector times a dense matrix.
  /*!
    \param In
           A - Multi-vector to be multiplied by dense matrix.
    \param In
           B - Dense matrix to multiply with A.
    \param Out
           \e this - Result of A times B.
           
    \return Integer error code, set to 0 if successful.
  */

  virtual int blockUpdate(scalarType scalarAB, const MultiVector<scalarType>& A, const DenseMatrix<scalarType>& B, 
			  scalarType scalar) {
    const Petra_RDP_MultiVector *A1 = dynamic_cast<const Petra_RDP_MultiVector *>(&A);
    if (A1==0) PETRA_CHK_ERR(-1);
    const Petra_RDP_MultiVector *B1 = dynamic_cast<const Petra_RDP_MultiVector *>(&B);
    if (B1==0) PETRA_CHK_ERR(-2);
    PETRA_CHK_ERR(Multiply('N', 'N', scalarAB, *A1, *B1, scalar));
  };


  //! Computes block dot product of two multivectors.
  /*!
    \param In
           A - Multi-vector to be used with the "\e this" multivector.
    \param Out
           result - result(i,j) will contain the dot product result of the ith column of \e this 
	   and the jth column of A.

    \return Integer error code, set to 0 if successful.
  */

  virtual int blockDotProduct(const MultiVector<scalarType>& A, 
			      DenseMatrix<scalarType>& result) const {
    const Petra_RDP_MultiVector *A1 = dynamic_cast<const Petra_RDP_MultiVector *>(&A);
    if (A1==0) PETRA_CHK_ERR(-1);
    // The this object will actually be an argument to the Petra_RDP_MultiVector Multiply method
    const Petra_RDP_MultiVector *B1 = dynamic_cast<const Petra_RDP_MultiVector *>(this);
    if (B1==0) PETRA_CHK_ERR(-2);
    // The "result" object will be the temporary this object
    Petra_RDP_MultiVector *tmpThis = dynamic_cast<Petra_RDP_MultiVector *>(&result);
    if (tmpThis==0) PETRA_CHK_ERR(-3);
    scalarType zero = 0.0;
    scalarType one = 1.0;
    PETRA_CHK_ERR(tmpThis->Multiply('T', 'N', one, *A1, *B1, zero));
  };
  //@}  

  // Random number utilities


  //@{ Random function utilities (used for setting up random() method).
  /*!
    \param In
           Seed - Should be an odd positive integer value (stored as scalarType).

    \return Integer error code, set to 0 if successful.
  */
  virtual int setSeed(scalarType seed) {return(Petra_RDP_MultiVector::SetSeed(seed));};

  //! Get seed from Random function.
  /*!
    \return Current random number seed.
  */
  virtual scalarType seed() {return(Petra_RDP_MultiVector::Seed());};
  //@}
  
  //@{ Attribute access functions
  
  //! Returns the number of vectors in the multi-vector.
  virtual int numVectors() const{return(Petra_RDP_MultiVector::NumVectors());};

  //! Returns the global vector length of vectors in the multi-vector.
  virtual int length() const {return(Petra_RDP_MultiVector::GlobalLength());};
  //@}

};

// Constructor implementations

//=============================================================================
// Petra_BlockMap Constructor

template<class scalarType>
MultiVector<scalarType>::MultiVector(const Petra_BlockMap& Map, int NumVectors)
  : Petra_RDP_MultiVector(Map, NumVectors){}

//==========================================================================
// Copy Constructor


template<class scalarType>
MultiVector<scalarType>::MultiVector(const Petra_TSF::MultiVector<scalarType>& Source)
  :Petra_RDP_MultiVector(Source) {}

//==========================================================================
// This constructor copies in or makes view of a standard Fortran array

template<>
MultiVector<double>::MultiVector(Petra_DataAccess CV, const Petra_BlockMap& Map, 
						     double *A, int MyLDA, int NumVectors)
  : Petra_RDP_MultiVector(CV, Map, A, MyLDA, NumVectors) {}

//==========================================================================
// This constructor copies in or makes view of a C/C++ array of pointer

template<>
MultiVector<double>::MultiVector(Petra_DataAccess CV, const Petra_BlockMap& Map, 
						     double **ArrayOfPointers, int NumVectors)
  : Petra_RDP_MultiVector(CV, Map, ArrayOfPointers, NumVectors) {}

//==========================================================================

// This constructor copies or makes view of selected vectors, specified in Indices, 
// from an existing MultiVector

template<>
MultiVector<double>::MultiVector(Petra_DataAccess CV, 
			 const Petra_RDP_MultiVector& Source, 
			 int *Indices, int NumVectors)
  : Petra_RDP_MultiVector(CV, Source, Indices, NumVectors) {}

//==========================================================================

// This interface copies or makes view of a range of vectors from an existing MultiVector

template<>
MultiVector<double>::MultiVector(Petra_DataAccess CV, 
				     const Petra_RDP_MultiVector<double>& Source, 
				     int StartIndex, int NumVectors)
  : Petra_RDP_MultiVector(CV, Source, StartIndex, NumVectors){}

} // namespace Petra_TSF
#endif /* _Petra_TSF_MULTIVECTOR_H_ */
