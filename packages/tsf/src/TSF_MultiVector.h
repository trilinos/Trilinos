#ifndef _TSF_MULTIVECTOR_H_
#define _TSF_MULTIVECTOR_H_

#include "TSF_Object.h"
#include "TSF_DenseMatrix.h"

namespace TSF {
//! TSF::Multivector:  The Trilinos Virtual MultiVector Class.
/*! The TSF::Multivector class is a pure virtual class that specifies the 
    interface that is generally required by the Trilinos Solver Framework solver classes, 
    especially block iterative methods.  
*/
template<class scalarType>
class MultiVector {
    
  public:

  //@{ \name Destructor 
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
	   
    \return Integer error code, set to 0 if successful.
  */

  virtual int cloneZero ( MultiVector<scalarType> * & clone, const int numVectors) = 0;

  //! Create a new multivector that is an exact full(deep) copy the \this multivector.
  /*! cloneCopy is essentially a call to the standard copy constructor of the multivector class of
      which the \this object is a member.  clone returns with a complete replication of \this.
    \param Out
           clone - The new multivector clone.
	   
    \return Integer error code, set to 0 if successful.
  */

  virtual int cloneCopy ( MultiVector<scalarType> * & clone) = 0;

  //! Create a new multivector that is a view (shallow copy) of all or some vectors in the \this multivector.
  /*! cloneView creates a multivector that shares data with the \this multivector.  Thus, changes to the
      values of one of these vectors also changes the values to the other.  
    \param Out
           clone - The new multivector.
	   
    \return Integer error code, set to 0 if successful.
  */

  virtual int cloneView ( MultiVector<scalarType> * & clone, int * indices, int numVectors) = 0;
  //@}

  //@{ \name Value initialization methods
  //! Set multi-vector values to random numbers.
  /*!
    \return Integer error code, set to 0 if successful.

  */
   virtual int random() = 0;
  //! Initialize all values in a multi-vector with constant value.
  /*!
    \param In
           Scalar - Value to use.

    \return Integer error code, set to 0 if successful.
  */
  virtual int putScalar (scalarType scalar) = 0;

  //! Copies a multivector into an existing multivector using = overloading.
  /*! Operator= copies the contents of one multivector to another existing multivector.  
    \param In
           source - Multivector to copy.
	   
    \return Contents of source in \e this.
  */

  virtual MultiVector<type> & operator =  ( const MultiVector<scalarType> & source) = 0;
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
  virtual int absoluteValue(const MultiVector<scalarType>& A) = 0;

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
  virtual int reciprocal(const MultiVector<scalarType>& A) = 0;

  //! Scale the current values of a multi-vector, \e this = scalar*\e this.
  /*!
    \param In
           scalar - Scale value.
    \param Out
           \e this - Multi-vector with scaled values.

    \return Integer error code, set to 0 if successful.
  */
  virtual int scale(scalarType scalar) = 0;

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
  virtual int scale(scalarType scalarA, const MultiVector<scalarType>& A) = 0;

  //! Multiply a Petra_MultiVector with another, element-by-element.
  /*! This function supports diagonal matrix multiply.  A is usually a single vector
      while B and \e this may have one or more columns.  Note that B and \e this must
      have the same shape.  A can be one vector or have the same shape as B.  The actual
      computation is \e this = scalar * \e this + scalarAB * B @ A where @ denotes element-wise
      multiplication.
  */
  virtual int hadamardProduct(scalarType scalarAB, const MultiVector<scalarType>& A, 
		       const MultiVector<scalarType>& B,
		       scalarType scalar ) = 0;


  //! Multiply a Petra_MultiVector by the reciprocal of another, element-by-element.
  /*! This function supports diagonal matrix scaling.  A is usually a single vector
      while B and \e this may have one or more columns.  Note that B and \e this must
      have the same shape.  A can be one vector or have the same shape as B. The actual
      computation is \e this = scalar * \e this + scalarAB * B @ A where @ denotes element-wise
      division.
  */
  virtual int hadamardReciprocalProduct(scalarType scalarAB, const MultiVector<scalarType>& A, 
				 const MultiVector<scalarType>& B,
				 scalarType scalar ) = 0;
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
  virtual int update(scalarType scalarA, const MultiVector<scalarType>& A, scalarType scalar) = 0;

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
		     scalarType scalarB, const MultiVector<scalarType>& B, scalarType scalar) = 0;

  //! Computes dot product of each corresponding pair of vectors.
  /*!
    \param In
           A - Multi-vector to be used with the "\e this" multivector.
    \param Out
           result - result[i] will contain the ith dot product result.

    \return Integer error code, set to 0 if successful.
  */

  virtual int dotProduct(const MultiVector<scalarType>& A, scalarType * & result) = 0;

  //! Compute 1-norm of each vector in multi-vector.
  /*!
    \param Out
           result - result[i] contains 1-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int norm1   (scalarType * & result) = 0;

  //! Compute 2-norm of each vector in multi-vector.
  /*!
    \param Out
           result - result[i] contains 2-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int norm2   (scalarType * & result) = 0;

  //! Compute Inf-norm of each vector in multi-vector.
  /*!
    \param Out
           result - result[i] contains Inf-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int normInf (scalarType * & result) = 0;

  //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
  /*!
    \param In
           weights - Multi-vector of weights.  If Weights contains a single vector,
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
  virtual int normWeighted   (const MultiVector<scalarType>& weights, scalarType * & result) = 0;

  //! Compute minimum value of each vector in multi-vector.
  /*!
    \param Out
           result - result[i] contains minimum value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int minValue  (scalarType * & result) = 0;

  //! Compute maximum value of each vector in multi-vector.
  /*!
    \param Out
           result - result[i] contains maximum value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int maxValue  (scalarType * & result) = 0;

  //! Compute mean (average) value of each vector in multi-vector.
  /*!
    \param Out
           result - result[i] contains mean value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int meanValue (scalarType * & result) = 0;
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
			  scalarType scalar) = 0;

  //! Computes block dot product of two multivectors.
  /*!
    \param In
           A - Multi-vector to be used with the "\e this" multivector.
    \param Out
           result - result(i,j) will contain the dot product result of the ith 
	   column of \e this and the jth column of A.

    \return Integer error code, set to 0 if successful.
  */

  virtual int blockDotProduct(const MultiVector<scalarType>& A, 
			      DenseMatrix<scalarType>& result) const = 0;
  //@}  

  // Random number utilities


  //@{ Random function utilities (used for setting up random() method).
  /*!
    \param In
           Seed - Should be an odd positive integer value (stored as scalarType).

    \return Integer error code, set to 0 if successful.
  */
  virtual int setSeed(scalarType seed) = 0;

  //! Get seed from Random function.
  /*!
    \return Current random number seed.
  */
  virtual scalarType seed() = 0;
  //@}
  
  //@{ Attribute access functions
  
  //! Returns the number of vectors in the multi-vector.
  virtual int numVectors() const = 0;

  //! Returns the global vector length of vectors in the multi-vector.
  virtual int length() const = 0;
  //@}
};

} // TSF namespace
#endif /* _TSF_MULTIVECTOR_H_ */
