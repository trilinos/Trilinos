#ifndef _TSF_RDP_MULTIVECTOR_H_
#define _TSF_RDP_MULTIVECTOR_H_

//! TSF_RDP_Multivector:  The Trilinos Virtual MultiVector Class.
/*! The TSF_RDP_Multivector class is a pure virtual class that specifies the 
    interface that is generally required by the Trilinos Solver Framework solver classes.  
*/
#include "TSF_Object.h"

class TSF_RDP_MultiVector {
    
  public:

  //! TSF_RDP_MultiVector Destructor.
  virtual ~TSF_RDP_MultiVector(void){};

  //! Create a new multivector of the same class as the \this multivector, fill with zeros.
  /*! CloneZero is essentially a call to the non-trivial default constructor of the multivector class of
      which the \this object is a member.  X returns with a new multivector in that class, 
      with NumVector columns.
    \param In
           NumVectors - The number of vectors that will be in X.
    \param Out
           X - The new multivector X.
	   
    \return Integer error code, set to 0 if successful.
  */

  virtual int CloneZero ( TSF_RDP_MultiVector * & X, const int NumVectors) = 0;

  //! Create a new multivector that is an exact full(deep) copy the \this multivector.
  /*! CloneCopy is essentially a call to the standard copy constructor of the multivector class of
      which the \this object is a member.  X returns with a complete replication of \this.
    \param Out
           X - The new multivector X.
	   
    \return Integer error code, set to 0 if successful.
  */

  virtual int CloneCopy ( TSF_RDP_MultiVector * & X) = 0;

  //! Create a new multivector that is a view (shallow copy) of all or some vectors in the \this multivector.
  /*! CloneView creates a multivector that shares data with the \this multivector.  Thus, changes to the
      values of one of these vectors also changes the values to the other.  
    \param Out
           X - The new multivector.
	   
    \return Integer error code, set to 0 if successful.
  */

  virtual int CloneView ( TSF_RDP_MultiVector * & X, int * Indices, int NumVectors) = 0;

  //! Set multi-vector values to random numbers.
  /*!
    \return Integer error code, set to 0 if successful.

  */
   virtual int Random() = 0;
  //! Initialize all values in a multi-vector with constant value.
  /*!
    \param In
           Scalar - Value to use.

    \return Integer error code, set to 0 if successful.
  */
  virtual int PutScalar (double Scalar) = 0;
  // Mathematical functions.

  //! Computes dot product of each corresponding pair of vectors.
  /*!
    \param In
           A - Multi-vector to be used with the "\e this" multivector.
    \param Out
           Result - Result[i] will contain the ith dot product result.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Dot(const TSF_RDP_MultiVector& A, double *Result) = 0;

  //! Puts element-wise absolute values of input Multi-vector in target.
  /*!
    \param In
           A - Input Multi-vector.
    \param Out
           \e this will contain the absolute values of the entries of A.

    \return Integer error code, set to 0 if successful.
    
    Note:  It is possible to use the same argument for A and \e this.
  */
  virtual int Abs(const TSF_RDP_MultiVector& A) = 0;

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
  virtual int Reciprocal(const TSF_RDP_MultiVector& A) = 0;

  //! Scale the current values of a multi-vector, \e this = Scalar*\e this.
  /*!
    \param In
           Scalar - Scale value.
    \param Out
           \e This - Multi-vector with scaled values.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Scale(double Scalar) = 0;

  //! Replace multi-vector values with scaled values of A, \e this = Scalar*A.
  /*!
    \param In
           ScalarA - Scale value.
    \param In
           A - Multi-vector to copy.
    \param Out
           \e This - Multi-vector with values overwritten by scaled values of A.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Scale(double ScalarA, const TSF_RDP_MultiVector& A) = 0;

  //! Update multi-vector values with scaled values of A, \e this = Scalar*\e this + ScalarA*A.
  /*!
    \param In
           ScalarA - Scale value for A.
    \param In
           A - Multi-vector to add.
    \param In
           Scalar - Scale value for \e this.
    \param Out
           \e This - Multi-vector with updatede values.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Update(double ScalarA, const TSF_RDP_MultiVector& A, double Scalar) = 0;

  //! Update multi-vector with scaled values of A and B, \e this = Scalar*\e this + ScalarA*A + ScalarB*B.
  /*!
    \param In
           ScalarA - Scale value for A.
    \param In
           A - Multi-vector to add.
    \param In
           ScalarB - Scale value for B.
    \param In
           B - Multi-vector to add.
    \param In
           Scalar - Scale value for \e this.
    \param Out
           \e This - Multi-vector with updatede values.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Update(double ScalarA, const TSF_RDP_MultiVector& A, 
		     double ScalarB, const TSF_RDP_MultiVector& B, double Scalar) = 0;

  //! Compute 1-norm of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains 1-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Norm1   (double * Result) = 0;

  //! Compute 2-norm of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains 2-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Norm2   (double * Result) = 0;

  //! Compute Inf-norm of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains Inf-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int NormInf (double * Result) = 0;

  //! Compute Weighted 2-norm (RMS Norm) of each vector in multi-vector.
  /*!
    \param In
           Weights - Multi-vector of weights.  If Weights contains a single vector,
           that vector will be used as the weights for all vectors of \e this.  Otherwise,
           Weights should have the same number of vectors as \e this.
    \param Out
           Result - Result[i] contains the weighted 2-norm of ith vector.  Specifically
           if we denote the ith vector in the multivector by \f$x\f$, and the ith weight
           vector by \f$w\f$ and let j represent the jth entry of each vector, on return
           Result[i] will contain the following result:
           \f[\sqrt{(1/n)\sum_{j=1}^n(w_jx_j)^2}\f],
           where \f$n\f$ is the global length of the vectors.

    \return Integer error code, set to 0 if successful.
  */
  virtual int NormWeighted   (const TSF_RDP_MultiVector& Weights, double * Result) = 0;

  //! Compute minimum value of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains minimum value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int MinValue  (double * Result) = 0;

  //! Compute maximum value of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains maximum value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int MaxValue  (double * Result) = 0;

  //! Compute mean (average) value of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains mean value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  virtual int MeanValue (double * Result) = 0;
  
  
  //! Matrix-Matrix multiplication, \e this = Scalar*\e this + ScalarAB*A*B.
  /*! This function performs a variety of matrix-matrix multiply operations, interpreting
      the Petra_MultiVectors (\e this-aka C , A and B) as 2D matrices.  Variations are due to
      the fact that A, B and C can be local replicated or global distributed
      Petra_MultiVectors and that we may or may not operate with the transpose of 
      A and B.  Possible cases are:
\verbatim

     Total of 32 case (2^5).
                                           Num
         OPERATIONS                        case  Notes
     1) C(local) = A^X(local) * B^X(local)  4   (X=Transpose or Not, No comm needed) 
     2) C(local) = A^T(distr) * B  (distr)  1   (2D dot product, replicate C)
     3) C(distr) = A  (distr) * B^X(local)  2   (2D vector update, no comm needed)

     Note that the following operations are not meaningful for 
     1D distributions:

     1) C(local) = A^T(distr) * B^T(distr)  1
     2) C(local) = A  (distr) * B^X(distr)  2
     3) C(distr) = A^X(local) * B^X(local)  4
     4) C(distr) = A^X(local) * B^X(distr)  4
     5) C(distr) = A^T(distr) * B^X(local)  2
     6) C(local) = A^X(distr) * B^X(local)  4
     7) C(distr) = A^X(distr) * B^X(local)  4
     8) C(local) = A^X(local) * B^X(distr)  4

\endverbatim

  \param In
         TransA - Operate with the transpose of A if = 'T', else no transpose if = 'N'.
  \param In
         TransB - Operate with the transpose of B if = 'T', else no transpose if = 'N'.

  \param In
         ScalarAB - Scalar to multiply with A*B.
  \param In
         A - Multi-vector.
  \param In
         B - Multi-vector.
  \param In
         Scalar - Scalar to multiply with \e this.

    \return Integer error code, set to 0 if successful.

\warning {Each multi-vector A, B and \e this is checked if it has constant stride using the
         ConstantStride() query function.  If it does not have constant stride, a temporary
	 copy is made and used for the computation.  This activity is transparent to the user,
	 except that there is memory and computation overhead.  All temporary space is deleted
	 prior to exit.}
	 
  */
  virtual int Multiply(char TransA, char TransB, double ScalarAB, 
		       const TSF_RDP_MultiVector& A, const TSF_RDP_MultiVector& B,
		       double Scalar ) = 0;
  


  //! Multiply a Petra_MultiVector with another, element-by-element.
  /*! This function supports diagonal matrix multiply.  A is usually a single vector
      while B and \e this may have one or more columns.  Note that B and \e this must
      have the same shape.  A can be one vector or have the same shape as B.  The actual
      computation is \e this = Scalar * \e this + ScalarAB * B @ A where @ denotes element-wise
      multiplication.
  */
  virtual int Multiply(double ScalarAB, const TSF_RDP_MultiVector& A, const TSF_RDP_MultiVector& B,
		       double Scalar ) = 0;


  //! Multiply a Petra_MultiVector by the reciprocal of another, element-by-element.
  /*! This function supports diagonal matrix scaling.  A is usually a single vector
      while B and \e this may have one or more columns.  Note that B and \e this must
      have the same shape.  A can be one vector or have the same shape as B. The actual
      computation is \e this = Scalar * \e this + ScalarAB * B @ A where @ denotes element-wise
      division.
  */
  virtual int ReciprocalMultiply(double ScalarAB, const TSF_RDP_MultiVector& A, const TSF_RDP_MultiVector& B,
		       double Scalar ) = 0;


  // Random number utilities


  //! Set seed for Random function.
  /*!
    \param In
           Seed - Should be an odd positive integer value (stored as double).

    \return Integer error code, set to 0 if successful.
  */
  virtual int SetSeed(double Seed) = 0;

  //! Get seed from Random function.
  /*!
    \return Current random number seed.
  */
  virtual double Seed() = 0;
  
  // Attribute access functions
  
  //! Returns the number of vectors in the multi-vector.
  virtual int NumVectors() const = 0;

  //! Returns the local vector length on the calling processor of vectors in the multi-vector.
  virtual int MyLength() const = 0;

  //! Returns the global vector length of vectors in the multi-vector.
  virtual int GlobalLength() const = 0;

  //! Returns the stride between  vectors in the multi-vector (only meaningful if ConstantStride() is true).
  virtual int Stride() const = 0;

};

#endif /* _TSF_RDP_MULTIVECTOR_H_ */
