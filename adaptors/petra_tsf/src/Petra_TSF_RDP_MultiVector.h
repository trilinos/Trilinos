#ifndef _PETRA_TSF_RDP_MULTIVECTOR_H_
#define _PETRA_TSF_RDP_MULTIVECTOR_H_

//! Petra_TSF_RDP_Multivector:  The Petra implementation of the Trilinos Virtual MultiVector Class.
/*! The Petra_TSF_RDP_Multivector class is a wrapper that uses the Petra_RDP_MultiVector class to
    implement the TSF_RDP_MultiVector pure virtual class.  
*/

#include "Petra_RDP_MultiVector.h"
#include "TSF_RDP_MultiVector.h"
class Petra_TSF_RDP_MultiVector: public TSF_RDP_MultiVector, public Petra_RDP_MultiVector {
    
  public:

  //! Basic Petra_TSF_RDP_MultiVector constuctor.
  /*! Creates a Petra_TSF_RDP_MultiVector object and fills with zero values.  

    \param In 
           Map - A Petra_LocalMap, Petra_Map or Petra_BlockMap.

	   \warning Note that, because Petra_LocalMap
	   derives from Petra_Map and Petra_Map derives from Petra_BlockMap, this constructor works
	   for all three types of Petra map classes.
    \param In 
           NumVectors - Number of vectors in multi-vector.

    \return Pointer to a Petra_TSF_RDP_MultiVector.

  */
  Petra_TSF_RDP_MultiVector(const Petra_BlockMap& Map, int NumVectors);

  //! Petra_TSF_RDP_MultiVector copy constructor.
  
  Petra_TSF_RDP_MultiVector(const Petra_RDP_MultiVector& Source);
  
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
  Petra_TSF_RDP_MultiVector(Petra_DataAccess CV, const Petra_BlockMap& Map, 
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
  Petra_TSF_RDP_MultiVector(Petra_DataAccess CV, const Petra_BlockMap& Map, 
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
  Petra_TSF_RDP_MultiVector(Petra_DataAccess CV,  
			const Petra_RDP_MultiVector& Source, int *Indices, int NumVectors);

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
  Petra_TSF_RDP_MultiVector(Petra_DataAccess CV, 
			const Petra_RDP_MultiVector& Source, int StartIndex, 
			int NumVectors);

  //! Petra_TSF_RDP_MultiVector Destructor.
  virtual ~Petra_TSF_RDP_MultiVector(void);

  //! Initialize all values in a multi-vector with constant value.
  /*!
    \param In
           Scalar - Value to use.

    \return Integer error code, set to 0 if successful.
  */
  int PutScalar (double Scalar){return(Petra_RDP_MultiVector::PutScalar(Scalar));};
  
  //! Set multi-vector values to random numbers.
  /*!
    \return Integer error code, set to 0 if successful.

  */
   int Random(){return(Petra_RDP_MultiVector::Random());};

  //! Scale the current values of a multi-vector, \e this = Scalar*\e this.
  /*!
    \param In
           Scalar - Scale value.
    \param Out
           \e This - Multi-vector with scaled values.

    \return Integer error code, set to 0 if successful.
  */
  int Scale(double Scalar){return(Petra_RDP_MultiVector::Scale(Scalar));};

  //! Create a new multivector of the same class as the \this multivector, fill with zeros.
  /*! CloneZero is essentially a call to the non-trivial default constructor of the multivector class of
      which the \this object is a member.  A returns with a new multivector in that class, 
      with NumVector columns.
    \param In
           NumVectors - The number of vectors that will be in A.
    \param Out
           A - The new multivector A.
	   
    \return Integer error code, set to 0 if successful.
  */

  virtual int CloneZero ( TSF_RDP_MultiVector * & A, const int NumVectors);

  //! Create a new multivector that is an exact full(deep) copy the \this multivector.
  /*! CloneCopy is essentially a call to the standard copy constructor of the multivector class of
      which the \this object is a member.  A returns with a complete replication of \this.
    \param Out
           A - The new multivector A.
	   
    \return Integer error code, set to 0 if successful.
  */

  virtual int CloneCopy ( TSF_RDP_MultiVector * & A);

  //! Create a new multivector that is a view (shallow copy) of all or some vectors in the \this multivector.
  /*! CloneView creates a multivector that shares data with the \this multivector.  Thus, changes to the
      values of one of these vectors also changes the values to the other.  
    \param Out
           A - The new multivector.
	   
    \return Integer error code, set to 0 if successful.
  */

  virtual int CloneView ( TSF_RDP_MultiVector * & A, int * Indices, int NumVectors);

  // Mathematical functions.

  //! Computes dot product of each corresponding pair of vectors.
  /*!
    \param In
           A - Multi-vector to be used with the "\e this" multivector.
    \param Out
           Result - Result[i] will contain the ith dot product result.

    \return Integer error code, set to 0 if successful.
  */
  virtual int Dot(const TSF_RDP_MultiVector& A, double *Result);

  //! Puts element-wise absolute values of input Multi-vector in target.
  /*!
    \param In
           A - Input Multi-vector.
    \param Out
           \e this will contain the absolute values of the entries of A.

    \return Integer error code, set to 0 if successful.
    
    Note:  It is possible to use the same argument for A and \e this.
  */
  virtual int Abs(const TSF_RDP_MultiVector& A);

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
  virtual int Reciprocal(const TSF_RDP_MultiVector& A);

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
  virtual int Scale(double ScalarA, const TSF_RDP_MultiVector& A);

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
  virtual int Update(double ScalarA, const TSF_RDP_MultiVector& A, double Scalar);

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
		     double ScalarB, const TSF_RDP_MultiVector& B, double Scalar);

  //! Compute 1-norm of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains 1-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  int Norm1   (double * Result){return(Petra_RDP_MultiVector::Norm1(Result));};

  //! Compute 2-norm of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains 2-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  int Norm2   (double * Result){return(Petra_RDP_MultiVector::Norm2(Result));};

  //! Compute Inf-norm of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains Inf-norm of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  int NormInf (double * Result){return(Petra_RDP_MultiVector::NormInf(Result));};

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
  virtual int NormWeighted   (const TSF_RDP_MultiVector& Weights, double * Result);
  
 
  //! Compute minimum value of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains minimum value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  int MinValue  (double * Result){return(Petra_RDP_MultiVector::MinValue(Result));};

  //! Compute maximum value of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains maximum value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  int MaxValue  (double * Result){return(Petra_RDP_MultiVector::MaxValue(Result));};

  //! Compute mean (average) value of each vector in multi-vector.
  /*!
    \param Out
           Result - Result[i] contains mean value of ith vector.

    \return Integer error code, set to 0 if successful.
  */
  int MeanValue (double * Result){return(Petra_RDP_MultiVector::MeanValue(Result));};
 
  //! Matrix-Matrix multiplication, \e this = Scalar*\e this + ScalarAB*A*B.
  /*! This function performs a variety of matrix-matrix multiply operations, interpreting
      the Petra_TSF_RDP_MultiVectors (\e this-aka C , A and B) as 2D matrices.  Variations are due to
      the fact that A, B and C can be local replicated or global distributed
      Petra_TSF_RDP_MultiVectors and that we may or may not operate with the transpose of 
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
		       double Scalar );
  


  //! Multiply a Petra_TSF_RDP_MultiVector with another, element-by-element.
  /*! This function supports diagonal matrix multiply.  A is usually a single vector
      while B and \e this may have one or more columns.  Note that B and \e this must
      have the same shape.  A can be one vector or have the same shape as B.  The actual
      computation is \e this = Scalar * \e this + ScalarAB * B @ A where @ denotes element-wise
      multiplication.
  */
  virtual int Multiply(double ScalarAB, const TSF_RDP_MultiVector& A, const TSF_RDP_MultiVector& B,
		       double Scalar );


  //! Multiply a Petra_TSF_RDP_MultiVector by the reciprocal of another, element-by-element.
  /*! This function supports diagonal matrix scaling.  A is usually a single vector
      while B and \e this may have one or more columns.  Note that B and \e this must
      have the same shape.  A can be one vector or have the same shape as B. The actual
      computation is \e this = Scalar * \e this + ScalarAB * B @ A where @ denotes element-wise
      division.
  */
  virtual int ReciprocalMultiply(double ScalarAB, const TSF_RDP_MultiVector& A, const TSF_RDP_MultiVector& B,
		       double Scalar );

  // Random number utilities


  //! Set seed for Random function.
  /*!
    \param In
           Seed - Should be an odd positive integer value (stored as double).

    \return Integer error code, set to 0 if successful.
  */
  int SetSeed(double Seed){return(Petra_RDP_MultiVector::SetSeed(Seed));};

  //! Get seed from Random function.
  /*!
    \return Current random number seed.
  */
  double Seed(){return(Petra_RDP_MultiVector::Seed());};

  // Attribute access functions
  
  //! Returns the number of vectors in the multi-vector.
  int NumVectors() const {return(Petra_RDP_MultiVector::NumVectors());};

  //! Returns the local vector length on the calling processor of vectors in the multi-vector.
  int MyLength() const {return(Petra_RDP_MultiVector::MyLength());};

  //! Returns the global vector length of vectors in the multi-vector.
  int GlobalLength()  const {return(Petra_RDP_MultiVector::GlobalLength());};

  //! Returns the stride between  vectors in the multi-vector (only meaningful if ConstantStride() is true).
  int Stride() const {return(Petra_RDP_MultiVector::Stride());};
};

#endif /* _Petra_TSF_RDP_MULTIVECTOR_H_ */
