#ifndef _Epetra_ESI_Vector_h_
#define _Epetra_ESI_Vector_h_

#include "Epetra_Vector.h"

//
//Epetra_ESI_IndexSpace.h includes everything "above" it,
//like Epetra_ESI_Object.h, Epetra_Comm.h, Epetra_SerialComm.h,
//Epetra_MpiComm.h, Epetra_Map.h and the ESI headers.
//
#include "Epetra_ESI_CHK_ERR.h"
#include "Epetra_ESI_IndexSpace.h"


namespace epetra_esi {

/** Petra's ESI Vector implementation.
This class implements these ESI interfaces:
<ul>
<li>   esi::Object
<li>   esi::Vector
</ul>
Note that although this class is templated, it may only be instantiated on
the type-pair double,int.
*/

template<class Scalar, class Ordinal>
class Vector : public virtual esi::Object,
               public virtual epetra_esi::Object,
                         public virtual esi::Vector<Scalar, Ordinal>,
                         public virtual Epetra_Vector
{
 public:
  /** Constructor. */
  Vector(epetra_esi::IndexSpace<Ordinal>& indexspace);

  /** Constructor that accepts a Epetra_Vector. No data copy is performed. */
  Vector(const Epetra_Vector& invec);

  /** Destructor. */
  virtual ~Vector();

  typedef TYPENAME esi::scalarTraits<Scalar>::magnitude_type magnitude_type;

  //
  //esi::Vector functions.
  //

  /** Produce a clone of 'this' vector.
    @param x Output. Newly allocated epetra_esi::Vector. Caller is responsible
       for destroying x.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode clone(esi::Vector<Scalar,Ordinal>*& x)
    { Epetra_BlockMap& pmap = (Epetra_Map&)Map();
      epetra_esi::IndexSpace<Ordinal>* pesimap =
           dynamic_cast<epetra_esi::IndexSpace<Ordinal>*>(&pmap);
      if (pesimap == NULL) return(-1);
      x = new epetra_esi::Vector<double,int>(*pesimap);
      if (x == NULL) return(-1); else return(0);
    }


  /** Get the mathematical global size of this vector.
    @param globalSize Output.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode getGlobalSize(Ordinal& globalSize)
    { globalSize = GlobalLength(); return(0); }

  /** Get the local size of this vector.
     @param localSize Output.
     @return error-code 0 if successful.
  */
  virtual esi::ErrorCode getLocalSize(Ordinal& localSize)
    { localSize = MyLength(); return(0); }

  /** Get the IndexSpace that this vector was constructed with.
    @param indexspace Output.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode getIndexSpace(esi::IndexSpace<Ordinal>*& indexspace)
    { Epetra_BlockMap& pmap = (Epetra_Map&)Map();
      indexspace = dynamic_cast<esi::IndexSpace<Ordinal>*>(&pmap);
      if (indexspace == NULL) EPETRA_ESI_ERR_BEHAVIOR(-1);
      return(0);
    }


  /** Set the coefficients of this vector.
    @param scalar Input. The value to be put throughout the vector.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode put(Scalar scalar)
    { EPETRA_ESI_ERR_BEHAVIOR(PutScalar(scalar)); }


  /** Scale this vector.
   @param scalar Input. The value by which, to multiply the coefficients of
               this vector.
   @return error-code, 0 if successful.
  */
  virtual esi::ErrorCode scale(Scalar scalar)
    { EPETRA_ESI_ERR_BEHAVIOR( Scale(scalar) ); }


  /** Obtain the 1-norm of this vector.
    @param norm Output. 1-norm of this vector.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode norm1(magnitude_type& norm)
    { EPETRA_ESI_ERR_BEHAVIOR( Norm1(&norm) ); }


  /** Obtain the 2-norm of this vector.
    @param norm Output. 2-norm of this vector.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode norm2(magnitude_type& norm)
    { EPETRA_ESI_ERR_BEHAVIOR( Norm2(&norm) ); }


  /** Obtain the square of the 2-norm of this vector. (Note from ABW: why is
    this an ESI function??)
    @param norm Output. Square of the 2-norm of this vector.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode norm2squared(magnitude_type& norm)
    { magnitude_type nrm;
      CHK_ERR( Norm2(&nrm) );
      norm = nrm*nrm;
      return(0);
    }


  /** Obtain the infinity-norm of this vector.
    @param norm Output. infinity-norm of this vector.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode normInfinity(magnitude_type& norm)
    { EPETRA_ESI_ERR_BEHAVIOR( NormInf(&norm) ); }



  /** Obtain the minimum absolute coefficient value. THIS FUNCTION IS NOT
   IMPLEMENTED.
   @param minCoef Output. 
   @param error-code 0 if successful.
  */
  virtual esi::ErrorCode minAbsCoef(magnitude_type& minCoef)
    { EPETRA_ESI_ERR_BEHAVIOR(-1); }


  /** Entry-by-entry multiplication: this[i] <- this[i] * x[i],
     for i = 0 to localSize */
  virtual esi::ErrorCode scaleDiagonal(esi::Vector<Scalar,Ordinal>& x)
    { epetra_esi::Vector<Scalar,Ordinal>* px =
          dynamic_cast<epetra_esi::Vector<Scalar,Ordinal>*>(&x);
      if (px == NULL) return(-1);
      EPETRA_ESI_ERR_BEHAVIOR( Multiply(1.0, *px, *this, 0.0) );
    }


  /** Copy the x vector into 'this' vector.
    @param x The run-time type of x is required to be epetra_esi::Vector.
    @return error-code 0 if successful
  */
  virtual esi::ErrorCode copy(esi::Vector<Scalar,Ordinal>& x)
    { epetra_esi::Vector<Scalar,Ordinal>* px =
             dynamic_cast<epetra_esi::Vector<Scalar,Ordinal>*>(&x);
      if (px == NULL) return(-1);
      EPETRA_ESI_ERR_BEHAVIOR(Update(1.0, *px, 0.0) );
    }


  /** Form product = x dot 'this'.
    @param x The run-time type of x is required to be epetra_esi::Vector.
    @param product Output. The result of the dot-product.
    @return error-code 0 if successful
  */
  virtual esi::ErrorCode dot(esi::Vector<Scalar,Ordinal>& x, Scalar & product)
    { epetra_esi::Vector<Scalar,Ordinal>* px =
             dynamic_cast<epetra_esi::Vector<Scalar,Ordinal>*>(&x);
      if (px == NULL) EPETRA_ESI_ERR_BEHAVIOR(-1);
      EPETRA_ESI_ERR_BEHAVIOR(Dot(*px, &product) );
    }


  /** Form this = this + scalar*x.
    @param x The run-time type of x is required to be epetra_esi::Vector.
    @param scalar Input.
    @return error-code 0 if successful
  */
  virtual esi::ErrorCode axpy(esi::Vector<Scalar,Ordinal>& x, Scalar scalar)
    { epetra_esi::Vector<Scalar,Ordinal>* px =
             dynamic_cast<epetra_esi::Vector<Scalar,Ordinal>*>(&x);
      if (px == NULL) EPETRA_ESI_ERR_BEHAVIOR(-1);
      EPETRA_ESI_ERR_BEHAVIOR(Update(scalar, *px, 1.0) );
    };


  /** Form this = scalar*this + x.
    @param x The run-time type of x is required to be epetra_esi::Vector.
    @param scalar Input.
    @return error-code 0 if successful
  */
  virtual esi::ErrorCode aypx(Scalar scalar, esi::Vector<Scalar,Ordinal>& x)
    { epetra_esi::Vector<Scalar,Ordinal>* px =
             dynamic_cast<epetra_esi::Vector<Scalar,Ordinal>*>(&x);
      if (px == NULL) return(-1);
      EPETRA_ESI_ERR_BEHAVIOR(Update(1.0, *px, scalar) );
    };


  /** Form this = a*x + b*y.
    @param a Input.
    @param x The run-time type of x is required to be epetra_esi::Vector.
    @param b Input.
    @param y The run-time type of x is required to be epetra_esi::Vector.
    @return error-code 0 if successful
  */
  virtual esi::ErrorCode axpby(Scalar a, esi::Vector<Scalar,Ordinal> &x, Scalar b, esi::Vector<Scalar,Ordinal> &y)
    { epetra_esi::Vector<Scalar,Ordinal>* px =
             dynamic_cast<epetra_esi::Vector<Scalar,Ordinal>*>(&x);
      if (px == NULL) return(-1);
      epetra_esi::Vector<Scalar,Ordinal>* py =
             dynamic_cast<epetra_esi::Vector<Scalar,Ordinal>*>(&y);
      if (py == NULL) EPETRA_ESI_ERR_BEHAVIOR(-1);
      EPETRA_ESI_ERR_BEHAVIOR(Update(a, *px, b, *py, 0.0) );
    };


  /** Get a pointer for read-access to the data array. This is intended to
    by read-only access, but we don't currently have a way to enforce that.
    Multiple calls may be made to getCoefPtrRead without intervening calls to
    releaseCoefPtrLock. However, all getCoefPtrRead calls must eventually be 
    matched with releaseCoefPtrLock calls.
    @param coef Output. This is the address of a pointer which will be set
       to point to this vector's internal data array.
    @return error-code 0 if successful
  */
  virtual esi::ErrorCode getCoefPtrReadLock(Scalar*& coef)
    {  if (currentReadWrites > 0) return(-1);
       currentReads++;
       int err = ExtractView(&coef);
       coefPtr = coef;
       EPETRA_ESI_ERR_BEHAVIOR(err);
    };


  /** Get a pointer for read-write-access to the data array. This call will
    fail if any calls to getCoefPtrRead or getCoefPtrReadWriteLock have been
    made without matching calls to releaseCoefPtrLock. In other words, only one
    read-write-lock may be 'out' at a time, and no read-locks may be out at
    the same time. Any calls to getArrayRead that are made after a call to this
    function, and before the matching call to releaseArray, will fail.
    @param coef Output. This is the address of a pointer which will be set
       to point to this vector's internal data array.
    @return error-code 0 if successful
  */
  virtual esi::ErrorCode getCoefPtrReadWriteLock(Scalar*& coef)
    { if (currentReads > 0 || currentReadWrites > 0) {
         EPETRA_ESI_ERR_BEHAVIOR(-1);
      }
      currentReadWrites++;
      int err = ExtractView(&coef);
      coefPtr = coef;
      EPETRA_ESI_ERR_BEHAVIOR(err);
    };
  

  /** Release the vector from a previous getCoefPtrRead* call.
    @param coef Input/Output. On input, must be a pointer obtained in a
      call to getCoefPtrRead*, and will be set to NULL on output.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode releaseCoefPtrLock(Scalar*& coef)
    { if (coef != coefPtr) EPETRA_ESI_ERR_BEHAVIOR(-1);
      if (currentReadWrites > 0) currentReadWrites = 0;
      if (currentReads > 0) currentReads--;
      coef = NULL;
      return(0);
    };


 private:
  epetra_esi::IndexSpace<Ordinal>* ispace_;
  int currentReads, currentReadWrites;
  Scalar* coefPtr;
  int whichConstructor;
};

}; //namespace epetra_esi

#ifdef EPETRA_ESI_INCLUDE_IMPLEMENTATION
#include "Epetra_ESI_Vector.cpp"
#endif

#endif

