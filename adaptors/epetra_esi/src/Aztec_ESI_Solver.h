#ifndef _Aztec_ESI_Solver_h_
#define _Aztec_ESI_Solver_h_

#include "AztecOO.h"

//
//Epetra_ESI_CrsMatrix.h includes everything "above" it,
//including Epetra_ESI_Object.h, Epetra_Comm.h, Epetra_SerialComm.h,
//Epetra_MpiComm.h, Epetra_Map.h,
//Epetra_ESI_Vector.h, Epetra_Vector.h, and the ESI headers.
//
#include "Epetra_ESI_CHK_ERR.h"
#include "Epetra_ESI_CrsMatrix.h"

namespace aztecoo_esi {

/** Aztec's ESI Solver implementation.
This class implements these ESI interfaces:
<ul>
<li> esi::Object
<li> esi::Operator
<li> esi::Solver
</ul>
Note that although this class is templated, it may only be instantiated on
the type-pair double,int.
*/

template<class Scalar, class Ordinal>
class Solver : public virtual epetra_esi::Object,
                      public virtual esi::Operator<Scalar, Ordinal>,
                      public virtual esi::Solver<Scalar, Ordinal>,
                      public virtual esi::SolverIterative<Scalar, Ordinal>,
                      public virtual AztecOO
{
 public:
  typedef TYPENAME esi::scalarTraits<Scalar>::magnitude_type magnitude_type;

  /** Constructor. */
  Solver(epetra_esi::CrsMatrix<Scalar,Ordinal>* A);

  /** Constructor that takes a esi::Operator. */
  Solver(esi::Operator<Scalar, Ordinal>* esi_op);

  /** Destructor. */
  virtual ~Solver()
    { if (whichConstructor_ == 1 || petraAalloced_) delete petraA_; }


  //
  //esi::Operator functions.
  //

  /** Function for signalling that the initialization is complete. 
   aztecoo_esi::Solver currently doesn't do anything in response to this function
   call.
   @return error-code 0 if successful.
  */
  virtual esi::ErrorCode setup( void )
    { return(0); }


  /** Function for applying this operator to a vector. (aztecoo_esi::Solver
   does a solve in response to this function call.) If the setup() function
   has not yet been called, it will be called internally, automatically at this
   point.
   @param b Input. The right-hand-side vector. The run-time type of this vector
       is required to be epetra_esi::Vector.
   @param x On exit, the contents of this will be the result of the solve.
       The run-time type of this vector must also be epetra_esi::Vector.
   @return error-code 0 if successful.
  */
  virtual esi::ErrorCode apply(esi::Vector<Scalar, Ordinal>& b, esi::Vector<Scalar, Ordinal>& x)
    { return( solve(b, x) ); }


  //
  //esi::Solver functions.
  //

  /** Function to attempt the solution of the linear system. If the setup()
   function has not yet been called, it will be called internally,
   automatically at this point.
   @param b Input. The right-hand-side vector. The run-time type of this vector
     is required to be epetra_esi::Vector.
   @param x On exit, the contents of this will be the result of the solve.
       The run-time type of this vector must also be epetra_esi::Vector.
   @return error-code 0 if successful.
  */
  virtual esi::ErrorCode solve(esi::Vector<Scalar, Ordinal>& b, esi::Vector<Scalar,
 Ordinal>& x)
    {

      Epetra_Vector* px = dynamic_cast<Epetra_Vector*>(&x);
      Epetra_Vector* pb = dynamic_cast<Epetra_Vector*>(&b);
      if (px != NULL || pb != NULL) {
        return( Iterate(petraA_, px, pb, maxIters_, tolerance_) ); 
      }

      //If the dynamic_cast operations failed (i.e., if the input vectors are
      //not Epetra_Vectors, we'll create some Epetra_Vectors...
      return( createPetraVectorsThenSolve(b, x) );
    }


  /** Function for setting control parameters. Each parameter string is
   usually a space-separated key-value pair, e.g., "AZ_tol 1.e-9".
    @param numParams Input, number of parameter strings.
    @param paramStrings Input, list of strings.
    @return error-code 0 if successful.
  */
  virtual esi::ErrorCode parameters(int numParams, char** paramStrings);


  //
  //esi::SolverIterative functions.
  //

  /** Query for the current Operator. */
  virtual esi::ErrorCode getOperator(esi::Operator<Scalar, Ordinal>*& A)
    {
      A = dynamic_cast<esi::Operator<Scalar, Ordinal>*>(petraA_);
      return(0);
    }

  /** Set the current Operator. */
  virtual esi::ErrorCode setOperator(esi::Operator<Scalar, Ordinal>& A)
    {
      petraA_ = NULL;
      petraA_ = dynamic_cast<epetra_esi::CrsMatrix<double,int>*>(&A);
      if (petraA_ != NULL) return(0);

      petraA_ = new epetra_esi::Operator<double,int>(A);
      if (petraA_ == NULL) {
        msgAbort("aztecoo_esi::Solver ctor failed to allocate epetra_esi::Operator.");
      }
      petraAalloced_ = true;
      return(0);
    }

  /** Query for the current Preconditioner. NOT IMPLEMENTED. */
  virtual esi::ErrorCode getPreconditioner(esi::Preconditioner<Scalar, Ordinal>*& A)
    { return(-1); }

  /** Set the current Preconditioner. NOT IMPLEMENTED. */
  virtual esi::ErrorCode setPreconditioner(esi::Preconditioner<Scalar, Ordinal>& A)
    { return(-1); }

  /** Get the convergence tolerance. */
  virtual esi::ErrorCode getTolerance( magnitude_type & tol )
    { tol = tolerance_; return(0); }

  /** Set the convergence tolerance. */
  virtual esi::ErrorCode setTolerance( magnitude_type tol )
    { tolerance_ = tol; return(0); }

  /** Get the maximum number of iterations. */
  virtual esi::ErrorCode getMaxIterations( Ordinal & maxIterations )
    { maxIterations = maxIters_; return(0); }

  /** Set the maximum number of iterations. */
  virtual esi::ErrorCode setMaxIterations( Ordinal maxIterations )
    { maxIters_ = maxIterations; return(0); }

  /** Query the number of iterations that were taken during the previous solve.
  */
  virtual esi::ErrorCode getNumIterationsTaken(Ordinal& itersTaken)
    { itersTaken = (Ordinal)(status_[AZ_its]); return(0); }

 private:
  int addInterfaces();

  int handleAzOption(int optionIndx,
                     const char* keyStr,
                     const char* valStr,
                     Epetra_Array<const char*>& azDefStrings);

  int handleAzParam(int paramIndx,
                    const char* keyStr,
                    const char* valStr,
                    Epetra_Array<const char*>& azDefStrings);

  int createPetraVectorsThenSolve(esi::Vector<Scalar, Ordinal>& b,
                                  esi::Vector<Scalar, Ordinal>& x);

  Epetra_RowMatrix* petraA_;
  Ordinal maxIters_;
  Scalar tolerance_;
  int whichConstructor_;
  bool petraAalloced_;
};

}; //namespace aztecoo_esi

#ifdef EPETRA_ESI_INCLUDE_IMPLEMENTATION
#include "Aztec_ESI_Solver.cpp"
#endif

#endif //_Aztec_ESI_Solver_h_

