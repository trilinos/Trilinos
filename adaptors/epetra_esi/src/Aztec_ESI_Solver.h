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
class Solver : public virtual esi::Object,
                      public virtual epetra_esi::Object,
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

  /** Non-ESI function for translating string parameters to Aztec options
      and params settings.
     @param numParams Input, number of parameter strings.
     @param paramStrings Input, list of strings.
     @param options Input/Output, user-allocated Aztec options array (of length
               AZ_OPTIONS_SIZE).
     @param params Input/Output, user-allocated Aztec params array (of length
               AZ_PARAMS_SIZE).
     @return error-code 0 if successful.
  */
  static int translateStringsToAztecSettings(int numParams, char** paramStrings,
                                             int* options, double* params);

//
//Note from ABW: The following is ugly, but for now I can't think of
//a clean way to translate string parameters into Aztec options/params....
//The string-lists need to be updated whenever changes are
//made in az_aztec_defs.h.
//
//az_def_map.h contains a list of symbol-value pairs gleaned directly
//from az_aztec_defs.h by a script that uses awk. This is what we'll use to
//figure out that, for example, "AZ_gmres" has the value 1.
//
  static Epetra_Array<const char*>& get_az_def_map()
    {
#include "az_def_map.h"
      static Epetra_Array<const char*> azDefMap(num_def_strs, num_def_strs,
                                                (const char**)az_def_map);
      return(azDefMap);
    }

//
//And here's the REALLY bad part: the order of the strings in az_option_strs
//and az_param_strs matters, because their positions are used for the
//index into the Aztec options_ and params_ arrays... In other words, don't
//mess with these string lists.
//Don't try this at home, folks.
//
  static Epetra_Array<const char*>& get_az_option_strs()
    {
      static const char* az_option_strs[] = {
      "AZ_solver",
      "AZ_scaling",
      "AZ_precond",
      "AZ_conv",
      "AZ_output",
      "AZ_pre_calc",
      "AZ_max_iter",
      "AZ_poly_ord",
      "AZ_overlap",
      "AZ_type_overlap",
      "AZ_kspace",
      "AZ_orthog",
      "AZ_aux_vec",
      "AZ_reorder",
      "AZ_keep_info",
      "AZ_recursion_level",
      "AZ_print_freq",
      "AZ_graph_fill",
      "AZ_subdomain_solve",
      "AZ_init_guess"
      };
       static int num_option_strs = 20;
      static Epetra_Array<const char*> azOptionStrs(num_option_strs,
                                                    num_option_strs,
                                                  (const char**)az_option_strs);
      return(azOptionStrs);
    }

  static Epetra_Array<const char*>& get_az_param_strs()
    {
      static const char* az_param_strs[] = {
        "AZ_tol",
        "AZ_drop",
        "AZ_ilut_fill",
        "AZ_omega",
        "AZ_rthresh",
        "AZ_athresh",
        "AZ_weights"
      };
      static int num_param_strs = 7;
      static Epetra_Array<const char*> azParamStrs(num_param_strs,
                                                   num_param_strs,
                                            (const char**)az_param_strs);
      return(azParamStrs);
    }

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

  static int azOptionValue(const char* valStr,
                           Epetra_Array<const char*>& azDefStrings);

  static float azParamValue(const char* valStr, 
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

