/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK_HYPRE_H
#define IFPACK_HYPRE_H

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_HYPRE

#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "krylov.h"
#include "_hypre_parcsr_mv.h"
#include "_hypre_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_Condest.h"
#include "Ifpack_ScalingType.h"
#include "Epetra_CompObject.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_BlockMap.h"
#include "Epetra_Map.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Time.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_MpiComm.h"

#ifndef HYPRE_ENUMS
#define HYPRE_ENUMS
//! This enumerated type defines the allowed solvers and preconditioners in Hypre. Some can be used as both solver and preconditioner.
enum Hypre_Solver{
    BoomerAMG,
    ParaSails,
    Euclid,
    AMS,
    Hybrid,
    PCG,
    GMRES,
    FlexGMRES,
    LGMRES,
    BiCGSTAB
};

//! This enumerated type defines the two options for applying inverse, either solve or apply the preconditioner.
enum Hypre_Chooser{
    Solver,
    Preconditioner
};
#endif //HYPRE_ENUMS

//! This class is used to help with passing parameters in the SetParameter() function. Use this class to call Hypre's internal parameters.
class FunctionParameter{
  public:
    //! Single int constructor.
    FunctionParameter(Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, int), int param1) :
      chooser_(chooser),
      option_(0),
      int_func_(funct_name),
      int_param1_(param1) {}

    //! Single double constructor.
    FunctionParameter(Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, double), double param1):
      chooser_(chooser),
      option_(1),
      double_func_(funct_name),
      double_param1_(param1) {}

    //! Single double, single int constructor.
    FunctionParameter(Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, double, int), double param1, int param2):
      chooser_(chooser),
      option_(2),
      double_int_func_(funct_name),
      int_param1_(param2),
      double_param1_(param1) {}

    //! Two ints constructor.
    FunctionParameter(Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, int, int), int param1, int param2):
      chooser_(chooser),
      option_(3),
      int_int_func_(funct_name),
      int_param1_(param1),
      int_param2_(param2) {}

    //! Int pointer constructor.
    FunctionParameter(Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, int*), int *param1):
      chooser_(chooser),
      option_(4),
      int_star_func_(funct_name),
      int_star_param_(param1) {}

    //! Double pointer constructor.
    FunctionParameter(Hypre_Chooser chooser, int (*funct_name)(HYPRE_Solver, double*), double* param1):
      chooser_(chooser),
      option_(5),
      double_star_func_(funct_name),
      double_star_param_(param1) {}

    //! Only method of this class. Calls the function pointer with the passed in HYPRE_Solver
    int CallFunction(HYPRE_Solver solver, HYPRE_Solver precond){
      if(chooser_ == Solver){
        if(option_ == 0){
          return int_func_(solver, int_param1_);
        } else if(option_ == 1){
          return double_func_(solver, double_param1_);
        } else if(option_ == 2){
          return double_int_func_(solver, double_param1_, int_param1_);
        } else if (option_ == 3){
          return int_int_func_(solver, int_param1_, int_param2_);
        } else if (option_ == 4){
          return int_star_func_(solver, int_star_param_);
        } else {
          return double_star_func_(solver, double_star_param_);
        }
      } else {
        if(option_ == 0){
          return int_func_(precond, int_param1_);
        } else if(option_ == 1){
          return double_func_(precond, double_param1_);
        } else if(option_ == 2){
          return double_int_func_(precond, double_param1_, int_param1_);
        } else if(option_ == 3) {
          return int_int_func_(precond, int_param1_, int_param2_);
        } else if(option_ == 4) {
          return int_star_func_(precond, int_star_param_);
        } else {
          return double_star_func_(precond, double_star_param_);
        }
      }
    }

  private:
    Hypre_Chooser chooser_;
    int option_;
    int (*int_func_)(HYPRE_Solver, int);
    int (*double_func_)(HYPRE_Solver, double);
    int (*double_int_func_)(HYPRE_Solver, double, int);
    int (*int_int_func_)(HYPRE_Solver, int, int);
    int (*int_star_func_)(HYPRE_Solver, int*);
    int (*double_star_func_)(HYPRE_Solver, double*);
    int int_param1_;
    int int_param2_;
    double double_param1_;
    int *int_star_param_;
    double *double_star_param_;
};

namespace Teuchos {
  class ParameterList;
}

//! Ifpack_Hypre: A class for constructing and using an ILU factorization of a given Epetra_RowMatrix, using the Hypre library by Lawrence Livermore National Laboratories.

/*!
Class Ifpack_Hypre: A class for using methods of Hypre with Epetra objects.
*/

class Ifpack_Hypre: public Ifpack_Preconditioner {

public:
  // @{ Constructors and destructors.
  //! Constructor
  Ifpack_Hypre(Epetra_RowMatrix* A);

  //! Destructor
  ~Ifpack_Hypre(){ Destroy();}

  // @}
  // @{ Construction methods

  //! Initialize the preconditioner, does not touch matrix values.
  int Initialize();

  //! Returns \c true if the preconditioner has been successfully initialized.
  bool IsInitialized() const{ return(IsInitialized_);}

  //! Compute ILU factors L and U using the specified graph, diagonal perturbation thresholds and relaxation parameters.
  /*! This function computes the ILU(k) factors.
   */
  int Compute();

  //! If factor is completed, this query returns true, otherwise it returns false.
  bool IsComputed() const{ return(IsComputed_);}


  //! Set parameters using a Teuchos::ParameterList object.
  /* This method is only available if the Teuchos package is enabled.
     This method recognizes six parameter names: Solver,
     Preconditioner, SolveOrPrecondition, SetPreconditioner, NumFunctions and Functions. These names are
     case sensitive. Solver requires an enumerated parameter of type Hypre_Solver. Preconditioner is similar
     except requires the type be a preconditioner. The options are listed below:
                       Solvers                            Preconditioners
                       BoomerAMG                          BoomerAMG
                       AMS                                ParaSails
                       Hybrid                             AMS
                       PCG (Default)                      Euclid (Default)
                       GMRES
                       FlexGMRES
                       LGMRES
                       BiCGSTAB
     SolveOrPrecondition takes enumerated type Hypre_Chooser, Solver will solve the system, Preconditioner will apply the preconditioner.
     SetPreconditioner takes a boolean, true means the solver will use the preconditioner.
     NumFunctions takes an int that describes how many parameters will be passed into Functions. (This needs to be correct.)
     Functions takes an array of Ref Counted Pointers to an object called FunctionParameter. This class is implemented in Ifpack_Hypre.h.
     The object takes whether it is Solver or Preconditioner that we are setting a parameter for.
     The function in Hypre that sets the parameter, and the parameters for that function. An example is below:

     RCP<FunctionParameter> functs[2];
     functs[0] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetMaxIter, 1000)); // max iterations
     functs[1] = rcp(new FunctionParameter(Solver, &HYPRE_PCGSetTol, 1e-7)); // conv. tolerance
     list.set("NumFunctions", 2);
     list.set<RCP<FunctionParameter>*>("Functions", functs);
     NOTE: SetParameters() must be called to use ApplyInverse(), the solvers will not be created otherwise. An empty list is acceptable to use defaults.
  */
  int SetParameters(Teuchos::ParameterList& parameterlist);

    //! Set a parameter that takes a single int.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set maximum iterations would be &HYPRE_BoomerAMGSetMaxIter
    \param parameter (In) -The integer parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int), int parameter);

    //! Set a parameter that takes a single double.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set tolerance would be &HYPRE_BoomerAMGSetTol
    \param parameter (In) -The double parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double), double parameter);

    //! Set a parameter that takes a double then an int.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set relaxation weight for a given level would be &HYPRE_BoomerAMGSetLevelRelaxWt
    \param parameter1 (In) -The double parameter being set.
    \param parameter2 (In) - The integer parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double, int), double parameter1, int parameter2);

    //! Set a parameter that takes two int parameters.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set relaxation type for a given level would be &HYPRE_BoomerAMGSetCycleRelaxType
    \param parameter1 (In) -The first integer parameter being set.
    \param parameter2 (In) - The second integer parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int, int), int parameter1, int parameter2);

    //! Set a parameter that takes a double*.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set relaxation weight would be &HYPRE_BoomerAMGSetRelaxWeight
    \param parameter (In) -The double* parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double*), double* parameter);

    //! Set a parameter that takes an int*.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type set to Solver or Preconditioner, whatever the parameter is setting for.
    \param *pt2Func (In) -The function that sets the parameter. It must set parameters for the type of solver or preconditioner that was created.
      An example is if the solver is BoomerAMG, the function to set grid relax type would be &HYPRE_BoomerAMGSetGridRelaxType
    \param parameter (In) -The int* parameter being set.

    \return Integer error code, set to 0 if successful.
   */
    int SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int*), int* parameter);

    //! Sets the solver that is used by the Solve() and ApplyInverse() methods. Until this is called, the default solver is PCG.
    /*!
    \param chooser (In) - A Hypre_Chooser enumerated type. If Solver, then we are selecting which solver, if Preconditioner, we are choosing which preconditioner to use.
    \param Solver (In) -A Hypre_Solver enumerated type to select the solver or preconditioner. Options for solver are:
    BoomerAMG, AMS, Hybrid, PCG, GMRES, FlexGMRES, LGMRES, and BiCGSTAB. See Hypre Ref Manual for more info on the solvers.
    Options for Preconditioner are: BoomerAMG, ParaSails, Euclid, and AMS.

    \return Integer error code, set to 0 if successful.
  */

    int SetParameter(Hypre_Chooser chooser, Hypre_Solver Solver);

    //! Sets the solver to use the selected preconditioner.
    /*!
    \param UsePreconditioner (In) -A boolean, true use preconditioner, false do not use the supplied preconditioner with the solver.
    The solver and preconditioner must have been selected and the solver must be one of the following solvers:
      Hybrid, PCG, GMRES, FlexGMRES, LGMRES, BiCGSTAB.

    \return Integer error code, set to 0 if successful.
  */

    int SetParameter(bool UsePreconditioner){ UsePreconditioner_ = UsePreconditioner; return 0;}

    //! Choose to solve the problem or apply the preconditioner.
    /*!
    \param chooser (In) -A Hypre_Chooser enumerated type, either Solver or Preconditioner.
    The chosen type must have been selected before this method is called.

    \return Integer error code, set to 0 if successful.
  */
    int SetParameter(Hypre_Chooser chooser) { SolveOrPrec_ = chooser; return 0;}

  //! Call all the function pointers stored in this object.
    int CallFunctions() const;

  //! If set true, transpose of this operator will be applied.
  /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
      affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface
      does not support transpose use, this method should return a value of -1.

      \param
       UseTranspose_in - (In) If true, multiply by the transpose of operator, otherwise just use operator.

      \return Always returns 0.
  */
  int SetUseTranspose(bool UseTranspose_in) {UseTranspose_ = UseTranspose_in; return(0);};

  // @}

  // @{ Mathematical functions.
  // Applies the matrix to X, returns the result in Y.
  int Apply(const Epetra_MultiVector& X,
               Epetra_MultiVector& Y) const{ return(Multiply(false,X,Y));}

  //! Returns the result of a Epetra_Operator multiplied with an Epetra_MultiVector X in Y.
  /*! In this implementation, we use the Hypre matrix to multiply with so that the map is the same
      as what is expected in solving methods.

    \param
    trans - (In) If true, use the transpose operation.
           X - (In) A Epetra_MultiVector of dimension NumVectors to mulitply with.
    \param Out
           Y - (Out) A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int Multiply(bool Trans, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
  /*! In this implementation, we use several existing attributes to determine how virtual
      method ApplyInverse() should call the concrete method Solve().  We pass in the UpperTriangular(),
      the Epetra_CrsMatrix::UseTranspose(), and NoDiagonal() methods. The most notable warning is that
      if a matrix has no diagonal values we assume that there is an implicit unit diagonal that should
      be accounted for when doing a triangular solve.

    \param
           X - (In) A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
           Y - (Out) A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Computes the estimated condition number and returns the value.
  double Condest(const Ifpack_CondestType CT = Ifpack_Cheap,
                 const int MaxIters = 1550,
                 const double Tol = 1e-9,
                 Epetra_RowMatrix* Matrix_in = 0);

  //! Returns the computed estimated condition number, or -1.0 if not computed.
  double Condest() const{ return(Condest_);}

  // @}
  // @{ Query methods

  //! Returns a character string describing the operator
  const char* Label() const {return(Label_);}

  //! Sets label for \c this object.
  int SetLabel(const char* Label_in)
  {
    strcpy(Label_,Label_in);
    return(0);
  }

  //! Returns a reference to the map that should be used for domain.
  const Epetra_Map& OperatorDomainMap() const{ return *GloballyContiguousRowMap_;}

  //! Returns a reference to the map that should be used for range.
  const Epetra_Map& OperatorRangeMap() const{ return *GloballyContiguousRowMap_;}

  //! Returns 0.0 because this class cannot compute Inf-norm.
  double NormInf() const {return(0.0);};

  //! Returns false because this class cannot compute an Inf-norm.
  bool HasNormInf() const {return(false);};

  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(UseTranspose_);};

  //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
  const Epetra_Comm & Comm() const{return(A_->Comm());};

  //! Returns a reference to the matrix to be preconditioned.
  const Epetra_RowMatrix& Matrix() const{ return(*A_);}

  //! Returns the Hypre matrix that was created upon construction.
  const HYPRE_IJMatrix& HypreMatrix()
  {
    if(IsInitialized() == false)
      Initialize();
    return(HypreA_);
  }

  //! Prints on stream basic information about \c this object.
  virtual std::ostream& Print(std::ostream& os) const;

  //! Returns the number of calls to Initialize().
  virtual int NumInitialize() const{ return(NumInitialize_);}

  //! Returns the number of calls to Compute().
  virtual int NumCompute() const{ return(NumCompute_);}

  //! Returns the number of calls to ApplyInverse().
  virtual int NumApplyInverse() const{ return(NumApplyInverse_);}

  //! Returns the time spent in Initialize().
  virtual double InitializeTime() const{ return(InitializeTime_);}

  //! Returns the time spent in Compute().
  virtual double ComputeTime() const{ return(ComputeTime_);}

  //! Returns the time spent in ApplyInverse().
  virtual double ApplyInverseTime() const{ return(ApplyInverseTime_);}

  //! Returns the number of flops in the initialization phase.
  virtual double InitializeFlops() const{ return(0.0);}

  //! Returns the number of flops in the compute phase.
  virtual double ComputeFlops() const{ return(ComputeFlops_);}

  //! Returns the number of flops in the apply inverse phase.
  virtual double ApplyInverseFlops() const{ return(ApplyInverseFlops_);}

private:

  // @}
  // @{ Private methods

  //! Copy constructor (should never be used)
  Ifpack_Hypre(const Ifpack_Hypre& RHS) : Time_(RHS.Comm()){}

  //! operator= (should never be used)
  Ifpack_Hypre& operator=(const Ifpack_Hypre& /*RHS*/){ return(*this);}

  //! Destroys all internal data
  void Destroy();

  //! Returns the MPI communicator used in the Epetra matrix
  MPI_Comm GetMpiComm() const
    { return (dynamic_cast<const Epetra_MpiComm*>(&A_->Comm()))->GetMpiComm();}

  //! Returns the result of a Ifpack_ILU forward/back solve on a Epetra_MultiVector X in Y.
  /*!
    \param In
    Trans -If true, solve transpose problem.
    \param
    X - (In) A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y - (Out) A Epetra_MultiVector of dimension NumVectorscontaining result.

    \return Integer error code, set to 0 if successful.
  */
  int Solve(bool Trans, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;


  //! Returns the number of global matrix rows.
  int NumGlobalRows() const {return(A_->NumGlobalRows());};

  //! Returns the number of global matrix columns.
  int NumGlobalCols() const {return(A_->NumGlobalCols());};

  //! Returns the number of local matrix rows.
  int NumMyRows() const {return(A_->NumMyRows());};

  //! Returns the number of local matrix columns.
  int NumMyCols() const {return(A_->NumMyCols());};

  //! Sets the solver type to be the passed in solver type.
  int SetSolverType(Hypre_Solver solver);

  //! Sets the preconditioner type to be the passed in type.
  int SetPrecondType(Hypre_Solver precond);

  //! Create the solver.
  int CreateSolver();

  //! Create the Preconditioner.
  int CreatePrecond();

  //! Copies matrix data from Epetra matrix to Hypre matrix.
  int CopyEpetraToHypre();

  //! Add a function to be called in Compute()
  int AddFunToList(Teuchos::RCP<FunctionParameter> NewFun);

  //! Create a BoomerAMG solver.
  int Hypre_BoomerAMGCreate(MPI_Comm /*comm*/, HYPRE_Solver *solver)
    { return HYPRE_BoomerAMGCreate(solver);}

  //! Create a ParaSails solver.
  int Hypre_ParaSailsCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParaSailsCreate(comm, solver);}

  //! Create a Euclid solver.
  int Hypre_EuclidCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_EuclidCreate(comm, solver);}

  //! Create an AMS solver.
  int Hypre_AMSCreate(MPI_Comm /*comm*/, HYPRE_Solver *solver)
    { return HYPRE_AMSCreate(solver);}

  //! Create a Hybrid solver.
  int Hypre_ParCSRHybridCreate(MPI_Comm /*comm*/, HYPRE_Solver *solver)
    { return HYPRE_ParCSRHybridCreate(solver);}

  //! Create a PCG solver.
  int Hypre_ParCSRPCGCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRPCGCreate(comm, solver);}

  //! Create a GMRES solver.
  int Hypre_ParCSRGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRGMRESCreate(comm, solver);}

  //! Create a FlexGMRES solver.
  int Hypre_ParCSRFlexGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRFlexGMRESCreate(comm, solver);}

  //! Create a LGMRES solver.
  int Hypre_ParCSRLGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRLGMRESCreate(comm, solver);}

  //! Create a BiCGSTAB solver.
  int Hypre_ParCSRBiCGSTABCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRBiCGSTABCreate(comm, solver);}

  // @}
  // @{ Internal data

  //! Pointer to the Epetra_RowMatrix to factorize
  Teuchos::RefCountPtr<Epetra_RowMatrix> A_;
  //! This objects copy of the ParameterList
  Teuchos::ParameterList List_;
  //! Needed to support Epetra_Operator abstract class
  bool UseTranspose_;
  //! A condition estimate for the preconditioner, -1 until Compute()
  double Condest_;
  //! If \c true, the preconditioner has been successfully initialized.
  bool IsInitialized_;
  //! If \c true, the preconditioner has been successfully computed.
  bool IsComputed_;
  //! Label of \c this object.
  char Label_[160];
  //! Contains the number of successful calls to Initialize().
  int NumInitialize_;
  //! Contains the number of successful call to Compute().
  int NumCompute_;
  //! Contains the number of successful call to ApplyInverse().
  mutable int NumApplyInverse_;
  //! Contains the time for all successful calls to Initialize().
  double InitializeTime_;
  //! Contains the time for all successful calls to Compute().
  double ComputeTime_;
  //! Contains the time for all successful calls to ApplyInverse().
  mutable double ApplyInverseTime_;
  //! Contains the number of flops for Compute().
  double ComputeFlops_;
  //! Contain sthe number of flops for ApplyInverse().
  mutable double ApplyInverseFlops_;
  //! Used for timing issues
  mutable Epetra_Time Time_;

  //! The Hypre matrix created in initialize()
  mutable HYPRE_IJMatrix HypreA_;
  //! Pointer to the CSR (same matrix)
  mutable HYPRE_ParCSRMatrix ParMatrix_;
  //! The Hypre Vector for input
  mutable HYPRE_IJVector XHypre_;
  //! The Hypre Vector for output
  mutable HYPRE_IJVector YHypre_;
  mutable HYPRE_ParVector ParX_;
  mutable HYPRE_ParVector ParY_;
  mutable hypre_ParVector *XVec_;
  mutable hypre_ParVector *YVec_;
  mutable hypre_Vector *XLocal_;
  mutable hypre_Vector *YLocal_;
  //! The Hypre Solver if doing a solve
  mutable HYPRE_Solver Solver_;
  //! The Hypre Solver if applying preconditioner
  mutable HYPRE_Solver Preconditioner_;
  //  The following are pointers to functions to use the solver and preconditioner.
  int (Ifpack_Hypre::*SolverCreatePtr_)(MPI_Comm, HYPRE_Solver*);
  int (*SolverDestroyPtr_)(HYPRE_Solver);
  int (*SolverSetupPtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  int (*SolverSolvePtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  int (*SolverPrecondPtr_)(HYPRE_Solver, HYPRE_PtrToParSolverFcn, HYPRE_PtrToParSolverFcn, HYPRE_Solver);
  int (Ifpack_Hypre::*PrecondCreatePtr_)(MPI_Comm, HYPRE_Solver*);
  int (*PrecondDestroyPtr_)(HYPRE_Solver);
  int (*PrecondSetupPtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);
  int (*PrecondSolvePtr_)(HYPRE_Solver, HYPRE_ParCSRMatrix, HYPRE_ParVector, HYPRE_ParVector);

  bool *IsSolverSetup_;
  bool *IsPrecondSetup_;
  //! Is the system to be solved or apply preconditioner
  Hypre_Chooser SolveOrPrec_;
  //! These are linear maps that meet the needs of Hypre
  Teuchos::RCP<const Epetra_Map> GloballyContiguousRowMap_;
  Teuchos::RCP<const Epetra_Map> GloballyContiguousColMap_;
  //! Counter of the number of parameters set
  int NumFunsToCall_;
  //! Which solver was chosen
  Hypre_Solver SolverType_;
  //! Which preconditioner was chosen
  Hypre_Solver PrecondType_;
  //! Should the preconditioner be used in the solver
  bool UsePreconditioner_;
  //! This contains a list of function pointers that will be called in compute
  std::vector<Teuchos::RCP<FunctionParameter> > FunsToCall_;
};

#endif // HAVE_HYPRE
#endif /* IFPACK_HYPRE_H */
