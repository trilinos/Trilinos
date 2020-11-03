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
#include "Ifpack_Hypre.h"
#if defined(HAVE_HYPRE) && defined(HAVE_MPI)
#include <stdexcept>

#include "Ifpack_Utils.h"
#include "Epetra_MpiComm.h"
#include "Epetra_IntVector.h"
#include "Epetra_Import.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "krylov.h"
#include "_hypre_parcsr_mv.h"
#include "_hypre_IJ_mv.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE.h"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcpFromRef;

// The Python script that generates the ParameterMap needs to be after these typedefs
typedef int (*int_func)(HYPRE_Solver, int);
typedef int (*double_func)(HYPRE_Solver, double);
typedef int (*double_int_func)(HYPRE_Solver, double, int);
typedef int (*int_int_func)(HYPRE_Solver, int, int);
typedef int (*int_star_func)(HYPRE_Solver, int*);
typedef int (*int_star_star_func)(HYPRE_Solver, int**);
typedef int (*double_star_func)(HYPRE_Solver, double*);
typedef int (*int_int_double_double_func)(HYPRE_Solver, int, int, double, double);
typedef int (*int_int_int_double_int_int_func)(HYPRE_Solver, int, int, int, double, int, int);
typedef int (*char_star_func)(HYPRE_Solver, char*);

//! This class is used to help with passing parameters in the SetParameter() function. Use this class to call Hypre's internal parameters.
class FunctionParameter{
  public:
    //! Single int constructor.
    FunctionParameter(Hypre_Chooser chooser, int_func funct, int param1) :
      chooser_(chooser),
      option_(0),
      int_func_(funct),
      int_param1_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, int param1) :
      chooser_(chooser),
      option_(0),
      int_func_(hypreMapIntFunc_.at(funct_name)),
      int_param1_(param1) {}

    //! Single double constructor.
    FunctionParameter(Hypre_Chooser chooser, double_func funct, double param1):
      chooser_(chooser),
      option_(1),
      double_func_(funct),
      double_param1_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, double param1):
      chooser_(chooser),
      option_(1),
      double_func_(hypreMapDoubleFunc_.at(funct_name)),
      double_param1_(param1) {}

    //! Single double, single int constructor.
    FunctionParameter(Hypre_Chooser chooser, double_int_func funct, double param1, int param2):
      chooser_(chooser),
      option_(2),
      double_int_func_(funct),
      int_param1_(param2),
      double_param1_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, double param1, int param2):
      chooser_(chooser),
      option_(2),
      double_int_func_(hypreMapDoubleIntFunc_.at(funct_name)),
      int_param1_(param2),
      double_param1_(param1) {}

    //! Two ints constructor.
    FunctionParameter(Hypre_Chooser chooser, int_int_func funct, int param1, int param2):
      chooser_(chooser),
      option_(3),
      int_int_func_(funct),
      int_param1_(param1),
      int_param2_(param2) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, int param1, int param2):
      chooser_(chooser),
      option_(3),
      int_int_func_(hypreMapIntIntFunc_.at(funct_name)),
      int_param1_(param1),
      int_param2_(param2) {}

    //! Int pointer constructor.
    FunctionParameter(Hypre_Chooser chooser, int_star_func funct, int *param1):
      chooser_(chooser),
      option_(4),
      int_star_func_(funct),
      int_star_param_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, int *param1):
      chooser_(chooser),
      option_(4),
      int_star_func_(hypreMapIntStarFunc_.at(funct_name)),
      int_star_param_(param1) {}

    //! Double pointer constructor.
    FunctionParameter(Hypre_Chooser chooser, double_star_func funct, double* param1):
      chooser_(chooser),
      option_(5),
      double_star_func_(funct),
      double_star_param_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, double* param1):
      chooser_(chooser),
      option_(5),
      double_star_func_(hypreMapDoubleStarFunc_.at(funct_name)),
      double_star_param_(param1) {}

    //! Two ints, two doubles constructor.
    FunctionParameter(Hypre_Chooser chooser, int_int_double_double_func funct, int param1, int param2, double param3, double param4):
      chooser_(chooser),
      option_(6),
      int_int_double_double_func_(funct),
      int_param1_(param1),
      int_param2_(param2),
      double_param1_(param3),
      double_param2_(param4) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, int param1, int param2, double param3, double param4):
      chooser_(chooser),
      option_(6),
      int_int_double_double_func_(hypreMapIntIntDoubleDoubleFunc_.at(funct_name)),
      int_param1_(param1),
      int_param2_(param2),
      double_param1_(param3),
      double_param2_(param4) {}

    //! Integer pointer to list of integer pointers
    FunctionParameter(Hypre_Chooser chooser, int_star_star_func funct, int ** param1):
      chooser_(chooser),
      option_(7),
      int_star_star_func_(funct),
      int_star_star_param_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, int** param1):
      chooser_(chooser),
      option_(7),
      int_star_star_func_(hypreMapIntStarStarFunc_.at(funct_name)),
      int_star_star_param_(param1) {}

    //! Five ints, one double constructor.
    FunctionParameter(Hypre_Chooser chooser, int_int_int_double_int_int_func funct, int param1, int param2, int param3, double param4, int param5, int param6):
      chooser_(chooser),
      option_(8),
      int_int_int_double_int_int_func_(funct),
      int_param1_(param1),
      int_param2_(param2),
      int_param3_(param3),
      int_param4_(param5),
      int_param5_(param6),
      double_param1_(param4) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, int param1, int param2, int param3, double param4, int param5, int param6):
      chooser_(chooser),
      option_(8),
      int_int_int_double_int_int_func_(hypreMapIntIntIntDoubleIntIntFunc_.at(funct_name)),
      int_param1_(param1),
      int_param2_(param2),
      int_param3_(param3),
      int_param4_(param5),
      int_param5_(param6),
      double_param1_(param4) {}

    //! Char pointer constructor.
    FunctionParameter(Hypre_Chooser chooser, char_star_func funct, char *param1):
      chooser_(chooser),
      option_(9),
      char_star_func_(funct),
      char_star_param_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, char *param1):
      chooser_(chooser),
      option_(9),
      char_star_func_(hypreMapCharStarFunc_.at(funct_name)),
      char_star_param_(param1) {}

  //! Only method of this class. Calls the function pointer with the passed in HYPRE_Solver
  int CallFunction(HYPRE_Solver solver, HYPRE_Solver precond) {
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
      } else if (option_ == 5){
        return double_star_func_(solver, double_star_param_);
      } else if (option_ == 6) {
        return int_int_double_double_func_(solver, int_param1_, int_param2_, double_param1_, double_param2_);
      } else if (option_ == 7) {
        return int_star_star_func_(solver, int_star_star_param_);
      } else if (option_ == 8) {
        return int_int_int_double_int_int_func_(solver, int_param1_, int_param2_, int_param3_, double_param1_, int_param4_, int_param5_);
      } else if (option_ == 9) {
        return char_star_func_(solver, char_star_param_);
      } else {
        IFPACK_CHK_ERR(-2);
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
      } else if(option_ == 5) {
        return double_star_func_(precond, double_star_param_);
      } else if (option_ == 6) {
        return int_int_double_double_func_(precond, int_param1_, int_param2_, double_param1_, double_param2_);
      } else if (option_ == 7) {
        return int_star_star_func_(precond, int_star_star_param_);
      } else if (option_ == 8) {
        return int_int_int_double_int_int_func_(precond, int_param1_, int_param2_, int_param3_, double_param1_, int_param4_, int_param5_);
      } else if (option_ == 9) {
        return char_star_func_(solver, char_star_param_);
      } else {
        IFPACK_CHK_ERR(-2);
      }
    }
  }
  
   static bool isFuncIntInt(std::string funct_name) {
    return (hypreMapIntIntFunc_.find(funct_name) != hypreMapIntIntFunc_.end());
  }
  
  static bool isFuncIntIntDoubleDouble(std::string funct_name) {
    return (hypreMapIntIntDoubleDoubleFunc_.find(funct_name) != hypreMapIntIntDoubleDoubleFunc_.end());
  }
  
  static bool isFuncIntIntIntDoubleIntInt(std::string funct_name) {
    return (hypreMapIntIntIntDoubleIntIntFunc_.find(funct_name) != hypreMapIntIntIntDoubleIntIntFunc_.end());
             }
  
  static bool isFuncIntStarStar(std::string funct_name) {
    return (hypreMapIntStarStarFunc_.find(funct_name) != hypreMapIntStarStarFunc_.end());
  }
  
  private:
    Hypre_Chooser chooser_;
    int option_;
    int_func int_func_;
    double_func double_func_;
    double_int_func double_int_func_;
    int_int_func int_int_func_;
    int_star_func int_star_func_;
    double_star_func double_star_func_;
    int_int_double_double_func int_int_double_double_func_;
    int_int_int_double_int_int_func int_int_int_double_int_int_func_;
    int_star_star_func int_star_star_func_;
    char_star_func char_star_func_;
    int int_param1_;
    int int_param2_;
    int int_param3_;
    int int_param4_;
    int int_param5_;
    double double_param1_;
    double double_param2_;
    int *int_star_param_;
    int **int_star_star_param_;
    double *double_star_param_;
    char *char_star_param_;

  static const std::map<std::string, int_func> hypreMapIntFunc_;
  static const std::map<std::string, double_func> hypreMapDoubleFunc_;
  static const std::map<std::string, double_int_func> hypreMapDoubleIntFunc_;
  static const std::map<std::string, int_int_func> hypreMapIntIntFunc_;
  static const std::map<std::string, int_star_func> hypreMapIntStarFunc_;
  static const std::map<std::string, double_star_func> hypreMapDoubleStarFunc_;
  static const std::map<std::string, int_int_double_double_func> hypreMapIntIntDoubleDoubleFunc_;
  static const std::map<std::string, int_int_int_double_int_int_func> hypreMapIntIntIntDoubleIntIntFunc_;
  static const std::map<std::string, int_star_star_func> hypreMapIntStarStarFunc_;
  static const std::map<std::string, char_star_func> hypreMapCharStarFunc_;

};

// NOTE: This really, really needs to be here and not up above, so please don't move it
#include "Ifpack_HypreParameterMap.h"

Ifpack_Hypre::Ifpack_Hypre(Epetra_RowMatrix* A):
  A_(rcp(A,false)),
  UseTranspose_(false),
  IsInitialized_(false),
  IsComputed_(false),
  Label_(),
  NumInitialize_(0),
  NumCompute_(0),
  NumApplyInverse_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyInverseTime_(0.0),
  ComputeFlops_(0.0),
  ApplyInverseFlops_(0.0),
  Time_(A_->Comm()),
  HypreA_(0),
  HypreG_(0),
  xHypre_(0),
  yHypre_(0),
  zHypre_(0),
  IsSolverCreated_(false),
  IsPrecondCreated_(false),
  SolveOrPrec_(Solver),
  NumFunsToCall_(0),
  SolverType_(PCG),
  PrecondType_(Euclid),
  UsePreconditioner_(false),
  Dump_(false)
{
  MPI_Comm comm = GetMpiComm();
  // Check that RowMap and RangeMap are the same.  While this could handle the
  // case where RowMap and RangeMap are permutations, other Ifpack PCs don't
  // handle this either.
  if (!A_->RowMatrixRowMap().SameAs(A_->OperatorRangeMap())) {
    IFPACK_CHK_ERRV(-1);
  }
  // Hypre expects the RowMap to be Linear.
  if (A_->RowMatrixRowMap().LinearMap()) {
    // note these are non-owning pointers, they are deleted by A_'s destructor
    GloballyContiguousRowMap_ = rcpFromRef(A_->RowMatrixRowMap());
    GloballyContiguousColMap_ = rcpFromRef(A_->RowMatrixColMap());
  } else {  
    // Must create GloballyContiguous Maps for Hypre
    if(A_->OperatorDomainMap().SameAs(A_->RowMatrixRowMap())) {
      Teuchos::RCP<const Epetra_RowMatrix> Aconst = A_;
      GloballyContiguousColMap_ = MakeContiguousColumnMap(Aconst);
      GloballyContiguousRowMap_ = rcp(new Epetra_Map(A_->RowMatrixRowMap().NumGlobalElements(),
                                                     A_->RowMatrixRowMap().NumMyElements(), 0, Comm()));
    }
    else {
      throw std::runtime_error("Ifpack_Hypre: Unsupported map configuration: Row/Domain maps do not match");
    }
  }
  // Next create vectors that will be used when ApplyInverse() is called
  int ilower = GloballyContiguousRowMap_->MinMyGID();
  int iupper = GloballyContiguousRowMap_->MaxMyGID();
  // X in AX = Y
  IFPACK_CHK_ERRV(HYPRE_IJVectorCreate(comm, ilower, iupper, &XHypre_));
  IFPACK_CHK_ERRV(HYPRE_IJVectorSetObjectType(XHypre_, HYPRE_PARCSR));
  IFPACK_CHK_ERRV(HYPRE_IJVectorInitialize(XHypre_));
  IFPACK_CHK_ERRV(HYPRE_IJVectorAssemble(XHypre_));
  IFPACK_CHK_ERRV(HYPRE_IJVectorGetObject(XHypre_, (void**) &ParX_));
  XVec_ = Teuchos::rcp((hypre_ParVector *) hypre_IJVectorObject(((hypre_IJVector *) XHypre_)),false);

  // Y in AX = Y
  IFPACK_CHK_ERRV(HYPRE_IJVectorCreate(comm, ilower, iupper, &YHypre_));
  IFPACK_CHK_ERRV(HYPRE_IJVectorSetObjectType(YHypre_, HYPRE_PARCSR));
  IFPACK_CHK_ERRV(HYPRE_IJVectorInitialize(YHypre_));
  IFPACK_CHK_ERRV(HYPRE_IJVectorAssemble(YHypre_));
  IFPACK_CHK_ERRV(HYPRE_IJVectorGetObject(YHypre_, (void**) &ParY_));
  YVec_ = Teuchos::rcp((hypre_ParVector *) hypre_IJVectorObject(((hypre_IJVector *) YHypre_)),false);

  // Cache
  VectorCache_.resize(A->RowMatrixRowMap().NumMyElements());
} //Constructor

//==============================================================================
void Ifpack_Hypre::Destroy(){
  if(IsInitialized()){
    IFPACK_CHK_ERRV(HYPRE_IJMatrixDestroy(HypreA_));
  }
  IFPACK_CHK_ERRV(HYPRE_IJVectorDestroy(XHypre_));
  IFPACK_CHK_ERRV(HYPRE_IJVectorDestroy(YHypre_));
  if(IsSolverCreated_){
    IFPACK_CHK_ERRV(SolverDestroyPtr_(Solver_));
  }
  if(IsPrecondCreated_){
    IFPACK_CHK_ERRV(PrecondDestroyPtr_(Preconditioner_));
  }

  // Maxwell
  if(HypreG_) {
    IFPACK_CHK_ERRV(HYPRE_IJMatrixDestroy(HypreG_));
  }
  if(xHypre_) {
    IFPACK_CHK_ERRV(HYPRE_IJVectorDestroy(xHypre_));
  }
  if(yHypre_) {
    IFPACK_CHK_ERRV(HYPRE_IJVectorDestroy(yHypre_));
  }
  if(zHypre_) {
    IFPACK_CHK_ERRV(HYPRE_IJVectorDestroy(zHypre_));
  }
} //Destroy()

//==============================================================================
int Ifpack_Hypre::Initialize(){
  Time_.ResetStartTime();
  if(IsInitialized_) return 0;
  // Create the Hypre matrix and copy values.  Note this uses values (which
  // Initialize() shouldn't do) but it doesn't care what they are (for
  // instance they can be uninitialized data even).  It should be possible to
  // set the Hypre structure without copying values, but this is the easiest
  // way to get the structure.
  MPI_Comm comm = GetMpiComm();
  int ilower = GloballyContiguousRowMap_->MinMyGID();
  int iupper = GloballyContiguousRowMap_->MaxMyGID();
  IFPACK_CHK_ERR(HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &HypreA_));
  IFPACK_CHK_ERR(HYPRE_IJMatrixSetObjectType(HypreA_, HYPRE_PARCSR));
  IFPACK_CHK_ERR(HYPRE_IJMatrixInitialize(HypreA_));
  CopyEpetraToHypre();
  if(SolveOrPrec_ == Solver) {
    IFPACK_CHK_ERR(SetSolverType(SolverType_));
    if (SolverPrecondPtr_ != NULL && UsePreconditioner_) {
      // both method allows a PC (first condition) and the user wants a PC (second)
      IFPACK_CHK_ERR(SetPrecondType(PrecondType_));
      CallFunctions();
      IFPACK_CHK_ERR(SolverPrecondPtr_(Solver_, PrecondSolvePtr_, PrecondSetupPtr_, Preconditioner_));
    } else {
      CallFunctions();
    }
  } else {
    IFPACK_CHK_ERR(SetPrecondType(PrecondType_));
    CallFunctions();
  }

  if (!G_.is_null()) {
    SetDiscreteGradient(G_);
  }

  if (!Coords_.is_null()) {
    SetCoordinates(Coords_);
  }

  // set flags
  IsInitialized_=true;
  NumInitialize_ = NumInitialize_ + 1;
  InitializeTime_ = InitializeTime_ + Time_.ElapsedTime();
  return 0;
} //Initialize()

//==============================================================================
int Ifpack_Hypre::SetParameters(Teuchos::ParameterList& list){

  std::map<std::string, Hypre_Solver> solverMap;
  solverMap["BoomerAMG"] = BoomerAMG;
  solverMap["ParaSails"] = ParaSails;
  solverMap["Euclid"] = Euclid;
  solverMap["AMS"] = AMS;
  solverMap["Hybrid"] = Hybrid;
  solverMap["PCG"] = PCG;
  solverMap["GMRES"] = GMRES;
  solverMap["FlexGMRES"] = FlexGMRES;
  solverMap["LGMRES"] = LGMRES;
  solverMap["BiCGSTAB"] = BiCGSTAB;

  std::map<std::string, Hypre_Chooser> chooserMap;
  chooserMap["Solver"] = Solver;
  chooserMap["Preconditioner"] = Preconditioner;

  List_ = list;
  Hypre_Solver solType;
  if (list.isType<std::string>("hypre: Solver"))
    solType = solverMap[list.get<std::string>("hypre: Solver")];
  else
    solType = list.get("hypre: Solver", PCG);
  SolverType_ = solType;
  Hypre_Solver precType;
  if (list.isType<std::string>("hypre: Preconditioner"))
    precType = solverMap[list.get<std::string>("hypre: Preconditioner")];
  else
    precType = list.get("hypre: Preconditioner", Euclid);
  PrecondType_ = precType;
  Hypre_Chooser chooser;
  if (list.isType<std::string>("hypre: SolveOrPrecondition"))
    chooser = chooserMap[list.get<std::string>("hypre: SolveOrPrecondition")];
  else
    chooser = list.get("hypre: SolveOrPrecondition", Solver);
  SolveOrPrec_ = chooser;
  bool SetPrecond = list.get("hypre: SetPreconditioner", false);
  IFPACK_CHK_ERR(SetParameter(SetPrecond));
  int NumFunctions = list.get("hypre: NumFunctions", 0);
  FunsToCall_.clear();
  NumFunsToCall_ = 0;
  if(NumFunctions > 0){
    RCP<FunctionParameter>* params = list.get<RCP<FunctionParameter>*>("hypre: Functions");
    for(int i = 0; i < NumFunctions; i++){
      IFPACK_CHK_ERR(AddFunToList(params[i]));
    }
  }

  if (list.isSublist("hypre: Solver functions")) {
    Teuchos::ParameterList solverList = list.sublist("hypre: Solver functions");
    for (auto it = solverList.begin(); it != solverList.end(); ++it) {
      std::string funct_name = it->first;
      if (it->second.isType<int>()) {
        IFPACK_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Solver, funct_name , Teuchos::getValue<int>(it->second)))));
      } else if (it->second.isType<double>()) {
        IFPACK_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Solver, funct_name , Teuchos::getValue<double>(it->second)))));
      } else {
        IFPACK_CHK_ERR(-1);
      }
    }
  }

  if (list.isSublist("hypre: Preconditioner functions")) {
    Teuchos::ParameterList precList = list.sublist("hypre: Preconditioner functions");
    for (auto it = precList.begin(); it != precList.end(); ++it) {
      std::string funct_name = it->first;
      if (it->second.isType<int>()) {
        IFPACK_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Preconditioner, funct_name , Teuchos::getValue<int>(it->second)))));
      } else if (it->second.isType<double>()) {
        IFPACK_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Preconditioner, funct_name , Teuchos::getValue<double>(it->second)))));
      } else if (it->second.isList()) {
        Teuchos::ParameterList pl = Teuchos::getValue<Teuchos::ParameterList>(it->second);
        if (FunctionParameter::isFuncIntInt(funct_name)) {
          int arg0 = pl.get<int>("arg 0");
          int arg1 = pl.get<int>("arg 1");
          IFPACK_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Preconditioner, funct_name , arg0, arg1))));
        } else if (FunctionParameter::isFuncIntIntDoubleDouble(funct_name)) {
          int arg0 = pl.get<int>("arg 0");
          int arg1 = pl.get<int>("arg 1");
          double arg2 = pl.get<double>("arg 2");
          double arg3 = pl.get<double>("arg 3");
          IFPACK_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Preconditioner, funct_name , arg0, arg1, arg2, arg3))));
        } else if (FunctionParameter::isFuncIntIntIntDoubleIntInt(funct_name)) {
          int arg0 = pl.get<int>("arg 0");
          int arg1 = pl.get<int>("arg 1");
          int arg2 = pl.get<int>("arg 2");
          double arg3 = pl.get<double>("arg 3");
          int arg4 = pl.get<int>("arg 4");
          int arg5 = pl.get<int>("arg 5");
          IFPACK_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Preconditioner, funct_name , arg0, arg1, arg2, arg3, arg4, arg5))));
        } else {
          IFPACK_CHK_ERR(-1);
        }
      }
    }
  }

  if (list.isSublist("Coordinates") && list.sublist("Coordinates").isType<Teuchos::RCP<Epetra_MultiVector> >("Coordinates"))
    Coords_ = list.sublist("Coordinates").get<Teuchos::RCP<Epetra_MultiVector> >("Coordinates");
  if (list.isSublist("Operators") && list.sublist("Operators").isType<Teuchos::RCP<const Epetra_CrsMatrix> >("G"))
    G_ = list.sublist("Operators").get<Teuchos::RCP<const Epetra_CrsMatrix> >("G");

  Dump_ = list.get("hypre: Dump", false);

  return 0;
} //SetParameters()

//==============================================================================
int Ifpack_Hypre::AddFunToList(RCP<FunctionParameter> NewFun){
  NumFunsToCall_ = NumFunsToCall_+1;
  FunsToCall_.resize(NumFunsToCall_);
  FunsToCall_[NumFunsToCall_-1] = NewFun;
  return 0;
} //AddFunToList()

//==============================================================================
int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int), int parameter){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter));
  IFPACK_CHK_ERR(AddFunToList(temp));
  return 0;
} //SetParameter() - int function pointer

//==============================================================================
int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double), double parameter){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter));
  IFPACK_CHK_ERR(AddFunToList(temp));
  return 0;
} //SetParameter() - double function pointer

//==============================================================================
int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double, int), double parameter1, int parameter2){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter1, parameter2));
  IFPACK_CHK_ERR(AddFunToList(temp));
  return 0;
} //SetParameter() - double,int function pointer

//==============================================================================
int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int, int), int parameter1, int parameter2){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter1, parameter2));
  IFPACK_CHK_ERR(AddFunToList(temp));
  return 0;
} //SetParameter() int,int function pointer

//==============================================================================
int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, double*), double* parameter){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter));
  IFPACK_CHK_ERR(AddFunToList(temp));
  return 0;
} //SetParameter() - double* function pointer

//==============================================================================
int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int*), int* parameter){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter));
  IFPACK_CHK_ERR(AddFunToList(temp));
  return 0;
} //SetParameter() - int* function pointer

//==============================================================================
int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, int (*pt2Func)(HYPRE_Solver, int**), int** parameter){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter));
  IFPACK_CHK_ERR(AddFunToList(temp));
  return 0;
} //SetParameter() - int** function pointer

//==============================================================================
int Ifpack_Hypre::SetParameter(Hypre_Chooser chooser, Hypre_Solver solver){
  if(chooser == Solver){
    SolverType_ = solver;
  } else {
    PrecondType_ = solver;
  }
  return 0;
} //SetParameter() - set type of solver

//==============================================================================
int Ifpack_Hypre::SetDiscreteGradient(Teuchos::RCP<const Epetra_CrsMatrix> G){

  // Sanity check
  if(!A_->RowMatrixRowMap().SameAs(G->RowMap()))
    throw std::runtime_error("Ifpack_Hypre: Edge map mismatch: A and discrete gradient");

  // Get the maps for the nodes (assuming the edge map from A is OK);
  GloballyContiguousNodeRowMap_ = rcp(new Epetra_Map(G->DomainMap().NumGlobalElements(),
                                                     G->DomainMap().NumMyElements(), 0, Comm()));
  Teuchos::RCP<const Epetra_RowMatrix> Grow = Teuchos::rcp_dynamic_cast<const Epetra_RowMatrix>(G);
  GloballyContiguousNodeColMap_ = MakeContiguousColumnMap(Grow);

  // Start building G
  MPI_Comm comm = GetMpiComm();
  int ilower = GloballyContiguousRowMap_->MinMyGID();
  int iupper = GloballyContiguousRowMap_->MaxMyGID();
  int jlower = GloballyContiguousNodeRowMap_->MinMyGID();
  int jupper = GloballyContiguousNodeRowMap_->MaxMyGID();
  IFPACK_CHK_ERR(HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, &HypreG_));
  IFPACK_CHK_ERR(HYPRE_IJMatrixSetObjectType(HypreG_, HYPRE_PARCSR));
  IFPACK_CHK_ERR(HYPRE_IJMatrixInitialize(HypreG_));

  std::vector<int> new_indices(G->MaxNumEntries());
  for(int i = 0; i < G->NumMyRows(); i++){
    int numEntries;
    double * values;
    int *indices;
    IFPACK_CHK_ERR(G->ExtractMyRowView(i, numEntries, values, indices));
    for(int j = 0; j < numEntries; j++){
      new_indices[j] = GloballyContiguousNodeColMap_->GID(indices[j]);
    }
    int GlobalRow[1];
    GlobalRow[0] = GloballyContiguousRowMap_->GID(i);
    IFPACK_CHK_ERR(HYPRE_IJMatrixSetValues(HypreG_, 1, &numEntries, GlobalRow, new_indices.data(), values));
  }
  IFPACK_CHK_ERR(HYPRE_IJMatrixAssemble(HypreG_));
  IFPACK_CHK_ERR(HYPRE_IJMatrixGetObject(HypreG_, (void**)&ParMatrixG_));

  if (Dump_)
    HYPRE_ParCSRMatrixPrint(ParMatrixG_,"G.mat");

  if(SolverType_ == AMS)
    HYPRE_AMSSetDiscreteGradient(Solver_, ParMatrixG_);
  if(PrecondType_ == AMS)
    HYPRE_AMSSetDiscreteGradient(Preconditioner_, ParMatrixG_);

  return 0;
} //SetDiscreteGradient()

//==============================================================================
int Ifpack_Hypre::SetCoordinates(Teuchos::RCP<Epetra_MultiVector> coords) {

  if(!G_.is_null() && !G_->DomainMap().SameAs(coords->Map()))
    throw std::runtime_error("Ifpack_Hypre: Node map mismatch: G->DomainMap() and coords");

  if(SolverType_ != AMS && PrecondType_ != AMS)
    return 0;

  double *xPtr;
  double *yPtr;
  double *zPtr;

  IFPACK_CHK_ERR(((*coords)(0))->ExtractView(&xPtr));
  IFPACK_CHK_ERR(((*coords)(1))->ExtractView(&yPtr));
  IFPACK_CHK_ERR(((*coords)(2))->ExtractView(&zPtr));

  MPI_Comm comm = GetMpiComm();
  int NumEntries = coords->MyLength();
  int * indices = GloballyContiguousNodeRowMap_->MyGlobalElements();

  int ilower = GloballyContiguousNodeRowMap_->MinMyGID();
  int iupper = GloballyContiguousNodeRowMap_->MaxMyGID();

  if( NumEntries != iupper-ilower+1) {
    std::cout<<"Ifpack_Hypre::SetCoordinates(): Error on rank "<<Comm().MyPID()<<": MyLength = "<<coords->MyLength()<<" GID range = ["<<ilower<<","<<iupper<<"]"<<std::endl;
    throw std::runtime_error("Ifpack_Hypre: SetCoordinates: Length mismatch");
  }

  IFPACK_CHK_ERR(HYPRE_IJVectorCreate(comm, ilower, iupper, &xHypre_));
  IFPACK_CHK_ERR(HYPRE_IJVectorSetObjectType(xHypre_, HYPRE_PARCSR));
  IFPACK_CHK_ERR(HYPRE_IJVectorInitialize(xHypre_));
  IFPACK_CHK_ERR(HYPRE_IJVectorSetValues(xHypre_,NumEntries,indices,xPtr));
  IFPACK_CHK_ERR(HYPRE_IJVectorAssemble(xHypre_));
  IFPACK_CHK_ERR(HYPRE_IJVectorGetObject(xHypre_, (void**) &xPar_));

  IFPACK_CHK_ERR(HYPRE_IJVectorCreate(comm, ilower, iupper, &yHypre_));
  IFPACK_CHK_ERR(HYPRE_IJVectorSetObjectType(yHypre_, HYPRE_PARCSR));
  IFPACK_CHK_ERR(HYPRE_IJVectorInitialize(yHypre_));
  IFPACK_CHK_ERR(HYPRE_IJVectorSetValues(yHypre_,NumEntries,indices,yPtr));
  IFPACK_CHK_ERR(HYPRE_IJVectorAssemble(yHypre_));
  IFPACK_CHK_ERR(HYPRE_IJVectorGetObject(yHypre_, (void**) &yPar_));

  IFPACK_CHK_ERR(HYPRE_IJVectorCreate(comm, ilower, iupper, &zHypre_));
  IFPACK_CHK_ERR(HYPRE_IJVectorSetObjectType(zHypre_, HYPRE_PARCSR));
  IFPACK_CHK_ERR(HYPRE_IJVectorInitialize(zHypre_));
  IFPACK_CHK_ERR(HYPRE_IJVectorSetValues(zHypre_,NumEntries,indices,zPtr));
  IFPACK_CHK_ERR(HYPRE_IJVectorAssemble(zHypre_));
  IFPACK_CHK_ERR(HYPRE_IJVectorGetObject(zHypre_, (void**) &zPar_));

  if (Dump_) {
    HYPRE_ParVectorPrint(xPar_,"coordX.dat");
    HYPRE_ParVectorPrint(yPar_,"coordY.dat");
    HYPRE_ParVectorPrint(zPar_,"coordZ.dat");
  }

  if(SolverType_ == AMS)
    HYPRE_AMSSetCoordinateVectors(Solver_, xPar_, yPar_, zPar_);
  if(PrecondType_ == AMS)
    HYPRE_AMSSetCoordinateVectors(Preconditioner_, xPar_, yPar_, zPar_);

  return 0;

} //SetCoordinates

//==============================================================================
int Ifpack_Hypre::Compute(){
  if(IsInitialized() == false){
    IFPACK_CHK_ERR(Initialize());
  }
  Time_.ResetStartTime();

  // Hypre Setup must be called after matrix has values
  if(SolveOrPrec_ == Solver){
    IFPACK_CHK_ERR(SolverSetupPtr_(Solver_, ParMatrix_, ParX_, ParY_));
  } else {
    IFPACK_CHK_ERR(PrecondSetupPtr_(Preconditioner_, ParMatrix_, ParX_, ParY_));
  }

  IsComputed_ = true;
  NumCompute_ = NumCompute_ + 1;
  ComputeTime_ = ComputeTime_ + Time_.ElapsedTime();
  return 0;
} //Compute()

//==============================================================================
int Ifpack_Hypre::CallFunctions() const{
  for(int i = 0; i < NumFunsToCall_; i++){
    IFPACK_CHK_ERR(FunsToCall_[i]->CallFunction(Solver_, Preconditioner_));
  }
  return 0;
} //CallFunctions()

//==============================================================================
int Ifpack_Hypre::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{
  if(IsComputed() == false){
    IFPACK_CHK_ERR(-1);
  }
  Time_.ResetStartTime();
  hypre_Vector *XLocal_ = hypre_ParVectorLocalVector(XVec_);
  hypre_Vector *YLocal_ = hypre_ParVectorLocalVector(YVec_);

  bool SameVectors = false;
  int NumVectors = X.NumVectors();
  if (NumVectors != Y.NumVectors()) IFPACK_CHK_ERR(-1);  // X and Y must have same number of vectors
  if(X.Pointers() == Y.Pointers() || (NumVectors == 1 && X[0] == Y[0])){
    SameVectors = true;
  }
  
  // NOTE: Here were assuming that the local ordering of Epetra's X/Y-vectors and 
  // Hypre's X/Y-vectors are the same.  Seeing as as this is more or less how we constructed
  // the Hypre matrices, this seems pretty reasoanble.

  for(int VecNum = 0; VecNum < NumVectors; VecNum++) {
    //Get values for current vector in multivector.
    // FIXME amk Nov 23, 2015: This will not work for funky data layouts
    double * XValues = const_cast<double*>(X[VecNum]);
    double * YValues;
    if(!SameVectors){
      YValues = const_cast<double*>(Y[VecNum]);
    } else {
      YValues = VectorCache_.getRawPtr();
    }
    // Temporarily make a pointer to data in Hypre for end
    double *XTemp = XLocal_->data;
    // Replace data in Hypre vectors with Epetra data
    XLocal_->data = XValues;
    double *YTemp = YLocal_->data;
    YLocal_->data = YValues;

    IFPACK_CHK_ERR(HYPRE_ParVectorSetConstantValues(ParY_, 0.0));
    if(SolveOrPrec_ == Solver){
      // Use the solver methods
      IFPACK_CHK_ERR(SolverSolvePtr_(Solver_, ParMatrix_, ParX_, ParY_));
    } else {
      // Apply the preconditioner
      IFPACK_CHK_ERR(PrecondSolvePtr_(Preconditioner_, ParMatrix_, ParX_, ParY_));
    }
    if(SameVectors){
      int NumEntries = Y.MyLength();
      for(int i = 0; i < NumEntries; i++)
	Y[VecNum][i] = YValues[i];      
    }
    XLocal_->data = XTemp;
    YLocal_->data = YTemp;
  }
  NumApplyInverse_ = NumApplyInverse_ + 1;
  ApplyInverseTime_ = ApplyInverseTime_ + Time_.ElapsedTime();
  return 0;
} //ApplyInverse()

//==============================================================================
int Ifpack_Hypre::Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const{
  if(IsInitialized() == false){
    IFPACK_CHK_ERR(-1);
  }
  hypre_Vector *XLocal_ = hypre_ParVectorLocalVector(XVec_);
  hypre_Vector *YLocal_ = hypre_ParVectorLocalVector(YVec_);
  bool SameVectors = false;
  int NumVectors = X.NumVectors();
  if (NumVectors != Y.NumVectors()) IFPACK_CHK_ERR(-1);  // X and Y must have same number of vectors
  if(X.Pointers() == Y.Pointers() || (NumVectors == 1 && X[0] == Y[0])){
    SameVectors = true;
  }

  // NOTE: Here were assuming that the local ordering of Epetra's X/Y-vectors and 
  // Hypre's X/Y-vectors are the same.  Seeing as as this is more or less how we constructed
  // the Hypre matrices, this seems pretty reasoanble.

  for(int VecNum = 0; VecNum < NumVectors; VecNum++) {
    //Get values for current vector in multivector.
    double * XValues=const_cast<double*>(X[VecNum]);
    double * YValues;
    double *XTemp = XLocal_->data;
    double *YTemp = YLocal_->data;
    if(!SameVectors){
      YValues = const_cast<double*>(Y[VecNum]);
    } else {
      YValues = VectorCache_.getRawPtr();
    }
    YLocal_->data = YValues;
    IFPACK_CHK_ERR(HYPRE_ParVectorSetConstantValues(ParY_,0.0));
    // Temporarily make a pointer to data in Hypre for end
    // Replace data in Hypre vectors with epetra values
    XLocal_->data = XValues;
    // Do actual computation.
    if(TransA) {
      // Use transpose of A in multiply
      IFPACK_CHK_ERR(HYPRE_ParCSRMatrixMatvecT(1.0, ParMatrix_, ParX_, 1.0, ParY_));
    } else {
      IFPACK_CHK_ERR(HYPRE_ParCSRMatrixMatvec(1.0, ParMatrix_, ParX_, 1.0, ParY_));
    }
    if(SameVectors){
      int NumEntries = Y.MyLength();
      for(int i = 0; i < NumEntries; i++)
	Y[VecNum][i] = YValues[i];      
    }
    XLocal_->data = XTemp;
    YLocal_->data = YTemp;
  }
  return 0;
} //Multiply()

//==============================================================================
std::ostream& Ifpack_Hypre::Print(std::ostream& os) const{
  using std::endl;
  if (!Comm().MyPID()) {
    os << endl;
    os << "================================================================================" << endl;
    os << "Ifpack_Hypre: " << Label() << endl << endl;
    os << "Using " << Comm().NumProc() << " processors." << endl;
    os << "Global number of rows            = " << A_->NumGlobalRows() << endl;
    os << "Global number of nonzeros        = " << A_->NumGlobalNonzeros() << endl;
    os << "Condition number estimate = " << Condest() << endl;
    os << endl;
    os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
    os << "-----           -------   --------------       ------------     --------" << endl;
    os << "Initialize()    "   << std::setw(5) << NumInitialize_
       << "  " << std::setw(15) << InitializeTime_
       << "              0.0              0.0" << endl;
    os << "Compute()       "   << std::setw(5) << NumCompute_
       << "  " << std::setw(15) << ComputeTime_
       << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_;
    if (ComputeTime_ != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ComputeFlops_ / ComputeTime_ << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "ApplyInverse()  "   << std::setw(5) << NumApplyInverse_
       << "  " << std::setw(15) << ApplyInverseTime_
       << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_;
    if (ApplyInverseTime_ != 0.0)
      os << "  " << std::setw(15) << 1.0e-6 * ApplyInverseFlops_ / ApplyInverseTime_ << endl;
    else
      os << "  " << std::setw(15) << 0.0 << endl;
    os << "================================================================================" << endl;
    os << endl;
  }
  return os;
} //Print()

//==============================================================================
double Ifpack_Hypre::Condest(const Ifpack_CondestType CT,
                             const int MaxIters,
                             const double Tol,
                             Epetra_RowMatrix* Matrix_in){
  if (!IsComputed()) // cannot compute right now
    return(-1.0);
  Condest_ = Ifpack_Condest(*this, CT, MaxIters, Tol, Matrix_in);
  return(Condest_);
} //Condest()

//==============================================================================
int Ifpack_Hypre::SetSolverType(Hypre_Solver Solver){
  switch(Solver) {
    case BoomerAMG:
      if(IsSolverCreated_){
        SolverDestroyPtr_(Solver_);
        IsSolverCreated_ = false;
      }
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_BoomerAMGCreate;
      SolverDestroyPtr_ = &HYPRE_BoomerAMGDestroy;
      SolverSetupPtr_ = &HYPRE_BoomerAMGSetup;
      SolverPrecondPtr_ = NULL;
      SolverSolvePtr_ = &HYPRE_BoomerAMGSolve;
      break;
    case AMS:
      if(IsSolverCreated_){
        SolverDestroyPtr_(Solver_);
        IsSolverCreated_ = false;
      }
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_AMSCreate;
      SolverDestroyPtr_ = &HYPRE_AMSDestroy;
      SolverSetupPtr_ = &HYPRE_AMSSetup;
      SolverSolvePtr_ = &HYPRE_AMSSolve;
      SolverPrecondPtr_ = NULL;
      break;
    case Hybrid:
      if(IsSolverCreated_){
        SolverDestroyPtr_(Solver_);
        IsSolverCreated_ = false;
      }
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_ParCSRHybridCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRHybridDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRHybridSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRHybridSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRHybridSetPrecond;
      break;
    case PCG:
      if(IsSolverCreated_){
        SolverDestroyPtr_(Solver_);
        IsSolverCreated_ = false;
      }
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_ParCSRPCGCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRPCGDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRPCGSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRPCGSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRPCGSetPrecond;
      break;
    case GMRES:
      if(IsSolverCreated_){
        SolverDestroyPtr_(Solver_);
        IsSolverCreated_ = false;
      }
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_ParCSRGMRESCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRGMRESDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRGMRESSetup;
      SolverPrecondPtr_ = &HYPRE_ParCSRGMRESSetPrecond;
      break;
    case FlexGMRES:
      if(IsSolverCreated_){
        SolverDestroyPtr_(Solver_);
        IsSolverCreated_ = false;
      }
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_ParCSRFlexGMRESCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRFlexGMRESDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRFlexGMRESSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRFlexGMRESSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRFlexGMRESSetPrecond;
      break;
    case LGMRES:
      if(IsSolverCreated_){
        SolverDestroyPtr_(Solver_);
        IsSolverCreated_ = false;
      }
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_ParCSRLGMRESCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRLGMRESDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRLGMRESSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRLGMRESSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRLGMRESSetPrecond;
      break;
    case BiCGSTAB:
      if(IsSolverCreated_){
        SolverDestroyPtr_(Solver_);
        IsSolverCreated_ = false;
      }
      SolverCreatePtr_ = &Ifpack_Hypre::Hypre_ParCSRBiCGSTABCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRBiCGSTABDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRBiCGSTABSetup;
      SolverSolvePtr_ = &HYPRE_ParCSRBiCGSTABSolve;
      SolverPrecondPtr_ = &HYPRE_ParCSRBiCGSTABSetPrecond;
      break;
    default:
      return -1;
    }
  CreateSolver();
  return 0;
} //SetSolverType()

//==============================================================================
int Ifpack_Hypre::SetPrecondType(Hypre_Solver Precond){
  switch(Precond) {
    case BoomerAMG:
      if(IsPrecondCreated_){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondCreated_ = false;
      }
      PrecondCreatePtr_ = &Ifpack_Hypre::Hypre_BoomerAMGCreate;
      PrecondDestroyPtr_ = &HYPRE_BoomerAMGDestroy;
      PrecondSetupPtr_ = &HYPRE_BoomerAMGSetup;
      PrecondSolvePtr_ = &HYPRE_BoomerAMGSolve;
      break;
    case ParaSails:
      if(IsPrecondCreated_){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondCreated_ = false;
      }
      PrecondCreatePtr_ = &Ifpack_Hypre::Hypre_ParaSailsCreate;
      PrecondDestroyPtr_ = &HYPRE_ParaSailsDestroy;
      PrecondSetupPtr_ = &HYPRE_ParaSailsSetup;
      PrecondSolvePtr_ = &HYPRE_ParaSailsSolve;
      break;
    case Euclid:
      if(IsPrecondCreated_){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondCreated_ = false;
      }
      PrecondCreatePtr_ = &Ifpack_Hypre::Hypre_EuclidCreate;
      PrecondDestroyPtr_ = &HYPRE_EuclidDestroy;
      PrecondSetupPtr_ = &HYPRE_EuclidSetup;
      PrecondSolvePtr_ = &HYPRE_EuclidSolve;
      break;
    case AMS:
      if(IsPrecondCreated_){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondCreated_ = false;
      }
      PrecondCreatePtr_ = &Ifpack_Hypre::Hypre_AMSCreate;
      PrecondDestroyPtr_ = &HYPRE_AMSDestroy;
      PrecondSetupPtr_ = &HYPRE_AMSSetup;
      PrecondSolvePtr_ = &HYPRE_AMSSolve;
      break;
    default:
      return -1;
    }
  CreatePrecond();
  return 0;

} //SetPrecondType()

//==============================================================================
int Ifpack_Hypre::CreateSolver(){
  MPI_Comm comm;
  HYPRE_ParCSRMatrixGetComm(ParMatrix_, &comm);
  int ierr = (this->*SolverCreatePtr_)(comm, &Solver_);
  IsSolverCreated_ = true;
  return ierr;
} //CreateSolver()

//==============================================================================
int Ifpack_Hypre::CreatePrecond(){
  MPI_Comm comm;
  HYPRE_ParCSRMatrixGetComm(ParMatrix_, &comm);
  int ierr = (this->*PrecondCreatePtr_)(comm, &Preconditioner_);
  IsPrecondCreated_ = true;
  return ierr;
} //CreatePrecond()


//==============================================================================
int Ifpack_Hypre::CopyEpetraToHypre(){
  Teuchos::RCP<const Epetra_CrsMatrix> Matrix = Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(A_);
  if(Matrix.is_null()) 
    throw std::runtime_error("Ifpack_Hypre: Unsupported matrix configuration: Epetra_CrsMatrix required");

  std::vector<int> new_indices(Matrix->MaxNumEntries());
  for(int i = 0; i < Matrix->NumMyRows(); i++){
    int numEntries;
    int *indices;
    double *values;
    IFPACK_CHK_ERR(Matrix->ExtractMyRowView(i, numEntries, values, indices));
    for(int j = 0; j < numEntries; j++){
      new_indices[j] = GloballyContiguousColMap_->GID(indices[j]);
    }
    int GlobalRow[1];
    GlobalRow[0] = GloballyContiguousRowMap_->GID(i);
    IFPACK_CHK_ERR(HYPRE_IJMatrixSetValues(HypreA_, 1, &numEntries, GlobalRow, new_indices.data(), values));
  }
  IFPACK_CHK_ERR(HYPRE_IJMatrixAssemble(HypreA_));
  IFPACK_CHK_ERR(HYPRE_IJMatrixGetObject(HypreA_, (void**)&ParMatrix_));
  if (Dump_)
    HYPRE_ParCSRMatrixPrint(ParMatrix_,"A.mat");
  return 0;
} //CopyEpetraToHypre()

//==============================================================================
int Ifpack_Hypre::Hypre_BoomerAMGCreate(MPI_Comm /*comm*/, HYPRE_Solver *solver)
    { return HYPRE_BoomerAMGCreate(solver);}

//==============================================================================
int Ifpack_Hypre::Hypre_ParaSailsCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParaSailsCreate(comm, solver);}

//==============================================================================
int Ifpack_Hypre::Hypre_EuclidCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_EuclidCreate(comm, solver);}

//==============================================================================
int Ifpack_Hypre::Hypre_AMSCreate(MPI_Comm /*comm*/, HYPRE_Solver *solver)
    { return HYPRE_AMSCreate(solver);}

//==============================================================================
int Ifpack_Hypre::Hypre_ParCSRHybridCreate(MPI_Comm /*comm*/, HYPRE_Solver *solver)
    { return HYPRE_ParCSRHybridCreate(solver);}

//==============================================================================
int Ifpack_Hypre::Hypre_ParCSRPCGCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRPCGCreate(comm, solver);}

//==============================================================================
int Ifpack_Hypre::Hypre_ParCSRGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRGMRESCreate(comm, solver);}

//==============================================================================
int Ifpack_Hypre::Hypre_ParCSRFlexGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRFlexGMRESCreate(comm, solver);}

//==============================================================================
int Ifpack_Hypre::Hypre_ParCSRLGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRLGMRESCreate(comm, solver);}

//==============================================================================
int Ifpack_Hypre::Hypre_ParCSRBiCGSTABCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRBiCGSTABCreate(comm, solver);}

//==============================================================================
Teuchos::RCP<const Epetra_Map> Ifpack_Hypre::MakeContiguousColumnMap(Teuchos::RCP<const Epetra_RowMatrix> &MatrixRow) const{
  // Must create GloballyContiguous DomainMap (which is a permutation of Matrix_'s
  // DomainMap) and the corresponding permuted ColumnMap.
  //   Epetra_GID  --------->   LID   ----------> HYPRE_GID
  //           via DomainMap.LID()       via GloballyContiguousDomainMap.GID()
  Teuchos::RCP<const Epetra_CrsMatrix> Matrix = Teuchos::rcp_dynamic_cast<const Epetra_CrsMatrix>(MatrixRow);
  if(Matrix.is_null()) 
    throw std::runtime_error("Ifpack_Hypre: Unsupported matrix configuration: Epetra_CrsMatrix required");
  const Epetra_Map & DomainMap = Matrix->DomainMap();
  const Epetra_Map & ColumnMap = Matrix->ColMap();
  const Epetra_Import * importer = Matrix->Importer();

  if(DomainMap.LinearMap() ) {
    // If the domain map is linear, then we can just use the column map as is.
    return rcpFromRef(ColumnMap);
  }
  else {
    // The domain map isn't linear, so we need a new domain map
    Teuchos::RCP<Epetra_Map> ContiguousDomainMap = rcp(new Epetra_Map(DomainMap.NumGlobalElements(),
                                                                      DomainMap.NumMyElements(), 0, Comm()));
    if(importer) {    
      // If there's an importer then we can use it to get a new column map
      Epetra_IntVector MyGIDsHYPRE(View,DomainMap,ContiguousDomainMap->MyGlobalElements());

      // import the HYPRE GIDs
      Epetra_IntVector ColGIDsHYPRE(ColumnMap);
      ColGIDsHYPRE.Import(MyGIDsHYPRE, *importer, Insert);
   
      // Make a HYPRE numbering-based column map.
      return Teuchos::rcp(new Epetra_Map(ColumnMap.NumGlobalElements(),ColGIDsHYPRE.MyLength(), &ColGIDsHYPRE[0], 0, Comm()));
    }
    else {
      // The problem has matching domain/column maps, and somehow the domain map isn't linear, so just use the new domain map
      return Teuchos::rcp(new Epetra_Map(ColumnMap.NumGlobalElements(),ColumnMap.NumMyElements(), ContiguousDomainMap->MyGlobalElements(), 0, Comm()));
    }
  }  
}



#endif // HAVE_HYPRE && HAVE_MPI
