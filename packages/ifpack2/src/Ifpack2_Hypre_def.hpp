/*@HEADER
// ***********************************************************************
//
//       Ifpack2:  Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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

#ifndef IFPACK2_HYPRE_DEF_HPP
#define IFPACK2_HYPRE_DEF_HPP

#include "Ifpack2_Hypre_decl.hpp"
#if defined(HAVE_IFPACK2_HYPRE) && defined(HAVE_IFPACK2_MPI)
#include <stdexcept>

#include "Tpetra_Import.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
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
typedef HYPRE_Int (*int_func)(HYPRE_Solver, HYPRE_Int);
typedef HYPRE_Int (*double_func)(HYPRE_Solver, double);
typedef HYPRE_Int (*double_int_func)(HYPRE_Solver, double, HYPRE_Int);
typedef HYPRE_Int (*int_int_func)(HYPRE_Solver, HYPRE_Int, HYPRE_Int);
typedef HYPRE_Int (*int_star_func)(HYPRE_Solver, HYPRE_Int*);
typedef HYPRE_Int (*int_star_star_func)(HYPRE_Solver, HYPRE_Int**);
typedef HYPRE_Int (*double_star_func)(HYPRE_Solver, double*);
typedef HYPRE_Int (*int_int_double_double_func)(HYPRE_Solver, HYPRE_Int, HYPRE_Int, double, double);
typedef HYPRE_Int (*int_int_int_double_int_int_func)(HYPRE_Solver, HYPRE_Int, HYPRE_Int, HYPRE_Int, double, HYPRE_Int, HYPRE_Int);
typedef HYPRE_Int (*char_star_func)(HYPRE_Solver, char*);

namespace Ifpack2 {

  void IFPACK2_CHK_ERRV(int code) {
    if(code<0) {
      std::ostringstream ofs;
      ofs << "Ifpack2::Hypre: Error with code "<<code<<std::endl;
      throw std::runtime_error(ofs.str());
    }
  }

  void IFPACK2_CHK_ERR(int code) {
    if(code<0) {
      std::ostringstream ofs;
      ofs << "Ifpack2::Hypre: Error with code "<<code<<std::endl;
      throw std::runtime_error(ofs.str());
    }
  }


//! This class is used to help with passing parameters in the SetParameter() function. Use this class to call Hypre's internal parameters.
class FunctionParameter{
  public:
    //! Single int constructor.
    FunctionParameter(Hypre_Chooser chooser, int_func funct, HYPRE_Int param1) :
      chooser_(chooser),
      option_(0),
      int_func_(funct),
      int_param1_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, HYPRE_Int param1) :
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
    FunctionParameter(Hypre_Chooser chooser, double_int_func funct, double param1, HYPRE_Int param2):
      chooser_(chooser),
      option_(2),
      double_int_func_(funct),
      int_param1_(param2),
      double_param1_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, double param1, HYPRE_Int param2):
      chooser_(chooser),
      option_(2),
      double_int_func_(hypreMapDoubleIntFunc_.at(funct_name)),
      int_param1_(param2),
      double_param1_(param1) {}

    //! Two ints constructor.
    FunctionParameter(Hypre_Chooser chooser, int_int_func funct, HYPRE_Int param1, HYPRE_Int param2):
      chooser_(chooser),
      option_(3),
      int_int_func_(funct),
      int_param1_(param1),
      int_param2_(param2) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, HYPRE_Int param1, HYPRE_Int param2):
      chooser_(chooser),
      option_(3),
      int_int_func_(hypreMapIntIntFunc_.at(funct_name)),
      int_param1_(param1),
      int_param2_(param2) {}

    //! Int pointer constructor.
    FunctionParameter(Hypre_Chooser chooser, int_star_func funct, HYPRE_Int *param1):
      chooser_(chooser),
      option_(4),
      int_star_func_(funct),
      int_star_param_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, HYPRE_Int *param1):
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
    FunctionParameter(Hypre_Chooser chooser, int_int_double_double_func funct, HYPRE_Int param1, HYPRE_Int param2, double param3, double param4):
      chooser_(chooser),
      option_(6),
      int_int_double_double_func_(funct),
      int_param1_(param1),
      int_param2_(param2),
      double_param1_(param3),
      double_param2_(param4) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, HYPRE_Int param1, HYPRE_Int param2, double param3, double param4):
      chooser_(chooser),
      option_(6),
      int_int_double_double_func_(hypreMapIntIntDoubleDoubleFunc_.at(funct_name)),
      int_param1_(param1),
      int_param2_(param2),
      double_param1_(param3),
      double_param2_(param4) {}

    //! Integer pointer to list of integer pointers
    FunctionParameter(Hypre_Chooser chooser, int_star_star_func funct, HYPRE_Int ** param1):
      chooser_(chooser),
      option_(7),
      int_star_star_func_(funct),
      int_star_star_param_(param1) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, HYPRE_Int** param1):
      chooser_(chooser),
      option_(7),
      int_star_star_func_(hypreMapIntStarStarFunc_.at(funct_name)),
      int_star_star_param_(param1) {}

    //! Five ints, one double constructor.
    FunctionParameter(Hypre_Chooser chooser, int_int_int_double_int_int_func funct, HYPRE_Int param1, HYPRE_Int param2, HYPRE_Int param3, double param4, HYPRE_Int param5, HYPRE_Int param6):
      chooser_(chooser),
      option_(8),
      int_int_int_double_int_int_func_(funct),
      int_param1_(param1),
      int_param2_(param2),
      int_param3_(param3),
      int_param4_(param5),
      int_param5_(param6),
      double_param1_(param4) {}

    FunctionParameter(Hypre_Chooser chooser, std::string funct_name, HYPRE_Int param1, HYPRE_Int param2, HYPRE_Int param3, double param4, HYPRE_Int param5, HYPRE_Int param6):
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
    if(chooser_ == Hypre_Is_Solver){
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
        IFPACK2_CHK_ERR(-2);
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
        IFPACK2_CHK_ERR(-2);
      }
    }
    return 0;
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
    HYPRE_Int int_param1_;
    HYPRE_Int int_param2_;
    HYPRE_Int int_param3_;
    HYPRE_Int int_param4_;
    HYPRE_Int int_param5_;
    double double_param1_;
    double double_param2_;
    HYPRE_Int *int_star_param_;
    HYPRE_Int **int_star_star_param_;
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
#include "Ifpack2_HypreParameterMap.hpp"


template<class MatrixType>
Hypre<MatrixType>::
Hypre(const Teuchos::RCP<const row_matrix_type>& A):
  A_(A),
  IsInitialized_(false),
  IsComputed_(false),  NumInitialize_(0),
  NumCompute_(0),
  NumApply_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyTime_(0.0),
  HypreA_(0),
  HypreG_(0),
  xHypre_(0),
  yHypre_(0),
  zHypre_(0),
  IsSolverCreated_(false),
  IsPrecondCreated_(false),
  SolveOrPrec_(Hypre_Is_Solver),
  NumFunsToCall_(0),
  SolverType_(PCG),
  PrecondType_(Euclid),
  UsePreconditioner_(false),
  Dump_(false)
{
  MPI_Comm comm = * (Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(A->getRowMap()->getComm())->getRawMpiComm());

  // Check that RowMap and RangeMap are the same.  While this could handle the
  // case where RowMap and RangeMap are permutations, other Ifpack PCs don't
  // handle this either.
  if (!A_->getRowMap()->isSameAs(*A_->getRangeMap())) {
    IFPACK2_CHK_ERRV(-1);
  }
  // Hypre expects the RowMap to be Linear.
  if (A_->getRowMap()->isContiguous()) {
    GloballyContiguousRowMap_ = A_->getRowMap();
    GloballyContiguousColMap_ = A_->getColMap();
  } else {  
    // Must create GloballyContiguous Maps for Hypre
    if(A_->getDomainMap()->isSameAs(*A_->getRowMap())) {
      Teuchos::RCP<const crs_matrix_type> Aconst = Teuchos::rcp_dynamic_cast<const crs_matrix_type>(A_);
      GloballyContiguousColMap_ = MakeContiguousColumnMap(Aconst);
      GloballyContiguousRowMap_ = rcp(new map_type(A_->getRowMap()->getGlobalNumElements(),
                                                   A_->getRowMap()->getNodeNumElements(), 0, A_->getRowMap()->getComm()));
    }
    else {
      throw std::runtime_error("Ifpack_Hypre: Unsupported map configuration: Row/Domain maps do not match");
    }
  }
  // Next create vectors that will be used when ApplyInverse() is called
  global_ordinal_type ilower = GloballyContiguousRowMap_->getMinGlobalIndex();
  global_ordinal_type iupper = GloballyContiguousRowMap_->getMaxGlobalIndex();
  // X in AX = Y
  IFPACK2_CHK_ERRV(HYPRE_IJVectorCreate(comm, ilower, iupper, &XHypre_));
  IFPACK2_CHK_ERRV(HYPRE_IJVectorSetObjectType(XHypre_, HYPRE_PARCSR));
  IFPACK2_CHK_ERRV(HYPRE_IJVectorInitialize(XHypre_));
  IFPACK2_CHK_ERRV(HYPRE_IJVectorAssemble(XHypre_));
  IFPACK2_CHK_ERRV(HYPRE_IJVectorGetObject(XHypre_, (void**) &ParX_));
  XVec_ = Teuchos::rcp((hypre_ParVector *) hypre_IJVectorObject(((hypre_IJVector *) XHypre_)),false);

  // Y in AX = Y
  IFPACK2_CHK_ERRV(HYPRE_IJVectorCreate(comm, ilower, iupper, &YHypre_));
  IFPACK2_CHK_ERRV(HYPRE_IJVectorSetObjectType(YHypre_, HYPRE_PARCSR));
  IFPACK2_CHK_ERRV(HYPRE_IJVectorInitialize(YHypre_));
  IFPACK2_CHK_ERRV(HYPRE_IJVectorAssemble(YHypre_));
  IFPACK2_CHK_ERRV(HYPRE_IJVectorGetObject(YHypre_, (void**) &ParY_));
  YVec_ = Teuchos::rcp((hypre_ParVector *) hypre_IJVectorObject(((hypre_IJVector *) YHypre_)),false);

  // Cache
  VectorCache_.resize(A->getRowMap()->getNodeNumElements());
} //Constructor

//==============================================================================
template<class MatrixType>
Hypre<MatrixType>::~Hypre() {
  Destroy();
}

//==============================================================================
template<class MatrixType>
void Hypre<MatrixType>::Destroy(){
  if(isInitialized()){
    IFPACK2_CHK_ERRV(HYPRE_IJMatrixDestroy(HypreA_));
  }
  IFPACK2_CHK_ERRV(HYPRE_IJVectorDestroy(XHypre_));
  IFPACK2_CHK_ERRV(HYPRE_IJVectorDestroy(YHypre_));
  if(IsSolverCreated_){
    IFPACK2_CHK_ERRV(SolverDestroyPtr_(Solver_));
  }
  if(IsPrecondCreated_){
    IFPACK2_CHK_ERRV(PrecondDestroyPtr_(Preconditioner_));
  }

  // Maxwell
  if(HypreG_) {
    IFPACK2_CHK_ERRV(HYPRE_IJMatrixDestroy(HypreG_));
  }
  if(xHypre_) {
    IFPACK2_CHK_ERRV(HYPRE_IJVectorDestroy(xHypre_));
  }
  if(yHypre_) {
    IFPACK2_CHK_ERRV(HYPRE_IJVectorDestroy(yHypre_));
  }
  if(zHypre_) {
    IFPACK2_CHK_ERRV(HYPRE_IJVectorDestroy(zHypre_));
  }
} //Destroy()

//==============================================================================
template<class MatrixType>
void Hypre<MatrixType>::initialize(){
  const std::string timerName ("Ifpack2::Hypre::initialize");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) timer = Teuchos::TimeMonitor::getNewCounter (timerName);

  if(IsInitialized_) return;  

  {
    Teuchos::TimeMonitor timeMon (*timer);
    
    // Create the Hypre matrix and copy values.  Note this uses values (which
    // Initialize() shouldn't do) but it doesn't care what they are (for
    // instance they can be uninitialized data even).  It should be possible to
    // set the Hypre structure without copying values, but this is the easiest
    // way to get the structure.
    MPI_Comm comm = * (Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(A_->getRowMap()->getComm())->getRawMpiComm());
    global_ordinal_type ilower = GloballyContiguousRowMap_->getMinGlobalIndex();
    global_ordinal_type iupper = GloballyContiguousRowMap_->getMaxGlobalIndex();
    IFPACK2_CHK_ERR(HYPRE_IJMatrixCreate(comm, ilower, iupper, ilower, iupper, &HypreA_));
    IFPACK2_CHK_ERR(HYPRE_IJMatrixSetObjectType(HypreA_, HYPRE_PARCSR));
    IFPACK2_CHK_ERR(HYPRE_IJMatrixInitialize(HypreA_));
    CopyTpetraToHypre();
    if(SolveOrPrec_ == Hypre_Is_Solver) {
      IFPACK2_CHK_ERR(SetSolverType(SolverType_));
      if (SolverPrecondPtr_ != NULL && UsePreconditioner_) {
        // both method allows a PC (first condition) and the user wants a PC (second)
        IFPACK2_CHK_ERR(SetPrecondType(PrecondType_));
        CallFunctions();
        IFPACK2_CHK_ERR(SolverPrecondPtr_(Solver_, PrecondSolvePtr_, PrecondSetupPtr_, Preconditioner_));
      } else {
        CallFunctions();
      }
    } else {
      IFPACK2_CHK_ERR(SetPrecondType(PrecondType_));
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
    NumInitialize_++;
  }
  InitializeTime_ = timer->totalElapsedTime ();
} //Initialize()

//==============================================================================
template<class MatrixType>
void Hypre<MatrixType>::setParameters(const Teuchos::ParameterList& list){

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
  chooserMap["Solver"] = Hypre_Is_Solver;
  chooserMap["Preconditioner"] = Hypre_Is_Preconditioner;

  List_ = list;
  Hypre_Solver solType;
  if (list.isType<std::string>("hypre: Solver"))
    solType = solverMap[list.get<std::string>("hypre: Solver")];
  else if(list.isParameter("hypre: Solver"))
    solType = (Hypre_Solver) list.get<int>("hypre: Solver");
  else 
    solType = PCG;
  SolverType_ = solType;
  Hypre_Solver precType;
  if (list.isType<std::string>("hypre: Preconditioner"))
    precType = solverMap[list.get<std::string>("hypre: Preconditioner")];
  else if(list.isParameter("hypre: Preconditioner"))
    precType = (Hypre_Solver) list.get<int>("hypre: Preconditioner");
  else 
    precType = Euclid;
  PrecondType_ = precType;
  Hypre_Chooser chooser;
  if (list.isType<std::string>("hypre: SolveOrPrecondition"))
    chooser = chooserMap[list.get<std::string>("hypre: SolveOrPrecondition")];
  else if(list.isParameter("hypre: SolveOrPrecondition"))
    chooser = (Hypre_Chooser) list.get<int>("hypre: SolveOrPrecondition");
  else 
    chooser = Hypre_Is_Solver;
  SolveOrPrec_ = chooser;
  bool SetPrecond = list.isParameter("hypre: SetPreconditioner") ? list.get<bool>("hypre: SetPreconditioner") : false;
  IFPACK2_CHK_ERR(SetParameter(SetPrecond));
  int NumFunctions = list.isParameter("hypre: NumFunctions") ? list.get<int>("hypre: NumFunctions") : 0;
  FunsToCall_.clear();
  NumFunsToCall_ = 0;
  if(NumFunctions > 0){
    RCP<FunctionParameter>* params = list.get<RCP<FunctionParameter>*>("hypre: Functions");
    for(int i = 0; i < NumFunctions; i++){
      IFPACK2_CHK_ERR(AddFunToList(params[i]));
    }
  }

  if (list.isSublist("hypre: Solver functions")) {
    Teuchos::ParameterList solverList = list.sublist("hypre: Solver functions");
    for (auto it = solverList.begin(); it != solverList.end(); ++it) {
      std::string funct_name = it->first;
      if (it->second.isType<HYPRE_Int>()) {
        IFPACK2_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Hypre_Is_Solver, funct_name , Teuchos::getValue<HYPRE_Int>(it->second)))));
      } else if (!std::is_same<HYPRE_Int,int>::value && it->second.isType<int>()) {
        IFPACK2_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Hypre_Is_Solver, funct_name , Teuchos::as<global_ordinal_type>(Teuchos::getValue<int>(it->second))))));
      } else if (it->second.isType<double>()) {
        IFPACK2_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Hypre_Is_Solver, funct_name , Teuchos::getValue<double>(it->second)))));
      } else {
        IFPACK2_CHK_ERR(-1);
      }
    }
  }

  if (list.isSublist("hypre: Preconditioner functions")) {
    Teuchos::ParameterList precList = list.sublist("hypre: Preconditioner functions");
    for (auto it = precList.begin(); it != precList.end(); ++it) {
      std::string funct_name = it->first;
      if (it->second.isType<HYPRE_Int>()) {
        IFPACK2_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Hypre_Is_Preconditioner, funct_name , Teuchos::getValue<HYPRE_Int>(it->second)))));
      } else if (!std::is_same<HYPRE_Int,int>::value && it->second.isType<int>()) {
        IFPACK2_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Hypre_Is_Preconditioner, funct_name ,  Teuchos::as<global_ordinal_type>(Teuchos::getValue<int>(it->second)))))); 
      } else if (it->second.isType<double>()) {
        IFPACK2_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Hypre_Is_Preconditioner, funct_name , Teuchos::getValue<double>(it->second)))));
      } else if (it->second.isList()) {
        Teuchos::ParameterList pl = Teuchos::getValue<Teuchos::ParameterList>(it->second);
        if (FunctionParameter::isFuncIntInt(funct_name)) {
          HYPRE_Int arg0 = pl.get<int>("arg 0");
          HYPRE_Int arg1 = pl.get<int>("arg 1");
          IFPACK2_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Hypre_Is_Preconditioner, funct_name , arg0, arg1))));
        } else if (FunctionParameter::isFuncIntIntDoubleDouble(funct_name)) {
          HYPRE_Int arg0 = pl.get<int>("arg 0");
          HYPRE_Int arg1 = pl.get<int>("arg 1");
          double arg2 = pl.get<double>("arg 2");
          double arg3 = pl.get<double>("arg 3");
          IFPACK2_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Hypre_Is_Preconditioner, funct_name , arg0, arg1, arg2, arg3))));
        } else if (FunctionParameter::isFuncIntIntIntDoubleIntInt(funct_name)) {
          HYPRE_Int arg0 = pl.get<int>("arg 0");
          HYPRE_Int arg1 = pl.get<int>("arg 1");
          HYPRE_Int arg2 = pl.get<int>("arg 2");
          double arg3 = pl.get<double>("arg 3");
          HYPRE_Int arg4 = pl.get<int>("arg 4");
          HYPRE_Int arg5 = pl.get<int>("arg 5");
          IFPACK2_CHK_ERR(AddFunToList(rcp(new FunctionParameter(Hypre_Is_Preconditioner, funct_name , arg0, arg1, arg2, arg3, arg4, arg5))));
        } else {
          IFPACK2_CHK_ERR(-1);
        }
      }
    }
  }

  if (list.isSublist("Coordinates") && list.sublist("Coordinates").isType<Teuchos::RCP<multivector_type> >("Coordinates"))
    Coords_ = list.sublist("Coordinates").get<Teuchos::RCP<multivector_type> >("Coordinates");
  if (list.isSublist("Operators") && list.sublist("Operators").isType<Teuchos::RCP<const crs_matrix_type> >("G"))
    G_ = list.sublist("Operators").get<Teuchos::RCP<const crs_matrix_type> >("G");

  Dump_ = list.isParameter("hypre: Dump") ? list.get<bool>("hypre: Dump") : false;
} //setParameters()

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::AddFunToList(RCP<FunctionParameter> NewFun){
  NumFunsToCall_ = NumFunsToCall_+1;
  FunsToCall_.resize(NumFunsToCall_);
  FunsToCall_[NumFunsToCall_-1] = NewFun;
  return 0;
} //AddFunToList()

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::SetParameter(Hypre_Chooser chooser, global_ordinal_type (*pt2Func)(HYPRE_Solver, global_ordinal_type), global_ordinal_type parameter){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter));
  IFPACK2_CHK_ERR(AddFunToList(temp));
  return 0;
} //SetParameter() - int function pointer

//=============================================================================
template<class MatrixType>
int Hypre<MatrixType>::SetParameter(Hypre_Chooser chooser, global_ordinal_type  (*pt2Func)(HYPRE_Solver, double), double parameter){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter));
  IFPACK2_CHK_ERR(AddFunToList(temp));
  return 0;
} //SetParameter() - double function pointer

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::SetParameter(Hypre_Chooser chooser, global_ordinal_type (*pt2Func)(HYPRE_Solver, double, global_ordinal_type), double parameter1, global_ordinal_type parameter2){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter1, parameter2));
  IFPACK2_CHK_ERR(AddFunToList(temp));
  return 0;
} //SetParameter() - double,int function pointer

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::SetParameter(Hypre_Chooser chooser, global_ordinal_type (*pt2Func)(HYPRE_Solver, global_ordinal_type, global_ordinal_type), global_ordinal_type parameter1, global_ordinal_type parameter2){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter1, parameter2));
  IFPACK2_CHK_ERR(AddFunToList(temp));
  return 0;
} //SetParameter() int,int function pointer

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::SetParameter(Hypre_Chooser chooser, global_ordinal_type (*pt2Func)(HYPRE_Solver, double*), double* parameter){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter));
  IFPACK2_CHK_ERR(AddFunToList(temp));
  return 0;
} //SetParameter() - double* function pointer

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::SetParameter(Hypre_Chooser chooser, global_ordinal_type (*pt2Func)(HYPRE_Solver, global_ordinal_type*), global_ordinal_type* parameter){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter));
  IFPACK2_CHK_ERR(AddFunToList(temp));
  return 0;
} //SetParameter() - int* function pointer

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::SetParameter(Hypre_Chooser chooser,  global_ordinal_type (*pt2Func)(HYPRE_Solver, global_ordinal_type**), global_ordinal_type** parameter){
  RCP<FunctionParameter> temp = rcp(new FunctionParameter(chooser, pt2Func, parameter));
  IFPACK2_CHK_ERR(AddFunToList(temp));
  return 0;
} //SetParameter() - int** function pointer

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::SetParameter(Hypre_Chooser chooser, Hypre_Solver solver){
  if(chooser == Hypre_Is_Solver){
    SolverType_ = solver;
  } else {
    PrecondType_ = solver;
  }
  return 0;
} //SetParameter() - set type of solver

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::SetDiscreteGradient(Teuchos::RCP<const crs_matrix_type> G){
  using LO = local_ordinal_type;
  using GO = global_ordinal_type;
  using SC = scalar_type;

  // Sanity check
  if(!A_->getRowMap()->isSameAs(*G->getRowMap()))
    throw std::runtime_error("Hypre<MatrixType>: Edge map mismatch: A and discrete gradient");

  // Get the maps for the nodes (assuming the edge map from A is OK);
  GloballyContiguousNodeRowMap_ = rcp(new map_type(G->getDomainMap()->getGlobalNumElements(),
                                                   G->getDomainMap()->getNodeNumElements(), 0, A_->getRowMap()->getComm()));
  GloballyContiguousNodeColMap_ = MakeContiguousColumnMap(G);

  // Start building G
  MPI_Comm comm = * (Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(A_->getRowMap()->getComm())->getRawMpiComm());
  GO ilower = GloballyContiguousRowMap_->getMinGlobalIndex();
  GO iupper = GloballyContiguousRowMap_->getMaxGlobalIndex();
  GO jlower = GloballyContiguousNodeRowMap_->getMinGlobalIndex();
  GO jupper = GloballyContiguousNodeRowMap_->getMaxGlobalIndex();
  IFPACK2_CHK_ERR(HYPRE_IJMatrixCreate(comm, ilower, iupper, jlower, jupper, &HypreG_));
  IFPACK2_CHK_ERR(HYPRE_IJMatrixSetObjectType(HypreG_, HYPRE_PARCSR));
  IFPACK2_CHK_ERR(HYPRE_IJMatrixInitialize(HypreG_));

  std::vector<GO> new_indices(G->getNodeMaxNumRowEntries());
  for(LO i = 0; i < (LO)G->getNodeNumRows(); i++){
    Teuchos::ArrayView<const SC> values;
    Teuchos::ArrayView<const LO> indices;
    G->getLocalRowView(i, indices, values);
    for(LO j = 0; j < (LO) indices.size(); j++){
      new_indices[j] = GloballyContiguousNodeColMap_->getGlobalElement(indices[j]);
    }
    GO GlobalRow[1];
    GO numEntries = (GO) indices.size();
    GlobalRow[0] = GloballyContiguousRowMap_->getGlobalElement(i);
    IFPACK2_CHK_ERR(HYPRE_IJMatrixSetValues(HypreG_, 1, &numEntries, GlobalRow, new_indices.data(), values.getRawPtr()));
  }
  IFPACK2_CHK_ERR(HYPRE_IJMatrixAssemble(HypreG_));
  IFPACK2_CHK_ERR(HYPRE_IJMatrixGetObject(HypreG_, (void**)&ParMatrixG_));

  if (Dump_)
    HYPRE_ParCSRMatrixPrint(ParMatrixG_,"G.mat");

  if(SolverType_ == AMS)
    HYPRE_AMSSetDiscreteGradient(Solver_, ParMatrixG_);
  if(PrecondType_ == AMS)
    HYPRE_AMSSetDiscreteGradient(Preconditioner_, ParMatrixG_);
  return 0;
} //SetDiscreteGradient()

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::SetCoordinates(Teuchos::RCP<multivector_type> coords) {

  if(!G_.is_null() && !G_->getDomainMap()->isSameAs(*coords->getMap()))
    throw std::runtime_error("Hypre<MatrixType>: Node map mismatch: G->DomainMap() and coords");

  if(SolverType_ != AMS && PrecondType_ != AMS)
    return 0;

  scalar_type *xPtr = coords->getDataNonConst(0).getRawPtr();
  scalar_type *yPtr = coords->getDataNonConst(1).getRawPtr();
  scalar_type *zPtr = coords->getDataNonConst(2).getRawPtr();

  MPI_Comm comm = * (Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(A_->getRowMap()->getComm())->getRawMpiComm());
  local_ordinal_type NumEntries = coords->getLocalLength();
  global_ordinal_type * indices = const_cast<global_ordinal_type*>(GloballyContiguousNodeRowMap_->getNodeElementList().getRawPtr());

  global_ordinal_type ilower = GloballyContiguousNodeRowMap_->getMinGlobalIndex();
  global_ordinal_type iupper = GloballyContiguousNodeRowMap_->getMaxGlobalIndex();

  if( NumEntries != iupper-ilower+1) {
    std::cout<<"Ifpack2::Hypre::SetCoordinates(): Error on rank "<<A_->getRowMap()->getComm()->getRank()<<": MyLength = "<<coords->getLocalLength()<<" GID range = ["<<ilower<<","<<iupper<<"]"<<std::endl;
    throw std::runtime_error("Hypre<MatrixType>: SetCoordinates: Length mismatch");
  }

  IFPACK2_CHK_ERR(HYPRE_IJVectorCreate(comm, ilower, iupper, &xHypre_));
  IFPACK2_CHK_ERR(HYPRE_IJVectorSetObjectType(xHypre_, HYPRE_PARCSR));
  IFPACK2_CHK_ERR(HYPRE_IJVectorInitialize(xHypre_));
  IFPACK2_CHK_ERR(HYPRE_IJVectorSetValues(xHypre_,NumEntries,indices,xPtr));
  IFPACK2_CHK_ERR(HYPRE_IJVectorAssemble(xHypre_));
  IFPACK2_CHK_ERR(HYPRE_IJVectorGetObject(xHypre_, (void**) &xPar_));

  IFPACK2_CHK_ERR(HYPRE_IJVectorCreate(comm, ilower, iupper, &yHypre_));
  IFPACK2_CHK_ERR(HYPRE_IJVectorSetObjectType(yHypre_, HYPRE_PARCSR));
  IFPACK2_CHK_ERR(HYPRE_IJVectorInitialize(yHypre_));
  IFPACK2_CHK_ERR(HYPRE_IJVectorSetValues(yHypre_,NumEntries,indices,yPtr));
  IFPACK2_CHK_ERR(HYPRE_IJVectorAssemble(yHypre_));
  IFPACK2_CHK_ERR(HYPRE_IJVectorGetObject(yHypre_, (void**) &yPar_));

  IFPACK2_CHK_ERR(HYPRE_IJVectorCreate(comm, ilower, iupper, &zHypre_));
  IFPACK2_CHK_ERR(HYPRE_IJVectorSetObjectType(zHypre_, HYPRE_PARCSR));
  IFPACK2_CHK_ERR(HYPRE_IJVectorInitialize(zHypre_));
  IFPACK2_CHK_ERR(HYPRE_IJVectorSetValues(zHypre_,NumEntries,indices,zPtr));
  IFPACK2_CHK_ERR(HYPRE_IJVectorAssemble(zHypre_));
  IFPACK2_CHK_ERR(HYPRE_IJVectorGetObject(zHypre_, (void**) &zPar_));

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
template<class MatrixType>
void Hypre<MatrixType>::compute(){
  const std::string timerName ("Ifpack2::Hypre::compute");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) timer = Teuchos::TimeMonitor::getNewCounter (timerName);

  // Start timing here.
  {
    Teuchos::TimeMonitor timeMon (*timer);

    if(isInitialized() == false){
      initialize();
    }
           
    // Hypre Setup must be called after matrix has values
    if(SolveOrPrec_ == Hypre_Is_Solver){
      IFPACK2_CHK_ERR(SolverSetupPtr_(Solver_, ParMatrix_, ParX_, ParY_));
    } else {
      IFPACK2_CHK_ERR(PrecondSetupPtr_(Preconditioner_, ParMatrix_, ParX_, ParY_));
    }
    
    IsComputed_ = true;
    NumCompute_++;
  }
  
  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  ComputeTime_ = timer->totalElapsedTime ();
} //Compute()

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::CallFunctions() const{
  for(int i = 0; i < NumFunsToCall_; i++){
    IFPACK2_CHK_ERR(FunsToCall_[i]->CallFunction(Solver_, Preconditioner_));
  }
  return 0;
} //CallFunctions()

//==============================================================================
template<class MatrixType>
void Hypre<MatrixType>::apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                               Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
                               Teuchos::ETransp mode,
                               scalar_type alpha,
                               scalar_type beta) const {
  using LO = local_ordinal_type;
  using SC = scalar_type;
  const std::string timerName ("Ifpack2::Hypre::apply");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) timer = Teuchos::TimeMonitor::getNewCounter (timerName);

  // Start timing here.
  {
    Teuchos::TimeMonitor timeMon (*timer);
    
    if(isComputed() == false){
      IFPACK2_CHK_ERR(-1);
    }
    hypre_Vector *XLocal_ = hypre_ParVectorLocalVector(XVec_);
    hypre_Vector *YLocal_ = hypre_ParVectorLocalVector(YVec_);
    
    bool SameVectors = false;
    size_t NumVectors = X.getNumVectors();
    if (NumVectors != Y.getNumVectors()) IFPACK2_CHK_ERR(-1);  // X and Y must have same number of vectors
    if(&X == &Y) { //FIXME: Maybe not the right way to check this
      SameVectors = true;
    }
    
    // NOTE: Here were assuming that the local ordering of Epetra's X/Y-vectors and 
    // Hypre's X/Y-vectors are the same.  Seeing as as this is more or less how we constructed
    // the Hypre matrices, this seems pretty reasoanble.
    
    for(int VecNum = 0; VecNum < (int) NumVectors; VecNum++) {
      //Get values for current vector in multivector.
      // FIXME amk Nov 23, 2015: This will not work for funky data layouts
      SC * XValues = const_cast<SC*>(X.getData(VecNum).getRawPtr());
      SC * YValues;
      if(!SameVectors){
        YValues = const_cast<SC*>(Y.getData(VecNum).getRawPtr());
      } else {
        YValues = VectorCache_.getRawPtr();
      }
      // Temporarily make a pointer to data in Hypre for end
      SC *XTemp = XLocal_->data;
      // Replace data in Hypre vectors with Epetra data
      XLocal_->data = XValues;
      SC *YTemp = YLocal_->data;
      YLocal_->data = YValues;
      
      IFPACK2_CHK_ERR(HYPRE_ParVectorSetConstantValues(ParY_, 0.0));
      if(SolveOrPrec_ == Hypre_Is_Solver){
        // Use the solver methods
        IFPACK2_CHK_ERR(SolverSolvePtr_(Solver_, ParMatrix_, ParX_, ParY_));
      } else {
        // Apply the preconditioner
        IFPACK2_CHK_ERR(PrecondSolvePtr_(Preconditioner_, ParMatrix_, ParX_, ParY_));
      }

      if(SameVectors){
        Teuchos::ArrayView<SC> Yv =  Y.getDataNonConst(VecNum)();
        LO NumEntries = Y.getLocalLength();
        for(LO i = 0; i < NumEntries; i++)
          Yv[i] = YValues[i];      
      }
      XLocal_->data = XTemp;
      YLocal_->data = YTemp;
    }
    NumApply_++;
  }
  // timer->totalElapsedTime() returns the total time over all timer
  // calls.  Thus, we use = instead of +=.
  ApplyTime_ = timer->totalElapsedTime ();
} //apply()


//==============================================================================
template<class MatrixType>
void Hypre<MatrixType>::applyMat (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
                                  Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
                                  Teuchos::ETransp mode) const {
  A_->apply(X,Y,mode);
} //applyMat()

//==============================================================================
template<class MatrixType>
std::string Hypre<MatrixType>::description() const {
  std::ostringstream out;

  // Output is a valid YAML dictionary in flow style.  If you don't
  // like everything on a single line, you should call describe()
  // instead.
  out << "\"Ifpack2::Hypre\": {";
  out << "Initialized: " << (isInitialized () ? "true" : "false") << ", "
      << "Computed: " << (isComputed () ? "true" : "false") << ", ";

  if (A_.is_null ()) {
    out << "Matrix: null";
  }
  else {
    out << "Global matrix dimensions: ["
        << A_->getGlobalNumRows () << ", "
        << A_->getGlobalNumCols () << "]"
        << ", Global nnz: " << A_->getGlobalNumEntries();
  }

  out << "}";
  return out.str ();
} //description()

//==============================================================================
template<class MatrixType>
void Hypre<MatrixType>::describe(Teuchos::FancyOStream &os, const Teuchos::EVerbosityLevel verbLevel) const {
  using std::endl;
  os << endl;
  os << "================================================================================" << endl;
  os << "Ifpack2::Hypre: " << endl << endl;
  os << "Using " << A_->getComm()->getSize() << " processors." << endl;
  os << "Global number of rows            = " << A_->getGlobalNumRows() << endl;
  os << "Global number of nonzeros        = " << A_->getGlobalNumEntries() << endl;
  //    os << "Condition number estimate = " << Condest() << endl;
  os << endl;
  os << "Phase           # calls   Total Time (s)"<<endl;
  os << "-----           -------   --------------"<<endl;
  os << "Initialize()    "   << std::setw(5) << NumInitialize_
     << "  " << std::setw(15) << InitializeTime_<<endl;
  os << "Compute()       "   << std::setw(5) << NumCompute_
     << "  " << std::setw(15) << ComputeTime_ << endl;
  os << "ApplyInverse()  "   << std::setw(5) << NumApply_
     << "  " << std::setw(15) << ApplyTime_ <<endl;
  os << "================================================================================" << endl;
  os << endl;
} //description

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::SetSolverType(Hypre_Solver Solver){
  switch(Solver) {
    case BoomerAMG:
      if(IsSolverCreated_){
        SolverDestroyPtr_(Solver_);
        IsSolverCreated_ = false;
      }
      SolverCreatePtr_ = &Hypre<MatrixType>::Hypre_BoomerAMGCreate;
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
      SolverCreatePtr_ = &Hypre<MatrixType>::Hypre_AMSCreate;
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
      SolverCreatePtr_ = &Hypre<MatrixType>::Hypre_ParCSRHybridCreate;
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
      SolverCreatePtr_ = &Hypre<MatrixType>::Hypre_ParCSRPCGCreate;
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
      SolverCreatePtr_ = &Hypre<MatrixType>::Hypre_ParCSRGMRESCreate;
      SolverDestroyPtr_ = &HYPRE_ParCSRGMRESDestroy;
      SolverSetupPtr_ = &HYPRE_ParCSRGMRESSetup;
      SolverPrecondPtr_ = &HYPRE_ParCSRGMRESSetPrecond;
      break;
    case FlexGMRES:
      if(IsSolverCreated_){
        SolverDestroyPtr_(Solver_);
        IsSolverCreated_ = false;
      }
      SolverCreatePtr_ = &Hypre<MatrixType>::Hypre_ParCSRFlexGMRESCreate;
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
      SolverCreatePtr_ = &Hypre<MatrixType>::Hypre_ParCSRLGMRESCreate;
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
      SolverCreatePtr_ = &Hypre<MatrixType>::Hypre_ParCSRBiCGSTABCreate;
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
template<class MatrixType>
int Hypre<MatrixType>::SetPrecondType(Hypre_Solver Precond){
  switch(Precond) {
    case BoomerAMG:
      if(IsPrecondCreated_){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondCreated_ = false;
      }
      PrecondCreatePtr_ = &Hypre<MatrixType>::Hypre_BoomerAMGCreate;
      PrecondDestroyPtr_ = &HYPRE_BoomerAMGDestroy;
      PrecondSetupPtr_ = &HYPRE_BoomerAMGSetup;
      PrecondSolvePtr_ = &HYPRE_BoomerAMGSolve;
      break;
    case ParaSails:
      if(IsPrecondCreated_){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondCreated_ = false;
      }
      PrecondCreatePtr_ = &Hypre<MatrixType>::Hypre_ParaSailsCreate;
      PrecondDestroyPtr_ = &HYPRE_ParaSailsDestroy;
      PrecondSetupPtr_ = &HYPRE_ParaSailsSetup;
      PrecondSolvePtr_ = &HYPRE_ParaSailsSolve;
      break;
    case Euclid:
      if(IsPrecondCreated_){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondCreated_ = false;
      }
      PrecondCreatePtr_ = &Hypre<MatrixType>::Hypre_EuclidCreate;
      PrecondDestroyPtr_ = &HYPRE_EuclidDestroy;
      PrecondSetupPtr_ = &HYPRE_EuclidSetup;
      PrecondSolvePtr_ = &HYPRE_EuclidSolve;
      break;
    case AMS:
      if(IsPrecondCreated_){
        PrecondDestroyPtr_(Preconditioner_);
        IsPrecondCreated_ = false;
      }
      PrecondCreatePtr_ = &Hypre<MatrixType>::Hypre_AMSCreate;
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
template<class MatrixType>
int Hypre<MatrixType>::CreateSolver(){
  MPI_Comm comm;
  HYPRE_ParCSRMatrixGetComm(ParMatrix_, &comm);
  int ierr = (this->*SolverCreatePtr_)(comm, &Solver_);
  IsSolverCreated_ = true;
  return ierr;
} //CreateSolver()

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::CreatePrecond(){
  MPI_Comm comm;
  HYPRE_ParCSRMatrixGetComm(ParMatrix_, &comm);
  int ierr = (this->*PrecondCreatePtr_)(comm, &Preconditioner_);
  IsPrecondCreated_ = true;
  return ierr;
} //CreatePrecond()


//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::CopyTpetraToHypre(){
  using LO = local_ordinal_type;
  using GO = global_ordinal_type;
  using SC = scalar_type;

  Teuchos::RCP<const crs_matrix_type> Matrix = Teuchos::rcp_dynamic_cast<const crs_matrix_type>(A_);
  if(Matrix.is_null()) 
    throw std::runtime_error("Hypre<MatrixType>: Unsupported matrix configuration: Tpetra::CrsMatrix required");

  std::vector<GO> new_indices(Matrix->getNodeMaxNumRowEntries());
  for(LO i = 0; i < (LO) Matrix->getNodeNumRows(); i++){
    Teuchos::ArrayView<const SC> values;
    Teuchos::ArrayView<const LO> indices;
    Matrix->getLocalRowView(i, indices, values);
    for(LO j = 0; j < (LO)indices.size(); j++){
      new_indices[j] = GloballyContiguousColMap_->getGlobalElement(indices[j]);
    }
    GO GlobalRow[1];
    GO numEntries = (GO) indices.size();
    GlobalRow[0] = GloballyContiguousRowMap_->getGlobalElement(i);    
    IFPACK2_CHK_ERR(HYPRE_IJMatrixSetValues(HypreA_, 1, &numEntries, GlobalRow, new_indices.data(), values.getRawPtr()));
  }
  IFPACK2_CHK_ERR(HYPRE_IJMatrixAssemble(HypreA_));
  IFPACK2_CHK_ERR(HYPRE_IJMatrixGetObject(HypreA_, (void**)&ParMatrix_));
  if (Dump_)
    HYPRE_ParCSRMatrixPrint(ParMatrix_,"A.mat");
  return 0;
} //CopyTpetraToHypre()

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::Hypre_BoomerAMGCreate(MPI_Comm /*comm*/, HYPRE_Solver *solver)
    { return HYPRE_BoomerAMGCreate(solver);}

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::Hypre_ParaSailsCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParaSailsCreate(comm, solver);}

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::Hypre_EuclidCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_EuclidCreate(comm, solver);}

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::Hypre_AMSCreate(MPI_Comm /*comm*/, HYPRE_Solver *solver)
    { return HYPRE_AMSCreate(solver);}

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::Hypre_ParCSRHybridCreate(MPI_Comm /*comm*/, HYPRE_Solver *solver)
    { return HYPRE_ParCSRHybridCreate(solver);}

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::Hypre_ParCSRPCGCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRPCGCreate(comm, solver);}

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::Hypre_ParCSRGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRGMRESCreate(comm, solver);}

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::Hypre_ParCSRFlexGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRFlexGMRESCreate(comm, solver);}

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::Hypre_ParCSRLGMRESCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRLGMRESCreate(comm, solver);}

//==============================================================================
template<class MatrixType>
int Hypre<MatrixType>::Hypre_ParCSRBiCGSTABCreate(MPI_Comm comm, HYPRE_Solver *solver)
    { return HYPRE_ParCSRBiCGSTABCreate(comm, solver);}

//==============================================================================
template<class MatrixType>
  Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> > 
Hypre<MatrixType>::MakeContiguousColumnMap(Teuchos::RCP<const crs_matrix_type> &Matrix) const{
  using import_type     = Tpetra::Import<local_ordinal_type,global_ordinal_type,node_type>;
  using go_vector_type  = Tpetra::Vector<global_ordinal_type,local_ordinal_type,global_ordinal_type,node_type>;

  // Must create GloballyContiguous DomainMap (which is a permutation of Matrix_'s
  // DomainMap) and the corresponding permuted ColumnMap.
  //   Epetra_GID  --------->   LID   ----------> HYPRE_GID
  //           via DomainMap.LID()       via GloballyContiguousDomainMap.GID()
  if(Matrix.is_null()) 
    throw std::runtime_error("Hypre<MatrixType>: Unsupported matrix configuration: Tpetra::CrsMatrix required");
  RCP<const map_type> DomainMap = Matrix->getDomainMap();
  RCP<const map_type> ColumnMap = Matrix->getColMap();
  RCP<const import_type> importer = Matrix->getGraph()->getImporter();

  if(DomainMap->isContiguous() ) {
    // If the domain map is linear, then we can just use the column map as is.
    return ColumnMap;
  }
  else {
    // The domain map isn't linear, so we need a new domain map
    Teuchos::RCP<map_type> ContiguousDomainMap = rcp(new map_type(DomainMap->getGlobalNumElements(),
                                                                      DomainMap->getNodeNumElements(), 0, DomainMap->getComm()));
    if(importer) {    
      // If there's an importer then we can use it to get a new column map
      go_vector_type MyGIDsHYPRE(DomainMap,ContiguousDomainMap->getNodeElementList());

      // import the HYPRE GIDs
      go_vector_type ColGIDsHYPRE(ColumnMap);
      ColGIDsHYPRE.doImport(MyGIDsHYPRE, *importer, Tpetra::INSERT);
   
      // Make a HYPRE numbering-based column map.
      return Teuchos::rcp(new map_type(ColumnMap->getGlobalNumElements(),ColGIDsHYPRE.getDataNonConst()(),0, ColumnMap->getComm()));
    }
    else {
      // The problem has matching domain/column maps, and somehow the domain map isn't linear, so just use the new domain map
      return Teuchos::rcp(new map_type(ColumnMap->getGlobalNumElements(),ContiguousDomainMap->getNodeElementList(), 0, ColumnMap->getComm()));
    }
  }  
}



template<class MatrixType>
int Hypre<MatrixType>::getNumInitialize() const {
  return NumInitialize_;
}


template<class MatrixType>
int Hypre<MatrixType>::getNumCompute() const {
  return NumCompute_;
}


template<class MatrixType>
int Hypre<MatrixType>::getNumApply() const {
  return NumApply_;
}


template<class MatrixType>
double Hypre<MatrixType>::getInitializeTime() const {
  return InitializeTime_;
}


template<class MatrixType>
double Hypre<MatrixType>::getComputeTime() const {
  return ComputeTime_;
}


template<class MatrixType>
double Hypre<MatrixType>::getApplyTime() const {
  return ApplyTime_;
}

template<class MatrixType>
Teuchos::RCP<const typename Hypre<MatrixType>::map_type>
Hypre<MatrixType>::
getDomainMap () const
{
  Teuchos::RCP<const row_matrix_type> A = getMatrix();
  TEUCHOS_TEST_FOR_EXCEPTION(
    A.is_null (), std::runtime_error, "Ifpack2::Hypre::getDomainMap: The "
    "input matrix A is null.  Please call setMatrix() with a nonnull input "
    "matrix before calling this method.");
  return A->getDomainMap ();
}


template<class MatrixType>
Teuchos::RCP<const typename Hypre<MatrixType>::map_type>
Hypre<MatrixType>::
getRangeMap () const
{
  Teuchos::RCP<const row_matrix_type> A = getMatrix();
  TEUCHOS_TEST_FOR_EXCEPTION(
    A.is_null (), std::runtime_error, "Ifpack2::Hypre::getRangeMap: The "
    "input matrix A is null.  Please call setMatrix() with a nonnull input "
    "matrix before calling this method.");
  return A->getRangeMap ();
}


template<class MatrixType>
void Hypre<MatrixType>::setMatrix (const Teuchos::RCP<const row_matrix_type>& A)
{
  if (A.getRawPtr () != getMatrix().getRawPtr ()) {
    IsInitialized_ = false;
    IsComputed_ = false;
    A_ = A;
  }
}


template<class MatrixType>
Teuchos::RCP<const typename Hypre<MatrixType>::row_matrix_type>
Hypre<MatrixType>::
getMatrix() const {
  return A_;
}


template<class MatrixType>
bool Hypre<MatrixType>::hasTransposeApply() const {
  return false;
}

}// Ifpack2 namespace


#define IFPACK2_HYPRE_INSTANT(S,LO,GO,N)                            \
  template class Ifpack2::Hypre< Tpetra::RowMatrix<S, LO, GO, N> >;


#endif // HAVE_HYPRE && HAVE_MPI
#endif // IFPACK2_HYPRE_DEF_HPP
