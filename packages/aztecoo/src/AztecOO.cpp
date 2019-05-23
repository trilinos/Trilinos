/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
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

#include <climits>
#include "AztecOO.h"
#include "Epetra_ConfigDefs.h"
#ifdef AZTEC_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_Operator.h"
#include "Epetra_Import.h"

#include <AztecOO_string_maps.h>
#include <AztecOO_Scaling.h>
#include <AZOO_printf.h>

#ifdef HAVE_AZTECOO_TEUCHOS
#include <Teuchos_ParameterList.hpp>
#endif

#ifdef HAVE_TEUCHOS_EXTENDED
#include <Teuchos_StrUtils.hpp>
#else
//If Trilinos wasn't configured with Teuchos-extended enabled, then we need
//to write our own function for converting a std::string to upper-case, and to do
//that we need the toupper() prototype that's declared in ctype.h.
#include <ctype.h>
#endif

AztecOO::MatrixData::~MatrixData()
{
  if (X!=0) delete X;
  if (Y!=0) delete Y;
  if (SourceVec!=0) delete SourceVec;
  if (TargetVec!=0) delete TargetVec;
}

AztecOO::OperatorData::~OperatorData()
{
  if (X!=0) delete X;
  if (Y!=0) delete Y;
}

//=============================================================================
AztecOO::AztecOO(Epetra_Operator    * A, 
		 Epetra_MultiVector * X,
		 Epetra_MultiVector * B) {
  out_stream_ = 0;
  err_stream_ = 0;

  inConstructor_ = true;  // Shut down complaints about zero pointers for a while
  AllocAzArrays();
  SetAztecDefaults();


  Epetra_RowMatrix * UserMatrix = dynamic_cast<Epetra_RowMatrix *>(A); // Try to cast operator to matrix 
  if (UserMatrix!=0) 
    SetUserMatrix(UserMatrix);
  else
    SetUserOperator(A);

  SetLHS(X);
  SetRHS(B);
  inConstructor_ = false;
}

//=============================================================================
AztecOO::AztecOO(Epetra_RowMatrix * A, 
		 Epetra_MultiVector * X,
		 Epetra_MultiVector * B) {
  out_stream_ = 0;
  err_stream_ = 0;

  inConstructor_ = true;  // Shut down complaints about zero pointers for a while
  AllocAzArrays();
  SetAztecDefaults();

  SetUserMatrix(A);
  
  SetLHS(X);
  SetRHS(B);
  inConstructor_ = false;
}

//=============================================================================
AztecOO::AztecOO(const Epetra_LinearProblem& prob) {
  out_stream_ = 0;
  err_stream_ = 0;

  inConstructor_ = true;  // Shut down complaints about zero pointers for a while
  AllocAzArrays();
  SetAztecDefaults();
  SetProblem(prob);

  inConstructor_ = false;
}

//=============================================================================
AztecOO::AztecOO() {
  out_stream_ = 0;
  err_stream_ = 0;

  inConstructor_ = true;  // Shut down complaints about zero pointers for a while
  AllocAzArrays();
  SetAztecDefaults();
  inConstructor_ = false;
}

//=============================================================================
AztecOO::AztecOO(const AztecOO& source) {
  out_stream_ = source.out_stream_;
  err_stream_ = source.err_stream_;

  inConstructor_ = true;  // Shut down complaints about zero pointers for a while
  AllocAzArrays();
  SetAztecDefaults();

  Epetra_LinearProblem *problem = source.GetProblem();
  if (problem) SetProblem(*problem);
  SetUserMatrix(source.GetUserMatrix());
  SetUserOperator(source.GetUserOperator());
  SetPrecMatrix(source.GetPrecMatrix());   // Assume user want to base preconditioner on this matrix
  SetPrecOperator(source.GetPrecOperator());   // Assume user want to base preconditioner on this matrix
  SetLHS(source.GetLHS());
  SetRHS(source.GetRHS());
  SetAllAztecOptions(source.options_);
  SetAllAztecParams(source.params_);
  inConstructor_ = false;
}

//=============================================================================
AztecOO::~AztecOO(void) {
  DeleteMemory();
  DeleteAzArrays();
  AZOO_clear_streams();
}

//=============================================================================
int AztecOO::AllocAzArrays()
{
  proc_config_ = new int[AZ_PROC_SIZE];

  options_ = new int[AZ_OPTIONS_SIZE];
  params_ = new double[AZ_PARAMS_SIZE];
  status_ = new double[AZ_STATUS_SIZE];  

  if (proc_config_ ==0 ||options_ == 0 || params_ == 0 || status_ == 0) EPETRA_CHK_ERR(-1);
  for ( int i = 0; i < AZ_STATUS_SIZE; ++i ) status_[i] = 0.0;

  return(0);
}
//=============================================================================
void AztecOO::DeleteMemory() {
  if (options_[AZ_keep_info] == 1 && Scaling_ != 0 && Scaling_->scale != 0) {
    Scaling_->scale(AZ_DESTROY_SCALING_DATA, Amat_, options_, x_, b_,
		    proc_config_, Scaling_);
  }  
  
  if (Prec_!=0) {AZ_precond_destroy(&Prec_); Prec_ = 0;}
  if (Pmat_ != 0) {
    AZ_rm_context(options_,params_,Pmat_->data_org);
    AZ_matrix_destroy(&Pmat_); Pmat_ = 0;
  }
  if (Amat_ != 0) {AZ_matrix_destroy(&Amat_); Amat_ = 0;}

  if (UserOperatorData_!=0) {delete UserOperatorData_; UserOperatorData_ = 0;}
  if (UserMatrixData_!=0) {delete UserMatrixData_; UserMatrixData_ = 0;}
  if (PrecOperatorData_!=0) {delete PrecOperatorData_; PrecOperatorData_ = 0;}
  if (PrecMatrixData_!=0) {delete PrecMatrixData_; PrecMatrixData_ = 0;}
  if (ResidualVector_!=0) {delete ResidualVector_; ResidualVector_ = 0;}
  if (conv_info_!=0) {AZ_converge_destroy(&conv_info_); conv_info_ = 0;}
  if (Scaling_created_ == true) {
    AZ_scaling_destroy(&Scaling_);
    Scaling_ = 0;
  }
  if (Label_!=0) {
    delete [] Label_;
    Label_ = 0;
  }
}

//=============================================================================
std::string AztecOO_uppercase(const std::string& s)
{
  //convert incoming std::string to uppercase, and prepend 'AZ_' if the first
  //two characters aren't already 'AZ'.

#ifdef HAVE_TEUCHOS_EXTENDED
  std::string upp = Teuchos::StrUtils::allCaps(s);
#else
  std::string upp(s);
  for(unsigned i=0; i<upp.length(); ++i) {
    upp[i] = toupper(upp[i]);
  }
#endif

  if (upp[0] == 'A' && upp[1] == 'Z') {
    return(upp);
  }

  std::string az_("AZ_");
  return(az_+upp);
}

#ifdef HAVE_AZTECOO_TEUCHOS
//=============================================================================
bool AztecOO_SetOptionOrParam(int offset,
                              const Teuchos::ParameterEntry& entry,
                              AztecOO* azoo)
{
  //Use the given parameter-list entry to set azoo::options_[offset] or
  //azoo::params_[offset], as appropriate.
  //
  //Return true if entry is used, false if not.
  //
  //The options_ and params_ arrays are allocated elsewhere, and have size
  //AZ_OPTIONS_SIZE and AZ_PARAMS_SIZE respectively. However, the positions
  //in those arrays that are eligible to be set by this function are
  //azoo::options_[0 ... AZ_FIRST_USER_OPTION] and
  //azoo::params_[0 ... AZ_FIRST_USER_PARAM]
  //All AZ_ #defines are listed in az_aztec_defs.h.

  bool entry_used = false;

  if (offset < 0) return(entry_used);

  int dummy_int;
  double dummy_double;
  std::string dummy_string;

  if (entry.isType<int>() || entry.isType<unsigned>()) {
    if (offset < AZ_FIRST_USER_OPTION) {
      int ival = entry.getValue(&dummy_int);
      azoo->SetAztecOption(offset, ival);
      entry_used = true;
    }
  }
  else if (entry.isType<std::string>()) {
    if (offset < AZ_FIRST_USER_OPTION) {
      std::string sname = AztecOO_uppercase(entry.getValue(&dummy_string));
      Teuchos::map<std::string,int>& val_map = AztecOO_value_map();

      Teuchos::map<std::string,int>::iterator result = val_map.find(sname);
      if (result != val_map.end()) {
        azoo->SetAztecOption(offset, (*result).second);
        entry_used = true;
      }
    }
  }
  else if (entry.isType<double>()) {
    if (offset < AZ_FIRST_USER_PARAM) {
      double entry_value = entry.getValue(&dummy_double);
      azoo->SetAztecParam(offset, entry_value);
      entry_used = true;
    }
  }

  return(entry_used);
}

//=============================================================================
int AztecOO::SetParameters(Teuchos::ParameterList& parameterlist,
                           bool cerr_warning_if_unused)
{
  //cerr_warning_if_unused is an optional argument, default value is false.

  AztecOO_initialize_maps();

  Teuchos::map<std::string,int>& azoo_key_map = AztecOO_key_map();

  //Iterate the ParameterList, setting any options/parameters for which the
  //ParameterEntry's name matches a key-word that we recoginze.

  Teuchos::ParameterList::ConstIterator
    pl_iter = parameterlist.begin(),
    pl_end  = parameterlist.end();

  for(; pl_iter != pl_end; ++pl_iter) {
    //create an upper-case copy of the entry's name and prepend AZ_ if necessary
    std::string name = AztecOO_uppercase((*pl_iter).first);

    const Teuchos::ParameterEntry& entry = (*pl_iter).second;

    Teuchos::map<std::string,int>::iterator result = azoo_key_map.find(name);
    bool entry_used = false;

    if (result != azoo_key_map.end()) {
      entry_used = AztecOO_SetOptionOrParam((*result).second, entry, this);
    }

    if (cerr_warning_if_unused && !entry_used) {
      std::ostream& ostrm = err_stream_ ? *err_stream_ : std::cerr;
      ostrm << "AztecOO:SetParameters warning: '"<<name<<"' not used."<<std::endl;
    }

    name = " ";
  }

  return(0);
}
#endif //HAVE_AZTECOO_TEUCHOS

//=============================================================================
int AztecOO::SetAztecDefaults() {

  // If not in constructor, then we need to make sure any allocated memory is
  // deleted before we zero out pointers.
  if (inConstructor_) 
    Label_ = 0;
  else
    DeleteMemory();

  SetLabel("AztecOO Object");
  AZ_defaults(options_, params_);
  options_[AZ_poly_ord] = 1; // Redefine default value to be 1 (instead of 3).
  UserOperatorData_ = 0;
  UserMatrixData_ = 0;
  PrecOperatorData_ = 0;
  PrecMatrixData_ = 0;
  Problem_ = 0;
  X_ = 0;
  B_ = 0;
  ResidualVector_ = 0;
  
  N_local_ = 0;
  x_LDA_ = 0;
  x_ = 0;
  b_LDA_ = 0;
  b_ = 0;
  Amat_ = 0;
  Pmat_ = 0;
  Prec_ = 0;
  Scaling_ = 0;
  Scaling_created_ = false;
  StatusTest_ = 0;
  conv_info_ = 0;
  
  condest_ = (-1.0); 
  useAdaptiveDefaults_ = true;
  NumTrials_ = 0;
  maxFill_ = 0;
  maxKspace_ = 0;
  athresholds_ = 0;
  rthresholds_ = 0;
  condestThreshold_ = 0.0;
  procConfigSet_ = false;
  
  return(0);
  
}

//=============================================================================
int AztecOO::SetProblem(const Epetra_LinearProblem& prob,
                        bool call_SetPrecMatrix) {

  bool prevInConstructor = inConstructor_;
  inConstructor_ = true;  // Shut down complaints about zero pointers for a while
                          //  Although this routine is not a constructor, we treat it like one
  
  // Retain old problem (if any)
  Epetra_LinearProblem * OldProblem = Problem_;

  Problem_ = (Epetra_LinearProblem *) &prob; // Record this object for later access if needed

  // Try to cast operator to matrix 
  Epetra_RowMatrix * UserMatrix = dynamic_cast<Epetra_RowMatrix *>(prob.GetOperator());
  if (UserMatrix!=0) 
    SetUserMatrix(UserMatrix, call_SetPrecMatrix);
  else
    SetUserOperator(prob.GetOperator());
  
  SetLHS(prob.GetLHS());
  SetRHS(prob.GetRHS());

  // Set PDL and symmetry the first time the problem is set
  if (OldProblem==0) SetProblemOptions(prob.GetPDL(), prob.IsOperatorSymmetric());
  else if (OldProblem!=Problem_) {
    if (OldProblem->GetPDL()!=Problem_->GetPDL() || 
	OldProblem->IsOperatorSymmetric()!=Problem_->IsOperatorSymmetric()) {
      EPETRA_CHK_ERR(1); // Warn that the PDL and Symmetry will not be reset
    }
  }
    
  inConstructor_ = prevInConstructor;  // Revert to previous value

  return(0);

}


//=============================================================================
int AztecOO::SetUserOperator(Epetra_Operator * UserOperator) {

  if (UserOperator == 0 && inConstructor_ == true) return(0);
  if (UserOperator == 0) EPETRA_CHK_ERR(-1);

  if (Amat_ != 0) {
    AZ_matrix_destroy(&Amat_);
    Amat_ = 0;
  }

  if (UserOperatorData_!=0) delete UserOperatorData_;
  UserOperatorData_ = new OperatorData(UserOperator); // Initialize User Operator Data

  SetProcConfig(UserOperator->Comm());

  N_local_ =  UserOperator->OperatorRangeMap().NumMyPoints();

  Amat_ = AZ_matrix_create(N_local_);
  AZ_set_MATFREE(Amat_, (void *) UserOperatorData_, Epetra_Aztec_operatorvec);
  
  // Aztec needs upper bound for matrix norm if doing polynomial preconditioning
  if (UserOperator->HasNormInf()) {
    AZ_set_MATFREE_matrix_norm(Amat_, UserOperator->NormInf()); 
    AZ_set_MATNORMINF(Amat_, (void*)UserOperatorData_, Epetra_Aztec_operatornorminf);
  }

  // Set operator label
   const char * label = UserOperator->Label();
   if (label!=0)
     AZ_set_matrix_print_string(Amat_,label);

  return(0);
}

//=============================================================================
int AztecOO::SetUserMatrix(Epetra_RowMatrix * UserMatrix,
                           bool call_SetPrecMatrix)
{
  if (UserMatrix == 0 && inConstructor_ == true) return(0);
  if (UserMatrix == 0) EPETRA_CHK_ERR(-1);

  if (UserMatrixData_!=0) delete UserMatrixData_;
  UserMatrixData_ = new MatrixData(UserMatrix); // Initialize user matrix data

  SetProcConfig(UserMatrix->Comm());
  EPETRA_CHK_ERR(SetUserOperator(UserMatrix));
  AZ_set_MATFREE(Amat_, (void *) UserMatrixData_, Epetra_Aztec_matvec);
  AZ_set_MATNORMINF(Amat_, (void*)UserMatrixData_, Epetra_Aztec_matnorminf);
  //int N_ghost = UserMatrix->NumMyCols() - UserMatrix->NumMyRows();
  int N_ghost = 0;
  if (UserMatrix->RowMatrixImporter()!=0)
    N_ghost = UserMatrix->RowMatrixImporter()->NumRecv();
  
  AZ_set_MATFREE_getrow(Amat_, (void *) UserMatrixData_, Epetra_Aztec_getrow,
			Epetra_Aztec_comm_wrapper,N_ghost,proc_config_);

  // Set matrix label
   const char * label = UserMatrix->Label();
   if (label!=0)
     AZ_set_matrix_print_string(Amat_,label);

  // If preconditioner not defined, set up to possibly use native Aztec precons
  if (Prec_==0 || call_SetPrecMatrix) {
    EPETRA_CHK_ERR(SetPrecMatrix(UserMatrix));
  }
  return(0);
}



//=============================================================================
void AztecOO::SetOutputStream(std::ostream& ostrm)
{
  out_stream_ = &ostrm;
  AZOO_set_stream_out(ostrm);
}

//=============================================================================
void AztecOO::SetErrorStream(std::ostream& ostrm)
{
  err_stream_ = &ostrm;
  AZOO_set_stream_err(ostrm);
}

//=============================================================================
int AztecOO::SetPrecMatrix(Epetra_RowMatrix * PrecMatrix) {

  if (PrecMatrix == 0 && inConstructor_ == true) return(0);
  if (PrecMatrix == 0) EPETRA_CHK_ERR(-1);
  if (Prec_!=0) {
    AZ_precond_destroy(&Prec_); 
    Prec_ = 0;
  }
  if (Pmat_ != 0) {
    AZ_rm_context(options_,params_,Pmat_->data_org);
    AZ_matrix_destroy(&Pmat_);
    Pmat_ = 0;
  }

  if (PrecMatrixData_!=0) delete PrecMatrixData_;
  PrecMatrixData_ = new MatrixData(PrecMatrix); // Initialize preconditioner matrix data
  SetProcConfig(PrecMatrix->Comm());
  Pmat_ = AZ_matrix_create(N_local_);
  AZ_set_MATFREE(Pmat_, (void *) PrecMatrixData_, Epetra_Aztec_matvec);
  AZ_set_MATNORMINF(Amat_, (void*)UserMatrixData_, Epetra_Aztec_matnorminf);
  
  // Aztec needs upper bound for matrix norm if doing polynomial preconditioning
  if (PrecMatrix->HasNormInf())
		AZ_set_MATFREE_matrix_norm(Pmat_, PrecMatrix->NormInf()); 
  //int N_ghost = PrecMatrix->NumMyCols() - PrecMatrix->NumMyRows();
  int N_ghost = 0;
  if (PrecMatrix->RowMatrixImporter()!=0)
    N_ghost = PrecMatrix->RowMatrixImporter()->NumRecv();

  AZ_set_MATFREE_getrow(Pmat_, (void *) PrecMatrixData_, Epetra_Aztec_getrow,
			Epetra_Aztec_comm_wrapper,N_ghost,proc_config_);
    

    Prec_ = AZ_precond_create(Pmat_, AZ_precondition, NULL);


  return(0);
}


//=============================================================================
int AztecOO::SetPrecOperator(Epetra_Operator * PrecOperator) {

  if (PrecOperator == 0 && inConstructor_ == true) return(0);
  if (PrecOperator == 0) EPETRA_CHK_ERR(-1);

  // Get rid of any other preconditioner
  if (Prec_!=0) {
    AZ_precond_destroy(&Prec_); 
    Prec_ = 0;
  }
  if (Pmat_ != 0) {
    AZ_matrix_destroy(&Pmat_);
    Pmat_ = 0;
  }

	if (PrecOperatorData_!=0) delete PrecOperatorData_;
  PrecOperatorData_ = new OperatorData(PrecOperator); // Initialize preconditioner operator data
  SetProcConfig(PrecOperator->Comm());

  if (Amat_==0) EPETRA_CHK_ERR(-2); // UserOperator must be defined first

   Prec_ = AZ_precond_create(Amat_, Epetra_Aztec_precond, (void *) PrecOperatorData_);

   options_[AZ_precond] = AZ_user_precond;
   const char * label = PrecOperator->Label();
   if (label==0)
     AZ_set_precond_print_string(Prec_,"User-defined preconditioner");
   else
     AZ_set_precond_print_string(Prec_,label);
   return(0);
}
//=============================================================================
int AztecOO::SetStatusTest(AztecOO_StatusTest * StatusTest) {

  if (StatusTest == 0) EPETRA_CHK_ERR(-1);
  if (Amat_==0) EPETRA_CHK_ERR(-2); // UserOperator Amat_ must be defined first
  if (UserOperatorData_->A==0) EPETRA_CHK_ERR(-3); //  ust be defined so we know how to make residual vector


  // Create AZ_CONVERGE_STRUCT if needed, associate with Amat
  if (conv_info_==0) {
    double dummy; // Need something to pass in to vector constructor, will not be used
    ResidualVector_ = new Epetra_Vector(View, UserOperatorData_->A->OperatorRangeMap(), &dummy);
    conv_info_ = AZ_converge_create();
    conv_info_->scaling = Scaling_;
    conv_info_->res_vec_object = (void *) ResidualVector_;
    Amat_->conv_info = conv_info_;
  }
  StatusTest_ = StatusTest;
  options_[AZ_conv] = AZTECOO_conv_test;
  conv_info_->conv_object = (void *) StatusTest_;
  conv_info_->conv_function = AztecOO_StatusTest_wrapper; 

  return(0);
}
//=============================================================================
int AztecOO::SetLHS(Epetra_MultiVector * X) {

  if (X == 0 && inConstructor_ == true) return(0);
  if (X == 0) EPETRA_CHK_ERR(-1);
  X_ = X;
  X_->ExtractView(&x_, &x_LDA_);
  return(0);
}
//=============================================================================
int AztecOO::SetRHS(Epetra_MultiVector * B) {

  if (B == 0 && inConstructor_ == true) return(0);
  if (B == 0) EPETRA_CHK_ERR(-1);
  B_ = B;
  B_->ExtractView(&b_, &b_LDA_);

  return(0);
}
//=============================================================================
int AztecOO::UnsetLHSRHS() {
  X_ = 0;
  x_ = 0;
  x_LDA_ = -1;
  B_ = 0;
  b_ = 0;
  b_LDA_ = -1;
  return(0);
}
//=============================================================================
int AztecOO::SetProcConfig(const Epetra_Comm & Comm) {

  if (!procConfigSet_) {
#ifdef AZTEC_MPI
    const Epetra_MpiComm* comm1 = dynamic_cast<const Epetra_MpiComm*> (&Comm);
    if (comm1 == 0) {
      std::ostream& ostrm = err_stream_ ? *err_stream_ : std::cerr;
      ostrm << "AztecOO::SetProcConfig ERROR, failed to dynamic_cast "
        << "Comm to Epetra_MpiComm."<<std::endl;
      return(-1);
    }
    AZ_set_proc_config(proc_config_, comm1->Comm());
#else
    AZ_set_proc_config(proc_config_, 0);
#endif
    procConfigSet_ = true;
  }
  return(0);
}
//=============================================================================
void AztecOO::DeleteAzArrays() {
  if (proc_config_!=0) {delete [] proc_config_; proc_config_ = 0;}
  if (options_!=0) {delete [] options_; options_ = 0;}
  if (params_!=0) {delete [] params_; params_ = 0;}
  if (status_!=0) {delete [] status_; status_ = 0;}
  if (athresholds_!=0) {delete [] athresholds_; athresholds_ = 0;}
  if (rthresholds_!=0) {delete [] rthresholds_; rthresholds_ = 0;}
}

//=============================================================================
int AztecOO::SetProblemOptions(ProblemDifficultyLevel PDL,
				bool ProblemSymmetric)
{
  if (ProblemSymmetric)
    {
      SetAztecOption(AZ_solver, AZ_cg);
      switch (PDL)
      {
      case easy:
        SetAztecOption(AZ_poly_ord, 1);
        SetAztecOption(AZ_precond, AZ_Jacobi);
        break;

      case moderate:

        SetAztecOption(AZ_precond, AZ_dom_decomp);
        SetAztecOption(AZ_subdomain_solve, AZ_icc);
        break;

      case hard:
      case unsure:
      default:
        SetAztecOption(AZ_precond, AZ_dom_decomp);
        SetAztecOption(AZ_subdomain_solve, AZ_icc);
        SetAztecParam(AZ_omega, 1.2);
      }
    }
  else
    {
      switch (PDL)
      {
      case easy:
        SetAztecOption(AZ_poly_ord, 1);
        SetAztecOption(AZ_precond, AZ_Jacobi);
        SetAztecOption(AZ_solver, AZ_bicgstab);
        break;

      case moderate:

        SetAztecOption(AZ_precond, AZ_dom_decomp);
        SetAztecOption(AZ_subdomain_solve, AZ_ilu);
        SetAztecOption(AZ_solver, AZ_gmres);
        break;

      case hard:
      case unsure:
      default:
        SetAztecOption(AZ_precond, AZ_dom_decomp);
        SetAztecOption(AZ_subdomain_solve, AZ_ilut);
        SetAztecOption(AZ_overlap, 1);
        SetAztecParam(AZ_ilut_fill, 3.0);
        //SetAztecParam(AZ_drop, 0.01);
        //SetAztecOption(AZ_kspace, 1000);
      }
    }

  return(0);
}

//=============================================================================
int AztecOO::SetPreconditioner(AZ_PREC_FUN  prec_function,
			       void *p_data)
{
  if (Pmat_==0) EPETRA_CHK_ERR(-1); // No matrix yet
  EPETRA_CHK_ERR(DestroyPreconditioner()); // Delete existing preconditioner if one exists
  if (Pmat_==0) EPETRA_CHK_ERR(-1); // No Pmat defined.
  Prec_ = AZ_precond_create(Pmat_,prec_function, p_data);
  options_[AZ_precond] = AZ_user_precond;

  return(0);
}
//
//=============================================================================
int AztecOO::ConstructPreconditioner(double & condest)
{
  if (PrecMatrixData_==0) EPETRA_CHK_ERR(-1); // No matrix yet

  Epetra_RowMatrix * PrecMatrix = PrecMatrixData_->A; // Extract Preconditioner matrix

  int precond_flag = options_[AZ_precond];

  if (precond_flag) {

    // Create default Aztec preconditioner if one not defined
    if (Prec_==0) {
      if (Pmat_==0)  EPETRA_CHK_ERR(-2); // No Pmat to use for building preconditioner
      Prec_ = AZ_precond_create(Pmat_, AZ_precondition, NULL);
    }
  
    AZ_mk_context(options_, params_, Pmat_->data_org, Prec_, proc_config_);

    int NN = PrecMatrix->NumMyCols();
    double * condvec = new double[NN];
    for (int i = 0 ; i < N_local_ ; i++ ) condvec[i] = 1.0;
    Prec_->prec_function(condvec,options_,proc_config_,params_,Pmat_,Prec_);
    condest_ = 0.0;
    {for (int i=0; i<N_local_; i++)
      if (fabs(condvec[i]) > condest_)
	condest_ = fabs(condvec[i]);}
    delete [] condvec;
    options_[AZ_pre_calc] = AZ_reuse;
    double tmp_condest = condest_;
    // if any processor has a tmp_condest==0.0, then it has a singular preconditioner, check for that first
    PrecMatrix->Comm().MinAll(&tmp_condest, &condest_, 1); // Get the min of all condition estimates
    if (condest_!=0.0)
      PrecMatrix->Comm().MaxAll(&tmp_condest, &condest_, 1); // Get the worst of all condition estimates

    condest = condest_;
  }

  return(0);
}
//
//=============================================================================
int AztecOO::DestroyPreconditioner() {

  if (Prec_!=0) {
    AZ_precond_destroy(&Prec_);
    Prec_ = 0;
    options_[AZ_pre_calc] = AZ_calc;
  }
  return(0);
}
//=============================================================================
int AztecOO::SetMatrixName(int label)

{
  Amat_->data_org[AZ_name] = label;

  return(0);
}

//=============================================================================
int AztecOO::recursiveIterate(int MaxIters, double Tolerance)
{
  SetAztecOption(AZ_max_iter, MaxIters);
  SetAztecParam(AZ_tol, Tolerance);

  int prec_allocated = 0;
  if (Prec_ == 0) {
    if (options_[AZ_precond] == AZ_user_precond) EPETRA_CHK_ERR(-10); // Cannot have user prec==0
    if (Pmat_!=0) {
      Prec_ = AZ_precond_create(Pmat_, AZ_precondition, NULL);
      prec_allocated = 1;
    }
  }

  for (int i=0; i< X_->NumVectors(); i++) {
    options_[AZ_recursion_level]++;
    AZ_oldsolve(x_ + i*x_LDA_, b_ + i*b_LDA_, options_, params_, status_, proc_config_,
		Amat_, Prec_, Scaling_);
    options_[AZ_recursion_level]--;
    if (prec_allocated == 1) {
      AZ_precond_destroy(&Prec_);
      Prec_ = 0;
      prec_allocated = 0;
    }
  }          


  // Determine end status

  //If status is AZ_normal or AZ_maxits, then make a quick exit:
  if (status_[AZ_why]==AZ_normal) {
    return(0);
  }
  else if (status_[AZ_why]==AZ_maxits) {
    return(1);
  }

  //If status is not AZ_normal or AZ_maxits, then we have to decide
  //which error-code to return, and also decide whether to print an
  //error-message to std::cerr, and what to print.

  int MyPID = X_->Map().Comm().MyPID();

  bool print_msg = (MyPID == 0) &&
                   (Epetra_Object::GetTracebackMode()>0) &&
                   (options_[AZ_diagnostics]!=AZ_none);

  std::ostream* ostrm = 0;
  if (print_msg) {
    ostrm = err_stream_ ? err_stream_ : &std::cerr;
  }

  int ierr = 0;
  if (status_[AZ_why] == AZ_param) {
    ierr = -1;
    if (print_msg) {
      *ostrm << "Aztec status AZ_param: option not implemented" << std::endl;
    }
  }
  else if (status_[AZ_why] == AZ_breakdown) {
    ierr = -2;
    if (print_msg) {
      *ostrm << "Aztec status AZ_breakdown: numerical breakdown" << std::endl;
    }
  }
  else if (status_[AZ_why] == AZ_loss) {
    ierr = -3;
    if (print_msg) {
      *ostrm << "Aztec status AZ_loss: loss of precision" << std::endl;
    }
  }
  else if (status_[AZ_why] == AZ_ill_cond) {
    ierr = -4;
    if (print_msg) {
      *ostrm << "Aztec status AZ_ill_cond: GMRES hessenberg ill-conditioned"
           << std::endl;
    }
  }
  else {
    throw B_->ReportError("Internal AztecOO Error", -5);
  }

  return(ierr);
}
//=============================================================================
int AztecOO::Iterate(long long MaxIters, double Tolerance)
{
  // jhurani@txcorp.com -- The "units/type" of MaxIters should be same as
  // the "units/type" of matrix size, and it can be long long after Epetra
  // has support for long long.  Hence the option should be set using
  // long long, and if called using int, the int value will be automatically
  // upgraded.  However, currently we can't store MaxIters value that is
  // greater than maximum int, so we error out.  This seems like a reasonable
  // design assumption since one would rarely want O(10^9) billion iterations.
  // In future, long long is stored in the options, the interface will not
  // have to change.

  if (X_ == 0 || B_ == 0 || UserOperatorData_ == 0 || MaxIters > INT_MAX) EPETRA_CHK_ERR(-11);

  if (UserMatrixData_!=0)
    if (GetUserMatrix()!=0) {

      long long nnz = GetUserMatrix()->NumGlobalNonzeros64();
      if (nnz==0) {
	EPETRA_CHK_ERR(-12); // Matrix has no entries.
      }
    }

  // MaxIters not greater than INT_MAX because of check above.
  int MaxIters_int = static_cast<int>(MaxIters);
  SetAztecOption(AZ_max_iter, MaxIters_int);
  SetAztecParam(AZ_tol, Tolerance);

  
  int prec_allocated = 0;
  if (Prec_ == 0) {
    if (options_[AZ_precond] == AZ_user_precond) EPETRA_CHK_ERR(-10); // Cannot have user prec==0
    if (Pmat_!=0) {
      Prec_ = AZ_precond_create(Pmat_, AZ_precondition, NULL);
      prec_allocated = 1;
    }
  }

  if (Scaling_ == 0 && options_[AZ_scaling] != AZ_none) {
    Scaling_ = AZ_scaling_create();
    Scaling_->scale = AztecOO_scale_epetra;
    Scaling_created_ = true;
  }

  AZ_iterate(x_, b_, options_, params_, status_, proc_config_,
	     Amat_, Prec_, Scaling_);

  if (prec_allocated == 1) {
    AZ_precond_destroy(&Prec_);
    Prec_ = 0;
    prec_allocated = 0;
  }
  
  if (options_[AZ_keep_info] != 1 &&
      Scaling_ != 0 && Scaling_->scale != 0) {
    Scaling_->scale(AZ_DESTROY_SCALING_DATA, Amat_, options_, x_, b_,
                    proc_config_, Scaling_);
  }

  // Determine end status

  //If status is AZ_normal or AZ_maxits, then make a quick exit:
  if (status_[AZ_why]==AZ_normal) {
    return(0);
  }
  else if (status_[AZ_why]==AZ_maxits) {
    return(1);
  }

  //If status is not AZ_normal or AZ_maxits, then we have to decide
  //which error-code to return, and also decide whether to print an
  //error-message to std::cerr, and what to print.

  int MyPID = X_->Map().Comm().MyPID();

  bool print_msg = (MyPID == 0) &&
                   (Epetra_Object::GetTracebackMode()>0) &&
                   (options_[AZ_diagnostics]!=AZ_none);

  std::ostream* ostrm = 0;
  if (print_msg) {
    ostrm = err_stream_ ? err_stream_ : &std::cerr;
  }

  int ierr = 0;
  if (status_[AZ_why] == AZ_param) {
    ierr = -1;
    if (print_msg) {
      *ostrm << "Aztec status AZ_param: option not implemented" << std::endl;
    }
  }
  else if (status_[AZ_why] == AZ_breakdown) {
    ierr = -2;
    if (print_msg) {
      *ostrm << "Aztec status AZ_breakdown: numerical breakdown" << std::endl;
    }
  }
  else if (status_[AZ_why] == AZ_loss) {
    ierr = -3;
    if (print_msg) {
      *ostrm << "Aztec status AZ_loss: loss of precision" << std::endl;
    }
  }
  else if (status_[AZ_why] == AZ_ill_cond) {
    ierr = -4;
    if (print_msg) {
      *ostrm << "Aztec status AZ_ill_cond: GMRES hessenberg ill-conditioned"
           << std::endl;
    }
  }
  else {
    throw B_->ReportError("Internal AztecOO Error", -5);
  }

  return(ierr);
}

//=============================================================================
int AztecOO::Iterate(Epetra_RowMatrix * A,
                      Epetra_MultiVector * X,
                      Epetra_MultiVector * B,
                      int MaxIters, double Tolerance)
{
  SetUserMatrix(A);
  
  SetLHS(X);
  SetRHS(B);

  EPETRA_CHK_ERR( Iterate(MaxIters, Tolerance) );
  return(0);
}

//=============================================================================
int AztecOO::SetAdaptiveParams(int NumTrials, double * athresholds, double * rthresholds,
				double condestThreshold, double maxFill, int maxKspace) {

  if (athresholds_!=0) delete [] athresholds_;
  if (rthresholds_!=0) delete [] rthresholds_;

  NumTrials_ = NumTrials;
  maxFill_ = maxFill;
  maxKspace_ = maxKspace;
  athresholds_ = new double[NumTrials];
  rthresholds_ = new double[NumTrials];
  for (int i=0; i<NumTrials; i++) {
    athresholds_[i] = athresholds[i];
    rthresholds_[i] = rthresholds[i];
  }
  if (condestThreshold>0) condestThreshold_ = condestThreshold;
  useAdaptiveDefaults_ = false;
  return(0);
}
//=============================================================================
int AztecOO::AdaptiveIterate(int MaxIters, int MaxSolveAttempts, double Tolerance) {

  // Check if adaptive strategy is appropriate (only works for domain decomp and if subdomain
  // solve is not AZ_lu.  If not, call standard Iterate() method

  if (options_[AZ_precond] != AZ_dom_decomp) {
    EPETRA_CHK_ERR(Iterate(MaxIters,Tolerance));
    return(0);
  }
  if (options_[AZ_subdomain_solve] == AZ_lu) {
    EPETRA_CHK_ERR(Iterate(MaxIters,Tolerance));
    return(0);
  }

  int NumSolveAttempts = 0;
  SetAztecOption(AZ_max_iter, MaxIters);
  SetAztecParam(AZ_tol, Tolerance);

  // Construct adaptive strategy if necessary
  if (useAdaptiveDefaults_) {

    if (options_[AZ_subdomain_solve] == AZ_bilu) {
      int NumTrials = 3;
      double athresholds[] = {0.0, 1.0E-14, 1.0E-3};
      double rthresholds[] = {0.0, 1.0E-14, 1.0E-3};
      double condestThreshold = 1.0E16;
      double maxFill = 4.0;
      int maxKspace = 4*options_[AZ_kspace];
      SetAdaptiveParams(NumTrials, athresholds, rthresholds, condestThreshold, maxFill, maxKspace);
    }
    else {
      int NumTrials = 7;
      double athresholds[] = {0.0, 1.0E-12, 1.0E-12, 1.0E-5, 1.0E-5, 1.0E-2, 1.0E-2};
      double rthresholds[] = {1.0, 1.0,     1.01,    1.0,    1.01,   1.01,   1.1   };
      double condestThreshold = 1.0E16;
      double maxFill = 4.0;
      int maxKspace = 4*options_[AZ_kspace];
      SetAdaptiveParams(NumTrials, athresholds, rthresholds, condestThreshold, maxFill, maxKspace);
    }
  }

  // If no trials were defined, just call regular Iterate method

  if (NumTrials_<1) {
    EPETRA_CHK_ERR(Iterate(MaxIters,Tolerance));
    return(0);
  }


  bool FirstCallToSolver = true;
  double oldResid=0.0;

  // Keep copy of best old solution
  Epetra_MultiVector Xold(*X_);
  
  // ********************************************************************************
  //  Phase:  Tweak fill level and drop tolerances
  // ********************************************************************************
  
  
  double fill;
  if (options_[AZ_subdomain_solve] == AZ_ilut) fill = params_[AZ_ilut_fill];
  else fill = (double) options_[AZ_graph_fill];

  double curMaxFill = EPETRA_MAX(fill, maxFill_);
  int curMaxKspace = EPETRA_MAX(options_[AZ_kspace], maxKspace_);

  // GMRES is only solver sensitive to kspace
  if (options_[AZ_solver]!=AZ_gmres) curMaxKspace = options_[AZ_kspace]; 
  int kspace =options_[AZ_kspace];

  int maxLossTries = 6;
  int numLossTries = 0;
 

  // *****************************************************************************
  // Master "while" loop 
  // *****************************************************************************

  while ((FirstCallToSolver || status_[AZ_why]!=AZ_normal) && NumSolveAttempts<MaxSolveAttempts) {

    SetAztecOption(AZ_kspace,kspace);

    if (options_[AZ_subdomain_solve] == AZ_ilut) params_[AZ_ilut_fill] = fill;
    else options_[AZ_graph_fill] = (int) fill;

    // ********************************************************************************
    //  Phase:  Find a preconditioner whose condest is below the condest threshold
    // ********************************************************************************

    // Start with first trial
    int curTrial = 0;
    // Initial guess for conditioner number
    condest_ = condestThreshold_; 
    
    // While condest is too large or iterative solver says the preconditioned subspace is too
    // ill-conditioned, try to form a better-conditioned preconditioner
    
    while ( curTrial<NumTrials_ && (condest_>=condestThreshold_ || condest_==0.0)) {
      if (Prec_!=0) DestroyPreconditioner(); // Get rid of any existing preconditioner
      
      SetAztecParam(AZ_athresh, athresholds_[curTrial]); // Set threshold values
      SetAztecParam(AZ_rthresh, rthresholds_[curTrial]);
      
      // Preconstruct preconditioner and get a condition number estimate
      ConstructPreconditioner(condest_);
      curTrial++;
    }

    // ********************************************************************************
    // Try solver now, but test if solver says  preconditioned subspace is too
    // ill-conditioned, if so try to form a better-conditioned preconditioner
    // ********************************************************************************

    if (!FirstCallToSolver) {
      // keep old residual and solution to see if progress is made (will keep old solution if not)
      oldResid = status_[AZ_r];
      Xold = *X_;
      // If a relative residual is being used as a stopping criterion,
      // adjust the tolerance to account for any progress made in the solution
      // The test for "status_[AZ_scaled_r]>0.0" insures that NaN will be ignored
      if (options_[AZ_conv]==AZ_r0 && status_[AZ_scaled_r]>0.0) {
	Tolerance = Tolerance/status_[AZ_scaled_r];
	SetAztecParam(AZ_tol, Tolerance);
      }
    }
    AZ_iterate(x_, b_, options_, params_, status_, proc_config_,
	       Amat_, Prec_, Scaling_);

    if (!FirstCallToSolver) {
      if (!(oldResid>status_[AZ_r])) *X_ = Xold;  // keep old solution
    }

    NumSolveAttempts++;
    FirstCallToSolver = false;
    
    while ( curTrial<NumTrials_ && 
	    (status_[AZ_why]==AZ_ill_cond || status_[AZ_why]==AZ_breakdown)) {
      
      if (Prec_!=0) DestroyPreconditioner(); // Get rid of any existing preconditioner
      
      SetAztecParam(AZ_athresh, athresholds_[curTrial]); // Set next threshold values
      SetAztecParam(AZ_rthresh, rthresholds_[curTrial]);
      
      // Preconstruct preconditioner and get a condition number estimate
      ConstructPreconditioner(condest_);
      curTrial++;

      // keep old residual and solution to see if progress is made (will keep old solution if not)
      oldResid = status_[AZ_r];
      Xold = *X_;
      // If a relative residual is being used as a stopping criterion,
      // adjust the tolerance to account for any progress made in the solution
      // The test for "status_[AZ_scaled_r]>0.0" insures that NaN will be ignored
      if (options_[AZ_conv]==AZ_r0 && status_[AZ_scaled_r]>0.0) {
	Tolerance = Tolerance/status_[AZ_scaled_r];
	SetAztecParam(AZ_tol, Tolerance);
      }

      AZ_iterate(x_, b_, options_, params_, status_, proc_config_,
		 Amat_, Prec_, Scaling_);
      NumSolveAttempts++;
      
      if (!(oldResid>status_[AZ_r])) *X_ = Xold;  // keep old solution
    }

    if (status_[AZ_why]==AZ_loss) { // This means we should try to solve one more time only
      if (numLossTries<maxLossTries) {
	numLossTries++;
      }
    }
    else if (status_[AZ_why]==AZ_maxits) { // This means we need more robust preconditioner
      if (fill<curMaxFill) {
	fill = EPETRA_MIN(2*fill, curMaxFill); // double fill
	curTrial = 0; // force new pass through preconditioner setup
      }
      else if (options_[AZ_subdomain_solve] == AZ_ilut && params_[AZ_drop]>0) {
	params_[AZ_drop] = 0.0; // if nonzero drop used with ILUT, try one more time with drop = 0
	curTrial = 0; // force new pass through preconditioner setup
      }
      else if (kspace<curMaxKspace && kspace<MaxIters) {
	kspace = EPETRA_MIN(curMaxKspace,2*kspace); // double kspace and try again
      }
    }
    else break;
  }
  if (Prec_!=0) DestroyPreconditioner(); // Delete preconditioner
  // Determine end status

  //If status is AZ_normal or AZ_maxits, then make a quick exit:
  if (status_[AZ_why]==AZ_normal) {
    return(0);
  }
  else if (status_[AZ_why]==AZ_maxits) {
    return(1);
  }

  //If status is not AZ_normal or AZ_maxits, then we have to decide
  //which error-code to return, and also decide whether to print an
  //error-message to std::cerr, and what to print.

  int MyPID = X_->Map().Comm().MyPID();

  bool print_msg = (MyPID == 0) &&
                   (Epetra_Object::GetTracebackMode()>0) &&
                   (options_[AZ_diagnostics]!=AZ_none);

  std::ostream* ostrm = 0;
  if (print_msg) {
    ostrm = err_stream_ ? err_stream_ : &std::cerr;
  }

  int ierr = 0;
  if (status_[AZ_why] == AZ_param) {
    ierr = -1;
    if (print_msg) {
      *ostrm << "Aztec status AZ_param: option not implemented" << std::endl;
    }
  }
  else if (status_[AZ_why] == AZ_breakdown) {
    ierr = -2;
    if (print_msg) {
      *ostrm << "Aztec status AZ_breakdown: numerical breakdown" << std::endl;
    }
  }
  else if (status_[AZ_why] == AZ_loss) {
    ierr = -3;
    if (print_msg) {
      *ostrm << "Aztec status AZ_loss: loss of precision" << std::endl;
    }
  }
  else if (status_[AZ_why] == AZ_ill_cond) {
    ierr = -4;
    if (print_msg) {
      *ostrm << "Aztec status AZ_ill_cond: GMRES hessenberg ill-conditioned"
           << std::endl;
    }
  }
  else {
    throw B_->ReportError("Internal AztecOO Error", -5);
  }

  return(ierr);
}

//=============================================================================
void AztecOO::PrintLinearSystem(const char* name)
{
  std::ostringstream aostr;
  std::ostringstream xostr;
  std::ostringstream bostr;

  aostr << "A_";
  xostr << "X_";
  bostr << "B_";

  if (name != 0) {
    aostr << name;
    xostr << name;
    bostr << name;
  }

  int MyPID = 0;
  int numProcs = 1;

  const Epetra_Comm* comm = 0;

  if (X_ != 0) comm = &(X_->Map().Comm());
  else {
    if (B_ != 0) comm = &(B_->Map().Comm());
    else if (GetUserMatrix() != 0) {
      comm = &(GetUserMatrix()->RowMatrixRowMap().Comm());
    }
  }

  std::ostream& ostrm = err_stream_ ? *err_stream_ : std::cerr;

  if (comm == 0) {

    ostrm << "Aztec::PrintLinearSystem: matrix and vectors seem to all"
        << " be NULL." << std::endl;
    return;
  }

  MyPID = comm->MyPID();
  numProcs = comm->NumProc();

  std::ios::openmode mode = std::ios::app;
  if (MyPID == 0) mode = std::ios::out;

  for(int p=0; p<numProcs; ++p) {
    comm->Barrier();

    if (p != MyPID) continue;

    if (UserMatrixData_ != 0) {
      std::string afilename = aostr.str();

      if (GetUserMatrix() != 0) {
        Epetra_CrsMatrix* crs =
          dynamic_cast<Epetra_CrsMatrix*>(GetUserMatrix());
        if (crs != 0) {
          std::ofstream afile(afilename.c_str(), mode);
          if (afile && !afile.bad()) {
            crs->Print(afile);
          }
          else {
            ostrm << "AztecOO::PrintLinearSystem: unable to open file '"
                 << afilename << "'."<<std::endl;
          }
        }
        else {

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES // CJ: TODO FIXME
          // ifdef-out dynamic_cast to Epetra_VbrMatrix to avoid linker
          // error of missing typeinfo.  Issue in Epetra 64 bit.

          Epetra_VbrMatrix* vbr =
            dynamic_cast<Epetra_VbrMatrix*>(GetUserMatrix());
          if (vbr != 0) {
            std::ofstream afile(afilename.c_str(), mode);
            if (afile && !afile.bad()) {
              crs->Print(afile);
            }
            else {
              ostrm << "AztecOO::PrintLinearSystem: unable to open file '"
                   << afilename << "'."<<std::endl;
            }
          }
#endif
        }
      }
    }

    if (X_ != 0) {
      std::string xfilename = xostr.str();

      std::ofstream xfile(xfilename.c_str(), mode);
      if (xfile && !xfile.bad()) {
        X_->Print(xfile);
      }
      else {
        ostrm << "AztecOO::PrintLinearSystem: unable to open file '"
             << xfilename << "'."<<std::endl;
      }
    }

    if (B_ != 0) {
      std::string bfilename = bostr.str();

      std::ofstream bfile(bfilename.c_str(), mode);
      if (bfile && !bfile.bad()) {
        B_->Print(bfile);
      }
      else {
        ostrm << "AztecOO::PrintLinearSystem: unable to open file '"
             << bfilename << "'."<<std::endl;
      }
    }
  }
}

//=============================================================================
const char * AztecOO::GetLabel() const {
  return(Label_);
}
//=============================================================================
void AztecOO::SetLabel(const char * const Label_in)
{ 
  if (Label_!=0) {
    delete [] Label_;
    Label_ = 0;
  }
  if (Label_in==0) return;
  Label_ = new char[strlen(Label_in)+1];
  strcpy(Label_,Label_in);
  return;
}
//=============================================================================
double Epetra_Aztec_matnorminf(AZ_MATRIX* Amat)
{
  AztecOO::MatrixData* Data = (AztecOO::MatrixData*)AZ_get_matvec_data(Amat);
  Epetra_RowMatrix* A = (Epetra_RowMatrix*)Data->A;
  return( A->NormInf() );
}

//=============================================================================
void Epetra_Aztec_matvec(double x[], double y[],
                         AZ_MATRIX *Amat, int proc_config[])
{
  (void)proc_config;

  AztecOO::MatrixData * Data = (AztecOO::MatrixData *) AZ_get_matvec_data(Amat);
  Epetra_RowMatrix * A = (Epetra_RowMatrix *) Data->A;
	Epetra_Vector * X = (Epetra_Vector *) Data->X;
	Epetra_Vector * Y = (Epetra_Vector *) Data->Y;

  if (X==0) {
		X = new Epetra_Vector(View, A->OperatorDomainMap(), x);
		X->SetLabel("Epetra_Aztec_matvec X Vector");
		Data->X = X;
		Y = new Epetra_Vector(View, A->OperatorRangeMap(), y);
		Y->SetLabel("Epetra_Aztec_matvec Y Vector");
		Data->Y = Y;
	}
	else {
		X->ResetView(x);
		Y->ResetView(y);
	}

  int ierr = A->Apply(*X, *Y);
  if (ierr!=0) throw X->ReportError("Error in call to Epetra_Operator for preconditioner", ierr);

}

//=============================================================================
double Epetra_Aztec_operatornorminf(AZ_MATRIX* Amat)
{
  AztecOO::OperatorData * Data = (AztecOO::OperatorData *) AZ_get_matvec_data(Amat);
  Epetra_Operator * A = (Epetra_Operator *) Data->A;

  return( A->HasNormInf() ? A->NormInf() : -1.0);
}

//=============================================================================
void Epetra_Aztec_operatorvec(double x[], double y[],
                              AZ_MATRIX *Amat, int proc_config[])
{
  (void)proc_config;
  AztecOO::OperatorData * Data = (AztecOO::OperatorData *) AZ_get_matvec_data(Amat);
  Epetra_Operator * A = (Epetra_Operator *) Data->A;
	Epetra_Vector * X = (Epetra_Vector *) Data->X;
	Epetra_Vector * Y = (Epetra_Vector *) Data->Y;

  if (X==0) {
		X = new Epetra_Vector(View, A->OperatorDomainMap(), x);
		X->SetLabel("Epetra_Aztec_operatorvec X Vector");
		Data->X = X;
		Y = new Epetra_Vector(View, A->OperatorRangeMap(), y);
		Y->SetLabel("Epetra_Aztec_operatorvec Y Vector");
		Data->Y = Y;
	}
	else {
		X->ResetView(x);
		Y->ResetView(y);
	}

  int ierr = A->Apply(*X, *Y);
  if (ierr!=0) throw X->ReportError("Error in call to Epetra_Operator for preconditioner", ierr);

}

//=============================================================================
void Epetra_Aztec_precond(double x[], int input_options[],
			  int proc_config[], double input_params[],
                          AZ_MATRIX *Amat, AZ_PRECOND *prec)
{
  (void)input_options;
  (void)proc_config;
  (void)input_params;
  (void)Amat;
  AztecOO::OperatorData * Data = (AztecOO::OperatorData *) AZ_get_precond_data(prec);
  Epetra_Operator * A = (Epetra_Operator *) Data->A;
	Epetra_Vector * X = (Epetra_Vector *) Data->X;
	Epetra_Vector * Y = (Epetra_Vector *) Data->Y;

  if (X==0) {
		X = new Epetra_Vector(View, A->OperatorDomainMap(), x);
		X->SetLabel("Epetra_Aztec_precond X Vector");
		Data->X = X;
		Y = new Epetra_Vector(View, A->OperatorRangeMap(), x);
		Y->SetLabel("Epetra_Aztec_precond Y Vector");
		Data->Y = Y;
	}
	else {
		X->ResetView(x);
		Y->ResetView(x);
	}

  int ierr = A->ApplyInverse(*X, *Y);
  if (ierr!=0) throw X->ReportError("Error in call to Epetra_Operator for preconditioner", ierr);

}

//=============================================================================
int Epetra_Aztec_getrow(int columns[], double values[],
     int row_lengths[], AZ_MATRIX *Amat, int N_requested_rows,
		       int requested_rows[], int allocated_space){
/*
 * Supply local matrix (without ghost node columns) for rows given by
 * requested_rows[0 ... N_requested_rows-1].  Return this information in
 * 'row_lengths, columns, values'.  If there is not enough space to complete
 * this operation, return 0. Otherwise, return 1.
 *
 * Parameters
 * ==========
 * Amat             On input, points to user's data containing matrix values.
 * N_requested_rows On input, number of rows for which nonzero are to be
 *                  returned.
 * requested_rows   On input, requested_rows[0...N_requested_rows-1] give the
 *                  row indices of the rows for which nonzero values are
 *                  returned.
 * row_lengths      On output, row_lengths[i] is the number of nonzeros in the
 *                  row 'requested_rows[i]'
 * columns,values   On output, columns[k] and values[k] contains the column
 *                  number and value of a matrix nonzero where all nonzeros for
 *                  requested_rows[i] appear before requested_rows[i+1]'s
 *                  nonzeros.  NOTE: Arrays are of size 'allocated_space'.
 * allocated_space  On input, indicates the space available in 'columns' and
 *                  'values' for storing nonzeros. If more space is needed,
 *                  return 0.
 */
  AztecOO::MatrixData * Data = (AztecOO::MatrixData *) AZ_get_matvec_data(Amat);
  Epetra_RowMatrix * A = (Epetra_RowMatrix *) Data->A;
  int NumEntries;

  double *Values = values;
  int * Indices = columns;
  int Length = allocated_space;

  for (int i = 0; i < N_requested_rows; i++) {
    int LocalRow = requested_rows[i];
    // Do copy, return if copy failed
    A->NumMyRowEntries(LocalRow, NumEntries);
    if (NumEntries>Length) return(0);  // This check should eliminate error reports from next line
    if (A->ExtractMyRowCopy(LocalRow, Length, NumEntries, Values, Indices)!=0) return(0);
    // update row lengths and pointers 
    row_lengths[i] = NumEntries;
    Values += NumEntries;
    Indices+= NumEntries;
    Length -= NumEntries;
  }
  return(1);
}
//=============================================================================
int Epetra_Aztec_comm_wrapper(double vec[], AZ_MATRIX *Amat) {
/*
 * Update vec's ghost node via communication. Note: the length of vec is
 * given by N_local + N_ghost where Amat was created via
 *                 AZ_matrix_create(N_local);
 * and a 'getrow' function was supplied via
 *                 AZ_set_MATFREE_getrow(Amat, , , , N_ghost, );
 *
 * Parameters
 * ==========
 * vec              On input, vec contains data. On output, ghost values
 *                  are updated.
 *
 * Amat             On input, points to user's data containing matrix values.
 *                  and communication information.
 */


  AztecOO::MatrixData * Data = (AztecOO::MatrixData *) AZ_get_matvec_data(Amat);
  Epetra_RowMatrix * A = (Epetra_RowMatrix *) Data->A;
  if (A->Comm().NumProc()==1) return(1); // Nothing to do in serial mode
  if (A->RowMatrixImporter()==0) return(1); // Nothing to do if no importer.

	Epetra_Vector * SourceVec = (Epetra_Vector *) Data->SourceVec;
	Epetra_Vector * TargetVec = (Epetra_Vector *) Data->TargetVec;

  if (SourceVec==0) {
		SourceVec = new Epetra_Vector(View, A->RowMatrixImporter()->SourceMap(), vec);
		SourceVec->SetLabel("Epetra_Aztec_comm_wrapper X source");
		Data->SourceVec = SourceVec;
		TargetVec = new Epetra_Vector(View, A->RowMatrixImporter()->TargetMap(), vec);
		TargetVec->SetLabel("Epetra_Aztec_comm_wrapper X target");
		Data->TargetVec = TargetVec;
	}
	else {
		SourceVec->ResetView(vec);
		TargetVec->ResetView(vec);
	}

#ifndef NDEBUG
  int retval =
#endif
  TargetVec->Import(*SourceVec, *(A->RowMatrixImporter()),Insert);
  assert(retval==0);

  return(1);
}
//=============================================================================
void AztecOO_StatusTest_wrapper(void * conv_test_obj,/* pointer to AztecOO_StatusTest object */
				void * res_vector_obj, /* pointer to Epetra_Vector res_vector */
				int iteration,       /* current iteration */
				double * res_vector, /* current natural residual vector */
				int print_info,      /* no info print if 0, else  print info */
				int sol_updated,      /* solution not updated if = 0, else it is
							 and is consistent with res_vector */
				int * converged,     /* = 0 on return if not converged, otherwise converged */
				int * isnan,         /* = 0 on return if not NaN, otherwise NaNs detected */
				double * rnorm,     /* = current norm on return */ 
				int * r_avail)     /* If set to AZ_TRUE on return, the residual 
						       vector is needed
						       by this convergence on subsequent calls 
						       and it should be 
						       supplied by the calling routine */
{


  AztecOO_StatusTest * StatusTest = (AztecOO_StatusTest *) conv_test_obj;

  Epetra_Vector * ResidualVector;
  if (res_vector==0) 
    ResidualVector = 0; // No residual vector passed in, so make the Epetra_Vector 0 too
  else {
    ResidualVector = (Epetra_Vector *) res_vector_obj; // Otherwise cast the pointer passed in
    ResidualVector->ResetView(res_vector);
  }

  bool SolutionUpdated = false;
  if (sol_updated) SolutionUpdated = true;
  AztecOO_StatusType Status = StatusTest->CheckStatus(iteration, ResidualVector, 
								  *rnorm, SolutionUpdated);
  
  if ((Status==Converged && print_info==0) || print_info==1) {
    StatusTest->Print(std::cout);
  }
  if (StatusTest->ResidualVectorRequired())
    *r_avail = AZ_TRUE;
  else
    *r_avail = AZ_FALSE;

  if (Status==Unconverged)
    *converged = 0;
  else if(Status==Converged)
    *converged = 1;
  else if (Status==NaN || Status==PartialFailed)
    *isnan = 1;
  else
    *isnan = 1; // Failed, treat same as isnan

  return;
}

