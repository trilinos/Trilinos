/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_config.h"
#include "ml_include.h"
#if defined(HAVE_ML_IFPACK) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRA)
#include "ml_epetra_utils.h"
#include "Epetra_Map.h" 
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h" 
#include "Epetra_VbrMatrix.h" 
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "ml_ifpack.h"
#include "ml_ifpack_wrap.h"
#include "Ifpack.h"
#include "Ifpack_Chebyshev.h"
#ifdef rst_dump
#include "ml_Ifpack_ML.h"
#endif
// converter from ML_Operator to Epetra_RowMatrix (only wraps)
#include "ml_RowMatrix.h"
// IFPACK factory class
#include "Ifpack.h"
#include "Ifpack_Chebyshev.h"

using namespace ML_Epetra;

namespace ML_Epetra{ 

  Epetra_Operator* ML_Gen_Smoother_Ifpack_Epetra(const Epetra_Operator *A,const Epetra_Vector *InvDiagonal, Teuchos::ParameterList & List,string printMsg,bool verbose){
  /* Variables */
  double lambda_min = 0.0;
  double lambda_max = 0.0;
  Teuchos::ParameterList IFPACKList=List.sublist("smoother: ifpack list");

  /* Parameter-list Options */
  string SmooType = List.get("smoother: type", "Chebyshev");
  if(SmooType=="IFPACK") SmooType=List.get("smoother: ifpack type","Chebyshev");
  int Sweeps = List.get("smoother: sweeps", 3);
  int MaximumIterations = List.get("eigen-analysis: max iters", 10);
  string EigenType_ = List.get("eigen-analysis: type", "cg");
  double boost = List.get("eigen-analysis: boost for lambda max", 1.0);
  double alpha = List.get("chebyshev: alpha",30.0001);  
  double omega = List.get("smoother: damping factor",1.0);
  Ifpack_Chebyshev* SmootherC_=0;
  Ifpack_Preconditioner* SmootherP_=0;

  /* Sanity Check*/
  if(Sweeps==0) return 0;

  /* Early Output*/
  int Nrows=A->OperatorDomainMap().NumGlobalElements();
  int Nnz=-1;
  const Epetra_RowMatrix *Arow=dynamic_cast<const Epetra_RowMatrix*>(A);
  if(Arow) Nnz=Arow->NumGlobalNonzeros();

  if(verbose && !A->Comm().MyPID())
    cout <<printMsg<<" # global rows = "<<Nrows<<" # estim. global nnz = "<<Nnz<<endl;


  /**********************************************/
  /***               Chebyshev                ***/
  /**********************************************/
  if(SmooType=="Chebyshev"){
    /* Grab Diagonal & invert if not provided */
    Epetra_Vector *InvDiagonal_;
    if(InvDiagonal) InvDiagonal_=const_cast<Epetra_Vector *>(InvDiagonal);
    else{
      const Epetra_CrsMatrix* Acrs=dynamic_cast<const Epetra_CrsMatrix*>(A);
      if(!Acrs) return 0;
      InvDiagonal_ = new Epetra_Vector(Acrs->RowMap());  
      Acrs->ExtractDiagonalCopy(*InvDiagonal_);
      for (int i = 0; i < InvDiagonal_->MyLength(); ++i)
	if ((*InvDiagonal_)[i] != 0.0)
	  (*InvDiagonal_)[i] = 1.0 / (*InvDiagonal_)[i];   
    }

    /* Do the eigenvalue estimation*/
    if (EigenType_ == "power-method") Ifpack_Chebyshev::PowerMethod(*A,*InvDiagonal_,MaximumIterations,lambda_max);
    else if(EigenType_ == "cg") Ifpack_Chebyshev::CG(*A,*InvDiagonal_,MaximumIterations,lambda_min,lambda_max);
    else ML_CHK_ERR(0); // not recognized
    
    lambda_min=lambda_max / alpha;
    
    /* Setup the Smoother's List*/
    IFPACKList.set("chebyshev: min eigenvalue", lambda_min);
    IFPACKList.set("chebyshev: max eigenvalue", boost * lambda_max);
    IFPACKList.set("chebyshev: ratio eigenvalue",alpha);
    IFPACKList.set("chebyshev: operator inv diagonal", InvDiagonal_);
    IFPACKList.set("chebyshev: degree", Sweeps);
    IFPACKList.set("chebyshev: zero starting solution",false);

    if(verbose && !A->Comm().MyPID()){
      cout << printMsg << "MLS/Chebyshev, polynomial order = "
	   <<  Sweeps
	   << ", alpha = " << alpha << endl;
      cout << printMsg << "lambda_min = " << lambda_min
             << ", lambda_max = " << boost*lambda_max << endl;
    }


    SmootherC_= new Ifpack_Chebyshev(A);
    if (SmootherC_ == 0) return 0;
    SmootherC_->SetParameters(IFPACKList);
    SmootherC_->Initialize();
    SmootherC_->Compute();
    return SmootherC_;
  }
  /**********************************************/
  /***               Gauss-Seidel             ***/
  /**********************************************/
  else if(SmooType=="Gauss-Seidel" || SmooType=="symmetric Gauss-Seidel" || SmooType=="point relaxation stand-alone"){      
    bool gs_type = List.get("smoother: Gauss-Seidel efficient symmetric",false);
    const Epetra_CrsMatrix* Acrs=dynamic_cast<const Epetra_CrsMatrix*>(A);
    if(!Acrs) return 0;
    string MyIfpackType="point relaxation stand-alone";
    string MyRelaxType="symmetric Gauss-Seidel";
    if(SmooType=="Gauss-Seidel" || SmooType=="symmetric Gauss-Seidel") MyRelaxType=SmooType;
    IFPACKList.set("relaxation: type", IFPACKList.get("relaxation: type",SmooType));
    IFPACKList.set("relaxation: sweeps", Sweeps);
    IFPACKList.set("relaxation: damping factor", omega);

    if(verbose && !A->Comm().MyPID()){
      cout << printMsg << IFPACKList.get("relaxation: type",SmooType).c_str()<<" (sweeps="
	   << Sweeps << ",omega=" << omega << endl;
    }
    	
    //NTS: Finish
    Ifpack Factory;
    SmootherP_ = Factory.Create(MyIfpackType,const_cast<Epetra_CrsMatrix*>(Acrs),0);
    if (SmootherP_ == 0) return 0;
    SmootherP_->SetParameters(IFPACKList);
    SmootherP_->Initialize();
    SmootherP_->Compute();
    return SmootherP_;


  }
  else{
    printf("ML_Gen_Smoother_Ifpack_New: Unknown preconditioner\n");
    ML_CHK_ERR(0);
  }
  return 0;

}/*ML_Gen_Smoother_Ifpack_New*/

}

#endif
