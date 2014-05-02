/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_config.h"
#include "ml_include.h"
#if defined(HAVE_ML_IFPACK) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRA)
#include "Ifpack_config.h"
#include "ml_epetra_utils.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Time.h"
#include "ml_ifpack.h"
#include "ml_ifpack_wrap.h"
#include "Ifpack_Chebyshev.h"
#ifdef rst_dump
#include "ml_Ifpack_ML.h"
#endif
// converter from ML_Operator to Epetra_RowMatrix (only wraps)
#include "ml_RowMatrix.h"
// IFPACK factory class
#ifdef HAVE_IFPACK_DYNAMIC_FACTORY
#include "Ifpack_DynamicFactory.h"
#else
#include "Ifpack.h"
#endif
#include "Ifpack_Chebyshev.h"
#include "Ifpack_BlockRelaxation.h"

using namespace ML_Epetra;

namespace ML_Epetra{

  Epetra_Operator* ML_Gen_Smoother_Ifpack_Epetra(const Epetra_Operator *A,const Epetra_Vector *InvDiagonal, Teuchos::ParameterList & List,std::string printMsg,bool verbose){
  /* Variables */
  double lambda_min = 0.0;
  double lambda_max = 0.0;
  Teuchos::ParameterList IFPACKList=List.sublist("smoother: ifpack list");

  /* Parameter-list Options */
  std::string PreOrPostSmoother = List.get("smoother: pre or post","both");
  std::string SmooType = List.get("smoother: type", "Chebyshev");
  if(SmooType=="IFPACK") SmooType=List.get("smoother: ifpack type","Chebyshev");
  int Sweeps = List.get("smoother: sweeps", 3);
  int IfpackOverlap = List.get("smoother: ifpack overlap",0);
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
    std::cout <<printMsg<<"# global rows = "<<Nrows<<" # estim. global nnz = "<<Nnz<<std::endl;


  /**********************************************/
  /***      Chebyshev (Including Block)       ***/
  /**********************************************/
  if(SmooType=="Chebyshev" || SmooType=="MLS" || SmooType=="IFPACK-Chebyshev" || SmooType=="IFPACK-Block Chebyshev"){
    bool allocated_inv_diagonal=false;
    int MaximumIterations = List.get("eigen-analysis: max iters", 10);
    std::string EigenType_ = List.get("eigen-analysis: type", "cg");
    double boost = List.get("eigen-analysis: boost for lambda max", 1.0);
    double alpha = List.get("chebyshev: alpha",30.0001);
    Epetra_Vector *InvDiagonal_=0;

    /* Block Chebyshev stuff if needed */
    int  MyCheby_nBlocks=     List.get("smoother: Block Chebyshev number of blocks",0);
    int* MyCheby_blockIndices=List.get("smoother: Block Chebyshev block list",(int*)0);
    int* MyCheby_blockStarts= List.get("smoother: Block Chebyshev block starts",(int*)0);
    bool MyCheby_NE=          List.get("smoother: chebyshev solve normal equations",false);

    if(SmooType == "IFPACK-Block Chebyshev" && MyCheby_blockIndices && MyCheby_blockStarts){
      // If we're using Block Chebyshev, it can compute it's own eigenvalue estimate..
      Teuchos::ParameterList PermuteList,BlockList;
      BlockList.set("apply mode","invert");
      PermuteList.set("number of local blocks",MyCheby_nBlocks);
      PermuteList.set("block start index",MyCheby_blockStarts);
      //        if(is_lid) PermuteList.set("block entry lids",Blockids_);
      //NTS: Add LID support
      PermuteList.set("block entry gids",MyCheby_blockIndices);
      PermuteList.set("blockdiagmatrix: list",BlockList);

      IFPACKList.set("chebyshev: use block mode",true);
      IFPACKList.set("chebyshev: block list",PermuteList);
      IFPACKList.set("chebyshev: eigenvalue max iterations",10);

      // EXPERIMENTAL: Cheby-NE
      IFPACKList.set("chebyshev: solve normal equations",MyCheby_NE);
    }
    else {
      /* Non-Blocked Chebyshev */
      /* Grab Diagonal & invert if not provided */
      if(InvDiagonal) InvDiagonal_=const_cast<Epetra_Vector *>(InvDiagonal);
      else{
	const Epetra_CrsMatrix* Acrs=dynamic_cast<const Epetra_CrsMatrix*>(A);
	if(!Acrs) return 0;
	allocated_inv_diagonal=true;
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
    }

    /* Setup the Smoother's List*/
    IFPACKList.set("chebyshev: ratio eigenvalue", alpha);
    IFPACKList.set("chebyshev: degree", Sweeps);
    IFPACKList.set("chebyshev: zero starting solution",false);

    // Setup
    SmootherC_= new Ifpack_Chebyshev(A);
    if (SmootherC_ == 0) return 0;
    SmootherC_->SetParameters(IFPACKList);
    SmootherC_->Initialize();
    SmootherC_->Compute();

    // Grab the lambda's if needed
    if(SmooType=="IFPACK-Block Chebyshev"){
      lambda_min=SmootherC_->GetLambdaMin();
      lambda_max=SmootherC_->GetLambdaMax();
    }

    // Smoother Info Output
    if(verbose && !A->Comm().MyPID()){
      if(SmooType=="IFPACK-Block Chebyshev") {
	std::cout << printMsg << "MLS/Block-Chebyshev, polynomial order = "
	     <<  Sweeps
	     << ", alpha = " << alpha << std::endl;
	std::cout << printMsg << "lambda_min = " << lambda_min
	     << ", lambda_max = " << lambda_max << std::endl;
      }
      else {
	std::cout << printMsg << "MLS/Chebyshev, polynomial order = "
	     <<  Sweeps
	     << ", alpha = " << alpha << std::endl;
	std::cout << printMsg << "lambda_min = " << lambda_min
             << ", lambda_max = " << boost*lambda_max << std::endl;
      }
    }


    // Cleanup:  Since Chebyshev will keep it's own copy of the Inverse Diagonal...
    if (allocated_inv_diagonal) delete InvDiagonal_;

    return SmootherC_;
  }
  /**********************************************/
  /***           Point Relaxation             ***/
  /**********************************************/
  else if(SmooType=="Gauss-Seidel" || SmooType=="symmetric Gauss-Seidel" || SmooType=="Jacobi"
	  || SmooType=="point relaxation stand-alone" || SmooType=="point relaxation" ){
    const Epetra_CrsMatrix* Acrs=dynamic_cast<const Epetra_CrsMatrix*>(A);
    if(!Acrs) return 0;
    std::string MyIfpackType="point relaxation stand-alone";
    if(IfpackOverlap > 0) MyIfpackType="point relaxation";
    std::string MyRelaxType="symmetric Gauss-Seidel";
    if(SmooType=="symmetric Gauss-Seidel") MyRelaxType=SmooType;
    else if(SmooType=="Jacobi") MyRelaxType=SmooType;

    IFPACKList.set("relaxation: type", IFPACKList.get("relaxation: type",MyRelaxType));
    IFPACKList.set("relaxation: sweeps", Sweeps);
    IFPACKList.set("relaxation: damping factor", omega);
    IFPACKList.set("relaxation: zero starting solution",false);

    if(verbose && !A->Comm().MyPID()){
      std::cout << printMsg << IFPACKList.get("relaxation: type",MyRelaxType).c_str()<<" (sweeps="
	   << Sweeps << ",omega=" << omega <<  ")" <<std::endl;
    }

#ifdef HAVE_IFPACK_DYNAMIC_FACTORY
    Ifpack_DynamicFactory Factory;
#else
    Ifpack Factory;
#endif
    SmootherP_ = Factory.Create(MyIfpackType,const_cast<Epetra_CrsMatrix*>(Acrs),IfpackOverlap);
    if (SmootherP_ == 0) return 0;
    SmootherP_->SetParameters(IFPACKList);
    SmootherP_->Initialize();
    SmootherP_->Compute();
    return SmootherP_;
  }
  /**********************************************/
  /***           Block Relaxation             ***/
  /**********************************************/
  else if(SmooType=="block Gauss-Seidel" || SmooType=="symmetric block Gauss-Seidel" || SmooType=="block Jacobi"
	  || SmooType=="block relaxation stand-alone" || SmooType=="block relaxation") {
    const Epetra_CrsMatrix* Acrs=dynamic_cast<const Epetra_CrsMatrix*>(A);
    if(!Acrs) return 0;
    std::string MyIfpackType="block relaxation stand-alone";
    if(IfpackOverlap > 0) MyIfpackType="block relaxation";
    std::string MyRelaxType="symmetric Gauss-Seidel";
    if(SmooType=="block Gauss-Seidel") MyRelaxType="Gauss-Seidel";
    else if(SmooType=="block Jacobi") MyRelaxType="Jacobi";

    IFPACKList.set("relaxation: type", IFPACKList.get("relaxation: type",MyRelaxType));
    IFPACKList.set("relaxation: sweeps", Sweeps);
    IFPACKList.set("relaxation: damping factor", omega);
    IFPACKList.set("relaxation: zero starting solution",false);
    bool use_line=false;
    if(List.isParameter("smoother: line detection threshold")) {
      use_line = true;
      IFPACKList.set("partitioner: line detection threshold",List.get("smoother: line detection threshold",-1.0));
      IFPACKList.set("partitioner: type","line");
      IFPACKList.set("partitioner: x-coordinates",List.get("x-coordinates",(double*)0));
      IFPACKList.set("partitioner: y-coordinates",List.get("y-coordinates",(double*)0));
      IFPACKList.set("partitioner: z-coordinates",List.get("z-coordinates",(double*)0));

    }

    if(verbose && !A->Comm().MyPID()){
      std::cout << printMsg << "block " << IFPACKList.get("relaxation: type",MyRelaxType).c_str()<<" (sweeps="
		<< Sweeps << ",omega=" << omega;
      if(use_line) std::cout << ",auto-line";	
    }
    
#ifdef HAVE_IFPACK_DYNAMIC_FACTORY
    Ifpack_DynamicFactory Factory;
#else
    Ifpack Factory;
#endif
    SmootherP_ = Factory.Create(MyIfpackType,const_cast<Epetra_CrsMatrix*>(Acrs),IfpackOverlap);
    if (SmootherP_ == 0) return 0;
    SmootherP_->SetParameters(IFPACKList);
    SmootherP_->Initialize();
    SmootherP_->Compute();

    int global_blocks=0,local_blocks=0;
    Ifpack_BlockRelaxation<Ifpack_DenseContainer> * temp0 = dynamic_cast<Ifpack_BlockRelaxation<Ifpack_DenseContainer>*>(SmootherP_);
    if(temp0) local_blocks = temp0->NumLocalBlocks();  
    A->Comm().SumAll(&local_blocks,&global_blocks,1);
    if(verbose && !A->Comm().MyPID())
      std::cout<<","<<global_blocks<<" blocks found)"<<std::endl;


    return SmootherP_;
  }
  /**********************************************/
  /***        Incomplete Factorization        ***/
  /**********************************************/
  else if(SmooType == "ILU" || SmooType == "IC" || SmooType == "ILUT"   ||
	  SmooType == "ICT" || SmooType == "SILU") {
    const Epetra_RowMatrix* Arow=dynamic_cast<const Epetra_RowMatrix*>(A);
    double MyLOF=0.0;
    if(SmooType=="ILUT" || SmooType=="ICT") MyLOF=List.get("smoother: ifpack level-of-fill",1.0);
    else MyLOF=List.get("smoother: ifpack level-of-fill",0.0);

    int MyIfpackOverlap = List.get("smoother: ifpack overlap", 0);
    double MyIfpackRT = List.get("smoother: ifpack relative threshold", 1.0);
    double MyIfpackAT = List.get("smoother: ifpack absolute threshold", 0.0);
    IFPACKList.set("ILU: sweeps",Sweeps);

    // Set the fact: LOF options, but only if they're not set already... All this sorcery is because level-of-fill
    // is an int for ILU and a double for ILUT.  Lovely.
    if(SmooType=="ILUT" || SmooType=="ICT"){
	IFPACKList.set("fact: level-of-fill", IFPACKList.get("fact: level-of-fill",MyLOF));
	IFPACKList.set("fact: ilut level-of-fill", IFPACKList.get("fact: ilut level-of-fill",MyLOF));
	IFPACKList.set("fact: ict level-of-fill", IFPACKList.get("fact: ict level-of-fill",MyLOF));
	MyLOF=IFPACKList.get("fact: level-of-fill",MyLOF);
    }
    else{
      IFPACKList.set("fact: level-of-fill", (int) IFPACKList.get("fact: level-of-fill",(int)MyLOF));
      MyLOF=IFPACKList.get("fact: level-of-fill",(int)MyLOF);
    }

    IFPACKList.set("fact: relative threshold", MyIfpackRT);
    IFPACKList.set("fact: absolute threshold", MyIfpackAT);

    if(verbose && !A->Comm().MyPID()){
      std::cout << printMsg << "IFPACK, type=`" << SmooType << "'," << std::endl
	   << printMsg << PreOrPostSmoother  << ",overlap=" << MyIfpackOverlap << std::endl;
      std::cout << printMsg << "level-of-fill=" << MyLOF;
      std::cout << ",rel. threshold=" << MyIfpackRT
	   << ",abs. threshold=" << MyIfpackAT << std::endl;
    }

#ifdef HAVE_IFPACK_DYNAMIC_FACTORY
    Ifpack_DynamicFactory Factory;
#else
    Ifpack Factory;
#endif
    SmootherP_ = Factory.Create(SmooType,const_cast<Epetra_RowMatrix*>(Arow),IfpackOverlap);
    if (SmootherP_ == 0) return 0;
    SmootherP_->SetParameters(IFPACKList);
    SmootherP_->Initialize();
    SmootherP_->Compute();
    return SmootherP_;
  }
  /**********************************************/
  /***                  SORa                  ***/
  /**********************************************/
  else if(SmooType=="SORa"){
    const Epetra_RowMatrix* Arow=dynamic_cast<const Epetra_RowMatrix*>(A);
    if(verbose && !A->Comm().MyPID()){
      std::cout << printMsg << "IFPACK/SORa("<<IFPACKList.get("sora: alpha",1.5)<<","<<IFPACKList.get("sora: gamma",1.0)<<")"
	   << ", sweeps = " <<IFPACKList.get("sora: sweeps",1)<<std::endl;
      if(IFPACKList.get("sora: oaz boundaries",false))
	std::cout << printMsg << "oaz boundary handling enabled"<<std::endl;
      if(IFPACKList.get("sora: use interproc damping",false))
	std::cout << printMsg << "interproc damping enabled"<<std::endl;
      if(IFPACKList.get("sora: use global damping",false))
	std::cout << printMsg << "global damping enabled"<<std::endl;
    }
#ifdef HAVE_IFPACK_DYNAMIC_FACTORY
    Ifpack_DynamicFactory Factory;
#else
    Ifpack Factory;
#endif
    SmootherP_ = Factory.Create(SmooType,const_cast<Epetra_RowMatrix*>(Arow),IfpackOverlap);
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
