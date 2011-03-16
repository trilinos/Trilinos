#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Stokhos_InterlacedTestSupport.hpp"

#include "Stokhos_Epetra.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_LocalMap.h"

// FEApp is defined in Trilinos/packages/sacado/example/FEApp
// #include "FEApp_ModelEvaluator.hpp"
  
Teuchos::RCP<Teuchos::ParameterList> buildAppParams(int num_KL,bool full_expansion)
{
   Teuchos::RCP<Teuchos::ParameterList> appParams = Teuchos::rcp(new Teuchos::ParameterList);

   Teuchos::ParameterList& problemParams = 
      appParams->sublist("Problem");
   problemParams.set("Name", "Heat Nonlinear Source");

   // Boundary conditions
   problemParams.set("Left BC", 0.0);
   problemParams.set("Right BC", 0.0);
 
   // Source function
   Teuchos::ParameterList& sourceParams = 
     problemParams.sublist("Source Function");
   sourceParams.set("Name", "Constant");
   sourceParams.set("Constant Value", 1.0);

   // Material
   Teuchos::ParameterList& matParams = 
     problemParams.sublist("Material Function");
   matParams.set("Name", "KL Exponential Random Field");
   matParams.set("Mean", 1.0);
   matParams.set("Standard Deviation", 0.5);
   matParams.set("Number of KL Terms", num_KL);
   Teuchos::Array<double> a(1), b(1), L(1);
   a[0] = 0.0; b[0] = 1.0; L[0] = 1.0;
   matParams.set("Domain Lower Bounds", a);
   matParams.set("Domain Upper Bounds", b);
   matParams.set("Correlation Lengths", L);

   // Response functions
   Teuchos::ParameterList& responseParams =
     problemParams.sublist("Response Functions");
   responseParams.set("Number", 1);
   responseParams.set("Response 0", "Solution Average");

    // Setup stochastic Galerkin algorithmic parameters
    Teuchos::RCP<Teuchos::ParameterList> sgParams = 
      Teuchos::rcp(&(appParams->sublist("SG Parameters")),false);
    if (!full_expansion) {
      sgParams->set("Parameter Expansion Type", "Linear");
      sgParams->set("Jacobian Expansion Type", "Linear");
    }
    Teuchos::ParameterList& sgOpParams = 
      sgParams->sublist("SG Operator");
    sgOpParams.set("Operator Method", "Matrix Free");
    Teuchos::ParameterList& sgPrecParams = 
      sgParams->sublist("SG Preconditioner");
    sgPrecParams.set("Preconditioner Method", "Mean-based");
    sgPrecParams.set("Mean Preconditioner Type", "ML");
    Teuchos::ParameterList& precParams = 
      sgPrecParams.sublist("Mean Preconditioner Parameters");
    precParams.set("default values", "SA");

   return appParams;
}

Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > buildBasis(int num_KL,int porder)
{
   // Create Stochastic Galerkin basis and expansion
   Teuchos::Array< Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(num_KL); 
   for(int i=0; i<num_KL; i++)
      bases[i] = Teuchos::rcp(new Stokhos::LegendreBasis<int,double>(porder));

   Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis
         = Teuchos::rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));

   return basis;
}

/*
Teuchos::RCP<EpetraExt::ModelEvaluator> buildFEAppME(int nelem,int num_KL,
                                                     const Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > & basis,
                                                     const Teuchos::RCP<const Epetra_Comm> & app_comm,
                                                     const Teuchos::RCP<Teuchos::ParameterList> & appParams)
{
   Teuchos::RCP<const Epetra_Comm> sg_comm = app_comm; // sg comm and application comm are
                                                       // the same...lost parallelism?!?!

   // build application
   Teuchos::RCP<FEApp::Application> app;
   {
      double h = 1.0/nelem;

      // Create mesh
      std::vector<double> x(nelem+1);
      for(int i=0; i<=nelem; i++)
         x[i] = h*i;

      app = Teuchos::rcp(new FEApp::Application(x, app_comm, appParams, false));
   }

   // Free parameters (determinisic, e.g., for sensitivities)
   Teuchos::RefCountPtr< Teuchos::Array<std::string> > free_param_names =
 	Teuchos::rcp(new Teuchos::Array<std::string>);
   free_param_names->push_back("Constant Source Function Value");

   // setup some named strings
   Teuchos::RefCountPtr< Teuchos::Array<std::string> > sg_param_names
         = Teuchos::rcp(new Teuchos::Array<std::string>);
   for(int i=0; i<num_KL; i++) {
      std::stringstream ss;
      ss << "KL Exponential Function Random Variable " << i;
      sg_param_names->push_back(ss.str());
   }

   Teuchos::RCP<EpetraExt::ModelEvaluator> model = 
        Teuchos::rcp(new FEApp::ModelEvaluator(app, free_param_names,
                                               sg_param_names));

   return model;
}
*/

Teuchos::RCP<Stokhos::ParallelData> buildParallelData(bool full_expansion,int num_KL,
                                                      const Teuchos::RCP<const Epetra_Comm> & globalComm,
                                                      const Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > & basis)
{
    Teuchos::ParameterList parallelParams;
    parallelParams.set("Number of Spatial Processors", globalComm->NumProc());

    int sz = basis->size();
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
    if (full_expansion)
      Cijk = basis->computeTripleProductTensor(sz);
    else
      Cijk = basis->computeTripleProductTensor(num_KL+1);

   Teuchos::RCP<Stokhos::ParallelData> sg_parallel_data 
      = Teuchos::rcp(new Stokhos::ParallelData(basis, Cijk, globalComm,
                                                parallelParams));

   return sg_parallel_data;
}
