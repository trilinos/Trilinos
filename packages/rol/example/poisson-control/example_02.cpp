
// Burgers includes
#include "example_02.hpp"
// ROL includes
#include "ROL_Algorithm.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_StdTeuchosBatchManager.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_Reduced_ParametrizedObjective_SimOpt.hpp"
#include "ROL_RiskNeutralObjective.hpp"
#include "ROL_Vector_SimOpt.hpp"
// Teuchos includes
#include "Teuchos_Time.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"

int main( int argc, char *argv[] ) {  

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  /***************************************************************************/
  /***************** GRAB INPUTS *********************************************/
  /***************************************************************************/
  // Get finite element parameter list
  std::string filename = "example_02.xml";
  Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
  Teuchos::updateParametersFromXmlFile( filename, Teuchos::Ptr<Teuchos::ParameterList>(&*parlist) );
  if ( parlist->get("Display Option",0) && (comm->getRank() > 0) ) {
    parlist->set("Display Option",0);
  }
  // Get ROL parameter list
  filename = "input.xml";
  Teuchos::RCP<Teuchos::ParameterList> ROL_parlist = Teuchos::rcp( new Teuchos::ParameterList() );
  Teuchos::updateParametersFromXmlFile( filename, Teuchos::Ptr<Teuchos::ParameterList>(&*ROL_parlist) );

  /***************************************************************************/
  /***************** INITIALIZE SAMPLERS *************************************/
  /***************************************************************************/
  int dim    = 2;
  bool useSA = parlist->get("Use Stochastic Approximation",false);
  int nSamp  = 1;
  if ( !useSA ) {
    nSamp  = parlist->get("Number of Monte Carlo Samples",1000);
  }
  std::vector<double> tmp(2); tmp[0] = -1.0; tmp[1] = 1.0;
  std::vector<std::vector<double> > bounds(dim,tmp);
  Teuchos::RCP<ROL::BatchManager<double> > bman
    = Teuchos::rcp(new ROL::StdTeuchosBatchManager<double,int>(comm));
  Teuchos::RCP<ROL::SampleGenerator<double> > sampler
    = Teuchos::rcp(new ROL::MonteCarloGenerator<double>(nSamp,bounds,bman,useSA));

  /***************************************************************************/
  /***************** INITIALIZE CONTROL VECTOR *******************************/
  /***************************************************************************/
  int nx = parlist->get("Number of Elements", 128);
  Teuchos::RCP<std::vector<double> > z_rcp = Teuchos::rcp( new std::vector<double>(nx+1, 0.0) );
  Teuchos::RCP<ROL::Vector<double> > z = Teuchos::rcp(new ROL::StdVector<double>(z_rcp));
  Teuchos::RCP<std::vector<double> > u_rcp = Teuchos::rcp( new std::vector<double>(nx-1, 0.0) );
  Teuchos::RCP<ROL::Vector<double> > u = Teuchos::rcp(new ROL::StdVector<double>(u_rcp));
  ROL::Vector_SimOpt<double> x(u,z);
  Teuchos::RCP<std::vector<double> > p_rcp = Teuchos::rcp( new std::vector<double>(nx-1, 0.0) );
  Teuchos::RCP<ROL::Vector<double> > p = Teuchos::rcp(new ROL::StdVector<double>(p_rcp));
  Teuchos::RCP<std::vector<double> > U_rcp = Teuchos::rcp( new std::vector<double>(nx+1, 35.0) );
  Teuchos::RCP<ROL::Vector<double> > U = Teuchos::rcp(new ROL::StdVector<double>(U_rcp));
  Teuchos::RCP<std::vector<double> > L_rcp = Teuchos::rcp( new std::vector<double>(nx+1, -5.0) );
  Teuchos::RCP<ROL::Vector<double> > L = Teuchos::rcp(new ROL::StdVector<double>(L_rcp));
  ROL::BoundConstraint<double> bnd(L,U);

  /***************************************************************************/
  /***************** INITIALIZE OBJECTIVE FUNCTION ***************************/
  /***************************************************************************/
  double alpha = parlist->get("Penalty Parameter", 1.e-4);
  Teuchos::RCP<FEM<double> > fem = Teuchos::rcp(new FEM<double>(nx));
  Teuchos::RCP<ROL::ParametrizedObjective_SimOpt<double> > pObj
    = Teuchos::rcp(new DiffusionObjective<double>(fem, alpha));
  Teuchos::RCP<ROL::ParametrizedEqualityConstraint_SimOpt<double> > pCon
    = Teuchos::rcp(new DiffusionEqualityConstraint<double>(fem));
  Teuchos::RCP<ROL::ParametrizedObjective<double> > robj
    = Teuchos::rcp(new ROL::Reduced_ParametrizedObjective_SimOpt<double>(pObj,pCon,u,p));
  ROL::RiskNeutralObjective<double> obj(robj,sampler);

  /***************************************************************************/
  /***************** RUN DERIVATIVE CHECK ************************************/
  /***************************************************************************/
  if (parlist->get("Run Derivative Check",false)) {
    // Direction to test finite differences
    Teuchos::RCP<std::vector<double> > dz_rcp = Teuchos::rcp( new std::vector<double>(nx+1, 0.0) );
    Teuchos::RCP<ROL::Vector<double> > dz = Teuchos::rcp(new ROL::StdVector<double>(dz_rcp));
    Teuchos::RCP<std::vector<double> > du_rcp = Teuchos::rcp( new std::vector<double>(nx-1, 0.0) );
    Teuchos::RCP<ROL::Vector<double> > du = Teuchos::rcp(new ROL::StdVector<double>(du_rcp));
    ROL::Vector_SimOpt<double> d(du,dz);
    // Set to random vectors
    srand(12345);
    for (int i=0; i<nx+1; i++) {
      (*dz_rcp)[i] = 2.0*(double)rand()/(double)RAND_MAX - 1.0;
      (*z_rcp)[i] = 2.0*(double)rand()/(double)RAND_MAX - 1.0;
    }
    for (int i=0; i<nx-1; i++) {
      (*du_rcp)[i] = 2.0*(double)rand()/(double)RAND_MAX - 1.0;
      (*u_rcp)[i] = 2.0*(double)rand()/(double)RAND_MAX - 1.0;
    }
    // Run derivative checks
    std::vector<double> param(dim,0.0);
    robj->setParameter(param);
    if ( comm->getRank() == 0 ) {
      std::cout << "\nRUN DERIVATIVE CHECK FOR PARAMETRIZED OBJECTIVE FUNCTION SIMOPT\n";
    }
    pObj->checkGradient(x,d,(comm->getRank()==0));
    pObj->checkHessVec(x,d,(comm->getRank()==0));
    if ( comm->getRank() == 0 ) {
      std::cout << "\nRUN DERIVATIVE CHECK FOR PARAMETRIZED EQUALITY CONSTRAINT SIMOPT\n";
    }
    pCon->checkApplyJacobian(x,d,*p,(comm->getRank()==0));
    pCon->checkApplyAdjointJacobian(x,*du,*p,x,(comm->getRank()==0));
    pCon->checkApplyAdjointHessian(x,*du,d,x,(comm->getRank()==0));
    if ( comm->getRank() == 0 ) {
      std::cout << "\nRUN DERIVATIVE CHECK FOR PARAMETRIZED OBJECTIVE FUNCTION\n";
    }
    robj->checkGradient(*z,*dz,(comm->getRank()==0));
    robj->checkHessVec(*z,*dz,(comm->getRank()==0));
    // Run derivative checks
    if ( comm->getRank() == 0 ) {
      std::cout << "\nRUN DERIVATIVE CHECK FOR RISK-NEUTRAL OBJECTIVE FUNCTION\n";
    }
    obj.checkGradient(*z,*dz,(comm->getRank()==0));
    obj.checkHessVec(*z,*dz,(comm->getRank()==0));
  }

  /***************************************************************************/
  /***************** INITIALIZE ROL ALGORITHM ********************************/
  /***************************************************************************/
  Teuchos::RCP<ROL::Algorithm<double> > algo; 
  if ( useSA ) {
    ROL_parlist->sublist("General").set("Recompute Objective Function",false);
    ROL_parlist->sublist("Step").sublist("Line Search").set("Initial Step Size",0.1/alpha);
    ROL_parlist->sublist("Step").sublist("Line Search").set("User Defined Initial Step Size",true);
    ROL_parlist->sublist("Step").sublist("Line Search").sublist("Line-Search Method").set("Type","Iteration Scaling");
    ROL_parlist->sublist("Step").sublist("Line Search").sublist("Descent Method").set("Type","Steepest Descent");
    ROL_parlist->sublist("Step").sublist("Line Search").sublist("Curvature Condition").set("Type","Null Curvature Condition");
    algo = Teuchos::rcp(new ROL::Algorithm<double>("Line Search",*ROL_parlist,false));
  } 
  else {
    algo = Teuchos::rcp(new ROL::Algorithm<double>("Trust Region",*ROL_parlist,false));
  }

  /***************************************************************************/
  /***************** PERFORM OPTIMIZATION ************************************/
  /***************************************************************************/
  Teuchos::Time timer("Optimization Time",true);
  z->zero();
  algo->run(*z,obj,bnd,(comm->getRank()==0));
  double optTime = timer.stop();

  /***************************************************************************/
  /***************** PRINT RESULTS *******************************************/
  /***************************************************************************/
  int my_number_samples = sampler->numMySamples(), number_samples = 0;
  Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_SUM,1,&my_number_samples,&number_samples);
  int my_number_solves  = Teuchos::rcp_dynamic_cast<DiffusionEqualityConstraint<double> >(pCon)->getNumSolves(), number_solves = 0;
  Teuchos::reduceAll<int,int>(*comm,Teuchos::REDUCE_SUM,1,&my_number_solves,&number_solves);
  if (comm->getRank() == 0) {
    std::cout << "Number of Samples    = " << number_samples << "\n";
    std::cout << "Number of Solves     = " << number_solves  << "\n";
    std::cout << "Optimization Time    = " << optTime        << "\n\n";
  }

  if ( comm->getRank() == 0 ) {
    std::ofstream file;
    if (useSA) {
      file.open("control_SA.txt");
    }
    else {
      file.open("control_SAA.txt");
    }
    std::vector<double> xmesh(fem->nz(),0.0);
    fem->build_mesh(xmesh);
    for (int i = 0; i < fem->nz(); i++ ) {
      file << std::setprecision(std::numeric_limits<double>::digits10) << std::scientific << xmesh[i] << "  "  
           << std::setprecision(std::numeric_limits<double>::digits10) << std::scientific << (*z_rcp)[i] 
           << "\n";
    }
    file.close();
  }

  return 0;
}




