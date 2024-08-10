// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the optimal control of Helmholtz problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_Time.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Stream.hpp"
#include "ROL_Solver.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_ScalarLinearConstraint.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_StdTeuchosBatchManager.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_PrimalDualRisk.hpp"

#include "../../../../TOOLS/linearpdeconstraint.hpp"
#include "../../../../TOOLS/pdeobjective.hpp"
#include "../../../../TOOLS/pdevector.hpp"
#include "../../../../TOOLS/meshreader.hpp"
#include "../../../../TOOLS/batchmanager.hpp"

#include "ROL_OED_Factory.hpp"
#include "ROL_OED_StdMomentOperator.hpp"
#include "../../../utilities/OED_SplitComm.hpp"

#include "../src/pde_poisson.hpp"
#include "../src/obj_poisson.hpp"
#include "../src/noise.hpp"
#include "../src/mesh.hpp"

template<typename Real>
void solve(const ROL::Ptr<ROL::OED::Factory<Real>>    &factory,
           const ROL::Ptr<ROL::SampleGenerator<Real>> &dsampler,
           const ROL::Ptr<ROL::SampleGenerator<Real>> &osampler,
           ROL::ParameterList                         &list,
           const std::string                          &dtype,
           std::ostream                               &stream     = std::cout,
           const bool                                  useUIG     = true,
           const bool                                  checkDeriv = false) {
  int nsampd = dsampler->numGlobalSamples();
  if (useUIG) {
    factory->setDesign(1.0/static_cast<Real>(nsampd));
  }
  else {
    std::stringstream dfile;
    dfile << dtype << "_optimal_design.txt";
    int err = factory->loadDesign(dfile.str(),2,nsampd);
    stream << "  Factory::loadDesign exited with error " << err << std::endl;
  }
  stream << "  Solve " << dtype << "-optimal design problem." << std::endl;
  std::clock_t timer = std::clock();
  ROL::Ptr<ROL::Problem<Real>> problem = factory->get(list,osampler);
  problem->setProjectionAlgorithm(list);
  problem->finalize(false,true,stream);
  std::string type = list.sublist("OED").get("Optimality Type","I");
  bool usePD = list.sublist("OED").sublist("R-Optimality").get("Use Primal-Dual Algorithm",false);
  if (type=="R" && usePD) {
    ROL::PrimalDualRisk<Real> solver(problem,osampler,list);
    // Commented out because check uses random vectors that do not respect bounds
    //if (checkDeriv) {
    //  factory->check(stream);
    //  //solver.check(stream);
    //  factory->reset();
    //}
    solver.run(stream);
  }
  else {
    ROL::Solver<Real> solver(problem,list);
    // Commented out because check uses random vectors that do not respect bounds
    //if (checkDeriv) {
    //  factory->check(stream);
    //  problem->check(true,stream);
    //  factory->reset();
    //}
    solver.solve(stream);
  }
  stream << "  " << dtype << "-optimal design time:      "
    << static_cast<Real>(std::clock()-timer)/static_cast<Real>(CLOCKS_PER_SEC)
    << " seconds" << std::endl;
  factory->profile(stream);
  factory->reset();
  std::stringstream dname, pname;
  dname << dtype << "_optimal_design";
  factory->printDesign(dname.str());
  pname << dtype << "_optimal_prediction_variance_samples";
  factory->printPredictionVariance(osampler,pname.str());
}

int main(int argc, char *argv[]) {
//  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  using RealT = double;

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  ROL::Ptr<const Teuchos::Comm<int>> dcomm, ocomm;
#ifdef HAVE_MPI
  OED_SplitComm<int,Teuchos::MpiComm<int>> splitcomm(1);
  dcomm = splitcomm.getDesignComm();
  ocomm = splitcomm.getSampleComm();
#else
  dcomm = ROL::makePtr<Teuchos::SerialComm<int>>();
  ocomm = Tpetra::getDefaultComm();
#endif
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  const int drank = dcomm->getRank();
  const int orank = ocomm->getRank();
  bool iprint = (argc > 1) && (drank == 0) && (orank == 0);
  iprint = ((argc > 2) ? true : iprint);
  ROL::Ptr<std::ostream> outStream
    = ROL::makeStreamPtr( std::cout, iprint );

  int errorFlag  = 0;

  // *** Example body.
  try {
    RealT tol(1e-8);// one(1);

    /*************************************************************************/
    /******* BUILD LINEAR REGRESSION MODEL BASED ON HELMHOLTZ ****************/
    /*************************************************************************/
    /*** Read in XML input ***/
    ROL::Ptr<ROL::ParameterList> parlist = ROL::getParametersFromXmlFile("input.xml");
    int      numSides = 4;
    bool   checkDeriv = parlist->sublist("Problem").get("Check Derivatives",  false);
    int     verbosity = parlist->sublist("General").get("Print Verbosity",0);
    bool      alPrint = parlist->sublist("Step").sublist("Augmented Lagrangian").get("Print Intermediate Optimization History",false);
    verbosity         = (iprint ? verbosity : 0);
    alPrint           = (iprint ? alPrint : false);
    parlist->sublist("General").set("Print Verbosity",verbosity);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Print Intermediate Optimization History",alPrint);

    // Initialize PDE describing Poisson equation.
    ROL::Ptr<MeshManager<RealT>>
      meshMgr = ROL::makePtr<MeshManager_OED_Poisson<RealT>>(*parlist);
    ROL::Ptr<PDE_OED_Poisson<RealT>>
      pde = ROL::makePtr<PDE_OED_Poisson<RealT>>(*parlist);
    ROL::Ptr<Linear_PDE_Constraint<RealT>>
      con = ROL::makePtr<Linear_PDE_Constraint<RealT>>(pde,meshMgr,dcomm,
                                                      *parlist,*outStream);
    ROL::Ptr<Assembler<RealT>> assembler = con->getAssembler();
    con->printMeshData(*outStream);

    // Create state vector.
    ROL::Ptr<Tpetra::MultiVector<>> u_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> r_ptr = assembler->createResidualVector();
    ROL::Ptr<ROL::Vector<RealT>> up, rp, zp;
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<ROL::StdVector<RealT>>(numSides);

    // Create observation function
    ROL::Ptr<QoI<RealT>>
      qoi = ROL::makePtr<QoI_Poisson_Observation<RealT>>(pde->getFE(),*parlist);
    ROL::Ptr<ROL::Objective_SimOpt<RealT>>
      obs = ROL::makePtr<PDE_Objective<RealT>>(qoi,assembler);

    // Run derivative checks
    if ( checkDeriv ) {
      ROL::Ptr<ROL::Vector<RealT>> dup, dzp;
      dup = up->clone(); dzp = zp->clone();
      up->randomize();  zp->randomize();
      dup->randomize(); dzp->randomize();
      ROL::Vector_SimOpt<RealT> uz(up,zp), duz(dup,dzp);
      std::vector<RealT> param(2,0);
      obs->setParameter(param);
      *outStream << "\n\nCheck Jacobian_1 of PDE Constraint\n";
      con->checkApplyJacobian_1(*up,*zp,*dup,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_2 of PDE Constraint\n";
      con->checkApplyJacobian_2(*up,*zp,*dzp,*rp,true,*outStream);
      *outStream << "\n\nCheck Gradient of Full Objective Function\n";
      obs->checkGradient(uz,duz,true,*outStream);
      *outStream << "\n\nCheck Hessian of Full Objective Function\n";
      obs->checkHessVec(uz,duz,true,*outStream);
    }

    // Compute states for factor decompostion and build regression model
    *outStream << std::endl << "Begin: PDE Solves" << std::endl;
    std::vector<ROL::Ptr<ROL::Vector<RealT>>> state(numSides);
    for (int i = 0; i < numSides; ++i) {
      std::clock_t time_r = std::clock();
      con->solve(*rp,*up,*zp->basis(i),tol);
      state[i] = up->clone(); state[i]->set(*up);
      std::stringstream real;
      real << "state_" << i << ".txt";
      con->outputTpetraVector(u_ptr,real.str());
      *outStream
        << "  Solve " << i << " time:      "
        << static_cast<RealT>(std::clock()-time_r)/static_cast<RealT>(CLOCKS_PER_SEC)
        << std::endl;
    }
    *outStream << "End: PDE Solves" << std::endl << std::endl;

    // Build nonlinear model and parameter vector
    ROL::Ptr<ROL::Vector<RealT>>
      theta = ROL::makePtr<ROL::StdVector<RealT>>(numSides,1);
    ROL::Ptr<ROL::Objective<RealT>>
      model = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obs,con,up,theta,up,false);

    /*************************************************************************/
    /******* BUILD EXPERIMENTAL DESIGN PROBLEM AND SOLVE *********************/
    /*************************************************************************/
    // Build samplers for experiment space
    int nsampd     = parlist->sublist("Problem").sublist("Design").get("Number of Samples",    5000);
    int nsampo     = parlist->sublist("Problem").sublist("Objective").get("Number of Samples", 200);
    std::vector<std::vector<RealT>> bounds(2);
    bounds[0] = {RealT(0), RealT(1)};
    bounds[1] = {RealT(0), RealT(1)};
    ROL::Ptr<ROL::BatchManager<RealT>> dbman, obman;
    dbman  = ROL::makePtr<ROL::StdTeuchosBatchManager<RealT,int>>(dcomm);
    obman  = ROL::makePtr<ROL::StdTeuchosBatchManager<RealT,int>>(ocomm);
    ROL::Ptr<ROL::SampleGenerator<RealT>> dsampler, osampler;
    dsampler  = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsampd,bounds,dbman,false,false,0,1);
    osampler  = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsampo,bounds,obman,false,false,0,nsampd+1);

    // Build OED problem factory.
    bool useUIG, homNoise, useBudget;
    std::string regType, dtype;
    ROL::Ptr<ROL::OED::Factory<RealT>>           factory;
    ROL::Ptr<ROL::OED::Noise<RealT>>             noise;
    ROL::Ptr<ROL::OED::StdMomentOperator<RealT>> cov;
    homNoise  = parlist->sublist("Problem").get("Homoscedastic Noise",true);
    regType   = parlist->sublist("Problem").get("Regression Type","Least Squares"); 
    useBudget = parlist->sublist("Problem").get("Use Budget Constraint",false);
    ROL::OED::RegressionType type = ROL::OED::StringToRegressionType(regType);
    noise    = ROL::makePtr<Poisson_Noise<RealT>>();
    cov      = ROL::makePtr<ROL::OED::StdMomentOperator<RealT>>(type,homNoise,noise);
    factory  = ROL::makePtr<ROL::OED::Factory<RealT>>(model,dsampler,theta,cov,*parlist);
    ROL::Ptr<ROL::Vector<RealT>> cost;
    if (useBudget) {
      cost = factory->getDesign()->dual().clone();
      const RealT half(0.5), c0(1), c1(8);
      std::vector<RealT> pt;
      for (int i = 0; i < dsampler->numMySamples(); ++i) {
        pt = dsampler->getMyPoint(i);
        RealT norm(0);
        for (const auto x : pt) norm = std::max(norm,std::abs(x-half));
        ROL::staticPtrCast<ROL::ProbabilityVector<RealT>>(cost)->setProbability(i,c0+c1*norm);
      }
      RealT budget(50);
      factory->setBudgetConstraint(cost,budget);
    }
    
    obman->barrier();

    // Print factors
    bool printFactors = parlist->sublist("Problem").get("Print Factors",false);
    if (printFactors) {
      std::stringstream factName_d, noisName_d, factName_o, noisName_o;
      std::ofstream factFile_d, noisFile_d, factFile_o, noisFile_o;
      ROL::Ptr<ROL::Vector<RealT>> Fp = zp->dual().clone();
      std::vector<RealT> pt;

      factName_d << "factors_des_" << dsampler->batchID() << ".txt";
      noisName_d << "noise_des_" << dsampler->batchID() << ".txt";
      factFile_d.open(factName_d.str());
      noisFile_d.open(noisName_d.str());
      factFile_d << std::scientific << std::setprecision(15);
      noisFile_d << std::scientific << std::setprecision(15);
      for (int i = 0; i < dsampler->numMySamples(); ++i) {
        pt = dsampler->getMyPoint(i);
        factory->getFactors()->evaluate(*Fp,pt);
        for (int j = 0; j < numSides; ++j) {
          factFile_d << std::right << std::setw(25)
                     << (*ROL::dynamicPtrCast<ROL::StdVector<RealT>>(Fp)->getVector())[j];
        }
        factFile_d << std::endl;
        for (int j = 0; j < 2; ++j) {
          noisFile_d << std::right << std::setw(25) << pt[j];
        }
        noisFile_d << std::right << std::setw(25) << noise->evaluate(pt) << std::endl;
      }
      factFile_d.close();
      noisFile_d.close();

      factName_o << "factors_opt_" << osampler->batchID() << ".txt";
      noisName_o << "noise_opt_" << osampler->batchID() << ".txt";
      factFile_o.open(factName_o.str());
      noisFile_o.open(noisName_o.str());
      factFile_o << std::scientific << std::setprecision(15);
      noisFile_o << std::scientific << std::setprecision(15);
      for (int i = 0; i < osampler->numMySamples(); ++i) {
        pt = osampler->getMyPoint(i);
        factory->getFactors()->evaluate(*Fp,pt);
        for (int j = 0; j < numSides; ++j) {
          factFile_o << std::right << std::setw(25)
                     << (*ROL::dynamicPtrCast<ROL::StdVector<RealT>>(Fp)->getVector())[j];
        }
        factFile_o << std::endl;
        for (int j = 0; j < 2; ++j) {
          noisFile_o << std::right << std::setw(25) << pt[j];
        }
        noisFile_o << std::right << std::setw(25) << noise->evaluate(pt) << std::endl;
      }
      factFile_o.close();
      noisFile_o.close();
    }

    // Solve OED problem
    dtype = parlist->sublist("OED").get("Optimality Type","A");
    useUIG = parlist->sublist("Problem").sublist("D-Optimal").get("Uniform Initial Guess", false);
    solve<RealT>(factory,dsampler,osampler,*parlist,dtype,*outStream,useUIG,checkDeriv);

    // Get a summary from the time monitor.
    Teuchos::TimeMonitor::summarize();
  }
  catch (std::logic_error &err) {
    std::cout << "Design MPI Rank = " << drank
              << "  Objective MPI Rank = " << orank
              << std::endl << err.what() << std::endl; 
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
