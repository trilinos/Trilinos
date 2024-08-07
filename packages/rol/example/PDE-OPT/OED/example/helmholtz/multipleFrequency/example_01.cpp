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
#include "ROL_UserInputGenerator.hpp"
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

#include "../src/pde_helmholtz.hpp"
#include "../src/obj_helmholtz.hpp"
#include "../src/noise.hpp"
#include "model.hpp"
#include "predfun.hpp"

template<typename Real>
void solve(const ROL::Ptr<ROL::OED::Factory<Real>>      &factory,
           const ROL::Ptr<ROL::Objective<Real>>         &predfun,
           const ROL::Ptr<ROL::SampleGenerator<Real>>   &dsampler,
           const ROL::Ptr<ROL::SampleGenerator<Real>>   &osampler,
           ROL::ParameterList                           &list,
           const std::string                            &dtype,
           std::ostream                                 &stream     = std::cout,
           const bool                                    useUIG     = true,
           const bool                                    checkDeriv = false) {
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
  ROL::Ptr<ROL::Problem<Real>> problem = factory->get(list,osampler,predfun);
  std::string type = list.sublist("OED").get("Optimality Type","I");
  bool usePD = list.sublist("OED").sublist("R-Optimality").get("Use Primal-Dual Algorithm",false);
  if (type=="R" && usePD) {
    ROL::PrimalDualRisk<Real> solver(problem,osampler,list);
    if (checkDeriv) {
      factory->check(stream);
      //solver.check(stream);
      factory->reset();
    }
    solver.run(stream);
  }
  else {
    ROL::Solver<Real> solver(problem,list);
    if (checkDeriv) {
      factory->check(stream);
      problem->check(true,stream);
      factory->reset();
    }
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
    //RealT tol(1e-8);// one(1);

    /*************************************************************************/
    /******* BUILD LINEAR REGRESSION MODEL BASED ON HELMHOLTZ ****************/
    /*************************************************************************/
    /*** Read in XML input ***/
    ROL::Ptr<ROL::ParameterList> parlist = ROL::getParametersFromXmlFile("input.xml");
    int   numSpeakers = 8;
    bool   checkDeriv = parlist->sublist("Problem").get("Check Derivatives",  false);
    int     verbosity = parlist->sublist("General").get("Print Verbosity",0);
    bool      alPrint = parlist->sublist("Step").sublist("Augmented Lagrangian").get("Print Intermediate Optimization History",false);
    verbosity         = (iprint ? verbosity : 0);
    alPrint           = (iprint ? alPrint : false);
    parlist->sublist("General").set("Print Verbosity",verbosity);
    parlist->sublist("Step").sublist("Augmented Lagrangian").set("Print Intermediate Optimization History",alPrint);

    // Initialize PDE describing Helmholtz equation.
    ROL::Ptr<MeshManager<RealT>> meshMgr;
    ROL::Ptr<PDE_Helmholtz_OCT<RealT>> pde;
    std::vector<ROL::Ptr<Linear_PDE_Constraint<RealT>>> cons(5);
    meshMgr = ROL::makePtr<MeshReader<RealT>>(*parlist);
    for (int i = 0; i < 5; ++i) {
      parlist->sublist("Problem").set("Frequency",static_cast<RealT>((i+1)*100));
      pde = ROL::makePtr<PDE_Helmholtz_OCT<RealT>>(*parlist);
      cons[i] = ROL::makePtr<Linear_PDE_Constraint<RealT>>(pde,meshMgr,dcomm,
                                                           *parlist,*outStream);
    }
    ROL::Ptr<Assembler<RealT>> assembler = cons[0]->getAssembler();
    cons[0]->printMeshData(*outStream);

    // Create state vector.
    ROL::Ptr<Tpetra::MultiVector<>> u_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> r_ptr = assembler->createResidualVector();
    ROL::Ptr<ROL::Vector<RealT>> up, rp, zp, theta, ovec;
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<ROL::StdVector<RealT>>(numSpeakers);
    theta = ROL::makePtr<ROL::StdVector<RealT>>(numSpeakers,1);
    ovec  = ROL::makePtr<ROL::StdVector<RealT>>(5,0);

    // Create observation function
    ROL::Ptr<QoI<RealT>> qoi;
    ROL::Ptr<ROL::Objective_SimOpt<RealT>> obs;
    qoi = ROL::makePtr<QoI_Helmholtz_Observation<RealT>>(pde->getFE(),*parlist);
    obs = ROL::makePtr<PDE_Objective<RealT>>(qoi,assembler);

    // Build nonlinear model and parameter vector
    std::vector<ROL::Ptr<ROL::Objective<RealT>>> objs(5);
    ROL::Ptr<ROL::Constraint<RealT>> model;
    for (int i = 0; i < 5; ++i)
      objs[i] = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obs,cons[i],up,theta,up,false);
    model = ROL::makePtr<Model<RealT>>(objs);

    // Build prediction function
    ROL::Ptr<std::vector<RealT>> vdata = ROL::makePtr<std::vector<RealT>>(5);
    (*vdata)[0] = static_cast<RealT>(0);
    (*vdata)[1] = static_cast<RealT>(1);
    (*vdata)[2] = static_cast<RealT>(2);
    (*vdata)[3] = static_cast<RealT>(3);
    (*vdata)[4] = static_cast<RealT>(4);
    ROL::Ptr<ROL::StdVector<RealT>> vec = ROL::makePtr<ROL::StdVector<RealT>>(vdata);
    ROL::Ptr<PredFun<RealT>> predfun = ROL::makePtr<PredFun<RealT>>(model,vec);

    // Run derivative checks
    if ( checkDeriv ) {
      ROL::Ptr<ROL::Vector<RealT>> dzp, rop, rzp;
      dzp = zp->clone(); rop = ovec->clone(); rzp = zp->clone();
      dzp->randomize();  rop->randomize();    rzp->randomize();
      std::vector<RealT> param(2,0);
      model->setParameter(param);
      *outStream << "\n\nCheck Jacobian of model\n";
      model->checkApplyJacobian(*rzp,*dzp,*rop,true,*outStream);
      *outStream << "\n\nCheck gradient of prediction function\n";
      model->checkAdjointConsistencyJacobian(*rop,*rzp,*dzp,true,*outStream);
      predfun->checkGradient(*rzp,*dzp,true,*outStream);
    }

    /*************************************************************************/
    /******* BUILD EXPERIMENTAL DESIGN PROBLEM AND SOLVE *********************/
    /*************************************************************************/
    // Build samplers for experiment space
    std::string pd = parlist->sublist("Problem").sublist("Design").get("Points File",          "points_5000.txt");
    std::string wd = parlist->sublist("Problem").sublist("Design").get("Weights File",         "weights_5000.txt");
    int nsampd     = parlist->sublist("Problem").sublist("Design").get("Number of Samples",    5000);
    std::string po = parlist->sublist("Problem").sublist("Objective").get("Points File",       "points_200.txt");
    std::string wo = parlist->sublist("Problem").sublist("Objective").get("Weights File",      "weights_200.txt");
    int nsampo     = parlist->sublist("Problem").sublist("Objective").get("Number of Samples", 200);
    *outStream << std::endl;
    *outStream << "  Design Points File:     " << pd << std::endl;
    *outStream << "  Design Weights File:    " << wd << std::endl;
    *outStream << "  Objective Points File:  " << po << std::endl;
    *outStream << "  Objective Weights File: " << wo << std::endl;
    *outStream << std::endl;
    ROL::Ptr<ROL::BatchManager<RealT>> dbman, obman;
    dbman  = ROL::makePtr<ROL::StdTeuchosBatchManager<RealT,int>>(dcomm);
    obman  = ROL::makePtr<ROL::StdTeuchosBatchManager<RealT,int>>(ocomm);
    ROL::Ptr<ROL::SampleGenerator<RealT>> dsampler, osampler;
    dsampler  = ROL::makePtr<ROL::UserInputGenerator<RealT>>(pd,wd,nsampd,2,dbman);
    osampler  = ROL::makePtr<ROL::UserInputGenerator<RealT>>(po,wo,nsampo,2,obman);

    // Build OED problem factory.
    bool useUIG, homNoise;
    std::string regType, dtype;
    ROL::Ptr<ROL::OED::Factory<RealT>>           factory;
    ROL::Ptr<ROL::OED::Noise<RealT>>             noise;
    ROL::Ptr<ROL::OED::StdMomentOperator<RealT>> cov;
    homNoise = parlist->sublist("Problem").get("Homoscedastic Noise",true);
    regType  = parlist->sublist("Problem").get("Regression Type","Least Squares"); 
    ROL::OED::RegressionType type = ROL::OED::StringToRegressionType(regType);
    noise    = ROL::makePtr<Helmholtz_Noise<RealT>>();
    cov      = ROL::makePtr<ROL::OED::StdMomentOperator<RealT>>(type,homNoise,noise);
    factory  = ROL::makePtr<ROL::OED::Factory<RealT>>(model,dsampler,theta,ovec,cov,*parlist);
    obman->barrier();

    // Solve OED problem
    dtype = parlist->sublist("OED").get("Optimality Type","A");
    useUIG = parlist->sublist("Problem").get("Uniform Initial Guess", false);
    solve<RealT>(factory,predfun,dsampler,osampler,*parlist,dtype,*outStream,useUIG,checkDeriv);

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
