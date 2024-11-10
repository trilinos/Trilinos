// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Shows how to solve the mother problem of PDE-constrained optimization:
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include "ROL_Algorithm.hpp"
#include "ROL_AugmentedLagrangianStep.hpp"
//#include "ROL_MoreauYosidaPenaltyStep.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_ScaledTpetraMultiVector.hpp"

#include "ROL_ScaledStdVector.hpp"

#include "ROL_Bounds.hpp"
#include "ROL_AugmentedLagrangian.hpp"
#include "ROL_Algorithm.hpp"

#include <iostream>
#include <algorithm>

#include "data.hpp"
#include "filter.hpp"
#include "objective.hpp"
#include "constraint.hpp"
#include "volume_constraint.hpp"

//#include <fenv.h>

typedef double RealT;

ROL::Ptr<Tpetra::MultiVector<> > createTpetraVector(const ROL::Ptr<const Tpetra::Map<> > &map) {
  return ROL::makePtr<Tpetra::MultiVector<>>(map, 1, true);
}

int main(int argc, char *argv[]) {
  //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0)) {
    outStream = ROL::makePtrFromRef(std::cout);
  }
  else {
    outStream = ROL::makePtrFromRef(bhs);
  }

  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename = "input.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize main data structure. ***/
    ROL::Ptr<ElasticitySIMPOperators<RealT> > data
      = ROL::makePtr<ElasticitySIMPOperators<RealT>>(comm, parlist, outStream);
    /*** Initialize density filter. ***/
    ROL::Ptr<DensityFilter<RealT> > filter
      = ROL::makePtr<DensityFilter<RealT>>(comm, parlist, outStream);
    /*** Build vectors and dress them up as ROL vectors. ***/
    ROL::Ptr<const Tpetra::Map<> > vecmap_u = data->getDomainMapA();
    ROL::Ptr<const Tpetra::Map<> > vecmap_z = data->getCellMap();
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr      = createTpetraVector(vecmap_u);
    ROL::Ptr<Tpetra::MultiVector<> > z_ptr      = createTpetraVector(vecmap_z);
    ROL::Ptr<Tpetra::MultiVector<> > du_ptr     = createTpetraVector(vecmap_u);
    ROL::Ptr<Tpetra::MultiVector<> > dw_ptr     = createTpetraVector(vecmap_u);
    ROL::Ptr<Tpetra::MultiVector<> > dz_ptr     = createTpetraVector(vecmap_z);
    ROL::Ptr<Tpetra::MultiVector<> > dz2_ptr    = createTpetraVector(vecmap_z);
    ROL::Ptr<std::vector<RealT> >    vc_ptr     = ROL::makePtr<std::vector<RealT>>(1, 0);
    ROL::Ptr<std::vector<RealT> >    vc_lam_ptr = ROL::makePtr<std::vector<RealT>>(1, 0);
    ROL::Ptr<std::vector<RealT> >    vscale_ptr = ROL::makePtr<std::vector<RealT>>(1, 0);
    // Set all values to 1 in u, z.
    u_ptr->putScalar(1.0);
    // Set z to gray solution.
    RealT volFrac = parlist->sublist("ElasticityTopoOpt").get<RealT>("Volume Fraction");
    z_ptr->putScalar(volFrac);
    // Set scaling vector for constraint
    RealT W = parlist->sublist("Geometry").get<RealT>("Width");
    RealT H = parlist->sublist("Geometry").get<RealT>("Height");
    RealT one(1), two(2);
    (*vscale_ptr)[0] = one/std::pow(W*H*(one-volFrac),two);
    // Set Scaling vector for density
    bool  useZscale = parlist->sublist("Problem").get<bool>("Use Scaled Density Vectors");
    RealT densityScaling = parlist->sublist("Problem").get<RealT>("Density Scaling");
    ROL::Ptr<Tpetra::MultiVector<> > scaleVec = createTpetraVector(vecmap_z);
    scaleVec->putScalar(densityScaling);
    if ( !useZscale ) {
      scaleVec->putScalar(one);
    }
    ROL::Ptr<const Tpetra::Vector<> > zscale_ptr = scaleVec->getVector(0);

   //test     
   /*data->updateMaterialDensity (z_ptr);
    ROL::Ptr<Tpetra::MultiVector<RealT> > rhs
      = ROL::makePtr<Tpetra::MultiVector<>> (data->getVecF()->getMap(), 1, true);
    data->ApplyMatAToVec(rhs, u_ptr);
    data->outputTpetraVector(rhs, "KU0.txt");
    data->ApplyInverseJacobian1ToVec(u_ptr, rhs, false);
    data->outputTpetraVector(u_ptr, "KKU0.txt");
    
    data->ApplyJacobian1ToVec(rhs, u_ptr);
    data->outputTpetraVector(rhs, "KU1.txt");
    data->ApplyInverseJacobian1ToVec(u_ptr, rhs, false);
    data->outputTpetraVector(u_ptr, "KKU1.txt");
  */
    //u_ptr->putScalar(1.0);
    //z_ptr->putScalar(1.0);
    // Randomize d vectors.
    du_ptr->randomize(); //du_ptr->scale(0);
    dw_ptr->randomize();
    dz_ptr->randomize(); //dz_ptr->scale(0);
    dz2_ptr->randomize();
    // Create ROL::TpetraMultiVectors.
    ROL::Ptr<ROL::Vector<RealT> > up   = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(u_ptr);
    ROL::Ptr<ROL::Vector<RealT> > dup  = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(du_ptr);
    ROL::Ptr<ROL::Vector<RealT> > dwp  = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(dw_ptr);
    ROL::Ptr<ROL::Vector<RealT> > zp 
      = ROL::makePtr<ROL::PrimalScaledTpetraMultiVector<RealT>>(z_ptr,zscale_ptr);
    ROL::Ptr<ROL::Vector<RealT> > dzp
      = ROL::makePtr<ROL::PrimalScaledTpetraMultiVector<RealT>>(dz_ptr,zscale_ptr);
    ROL::Ptr<ROL::Vector<RealT> > dz2p
      = ROL::makePtr<ROL::PrimalScaledTpetraMultiVector<RealT>>(dz2_ptr,zscale_ptr);
    ROL::Ptr<ROL::Vector<RealT> > vcp
      = ROL::makePtr<ROL::PrimalScaledStdVector<RealT>>(vc_ptr,vscale_ptr);
    ROL::Ptr<ROL::Vector<RealT> > vc_lamp
      = ROL::makePtr<ROL::DualScaledStdVector<RealT>>(vc_lam_ptr,vscale_ptr);
    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);
    ROL::Vector_SimOpt<RealT> d2(dwp,dz2p);

    /*** Build objective function, constraint and reduced objective function. ***/
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj
       = ROL::makePtr<Objective_PDEOPT_ElasticitySIMP<RealT>>(data, filter, parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
       = ROL::makePtr<EqualityConstraint_PDEOPT_ElasticitySIMP<RealT>>(data, filter, parlist);
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > objReduced
       = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, dwp);
    ROL::Ptr<ROL::Constraint<RealT> > volcon
       = ROL::makePtr<EqualityConstraint_PDEOPT_ElasticitySIMP_Volume<RealT>>(data, parlist);

    /*** Build bound constraint ***/
    ROL::Ptr<Tpetra::MultiVector<> > lo_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > hi_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    lo_ptr->putScalar(0.0); hi_ptr->putScalar(1.0);
    ROL::Ptr<ROL::Vector<RealT> > lop
      = ROL::makePtr<ROL::PrimalScaledTpetraMultiVector<RealT>>(lo_ptr, zscale_ptr);
    ROL::Ptr<ROL::Vector<RealT> > hip
      = ROL::makePtr<ROL::PrimalScaledTpetraMultiVector<RealT>>(hi_ptr, zscale_ptr);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd = ROL::makePtr<ROL::Bounds<RealT>>(lop,hip);

    /*** Check functional interface. ***/
    *outStream << "Checking Objective:" << "\n";
    obj->checkGradient(x,d,true,*outStream);
    obj->checkHessVec(x,d,true,*outStream);
    obj->checkHessSym(x,d,d2,true,*outStream);
    *outStream << "Checking Constraint:" << "\n";
    con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
    con->checkAdjointConsistencyJacobian_1(*dwp, *dup, *up, *zp, true, *outStream);
    con->checkAdjointConsistencyJacobian_2(*dwp, *dzp, *up, *zp, true, *outStream);
    con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
    con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);
    con->checkApplyJacobian(x,d,*up,true,*outStream);
    con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
    *outStream << "Checking Reduced Objective:" << "\n";
    objReduced->checkGradient(*zp,*dzp,true,*outStream);
    objReduced->checkHessVec(*zp,*dzp,true,*outStream);
    *outStream << "Checking Volume Constraint:" << "\n";
    volcon->checkAdjointConsistencyJacobian(*vcp,*dzp,*zp,true,*outStream);
    volcon->checkApplyJacobian(*zp,*dzp,*vcp,true,*outStream);
    volcon->checkApplyAdjointHessian(*zp,*vcp,*dzp,*zp,true,*outStream);

    /*** Run optimization ***/
    ROL::Ptr<ROL::Step<RealT>>
      step = ROL::makePtr<ROL::AugmentedLagrangianStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      status = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algo(step,status,false);
    ROL::AugmentedLagrangian<RealT> augLag(objReduced,volcon,*vc_lamp,1.0,*zp,*vcp,*parlist);
    algo.run(*zp,*vc_lamp,augLag,*volcon,*bnd,true,*outStream);
    //ROL::Ptr<ROL::Step<RealT>>
    //  step = ROL::makePtr<ROL::MoreauYosidaPenaltyStep<RealT>>(*parlist);
    //ROL::Ptr<ROL::StatusTest<RealT>>
    //  status = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(*parlist);
    //ROL::MoreauYosidaPenalty<RealT> MYpen(objReduced,bnd,*zp,*parlist);
    //ROL::Algorithm<RealT> algo(step,status,false);
    //algo.run(*zp,*vc_lamp,MYpen,*volcon,*bnd,true,*outStream);

    // new filter, for testing
    /*parlist->sublist("Density Filter").set("Enable", true);
    ROL::Ptr<DensityFilter<RealT> > testfilter
      = ROL::makePtr<DensityFilter<RealT>>(comm, parlist, outStream);
    ROL::Ptr<Tpetra::MultiVector<> > z_filtered_ptr = ROL::makePtr<Tpetra::MultiVector<>>(*z_ptr, Teuchos::Copy);
    testfilter->apply(z_filtered_ptr, z_ptr);
    ROL::Ptr<Tpetra::MultiVector<> > cm_ptr = data->getCellAreas();
    ROL::Ptr<Tpetra::MultiVector<> > icm_ptr = ROL::makePtr<Tpetra::MultiVector<>>(*cm_ptr, Teuchos::Copy);
    ROL::Ptr<Tpetra::MultiVector<> > zf_scaled_ptr = ROL::makePtr<Tpetra::MultiVector<>>(*z_ptr, Teuchos::Copy);
    icm_ptr->reciprocal(*cm_ptr);
    zf_scaled_ptr->elementWiseMultiply(1.0, *(icm_ptr->getVector(0)), *z_filtered_ptr, 0.0);
    data->outputTpetraVector(zf_scaled_ptr, "density_filtered_scaled.txt");
    data->outputTpetraVector(z_filtered_ptr, "density_filtered.txt");*/
    
    data->outputTpetraVector(z_ptr, "density.txt");
    data->outputTpetraVector(u_ptr, "state.txt");
    data->outputTpetraVector(zscale_ptr, "weights.txt");

    // Get a summary from the time monitor.
    Teuchos::TimeMonitor::summarize();
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
