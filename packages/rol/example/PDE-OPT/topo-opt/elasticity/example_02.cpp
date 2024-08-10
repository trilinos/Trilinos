// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_02.cpp
    \brief Shows how to solve the stochastic stuctural topology optimization problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_Time.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Solver.hpp"
#include "ROL_StochasticProblem.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_LinearCombinationObjective.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"
#include "ROL_SingletonTeuchosBatchManager.hpp"

#include "../../TOOLS/pdeobjective.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "../../TOOLS/meshreader.hpp"

#include "src/mesh_topo-opt.hpp"
#include "src/pde_elasticity.hpp"
#include "src/sampler.hpp"

#include "../src/pde_filter.hpp"
#include "../src/filtered_compliance_robj.hpp"
#include "../src/volume_con.hpp"
#include "../src/volume_obj.hpp"
#include "../src/obj_volume.hpp"
#include "../src/obj_phasefield.hpp"
#include "../src/printCDF.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::getDefaultComm();
  ROL::Ptr<const Teuchos::Comm<int> > serial_comm
    = ROL::makePtr<Teuchos::SerialComm<int>>();
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
    RealT tol(1e-8);

    /*** Read in XML input ***/
    std::string filename = "input_ex02.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    // Retrieve parameters.
    const std::string example = parlist->sublist("Problem").get("Example", "Default");
    const RealT volFraction   = parlist->sublist("Problem").get("Volume Fraction", 0.4);
    const RealT initDens      = parlist->sublist("Problem").get("Initial Density",volFraction);
    const RealT cmpFactor     = parlist->sublist("Problem").get("Compliance Factor", 1.1);
    int probDim               = parlist->sublist("Problem").get("Problem Dimension", 2);
    const std::string minType = parlist->sublist("Problem").get("Minimization Type", "Volume");
    const bool usePhaseField  = parlist->sublist("Problem").get("Use Phase Field", false);
    bool useFilter            = parlist->sublist("Problem").get("Use Filter", true);
    const bool normalizeObj   = parlist->sublist("Problem").get("Normalize Compliance", true);
    const bool volEq          = parlist->sublist("Problem").get("Use Volume Equality Constraint",true);
    std::string hessAppr      = parlist->sublist("Problem").get("Hessian Approximation","None");
    hessAppr                  = (minType == "Compliance" ? hessAppr : "None");
    if (example == "2D Wheel"                   ||
        example == "2D Truss"                   ||
        example == "2D Cantilever with 1 Load"  ||
        example == "2D Cantilever with 3 Loads" ||
        example == "2D Beams"                   ||
        example == "2D Carrier Plate")
      probDim = 2;
    else if (example == "3D Cantilever" ||
             example == "3D L Beam")
      probDim = 3;
    if (usePhaseField) useFilter = false;

    /*** Initialize main data structure. ***/
    TEUCHOS_TEST_FOR_EXCEPTION(probDim<2 || probDim>3, std::invalid_argument,
      ">>> PDE-OPT/topo-opt/elasticity/example_01.cpp: Problem dim is not 2 or 3!");
    ROL::Ptr<MeshManager<RealT>> meshMgr;
    if (probDim == 2) meshMgr = ROL::makePtr<MeshManager_TopoOpt<RealT>>(*parlist);
    else              meshMgr = ROL::makePtr<MeshReader<RealT>>(*parlist);
    // Initialize PDE describing elasticity equations.
    ROL::Ptr<PDE_Elasticity<RealT>>
      pde = ROL::makePtr<PDE_Elasticity<RealT>>(*parlist);

    // Initialize the filter PDE.
    ROL::Ptr<PDE<RealT>> pdeFilter = pde;

    // Initialize reduced compliance objective function.
    ROL::Ptr<ROL::Objective<RealT>> robj_com;
    ROL::Ptr<Assembler<RealT>> assembler, assemblerFilter;
    if (useFilter && !usePhaseField) {
      pdeFilter = ROL::makePtr<PDE_Filter<RealT>>(*parlist);
      pde->setDensityFields(pdeFilter->getFields());
      robj_com = ROL::makePtr<TopOptFilteredComplianceObjective<RealT>>(
                 pde,pdeFilter,meshMgr,comm,*parlist,*outStream);
      assembler = ROL::dynamicPtrCast<TopOptFilteredComplianceObjective<RealT>>(robj_com)->getAssembler();
      assemblerFilter = ROL::dynamicPtrCast<TopOptFilteredComplianceObjective<RealT>>(robj_com)->getFilterAssembler();
    }
    else {
      robj_com = ROL::makePtr<TopOptComplianceObjective<RealT>>(
                 pde,meshMgr,comm,*parlist,*outStream);
      assembler = ROL::dynamicPtrCast<TopOptComplianceObjective<RealT>>(robj_com)->getAssembler();
      assemblerFilter = assembler;
    }

    // Create vectors.
    ROL::Ptr<Tpetra::MultiVector<>> u_ptr, p_ptr, r_ptr, z_ptr;
    u_ptr = assembler->createStateVector();         u_ptr->putScalar(0.0);
    p_ptr = assembler->createStateVector();         p_ptr->putScalar(0.0);
    r_ptr = assembler->createResidualVector();      r_ptr->putScalar(0.0);
    z_ptr = assemblerFilter->createControlVector(); z_ptr->putScalar(1.0);
    ROL::Ptr<ROL::Vector<RealT>> up, pp, rp, zp;
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pdeFilter,assemblerFilter,*parlist);

    // Build volume objective function.
    ROL::Ptr<QoI<RealT>> qoi_vol;
    if (minType == "Compliance" && volEq) {
      if (useFilter && !usePhaseField)
        qoi_vol = ROL::makePtr<QoI_Volume_TopoOpt<RealT>>(ROL::dynamicPtrCast<PDE_Filter<RealT>>(pdeFilter)->getDensityFE(),
                                                          volFraction);
      else
        qoi_vol = ROL::makePtr<QoI_Volume_TopoOpt<RealT>>(pde->getDensityFE(),
                                                          volFraction);
    }
    else {
      if (useFilter && !usePhaseField)
        qoi_vol = ROL::makePtr<QoI_Volume_TopoOpt<RealT>>(ROL::dynamicPtrCast<PDE_Filter<RealT>>(pdeFilter)->getDensityFE());
      else
        qoi_vol = ROL::makePtr<QoI_Volume_TopoOpt<RealT>>(pde->getDensityFE());
    }
    ROL::Ptr<TopOptVolumeObjective<RealT>>
      obj_vol = ROL::makePtr<TopOptVolumeObjective<RealT>>(qoi_vol,assemblerFilter,zp);
    ROL::Ptr<TopOptVolumeConstraint<RealT>>
      con_vol = ROL::makePtr<TopOptVolumeConstraint<RealT>>(qoi_vol,assemblerFilter,zp);

    // Compute volume
    zp->setScalar(1.0);
    RealT vol = obj_vol->value(*zp,tol);

    // Normalize compliance objective function
    RealT cs(1);
    zp->setScalar(initDens);
    if (normalizeObj) {
      if (useFilter && !usePhaseField)
        cs = ROL::dynamicPtrCast<TopOptFilteredComplianceObjective<RealT>>(robj_com)->normalize(*zp,tol);
      else
        cs = ROL::dynamicPtrCast<TopOptComplianceObjective<RealT>>(robj_com)->normalize(*zp,tol);
    }

    // Output problem details
    *outStream << std::endl;
    *outStream << "Problem Data"          << std::endl;
    *outStream << "  Example:           " << example << std::endl;
    *outStream << "  Dimension:         " << probDim << std::endl;
    *outStream << "  Minimize Type:     " << minType << std::endl;
    *outStream << "  Use Phase Field:   " << usePhaseField << std::endl;
    if (!usePhaseField) {
      *outStream << "  SIMP Power:        "
                 << parlist->sublist("Problem").get("SIMP Power",3.0) << std::endl;
      *outStream << "  Use Filter:        " << useFilter << std::endl;
    }
    if (minType == "Volume") {
      *outStream << "  Compliance Factor: " << cmpFactor << std::endl;
    }
    else if (minType == "Compliance") {
      *outStream << "  Volume Fraction:   " << volFraction << std::endl;
      *outStream << "  Initial Density:   " << initDens    << std::endl;
      *outStream << "  Use Equality:      " << volEq       << std::endl;
      *outStream << "  Hessian Approx:    " << hessAppr    << std::endl;
    }
    *outStream << "  Domain Volume:     " << vol << std::endl;
    *outStream << "  Compliance Scale:  " << cs  << std::endl;
    *outStream << std::endl;

    // Create objective, constraint, multiplier and bounds
    ROL::Ptr<ROL::Objective<RealT>>       obj;
    ROL::Ptr<ROL::Constraint<RealT>>      icon;
    ROL::Ptr<ROL::Vector<RealT>>          iup, ilp;
    ROL::Ptr<ROL::Vector<RealT>>          imul;
    ROL::Ptr<ROL::BoundConstraint<RealT>> ibnd;
    if (!usePhaseField) {
      if (minType == "Volume") {
        obj  = obj_vol;
        icon = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>(robj_com);
        // Set upper bound to compliance for solid beam.
        zp->setScalar(1.0);
        RealT comp = robj_com->value(*zp,tol);
        ilp  = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
        iup  = ROL::makePtr<ROL::SingletonVector<RealT>>(cmpFactor*comp);
        imul = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
        ibnd = ROL::makePtr<ROL::Bounds<RealT>>(ilp,iup);
      }
      else if (minType == "Compliance") {
        obj  = robj_com;
        icon = con_vol;
        imul = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
        if (volEq) {
          ilp  = ROL::nullPtr;
          iup  = ROL::nullPtr;
          ibnd = ROL::nullPtr;
        }
        else {
          // Set upper bound to fraction of total volume.
          ilp  = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
          iup  = ROL::makePtr<ROL::SingletonVector<RealT>>(volFraction*vol);
          ibnd = ROL::makePtr<ROL::Bounds<RealT>>(ilp,iup);
        }
      }
      else if (minType == "Total") {
        std::vector<ROL::Ptr<ROL::Objective<RealT>>> obj_vec(2,ROL::nullPtr);
        obj_vec[0] = robj_com;
        obj_vec[1] = obj_vol;
        std::vector<RealT> weights(2,0.0);
        weights[0] = static_cast<RealT>(1);
        weights[1] = parlist->sublist("Problem").get("Volume Objective Scale",0.04096);
        obj  = ROL::makePtr<ROL::LinearCombinationObjective<RealT>>(weights,obj_vec);
        icon = ROL::nullPtr;
        ilp  = ROL::nullPtr;
        iup  = ROL::nullPtr;
        imul = ROL::nullPtr;
        ibnd = ROL::nullPtr;
      }
      else {
        throw ROL::Exception::NotImplemented("Unknown minimization type!");
      }
    }
    else {
      // Build Modica-Mortola Energy objective function.
      RealT penParam = parlist->sublist("Problem").get("Phase Field Penalty Parameter",1e-1);
      ROL::Ptr<QoI<RealT>>
        qoi_pfe = ROL::makePtr<QoI_ModicaMortolaEnergy_PhaseField<RealT>>(pde->getDensityFE(),
                                                                          penParam);
      ROL::Ptr<IntegralOptObjective<RealT>>
        obj_pfe = ROL::makePtr<IntegralOptObjective<RealT>>(qoi_pfe,assembler);
      // Get weights for linear combination objective.
      std::vector<RealT> weights(1,0.0);
      weights[0] = parlist->sublist("Problem").get("Phase Field Energy Objective Scale",0.00064);
      std::vector<ROL::Ptr<ROL::Objective<RealT>>> obj_vec(1,ROL::nullPtr);
      obj_vec[0] = obj_pfe;
      if (minType == "Volume") {
        weights.push_back(parlist->sublist("Problem").get("Volume Objective Scale",0.04096));
        obj_vec.push_back(obj_vol);
        obj  = ROL::makePtr<ROL::LinearCombinationObjective<RealT>>(weights,obj_vec);
        icon = ROL::makePtr<ROL::ConstraintFromObjective<RealT>>(robj_com);
        // Set upper bound to compliance for solid beam.
        RealT comp = robj_com->value(*zp,tol);
        ilp  = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
        iup  = ROL::makePtr<ROL::SingletonVector<RealT>>(cmpFactor*comp);
        imul = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
        ibnd = ROL::makePtr<ROL::Bounds<RealT>>(ilp,iup);
      }
      else if (minType == "Compliance") {
        weights.push_back(static_cast<RealT>(1));
        obj_vec.push_back(robj_com);
        obj  = ROL::makePtr<ROL::LinearCombinationObjective<RealT>>(weights,obj_vec);
        icon = con_vol;
        imul = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
        if (volEq) {
          ilp  = ROL::nullPtr;
          iup  = ROL::nullPtr;
          ibnd = ROL::nullPtr;
        }
        else {
          // Set upper bound to fraction of total volume.
          ilp  = ROL::makePtr<ROL::SingletonVector<RealT>>(0);
          iup  = ROL::makePtr<ROL::SingletonVector<RealT>>(volFraction*vol);
          ibnd = ROL::makePtr<ROL::Bounds<RealT>>(ilp,iup);
        }
      }
      else if (minType == "Total") {
        weights.push_back(static_cast<RealT>(1));
        obj_vec.push_back(robj_com);
        weights.push_back(parlist->sublist("Problem").get("Volume Objective Scale",0.04096));
        obj_vec.push_back(obj_vol);
        obj  = ROL::makePtr<ROL::LinearCombinationObjective<RealT>>(weights,obj_vec);
        icon = ROL::nullPtr;
        ilp  = ROL::nullPtr;
        iup  = ROL::nullPtr;
        imul = ROL::nullPtr;
        ibnd = ROL::nullPtr;
      }
      else {
        throw ROL::Exception::NotImplemented("Unknown minimization type!");
      }
    }

    // Initialize bound constraints.
    RealT lval = (usePhaseField ? -1.0 : 0.0), uval = 1.0;
    ROL::Ptr<ROL::Vector<RealT>> lop = zp->clone(); lop->setScalar(lval);
    ROL::Ptr<ROL::Vector<RealT>> hip = zp->clone(); hip->setScalar(uval);
    ROL::Ptr<ROL::BoundConstraint<RealT>>
      bnd = ROL::makePtr<ROL::Bounds<RealT>>(lop,hip);
    if (usePhaseField) bnd->deactivate();

    // Build sampler.
    BuildSampler<RealT> buildSampler(parlist->sublist("Problem"),probDim,example);
    //int stochDim = buildSampler.getDimension();
    int nsamp = parlist->sublist("Problem").get("Number of samples", 1);
    ROL::Ptr<ROL::BatchManager<RealT>>
      bman = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT>>
      sampler = buildSampler.get(nsamp,bman);

    // Set up optimization problem.
    ROL::Ptr<ROL::StochasticProblem<RealT>>
      prob = ROL::makePtr<ROL::StochasticProblem<RealT>>(obj,zp);
    if (bnd->isActivated()) prob->addBoundConstraint(bnd);
    if ( minType == "Compliance" ) {
      if ( volEq ) prob->addLinearConstraint("Volume",icon,imul);
      else         prob->addLinearConstraint("Volume",icon,imul,ibnd);
    }
    else if ( minType == "Volume" )
      prob->addConstraint("Compliance",icon,imul,ibnd);
    if ( minType == "Volume" ) {
      ROL::Ptr<ROL::BatchManager<RealT>>
        cbman = ROL::makePtr<ROL::SingletonTeuchosBatchManager<RealT,int>>(comm);
      parlist->sublist("SOL").sublist("Compliance") = parlist->sublist("SOL");
      prob->makeConstraintStochastic("Compliance",*parlist,sampler,cbman);
    }
    else {
      parlist->sublist("SOL").sublist("Objective") = parlist->sublist("SOL");
      prob->makeObjectiveStochastic(*parlist,sampler);
    }
    prob->setProjectionAlgorithm(*parlist);
    prob->finalize(false,true,*outStream);

    // Check derivatives.
    bool derivCheck = parlist->sublist("Problem").get("Check derivatives",false);
    if (derivCheck) prob->check(true,*outStream);

    // Solve optimization problem.
    Teuchos::Time algoTimer("Algorithm Time", true);
    ROL::Solver<RealT> solver(prob,*parlist);
    solver.solve(*outStream);
    algoTimer.stop();
    *outStream << "Total optimization time = " << algoTimer.totalElapsedTime() << " seconds.\n";

    // Output.
    if (useFilter && !usePhaseField)
      ROL::dynamicPtrCast<TopOptFilteredComplianceObjective<RealT>>(robj_com)->summarize(*outStream);
    else
      ROL::dynamicPtrCast<TopOptComplianceObjective<RealT>>(robj_com)->summarize(*outStream);
    if (minType == "Compliance") con_vol->summarize(*outStream);
    else                         obj_vol->summarize(*outStream);
    assembler->printMeshData(*outStream);
    up->zero(); pp->zero();
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      robj_com->setParameter(sampler->getMyPoint(i));
      if (useFilter && !usePhaseField)
        ROL::dynamicPtrCast<TopOptFilteredComplianceObjective<RealT>>(robj_com)->solveState(*up,*zp);
      else
        ROL::dynamicPtrCast<TopOptComplianceObjective<RealT>>(robj_com)->solveState(*up,*zp);
      pp->axpy(sampler->getMyWeight(i),*up);
    }
    sampler->sumAll(*pp,*up);
    assembler->outputTpetraVector(u_ptr,"mean_state.txt");
    assembler->outputTpetraVector(z_ptr,"density.txt");
    int nsamp_dist = parlist->sublist("Problem").get("Number of output samples",1000);
    ROL::Ptr<ROL::SampleGenerator<RealT>>
      sampler_dist = buildSampler.get(nsamp_dist,bman);
    if (minType != "Volume")
      TopOptPrintCDF<RealT>(*obj,*zp,*sampler_dist,nsamp_dist,comm,"obj_samples.txt");
    else
      TopOptPrintCDF<RealT>(*robj_com,*zp,*sampler_dist,nsamp_dist,comm,"obj_samples.txt");

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
