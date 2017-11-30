// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

/*! \file  example_02.cpp
    \brief Shows how to solve the Poisson-Boltzmann problem.
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include <iostream>
#include <algorithm>

#include "ROL_Algorithm.hpp"
#include "ROL_UnaryFunctions.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_OptimizationProblem.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"
#include "ROL_CompositeConstraint_SimOpt.hpp"

#include "../TOOLS/linearpdeconstraint.hpp"
#include "../TOOLS/pdeconstraint.hpp"
#include "../TOOLS/pdeobjective.hpp"
#include "../TOOLS/pdevector.hpp"
#include "../TOOLS/batchmanager.hpp"
#include "pde_poisson_boltzmann_ex02.hpp"
#include "mesh_poisson_boltzmann_ex02.hpp"
#include "obj_poisson_boltzmann_ex02.hpp"

typedef double RealT;

template<class Real>
Real random(const Teuchos::Comm<int> &comm,
            const Real a = -1, const Real b = 1) {
  Real val(0), u(0);
  if ( Teuchos::rank<int>(comm)==0 ) {
    u   = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
    val = (b-a)*u + a;
  }
  Teuchos::broadcast<int,Real>(comm,0,1,&val);
  return val;
}

int main(int argc, char *argv[]) {
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::Ptr<const Teuchos::Comm<int> > comm
    = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
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

    /*** Read in XML input ***/
    std::string filename = "input_ex02.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );

    /*** Initialize main data structure. ***/
    ROL::Ptr<MeshManager<RealT> > meshMgr
      = ROL::makePtr<MeshManager_Example02<RealT>>(*parlist);
    // Initialize PDE describe Poisson's equation
    ROL::Ptr<PDE_Poisson_Boltzmann_ex02<RealT> > pde
      = ROL::makePtr<PDE_Poisson_Boltzmann_ex02<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con
      = ROL::makePtr<PDE_Constraint<RealT>>(pde,meshMgr,serial_comm,*parlist,*outStream);
    ROL::Ptr<PDE_Constraint<RealT> > pdeCon
      = ROL::dynamicPtrCast<PDE_Constraint<RealT> >(con);
    ROL::Ptr<PDE_Doping<RealT> > pdeDoping
      = ROL::makePtr<PDE_Doping<RealT>>(*parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > conDoping
      = ROL::makePtr<Linear_PDE_Constraint<RealT>>(pdeDoping,meshMgr,serial_comm,*parlist,*outStream,true);
    const ROL::Ptr<Assembler<RealT> > assembler = pdeCon->getAssembler();
    assembler->printMeshData(*outStream);
    con->setSolveParameters(*parlist);

    /*************************************************************************/
    /***************** BUILD VECTORS *****************************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<> >  u_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> >  p_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> >  r_ptr = assembler->createResidualVector();
    ROL::Ptr<Tpetra::MultiVector<> >  z_ptr = assembler->createControlVector();
    ROL::Ptr<ROL::Vector<RealT> > up, pp, rp, zp;
    u_ptr->randomize();  //u_ptr->putScalar(static_cast<RealT>(1));
    p_ptr->randomize();  //p_ptr->putScalar(static_cast<RealT>(1));
    r_ptr->randomize();  //r_ptr->putScalar(static_cast<RealT>(1));
    z_ptr->randomize();  //z_ptr->putScalar(static_cast<RealT>(1));
    up = ROL::makePtr<PDE_PrimalSimVector<RealT>>(u_ptr,pde,assembler,*parlist);
    pp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(p_ptr,pde,assembler,*parlist);
    rp = ROL::makePtr<PDE_DualSimVector<RealT>>(r_ptr,pde,assembler,*parlist);
    zp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(z_ptr,pde,assembler,*parlist);

    /*************************************************************************/
    /***************** BUILD SAMPLER *****************************************/
    /*************************************************************************/
    int stochDim = 3;
    std::vector<ROL::Ptr<ROL::Distribution<RealT> > > distVec(stochDim);
    // Build lambda2 distribution
    Teuchos::ParameterList LvolList;
    LvolList.sublist("Distribution").set("Name","Uniform");
    LvolList.sublist("Distribution").sublist("Uniform").set("Lower Bound", -2.0);
    LvolList.sublist("Distribution").sublist("Uniform").set("Upper Bound", -1.0);
    distVec[0] = ROL::DistributionFactory<RealT>(LvolList);
    // Build delta2 distribution
    Teuchos::ParameterList DvolList;
    DvolList.sublist("Distribution").set("Name","Uniform");
    DvolList.sublist("Distribution").sublist("Uniform").set("Lower Bound", -1.0);
    DvolList.sublist("Distribution").sublist("Uniform").set("Upper Bound",  0.0);
    distVec[1] = ROL::DistributionFactory<RealT>(DvolList);
    // Build doping diffusion distribution
    Teuchos::ParameterList dopeList;
    dopeList.sublist("Distribution").set("Name","Uniform");
    dopeList.sublist("Distribution").sublist("Uniform").set("Lower Bound", 2e-6);
    dopeList.sublist("Distribution").sublist("Uniform").set("Upper Bound", 2e-4);
    distVec[2] = ROL::DistributionFactory<RealT>(dopeList);
    // Build sampler
    int nsamp = parlist->sublist("Problem").get("Number of Samples",100);
    ROL::Ptr<ROL::BatchManager<RealT> > bman
      = ROL::makePtr<ROL::TpetraTeuchosBatchManager<RealT>>(comm);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp,distVec,bman);
    // Print samples
    std::vector<RealT> sample(stochDim), Lmean(stochDim), Gmean(stochDim);
    std::stringstream name_samp;
    name_samp << "samples_" << bman->batchID() << ".txt";
    std::ofstream file_samp;
    file_samp.open(name_samp.str());
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      sample = sampler->getMyPoint(i);
      file_samp << std::scientific << std::setprecision(15);
      for (int j = 0; j < stochDim; ++j) {
        file_samp << std::setw(25) << std::left << sample[j];
        Lmean[j] += sampler->getMyWeight(i)*sample[j];
      }
      file_samp << std::endl;
    }
    file_samp.close();
    bman->sumAll(&Lmean[0],&Gmean[0],stochDim);

    /*************************************************************************/
    /***************** BUILD REFERENCE DOPING AND POTENTIAL ******************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<> > ru_ptr = assembler->createStateVector();
    ROL::Ptr<ROL::Vector<RealT> > rup
      = ROL::makePtr<PDE_PrimalSimVector<RealT>>(ru_ptr,pde,assembler,*parlist);
    ROL::Ptr<Tpetra::MultiVector<> > rz_ptr = assembler->createControlVector();
    ROL::Ptr<ROL::Vector<RealT> > rzp
      = ROL::makePtr<PDE_PrimalOptVector<RealT>>(rz_ptr,pde,assembler,*parlist);
    ROL::Ptr<Doping<RealT> > dope
      = ROL::makePtr<Doping<RealT>>(pde->getFE(), pde->getCellNodes(),
                                       assembler->getDofManager()->getCellDofs(),
                                       assembler->getCellIds(),*parlist);
    // Initialize "filtered" of "unfiltered" constraint.
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > pdeWithDoping
      = ROL::makePtr<ROL::CompositeConstraint_SimOpt<RealT>>(con, conDoping,
        *rp, *rp, *up, *zp, *zp, true, true);
    pdeWithDoping->setSolveParameters(*parlist);
    dope->build(rz_ptr);
    RealT tol(1.e-8);
    pdeWithDoping->setParameter(Gmean);
    pdeWithDoping->solve(*rp,*rup,*rzp,tol);
    pdeCon->outputTpetraVector(ru_ptr,"reference_state.txt");
    pdeCon->outputTpetraVector(rz_ptr,"reference_control.txt");

    /*************************************************************************/
    /***************** BUILD COST FUNCTIONAL *********************************/
    /*************************************************************************/
    std::vector<ROL::Ptr<QoI<RealT> > > qoi_vec(3,ROL::nullPtr);
    // Current flow over drain
    qoi_vec[0] = ROL::makePtr<QoI_State_Cost_1_Poisson_Boltzmann<RealT>>(pde->getFE(),
                                    pde->getBdryFE(),pde->getBdryCellLocIds(),*parlist);
    ROL::Ptr<IntegralObjective<RealT> > stateObj
      = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[0],assembler);
    // Deviation from reference doping
    qoi_vec[1] = ROL::makePtr<QoI_Control_Cost_1_Poisson_Boltzmann<RealT>>(pde->getFE(),dope);
    ROL::Ptr<IntegralObjective<RealT> > ctrlObj1
      = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[1],assembler);
    // H1-Seminorm of doping
    qoi_vec[2] = ROL::makePtr<QoI_Control_Cost_2_Poisson_Boltzmann<RealT>>(pde->getFE());
    ROL::Ptr<IntegralObjective<RealT> > ctrlObj2
      = ROL::makePtr<IntegralObjective<RealT>>(qoi_vec[2],assembler);
    // Build standard vector objective function
    RealT currentWeight = parlist->sublist("Problem").get("Desired Current Scale",1.5);
    RealT J = stateObj->value(*rup,*rzp,tol); // Reference current flow over drain
    J *= currentWeight; // Increase current flow by 50%
    RealT w1 = parlist->sublist("Problem").get("State Cost Parameter",1e-3);
    RealT w2 = parlist->sublist("Problem").get("Control Misfit Parameter",1e-2);
    RealT w3 = parlist->sublist("Problem").get("Control Cost Parameter",1e-8);
    ROL::Ptr<ROL::StdObjective<RealT> > std_obj
      = ROL::makePtr<StdObjective_Poisson_Boltzmann<RealT>>(J,w1,w2,w3);
    // Build full-space objective
    ROL::Ptr<PDE_Objective<RealT> > obj
      = ROL::makePtr<PDE_Objective<RealT>>(qoi_vec,std_obj,assembler);
    // Build reduced-space objective
    bool storage = parlist->sublist("Problem").get("Use state storage",true);
    ROL::Ptr<ROL::SimController<RealT> > stateStore
      = ROL::makePtr<ROL::SimController<RealT>>();
    ROL::Ptr<ROL::Reduced_Objective_SimOpt<RealT> > objReduced
      = ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(
                       obj,pdeWithDoping,stateStore,up,zp,pp,storage);
 
    /*************************************************************************/
    /***************** BUILD BOUND CONSTRAINT ********************************/
    /*************************************************************************/
    ROL::Ptr<Tpetra::MultiVector<> > lo_ptr = assembler->createControlVector();
    ROL::Ptr<Tpetra::MultiVector<> > hi_ptr = assembler->createControlVector();
    ROL::Ptr<DopingBounds<RealT> > dopeBnd
      = ROL::makePtr<DopingBounds<RealT>>(pde->getFE(),pde->getCellNodes(),
                                             assembler->getDofManager()->getCellDofs(),
                                             assembler->getCellIds(),*parlist);
    dopeBnd->build(lo_ptr,hi_ptr);
    ROL::Ptr<ROL::Vector<RealT> > lop, hip;
    lop = ROL::makePtr<PDE_PrimalOptVector<RealT>>(lo_ptr,pde,assembler);
    hip = ROL::makePtr<PDE_PrimalOptVector<RealT>>(hi_ptr,pde,assembler);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd
      = ROL::makePtr<ROL::Bounds<RealT>>(lop,hip);
    bool deactivate = parlist->sublist("Problem").get("Deactivate Bound Constraints",false);
    if (deactivate) {
      bnd->deactivate();
    }

    /*************************************************************************/
    /***************** BUILD STOCHASTIC PROBLEM ******************************/
    /*************************************************************************/
    ROL::OptimizationProblem<RealT> opt(objReduced,zp,bnd);
    parlist->sublist("SOL").set("Initial Statistic", static_cast<RealT>(1));
    opt.setStochasticObjective(*parlist,sampler);

    /*************************************************************************/
    /***************** RUN VECTOR AND DERIVATIVE CHECKS **********************/
    /*************************************************************************/
    bool checkDeriv = parlist->sublist("Problem").get("Check Derivatives",false);
    if ( checkDeriv ) {
      ROL::Ptr<Tpetra::MultiVector<> > du_ptr = assembler->createStateVector();
      ROL::Ptr<Tpetra::MultiVector<> > dz_ptr = assembler->createControlVector();
      du_ptr->randomize(); //du_ptr->putScalar(static_cast<RealT>(0));
      dz_ptr->randomize(); //dz_ptr->putScalar(static_cast<RealT>(1));
      ROL::Ptr<ROL::Vector<RealT> > dup, dzp;
      dup = ROL::makePtr<PDE_PrimalSimVector<RealT>>(du_ptr,pde,assembler,*parlist);
      dzp = ROL::makePtr<PDE_PrimalOptVector<RealT>>(dz_ptr,pde,assembler,*parlist);
      // Create ROL SimOpt vectors
      ROL::Vector_SimOpt<RealT> x(up,zp);
      ROL::Vector_SimOpt<RealT> d(dup,dzp);

      objReduced->setParameter(Gmean);
      *outStream << "\n\nCheck Gradient of Full Objective Function\n";
      obj->checkGradient(x,d,true,*outStream);
      *outStream << "\n\nCheck Hessian of Full Objective Function\n";
      obj->checkHessVec(x,d,true,*outStream);
      *outStream << "\n\nCheck Full Jacobian of PDE Constraint\n";
      pdeWithDoping->checkApplyJacobian(x,d,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_1 of PDE Constraint\n";
      pdeWithDoping->checkApplyJacobian_1(*up,*zp,*dup,*rp,true,*outStream);
      *outStream << "\n\nCheck Jacobian_2 of PDE Constraint\n";
      pdeWithDoping->checkApplyJacobian_2(*up,*zp,*dzp,*rp,true,*outStream);
      *outStream << "\n\nCheck Full Hessian of PDE Constraint\n";
      pdeWithDoping->checkApplyAdjointHessian(x,*pp,d,x,true,*outStream);
      *outStream << "\n";
      pdeWithDoping->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
      *outStream << "\n";
      pdeWithDoping->checkInverseJacobian_1(*rp,*dup,*up,*zp,true,*outStream);
      *outStream << "\n";
      *outStream << "\n\nCheck Gradient of Reduced Objective Function\n";
      objReduced->checkGradient(*zp,*dzp,true,*outStream);
      *outStream << "\n\nCheck Hessian of Reduced Objective Function\n";
      objReduced->checkHessVec(*zp,*dzp,true,*outStream);

      opt.check(*outStream);
    }

    /*************************************************************************/
    /***************** SOLVE OPTIMIZATION PROBLEM ****************************/
    /*************************************************************************/
    ROL::Algorithm<RealT> algo("Trust Region",*parlist,false);
    zp->set(*rzp);
    std::clock_t timer = std::clock();
    algo.run(opt,true,*outStream);
    *outStream << "Optimization time: "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;

    /*************************************************************************/
    /***************** OUTPUT RESULTS ****************************************/
    /*************************************************************************/
    std::clock_t timer_print = std::clock();
    // Output control to file
    pdeCon->outputTpetraVector(z_ptr,"control.txt");
    // Output expected state and samples to file
    ROL::Ptr<Tpetra::MultiVector<> > Lu_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> > Lv_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> > Gu_ptr = assembler->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<> > Gv_ptr = assembler->createStateVector();
    ROL::Ptr<ROL::Vector<RealT> > Lup, Gup, Lvp, Gvp;
    Lup = ROL::makePtr<PDE_PrimalSimVector<RealT>>(Lu_ptr,pde,assembler,*parlist);
    Lvp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(Lv_ptr,pde,assembler,*parlist);
    Gup = ROL::makePtr<PDE_PrimalSimVector<RealT>>(Gu_ptr,pde,assembler,*parlist);
    Gvp = ROL::makePtr<PDE_PrimalSimVector<RealT>>(Gv_ptr,pde,assembler,*parlist);
    ROL::Elementwise::Power<RealT> sqr(2.0);
    for (int i = 0; i < sampler->numMySamples(); ++i) {
      sample = sampler->getMyPoint(i);
      con->setParameter(sample);
      con->solve(*rp,*Gup,*zp,tol);
      Gvp->set(*Gup); Gvp->applyUnary(sqr);
      Lup->axpy(sampler->getMyWeight(i),*Gup);
      Lvp->axpy(sampler->getMyWeight(i),*Gvp);
    }
    bman->sumAll(*Lup,*Gup);
    pdeCon->outputTpetraVector(Gu_ptr,"mean_state.txt");
    bman->sumAll(*Lvp,*Gvp);
    Gup->applyUnary(sqr);
    Gvp->axpy(-1.0,*Gup);
    pdeCon->outputTpetraVector(Gv_ptr,"variance_state.txt");
    // Build objective function distribution
    RealT val(0), val1(0), val2(0);
    int nsamp_dist = parlist->sublist("Problem").get("Number of Output Samples",100);
    ROL::Ptr<ROL::SampleGenerator<RealT> > sampler_dist
      = ROL::makePtr<ROL::MonteCarloGenerator<RealT>>(nsamp_dist,distVec,bman);
    std::stringstream name;
    name << "obj_samples_" << bman->batchID() << ".txt";
    std::ofstream file;
    file.open(name.str());
    for (int i = 0; i < sampler_dist->numMySamples(); ++i) {
      sample = sampler_dist->getMyPoint(i);
      objReduced->setParameter(sample);
      val  = objReduced->value(*zp,tol);
      val1 = ctrlObj1->value(*up,*zp,tol);
      val2 = ctrlObj2->value(*up,*zp,tol);
      file << std::scientific << std::setprecision(15);
      for (int j = 0; j < stochDim; ++j) {
        file << std::setw(25) << std::left << sample[j];
      }
      file << std::setw(25) << std::left << val;
      file << std::setw(25) << std::left << (val-w2*val1-w3*val2)/w1;
      file << std::setw(25) << std::left << val1;
      file << std::setw(25) << std::left << val2;
      file << std::endl;
    }
    file.close();
    *outStream << "Output time: "
               << static_cast<RealT>(std::clock()-timer_print)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds." << std::endl << std::endl;
  }
  catch (std::logic_error err) {
    *outStream << err.what() << "\n";
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED\n";
  else
    std::cout << "End Result: TEST PASSED\n";

  return 0;
}
