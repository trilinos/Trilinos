// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_01.cpp
    \brief Solves a source inversion problem governed by the
           advection-diffusion equation.  Performs optimal
           experimental design (OED).
*/

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"

#include "ROL_Algorithm.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_TrustRegionStep.hpp"
#include "ROL_CompositeStep.hpp"
#include "ROL_Reduced_Objective_SimOpt.hpp"
#include "ROL_TpetraMultiVector.hpp"

#include <iostream>
#include <algorithm>

#include "data.hpp"
#include "objective.hpp"
#include "constraint.hpp"
#include "ROL_ExperimentDesignObjective.hpp"

typedef double RealT;

class MyInterfaceOED : public ROL::ExperimentDesignInterface<RealT> {
public:
  MyInterfaceOED(const ROL::Ptr<ROL::Objective_SimOpt<RealT> > &obj,
                 const ROL::Ptr<ROL::Constraint_SimOpt<RealT> > &con,
                 const ROL::Ptr<ROL::Vector<RealT> > &state,
                 const ROL::Ptr<ROL::Vector<RealT> > &stateDual,
                 const ROL::Ptr<ROL::Vector<RealT> > &control,
                 const ROL::Ptr<ROL::Vector<RealT> > &controlDual,
                 const ROL::Ptr<ROL::Vector<RealT> > &constraint,
                 const ROL::Ptr<ROL::Vector<RealT> > &constraintDual,
                 const ROL::Ptr<ROL::Vector<RealT> > &observation,
                 const ROL::Ptr<ROL::Vector<RealT> > &observationDual,
                 const std::vector<ROL::Ptr<ROL::Vector<RealT> > > &randvecs,
                 const std::vector<ROL::Ptr<ROL::Vector<RealT> > > &training) :
    ExperimentDesignInterface<RealT>(obj, con, state, stateDual, control, controlDual, constraint, constraintDual, observation, observationDual, randvecs, training) {}

  // Override interface functions; in this case, they are the same.

  virtual void applyObserveOp(ROL::Vector<RealT> &obsv, const ROL::Vector<RealT> &v) const {
    obsv.set(v);
  }

  virtual void applyAdjointObserveOp(ROL::Vector<RealT> &aobsv, const ROL::Vector<RealT> &v) const {
    aobsv.set(v);
  }

  virtual void applyWeightOp(ROL::Vector<RealT> &weightv, const ROL::Vector<RealT> &v, const ROL::Vector<RealT> &w) const {
    weightv.set(v.dual());
    weightv.applyBinary(ROL::Elementwise::Multiply<RealT>(), w);
  }

}; // class MyInterfaceOED

// Used for elementwise is-greather-than?
template<class Real>
class IsGreaterThan : public ROL::Elementwise::UnaryFunction<Real> {
private:
  const Real num_, yes_, no_;
public:
  IsGreaterThan(const Real &num, const Real &yes, const Real &no) : num_(num), yes_(yes), no_(no) {}
  Real apply(const Real &x) const {
    if (x > num_) {
      return static_cast<Real>(yes_);
    }
    else {
      return static_cast<Real>(no_);
    }
  }
}; // class IsGreaterThan

int main(int argc, char *argv[]) {

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
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
    std::string filename    = "input.xml";
    std::string filenameOED = "inputOED.xml";
    Teuchos::RCP<Teuchos::ParameterList> parlist = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::RCP<Teuchos::ParameterList> parlistOED = Teuchos::rcp( new Teuchos::ParameterList() );
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    Teuchos::updateParametersFromXmlFile( filenameOED, parlistOED.ptr() );

    /*** Initialize main data structure. ***/
    ROL::Ptr<PoissonData<RealT> > data = ROL::makePtr<PoissonData<RealT>>(comm, parlist, outStream);

    // Get random weights parameter.
    RealT fnzw = parlist->sublist("Problem").get("Fraction of nonzero weights", 0.5);
    fnzw = 1.0 - 2.0*fnzw;

    /*** Build vectors and dress them up as ROL vectors. ***/
    ROL::Ptr<const Tpetra::Map<> > vecmap_u = data->getMatA()->getDomainMap();
    ROL::Ptr<const Tpetra::Map<> > vecmap_z = data->getMatB()->getDomainMap();
    ROL::Ptr<const Tpetra::Map<> > vecmap_c = data->getMatA()->getRangeMap();
    ROL::Ptr<Tpetra::MultiVector<> > u_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > p_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > w_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > wup_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > wlo_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > z_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > c_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_c, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > du_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > dw_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::Ptr<Tpetra::MultiVector<> > dz_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    // Set all values to 1 in u, z and c.
    u_ptr->putScalar(1.0);
    p_ptr->putScalar(1.0);
    z_ptr->putScalar(1.0);
    c_ptr->putScalar(1.0);
    w_ptr->randomize();
    wlo_ptr->putScalar(0.0);
    wup_ptr->putScalar(1.0);
    // Randomize d vectors.
    du_ptr->randomize();
    dw_ptr->randomize();
    dz_ptr->randomize();
    // Create ROL::TpetraMultiVectors.
    ROL::Ptr<ROL::Vector<RealT> > up = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(u_ptr);
    ROL::Ptr<ROL::Vector<RealT> > pp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(p_ptr);
    ROL::Ptr<ROL::Vector<RealT> > wp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(w_ptr);
    ROL::Ptr<ROL::Vector<RealT> > wlop = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(wlo_ptr);
    ROL::Ptr<ROL::Vector<RealT> > wupp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(wup_ptr);
    ROL::Ptr<ROL::Vector<RealT> > zp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(z_ptr);
    ROL::Ptr<ROL::Vector<RealT> > cp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(c_ptr);
    ROL::Ptr<ROL::Vector<RealT> > dup = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(du_ptr);
    ROL::Ptr<ROL::Vector<RealT> > dwp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(dw_ptr);
    ROL::Ptr<ROL::Vector<RealT> > dzp = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(dz_ptr);
    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    /*** Build objective function, constraint and reduced objective function. ***/
    wp->applyUnary(IsGreaterThan<RealT>(fnzw, 1.0, 0.0));
    ROL::Ptr<ROL::Objective_SimOpt<RealT> > obj =
      ROL::makePtr<Objective_PDEOPT_Poisson<RealT>>(data, w_ptr, parlist);
    ROL::Ptr<ROL::Constraint_SimOpt<RealT> > con =
      ROL::makePtr<EqualityConstraint_PDEOPT_Poisson<RealT>>(data, parlist);
    ROL::Ptr<ROL::Objective<RealT> > objReduced =
      ROL::makePtr<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp);

    /*** Check functional interface. ***/
    obj->checkGradient(x,d,true,*outStream);
    obj->checkHessVec(x,d,true,*outStream);
    con->checkApplyJacobian(x,d,*up,true,*outStream);
    con->checkApplyAdjointHessian(x,*dup,d,x,true,*outStream);
    con->checkAdjointConsistencyJacobian(*dup,d,x,true,*outStream);
    con->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
    con->checkInverseAdjointJacobian_1(*up,*up,*up,*zp,true,*outStream);
    objReduced->checkGradient(*zp,*dzp,true,*outStream);
    objReduced->checkHessVec(*zp,*dzp,true,*outStream);

    RealT tol = 1e-8;
    zp->zero();
    con->solve(*cp, *up, *zp, tol);
    data->outputTpetraVector(u_ptr, "data.txt");
    data->outputTpetraVector(data->getVecF(), "sources.txt");

    data->setVecUd(u_ptr);
    data->zeroRHS();

    /***
         Solve source inversion optimization problem with prescribed sensor locations.
    ***/
    ROL::Ptr<ROL::Step<RealT>> step = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>> status = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algo_tr(step,status,false);
    zp->zero(); // set zero initial guess
    algo_tr.run(*zp, *objReduced, true, *outStream);
    con->solve(*cp, *up, *zp, tol);

    //ROL::Algorithm<RealT> algo_cs("Composite Step",*parlist,false);
    //x.zero(); // set zero initial guess
    //algo_cs.run(x, *cp, *obj, *con, true, *outStream);

    data->outputTpetraVector(u_ptr, "state.txt");
    data->outputTpetraVector(z_ptr, "control.txt");
    data->outputTpetraVector(w_ptr, "weights.txt");
    //std::cout << std::endl << "Sum of random 0/1 entries: " << wp->reduce(ROL::Elementwise::ReductionSum<RealT>()) << std::endl;
    //data->outputTpetraData();

    /***
         Solve OED problem to obtain sparse sensor locations.
    ***/
    std::vector<ROL::Ptr<Tpetra::MultiVector<> > > randvecs_ptr;
    std::vector<ROL::Ptr<ROL::Vector<RealT> > > randvecs;
    std::vector<ROL::Ptr<ROL::Vector<RealT> > > training_models;
    // Get number of random vectors for trace estimation.
    int numRandVecs = parlistOED->sublist("Problem").get("OED Number of random vectors", 1);
    for (int i=0; i<numRandVecs; ++i) {
      randvecs_ptr.push_back( ROL::makePtr<Tpetra::MultiVector<>>(vecmap_z, 1, true));
      //ROL::Ptr<Tpetra::MultiVector<> > rand01_ptr = ROL::makePtr<Tpetra::MultiVector<>>(vecmap_z, 1, true);
      randvecs.push_back( ROL::makePtr<ROL::TpetraMultiVector<RealT>>(randvecs_ptr[i]));
      //ROL::Ptr<ROL::Vector<RealT> > rand01p = ROL::makePtr<ROL::TpetraMultiVector<RealT>>(rand01_ptr);
      randvecs_ptr[i]->randomize();
      randvecs[i]->applyUnary(IsGreaterThan<RealT>(0.0, 1.0, -1.0));
      std::string fname = "rand" + std::to_string(i) + ".txt";
      data->outputTpetraVector(randvecs_ptr[i], fname);
    }

    ROL::Ptr<MyInterfaceOED> oed =
      ROL::makePtr<MyInterfaceOED>(obj, con, up, up, zp, zp, cp, cp, up, up, randvecs, training_models);
    ROL::ExperimentDesignObjective<RealT> objOED(oed, parlistOED);
    ROL::Bounds<RealT> bconOED(wlop, wupp);
    w_ptr->putScalar(1e-2);
    *outStream << std::endl << "Checking OED objective gradient:" << std::endl;
    dwp->scale(1e-2);
    objOED.checkGradient(*wp,*dwp,true,*outStream);
    ROL::Ptr<ROL::Step<RealT>> step_oed = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlistOED);
    ROL::Ptr<ROL::StatusTest<RealT>> status_oed = ROL::makePtr<ROL::StatusTest<RealT>>(*parlistOED);
    ROL::Algorithm<RealT> algo_tr_oed(step_oed,status_oed,false);
    wp->zero(); // set zero initial guess
    algo_tr_oed.run(*wp, objOED, bconOED, true, *outStream);
    data->outputTpetraVector(w_ptr, "weightsOED.txt");

    /***
         Solve source inversion optimization problem with optimal sensor locations.
    ***/
    wp->applyUnary(IsGreaterThan<RealT>(1e-1, 1.0, 0.0));
    RealT numLocOED = wp->reduce(ROL::Elementwise::ReductionSum<RealT>());
    *outStream << std::endl << "Number of nonzero OED locations: " << numLocOED << std::endl;
    zp->zero(); // set zero initial guess
    ROL::Ptr<ROL::Step<RealT>> step_optimal = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>> status_optimal = ROL::makePtr<ROL::StatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algo_tr_optimal(step_optimal,status_optimal,false);
    algo_tr_optimal.run(*zp, *objReduced, true, *outStream);
    con->solve(*cp, *up, *zp, tol);
    data->outputTpetraVector(u_ptr, "stateOED.txt");
    data->outputTpetraVector(z_ptr, "controlOED.txt");
    data->outputTpetraVector(w_ptr, "weightsOED.txt");

    /***
         Solve source inversion optimization problem with random sensor locations.
    ***/
    w_ptr->randomize();
    RealT numLocTotal = w_ptr->getGlobalLength();
    wp->applyUnary(IsGreaterThan<RealT>((1-2*numLocOED/numLocTotal), 1.0, 0.0));
    *outStream << std::endl << "Number of nonzero random locations: " << wp->reduce(ROL::Elementwise::ReductionSum<RealT>()) << std::endl;
    zp->zero(); // set zero initial guess
    ROL::Ptr<ROL::Step<RealT>> step_random = ROL::makePtr<ROL::TrustRegionStep<RealT>>(*parlistOED);
    ROL::Ptr<ROL::StatusTest<RealT>> status_random = ROL::makePtr<ROL::StatusTest<RealT>>(*parlistOED);
    ROL::Algorithm<RealT> algo_tr_random(step_random,status_random,false);
    algo_tr_random.run(*zp, *objReduced, true, *outStream);
    con->solve(*cp, *up, *zp, tol);
    data->outputTpetraVector(u_ptr, "stateRandom.txt");
    data->outputTpetraVector(z_ptr, "controlRandom.txt");
    data->outputTpetraVector(w_ptr, "weightsRandom.txt");

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
