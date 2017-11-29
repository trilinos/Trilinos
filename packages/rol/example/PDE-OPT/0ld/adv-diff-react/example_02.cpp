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

/*! \file  example_01.cpp
    \brief Solves a source inversion problem governed by the
           advection-diffusion equation.  Performs optimal
           experimental design (OED).
*/

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"

#include "ROL_Algorithm.hpp"
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
  MyInterfaceOED(const ROL::SharedPointer<ROL::Objective_SimOpt<RealT> > &obj,
                 const ROL::SharedPointer<ROL::Constraint_SimOpt<RealT> > &con,
                 const ROL::SharedPointer<ROL::Vector<RealT> > &state,
                 const ROL::SharedPointer<ROL::Vector<RealT> > &stateDual,
                 const ROL::SharedPointer<ROL::Vector<RealT> > &control,
                 const ROL::SharedPointer<ROL::Vector<RealT> > &controlDual,
                 const ROL::SharedPointer<ROL::Vector<RealT> > &constraint,
                 const ROL::SharedPointer<ROL::Vector<RealT> > &constraintDual,
                 const ROL::SharedPointer<ROL::Vector<RealT> > &observation,
                 const ROL::SharedPointer<ROL::Vector<RealT> > &observationDual,
                 const std::vector<ROL::SharedPointer<ROL::Vector<RealT> > > &randvecs,
                 const std::vector<ROL::SharedPointer<ROL::Vector<RealT> > > &training) :
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
  ROL::SharedPointer<std::ostream> outStream;
  Teuchos::oblackholestream bhs; // outputs nothing

  /*** Initialize communicator. ***/
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  ROL::SharedPointer<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  const int myRank = comm->getRank();
  if ((iprint > 0) && (myRank == 0)) {
    outStream = ROL::makeSharedFromRef(std::cout);
  }
  else {
    outStream = ROL::makeSharedFromRef(bhs);
  }

  int errorFlag  = 0;

  // *** Example body.
  try {

    /*** Read in XML input ***/
    std::string filename    = "input.xml";
    std::string filenameOED = "inputOED.xml";
    ROL::SharedPointer<Teuchos::ParameterList> parlist = ROL::makeShared<Teuchos::ParameterList>();
    ROL::SharedPointer<Teuchos::ParameterList> parlistOED = ROL::makeShared<Teuchos::ParameterList>();
    Teuchos::updateParametersFromXmlFile( filename, parlist.ptr() );
    Teuchos::updateParametersFromXmlFile( filenameOED, parlistOED.ptr() );

    /*** Initialize main data structure. ***/
    ROL::SharedPointer<PoissonData<RealT> > data = ROL::makeShared<PoissonData<RealT>>(comm, parlist, outStream);

    // Get random weights parameter.
    RealT fnzw = parlist->sublist("Problem").get("Fraction of nonzero weights", 0.5);
    fnzw = 1.0 - 2.0*fnzw;

    /*** Build vectors and dress them up as ROL vectors. ***/
    ROL::SharedPointer<const Tpetra::Map<> > vecmap_u = data->getMatA()->getDomainMap();
    ROL::SharedPointer<const Tpetra::Map<> > vecmap_z = data->getMatB()->getDomainMap();
    ROL::SharedPointer<const Tpetra::Map<> > vecmap_c = data->getMatA()->getRangeMap();
    ROL::SharedPointer<Tpetra::MultiVector<> > u_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::SharedPointer<Tpetra::MultiVector<> > p_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::SharedPointer<Tpetra::MultiVector<> > w_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::SharedPointer<Tpetra::MultiVector<> > wup_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::SharedPointer<Tpetra::MultiVector<> > wlo_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::SharedPointer<Tpetra::MultiVector<> > z_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_z, 1, true);
    ROL::SharedPointer<Tpetra::MultiVector<> > c_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_c, 1, true);
    ROL::SharedPointer<Tpetra::MultiVector<> > du_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::SharedPointer<Tpetra::MultiVector<> > dw_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_u, 1, true);
    ROL::SharedPointer<Tpetra::MultiVector<> > dz_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_z, 1, true);
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
    ROL::SharedPointer<ROL::Vector<RealT> > up = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(u_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > pp = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(p_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > wp = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(w_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > wlop = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(wlo_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > wupp = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(wup_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > zp = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(z_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > cp = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(c_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > dup = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(du_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > dwp = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(dw_ptr);
    ROL::SharedPointer<ROL::Vector<RealT> > dzp = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(dz_ptr);
    // Create ROL SimOpt vectors.
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> d(dup,dzp);

    /*** Build objective function, constraint and reduced objective function. ***/
    wp->applyUnary(IsGreaterThan<RealT>(fnzw, 1.0, 0.0));
    ROL::SharedPointer<ROL::Objective_SimOpt<RealT> > obj =
      ROL::makeShared<Objective_PDEOPT_Poisson<RealT>>(data, w_ptr, parlist);
    ROL::SharedPointer<ROL::Constraint_SimOpt<RealT> > con =
      ROL::makeShared<EqualityConstraint_PDEOPT_Poisson<RealT>>(data, parlist);
    ROL::SharedPointer<ROL::Objective<RealT> > objReduced =
      ROL::makeShared<ROL::Reduced_Objective_SimOpt<RealT>>(obj, con, up, zp, pp);

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
    ROL::Algorithm<RealT> algo_tr("Trust Region",*parlist,false);
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
    std::vector<ROL::SharedPointer<Tpetra::MultiVector<> > > randvecs_ptr;
    std::vector<ROL::SharedPointer<ROL::Vector<RealT> > > randvecs;
    std::vector<ROL::SharedPointer<ROL::Vector<RealT> > > training_models;
    // Get number of random vectors for trace estimation.
    int numRandVecs = parlistOED->sublist("Problem").get("OED Number of random vectors", 1);
    for (int i=0; i<numRandVecs; ++i) {
      randvecs_ptr.push_back( ROL::makeShared<Tpetra::MultiVector<>>(vecmap_z, 1, true));
      //ROL::SharedPointer<Tpetra::MultiVector<> > rand01_ptr = ROL::makeShared<Tpetra::MultiVector<>>(vecmap_z, 1, true);
      randvecs.push_back( ROL::makeShared<ROL::TpetraMultiVector<RealT>>(randvecs_ptr[i]));
      //ROL::SharedPointer<ROL::Vector<RealT> > rand01p = ROL::makeShared<ROL::TpetraMultiVector<RealT>>(rand01_ptr);
      randvecs_ptr[i]->randomize();
      randvecs[i]->applyUnary(IsGreaterThan<RealT>(0.0, 1.0, -1.0));
      std::string fname = "rand" + std::to_string(i) + ".txt";
      data->outputTpetraVector(randvecs_ptr[i], fname);
    }

    ROL::SharedPointer<MyInterfaceOED> oed =
      ROL::makeShared<MyInterfaceOED>(obj, con, up, up, zp, zp, cp, cp, up, up, randvecs, training_models);
    ROL::ExperimentDesignObjective<RealT> objOED(oed, parlistOED);
    ROL::Bounds<RealT> bconOED(wlop, wupp);
    w_ptr->putScalar(1e-2);
    *outStream << std::endl << "Checking OED objective gradient:" << std::endl;
    dwp->scale(1e-2);
    objOED.checkGradient(*wp,*dwp,true,*outStream);
    ROL::Algorithm<RealT> algo_tr_oed("Trust Region",*parlistOED,false);
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
    ROL::Algorithm<RealT> algo_tr_optimal("Trust Region",*parlist,false);
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
    ROL::Algorithm<RealT> algo_tr_random("Trust Region",*parlist,false);
    algo_tr_random.run(*zp, *objReduced, true, *outStream);
    con->solve(*cp, *up, *zp, tol);
    data->outputTpetraVector(u_ptr, "stateRandom.txt");
    data->outputTpetraVector(z_ptr, "controlRandom.txt");
    data->outputTpetraVector(w_ptr, "weightsRandom.txt");

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
