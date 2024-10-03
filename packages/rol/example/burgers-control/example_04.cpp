// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  example_04.cpp
    \brief Shows how to solve a steady Burgers' optimal control problem using
           full-space methods.
*/

#include "ROL_Algorithm.hpp"
#include "ROL_MoreauYosidaPenaltyStep.hpp"
#include "ROL_BoundConstraint_SimOpt.hpp"
#include "ROL_Vector_SimOpt.hpp"
#include "ROL_ParameterList.hpp"

#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include <iostream>
#include <algorithm>

#include "example_04.hpp"

typedef double RealT;
typedef H1VectorPrimal<RealT> PrimalStateVector;
typedef H1VectorDual<RealT> DualStateVector;
typedef L2VectorPrimal<RealT> PrimalControlVector;
typedef L2VectorDual<RealT> DualControlVector;
typedef H1VectorDual<RealT> PrimalConstraintVector;
typedef H1VectorPrimal<RealT> DualConstraintVector;

int main(int argc, char *argv[]) {

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Example body.
  try {
    /*************************************************************************/
    /************* INITIALIZE BURGERS FEM CLASS ******************************/
    /*************************************************************************/
    int nx      = 128;   // Set spatial discretization.
    RealT alpha = 1.e-3; // Set penalty parameter.
    RealT nu    = 1e-2;  // Viscosity parameter.
    RealT nl    = 1.0;   // Nonlinearity parameter (1 = Burgers, 0 = linear).
    RealT u0    = 1.0;   // Dirichlet boundary condition at x=0.
    RealT u1    = 0.0;   // Dirichlet boundary condition at x=1.
    RealT f     = 0.0;   // Constant volumetric force.
    RealT cH1   = 1.0;   // Scale for derivative term in H1 norm.
    RealT cL2   = 0.0;   // Scale for mass term in H1 norm.
    ROL::Ptr<BurgersFEM<RealT> > fem
      = ROL::makePtr<BurgersFEM<RealT>>(nx,nu,nl,u0,u1,f,cH1,cL2);
    fem->test_inverse_mass(*outStream);
    fem->test_inverse_H1(*outStream);
    /*************************************************************************/
    /************* INITIALIZE SIMOPT OBJECTIVE FUNCTION **********************/
    /*************************************************************************/
    ROL::Ptr<std::vector<RealT> > ud_ptr
      = ROL::makePtr<std::vector<RealT>>(nx, 1.);
    ROL::Ptr<ROL::Vector<RealT> > ud
      = ROL::makePtr<L2VectorPrimal<RealT>>(ud_ptr,fem);
    Objective_BurgersControl<RealT> obj(fem,ud,alpha);
    /*************************************************************************/
    /************* INITIALIZE SIMOPT EQUALITY CONSTRAINT *********************/
    /*************************************************************************/
    bool useEChessian = true;
    Constraint_BurgersControl<RealT> con(fem, useEChessian);
    /*************************************************************************/
    /************* INITIALIZE BOUND CONSTRAINTS ******************************/
    /*************************************************************************/
    // INITIALIZE STATE CONSTRAINTS
    std::vector<RealT> Ulo(nx, 0.), Uhi(nx, 1.);
    //std::vector<RealT> Ulo(nx, -1.e8), Uhi(nx, 1.e8);
    ROL::Ptr<ROL::BoundConstraint<RealT> > Ubnd
       = ROL::makePtr<H1BoundConstraint<RealT>>(Ulo,Uhi,fem);
    //Ubnd->deactivate();
    // INITIALIZE CONTROL CONSTRAINTS
    //std::vector<RealT> Zlo(nx+2, -1.e8), Zhi(nx+2, 1.e8);
    std::vector<RealT> Zlo(nx+2,0.), Zhi(nx+2,2.);
    ROL::Ptr<ROL::BoundConstraint<RealT> > Zbnd
      = ROL::makePtr<L2BoundConstraint<RealT>>(Zlo,Zhi,fem);
    //Zbnd->deactivate();
    // INITIALIZE SIMOPT BOUND CONSTRAINTS
    ROL::BoundConstraint_SimOpt<RealT> bnd(Ubnd,Zbnd);
    bnd.deactivate();
    /*************************************************************************/
    /************* INITIALIZE VECTOR STORAGE *********************************/
    /*************************************************************************/
    // INITIALIZE CONTROL VECTORS
    ROL::Ptr<std::vector<RealT> > z_ptr
      = ROL::makePtr<std::vector<RealT>>(nx+2, 0.);
    ROL::Ptr<std::vector<RealT> > zrand_ptr
      = ROL::makePtr<std::vector<RealT>>(nx+2, 1.);
    ROL::Ptr<std::vector<RealT> > gz_ptr
      = ROL::makePtr<std::vector<RealT>>(nx+2, 1.);
    ROL::Ptr<std::vector<RealT> > yz_ptr
      = ROL::makePtr<std::vector<RealT>>(nx+2, 1.);
    for (int i=0; i<nx+2; i++) {
      (*zrand_ptr)[i] = 10.*(RealT)rand()/(RealT)RAND_MAX-5.;
      (*yz_ptr)[i] = 10.*(RealT)rand()/(RealT)RAND_MAX-5.;
    }
    ROL::Ptr<ROL::Vector<RealT> > zp
      = ROL::makePtr<PrimalControlVector>(z_ptr,fem);
    ROL::Ptr<ROL::Vector<RealT> > zrandp
      = ROL::makePtr<PrimalControlVector>(zrand_ptr,fem);
    ROL::Ptr<ROL::Vector<RealT> > gzp
      = ROL::makePtr<DualControlVector>(gz_ptr,fem);
    ROL::Ptr<ROL::Vector<RealT> > yzp
      = ROL::makePtr<PrimalControlVector>(yz_ptr,fem);
    // INITIALIZE STATE VECTORS
    ROL::Ptr<std::vector<RealT> > u_ptr
      = ROL::makePtr<std::vector<RealT>>(nx, 1.);
    ROL::Ptr<std::vector<RealT> > gu_ptr
      = ROL::makePtr<std::vector<RealT>>(nx, 1.);
    ROL::Ptr<std::vector<RealT> > yu_ptr
      = ROL::makePtr<std::vector<RealT>>(nx, 1.);
    for (int i=0; i<nx; i++) {
      (*yu_ptr)[i] = 10.*(RealT)rand()/(RealT)RAND_MAX-5.;
    }
    ROL::Ptr<ROL::Vector<RealT> > up
      = ROL::makePtr<PrimalStateVector>(u_ptr,fem);
    ROL::Ptr<ROL::Vector<RealT> > gup
      = ROL::makePtr<DualStateVector>(gu_ptr,fem);
    ROL::Ptr<ROL::Vector<RealT> > yup
      = ROL::makePtr<PrimalStateVector>(yu_ptr,fem);
    // INITIALIZE CONSTRAINT VECTORS
    ROL::Ptr<std::vector<RealT> > c_ptr
      = ROL::makePtr<std::vector<RealT>>(nx, 1.);
    ROL::Ptr<std::vector<RealT> > l_ptr
      = ROL::makePtr<std::vector<RealT>>(nx, 1.);
    for (int i=0; i<nx; i++) {
      (*l_ptr)[i] = (RealT)rand()/(RealT)RAND_MAX;
    }
    PrimalConstraintVector c(c_ptr,fem);
    DualConstraintVector l(l_ptr,fem);
    // INITIALIZE SIMOPT VECTORS
    ROL::Vector_SimOpt<RealT> x(up,zp);
    ROL::Vector_SimOpt<RealT> g(gup,gzp);
    ROL::Vector_SimOpt<RealT> y(yup,yzp);
    // READ IN XML INPUT
    std::string filename = "input.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );
    
    /*************************************************************************/
    /************* CHECK DERIVATIVES AND CONSISTENCY *************************/
    /*************************************************************************/
    zp->set(*zrandp);
    // CHECK OBJECTIVE DERIVATIVES
    obj.checkGradient(x,g,y,true,*outStream);
    obj.checkHessVec(x,g,y,true,*outStream);
    // CHECK EQUALITY CONSTRAINT DERIVATIVES
    con.checkApplyJacobian(x,y,c,true,*outStream);
    con.checkApplyAdjointHessian(x,*yup,y,g,true,*outStream);
    // CHECK EQUALITY CONSTRAINT CONSISTENCY
    con.checkSolve(*up,*zp,c,true,*outStream);
    con.checkAdjointConsistencyJacobian_1(l,*yup,*up,*zp,true,*outStream);
    con.checkAdjointConsistencyJacobian_2(l,*yzp,*up,*zp,true,*outStream);
    con.checkInverseJacobian_1(c,*yup,*up,*zp,true,*outStream);
    con.checkInverseAdjointJacobian_1(c,*yup,*up,*zp,true,*outStream);
    *outStream << "\n";
    // CHECK PENALTY OBJECTIVE DERIVATIVES
    ROL::Ptr<ROL::Objective<RealT> > obj_ptr = ROL::makePtrFromRef(obj);
    ROL::Ptr<ROL::Constraint<RealT> > con_ptr = ROL::makePtrFromRef(con);
    ROL::Ptr<ROL::BoundConstraint<RealT> > bnd_ptr = ROL::makePtrFromRef(bnd);
    ROL::MoreauYosidaPenalty<RealT> myPen(obj_ptr,bnd_ptr,x,*parlist);
    myPen.checkGradient(x, y, true, *outStream);
    myPen.checkHessVec(x, g, y, true, *outStream);
    ROL::AugmentedLagrangian<RealT> myAugLag(obj_ptr,con_ptr,l,1.,x,c,*parlist);
    myAugLag.checkGradient(x, y, true, *outStream);
    myAugLag.checkHessVec(x, g, y, true, *outStream);
    /*************************************************************************/
    /************* RUN OPTIMIZATION ******************************************/
    /*************************************************************************/
    // SOLVE USING MOREAU-YOSIDA PENALTY
    ROL::Ptr<ROL::Step<RealT>>
      stepMY = ROL::makePtr<ROL::MoreauYosidaPenaltyStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      statusMY = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algoMY(stepMY,statusMY,false);
    zp->set(*zrandp);
    RealT zerotol = std::sqrt(ROL::ROL_EPSILON<RealT>());
    con.solve(c,*up,*zp,zerotol);
    obj.gradient_1(*gup,*up,*zp,zerotol);
    gup->scale(-1.0);
    con.applyInverseAdjointJacobian_1(l,*gup,*up,*zp,zerotol);
    gup->zero(); c.zero();
    algoMY.run(x, g, l, c, myPen, con, bnd, true, *outStream);
    ROL::Ptr<ROL::Vector<RealT> > xMY = x.clone();
    xMY->set(x);
    // SOLVE USING AUGMENTED LAGRANGIAN
    ROL::Ptr<ROL::Step<RealT>>
      stepAL = ROL::makePtr<ROL::AugmentedLagrangianStep<RealT>>(*parlist);
    ROL::Ptr<ROL::StatusTest<RealT>>
      statusAL = ROL::makePtr<ROL::ConstraintStatusTest<RealT>>(*parlist);
    ROL::Algorithm<RealT> algoAL(stepAL,statusAL,false);
    zp->set(*zrandp);
    con.solve(c,*up,*zp,zerotol);
    obj.gradient_1(*gup,*up,*zp,zerotol);
    gup->scale(-1.0);
    con.applyInverseAdjointJacobian_1(l,*gup,*up,*zp,zerotol);
    gup->zero(); c.zero();
    algoAL.run(x, g, l, c, myAugLag, con, bnd, true, *outStream);
    // COMPARE SOLUTIONS
    ROL::Ptr<ROL::Vector<RealT> > err = x.clone();
    err->set(x); err->axpy(-1.,*xMY);
    errorFlag += ((err->norm() > 1.e-7*x.norm()) ? 1 : 0);
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
