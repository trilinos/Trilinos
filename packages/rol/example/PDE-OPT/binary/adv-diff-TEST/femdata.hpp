
#ifndef FEMDATA_H
#define FEMDATA_H

#include "ROL_Objective.hpp"
#include "ROL_SimController.hpp"
#include "../../TOOLS/assembler.hpp"
#include "../../TOOLS/solver.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "pde_adv_diff.hpp"
#include "qoi_adv_diff.hpp"
#include "mesh_adv_diff.hpp"

template <class Real>
class FEMdata {
private:
  // Advection Diffusion PDE
  ROL::Ptr<PDE_adv_diff<Real>> pde_;
  ROL::Ptr<Assembler<Real>> assembler_;
  ROL::Ptr<Solver<Real>> solver_;
  ROL::Ptr<Tpetra::CrsMatrix<>> matJ1_, matJ2_, matH_;
  ROL::Ptr<Tpetra::MultiVector<>> vecG_, vecR_, vecJ3_;
  Real c0_;

  bool usePC_;
  int psize_;

  mutable unsigned napJ1_, napJ2_, napH1_, nasse_;

public:
  FEMdata(const ROL::Ptr<const Teuchos::Comm<int>>     &comm,
          ROL::ParameterList                           &list,
          std::ostream                                 &stream = std::cout)
    : napJ1_(0), napJ2_(0), napH1_(0), nasse_(0) {
    usePC_  = list.sublist("Problem").get("Piecewise Constant Controls",false);
    int nx  = list.sublist("Problem").get("Number of X-Cells",1);
    int ny  = list.sublist("Problem").get("Number of Y-Cells",1);
    psize_  = nx*ny;

    assemble(comm,list,stream);
  }

  const ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return assembler_;
  }

  const ROL::Ptr<FE<Real>> getFE(void) const {
    return pde_->getFE();
  }

  const ROL::Ptr<FE<Real>> getFE2(void) const {
    return pde_->getFE2();
  }

  // PDE definitions
  void addPDErhs(ROL::Vector<Real> &rhs) const {
    const Real one(1);
    ROL::Ptr<Tpetra::MultiVector<>> rhsf = getField(rhs);
    rhsf->update(one,*vecR_,one);
  }

  void applyInversePDEJacobian(ROL::Vector<Real>       &x,
                               const ROL::Vector<Real> &b,
                               bool transpose) {
    ROL::Ptr<Tpetra::MultiVector<>>       xf = getField(x);
    ROL::Ptr<const Tpetra::MultiVector<>> bf = getConstField(b);
    solver_->solve(xf,bf,transpose);
    napJ1_++;
  }

  void applyControlJacobian(ROL::Vector<Real>       &Jv,
                            const ROL::Vector<Real> &v,
                            const bool transpose) {
    if (usePC_) {
      if (transpose) {
        ROL::Ptr<std::vector<Real>> Jvp = getParameter(Jv);
        ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
        Teuchos::Array<Real> val(1,0);
        const size_t size = Jvp->size();
        for (size_t i = 0; i < size; ++i) {
          Teuchos::ArrayView<const size_t> col(&i,1);
          vecJ3_->subView(col)->dot(*vf, val.view(0,1));
          (*Jvp)[i] = val[0];
        }        
      }
      else {
        Jv.zero();
        ROL::Ptr<Tpetra::MultiVector<>> Jvf = getField(Jv);
        ROL::Ptr<const std::vector<Real>> vp = getConstParameter(v);
        const Real one(1);
        const size_t size = vp->size();
        for (size_t i = 0; i < size; ++i) {
          Teuchos::ArrayView<const size_t> col(&i,1);
          Jvf->update((*vp)[i],*vecJ3_->subView(col),one);
        }
      }
    }
    else {
      ROL::Ptr<Tpetra::MultiVector<>> Jvf = getField(Jv);
      ROL::Ptr<const Tpetra::MultiVector<>> vf = getConstField(v);
      if (transpose) {
        matJ2_->apply(*vf,*Jvf,Teuchos::TRANS);
      }
      else {
        matJ2_->apply(*vf,*Jvf);
      }
    }
    napJ2_++;
  }

  // Objective definitions
  Real getObjectiveConstant(void) const {
    return c0_;
  }

  void addObjectiveGradient(ROL::Vector<Real> &g) const {
    const Real one(1);
    ROL::Ptr<Tpetra::MultiVector<>> gf = getField(g);
    gf->update(one,*vecG_,one);
  }

  void applyObjectiveHessian(ROL::Vector<Real>       &Hv,
                             const ROL::Vector<Real> &v) const {
    ROL::Ptr<Tpetra::MultiVector<>>       Hvf = getField(Hv);
    ROL::Ptr<const Tpetra::MultiVector<>> vf  = getConstField(v);
    matH_->apply(*vf,*Hvf);
    napH1_++;
  }

  ROL::Ptr<ROL::Vector<Real>> createStateVector(ROL::ParameterList &list) const {
    return ROL::makePtr<PDE_PrimalSimVector<Real>>(assembler_->createStateVector(),pde_,assembler_,list);
  }

  ROL::Ptr<ROL::Vector<Real>> createResidualVector(ROL::ParameterList &list) const {
    return ROL::makePtr<PDE_DualSimVector<Real>>(assembler_->createResidualVector(),pde_,assembler_,list);
  }

  ROL::Ptr<ROL::Vector<Real>> createControlVector(ROL::ParameterList &list) const {
    return ROL::makePtr<PDE_PrimalOptVector<Real>>(assembler_->createControlVector(),pde_,assembler_,list);
  }

  void summarize(std::ostream &stream) const {
    stream << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << "  FEMdata::summarize" << std::endl;
    stream << "    Number of assemblies:                 " << nasse_ << std::endl;
    stream << "    Number of Jacobian_1 sovles:          " << napJ1_ << std::endl;
    stream << "    Number of Jacobian_2 applies:         " << napJ2_ << std::endl;
    stream << "    Number of Hessian_11 applies:         " << napH1_ << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << std::endl;
  }

private:

  // Finite Element Assembly
  void assemble(const ROL::Ptr<const Teuchos::Comm<int>>     &comm,
                ROL::ParameterList                           &list,
                std::ostream                                 &stream) {
    nasse_++;

    pde_ = ROL::makePtr<PDE_adv_diff<Real>>(list);
    ROL::Ptr<MeshManager<Real>>  mesh = ROL::makePtr<MeshManager_adv_diff<Real>>(list);
    assembler_ = ROL::makePtr<Assembler<Real>>(pde_->getFields(),pde_->getFields2(),
                                               mesh,comm,list,stream);
    assembler_->setCellNodes(*pde_);
    ROL::Ptr<QoI<Real>> qoi = ROL::makePtr<QoI_State_Cost_adv_diff<Real>>(pde_->getFE(),list);
    solver_ = ROL::makePtr<Solver<Real>>(list.sublist("Solver"));

    // Create zero vectors
    const Real zero(0);
    ROL::Ptr<Tpetra::MultiVector<>> stateZero   = assembler_->createStateVector();
    ROL::Ptr<Tpetra::MultiVector<>> controlZero = ROL::nullPtr;
    ROL::Ptr<std::vector<Real>>     paramZero   = ROL::nullPtr;
    stateZero->putScalar(zero);
    if (usePC_) {
      paramZero   = ROL::makePtr<std::vector<Real>>(psize_,zero);
    }
    else {
      controlZero = assembler_->createControlVector();
      controlZero->putScalar(zero);
    }

    // Assemble components of the linear PDE
    assembler_->assemblePDEResidual(vecR_,pde_,stateZero,controlZero,paramZero);
    assembler_->assemblePDEJacobian1(matJ1_,pde_,stateZero,controlZero,paramZero);
    if (usePC_) {
      assembler_->assemblePDEJacobian3(vecJ3_,pde_,stateZero,controlZero,paramZero);
    }
    else {
      assembler_->assemblePDEJacobian2(matJ2_,pde_,stateZero,controlZero,paramZero);
    }
    solver_->setA(matJ1_);

    // Assemble components of the quadratic objective function
    c0_ = assembler_->assembleQoIValue(qoi,stateZero);
    assembler_->assembleQoIGradient1(vecG_,qoi,stateZero);
    assembler_->assembleQoIHessian11(matH_,qoi,stateZero);
  }

  ROL::Ptr<const Tpetra::MultiVector<>> getConstField(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const Tpetra::MultiVector<>> xp;
    try {
      xp = dynamic_cast<const ROL::TpetraMultiVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<const ROL::TpetraMultiVector<Real>> xvec
          = dynamic_cast<const PDE_OptVector<Real>&>(x).getField();
        if (xvec == ROL::nullPtr) {
          xp = ROL::nullPtr;
        }
        else {
          xp = xvec->getVector();
        }
      }
      catch (std::exception &ee) {
        xp = ROL::nullPtr;
      }
    }
    return xp;
  }

  ROL::Ptr<Tpetra::MultiVector<>> getField(ROL::Vector<Real> &x) const {
    ROL::Ptr<Tpetra::MultiVector<>> xp;
    try {
      xp = dynamic_cast<ROL::TpetraMultiVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<ROL::TpetraMultiVector<Real>> xvec
          = dynamic_cast<PDE_OptVector<Real>&>(x).getField();
        if ( xvec == ROL::nullPtr ) {
          xp = ROL::nullPtr;
        }
        else {
          xp = xvec->getVector();
        }
      }
      catch (std::exception &ee) {
        xp = ROL::nullPtr;
      }
    }
    return xp;
  }

  ROL::Ptr<const std::vector<Real> > getConstParameter(const ROL::Vector<Real> &x) const {
    ROL::Ptr<const std::vector<Real> > xp;
    try {
      xp = dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<const ROL::StdVector<Real> > xvec
          = dynamic_cast<const PDE_OptVector<Real>&>(x).getParameter();
        if ( xvec == ROL::nullPtr ) {
          xp = ROL::nullPtr;
        }
        else {
          xp = xvec->getVector();
        }
      }
      catch (std::exception &ee) {
        xp = ROL::nullPtr;
      }
    }
    return xp;
  }

  ROL::Ptr<std::vector<Real> > getParameter(ROL::Vector<Real> &x) const {
    ROL::Ptr<std::vector<Real> > xp;
    try {
      xp = dynamic_cast<ROL::StdVector<Real>&>(x).getVector();
    }
    catch (std::exception &e) {
      try {
        ROL::Ptr<ROL::StdVector<Real> > xvec
          = dynamic_cast<PDE_OptVector<Real>&>(x).getParameter();
        if ( xvec == ROL::nullPtr ) {
          xp = ROL::nullPtr;
        }
        else {
          xp = xvec->getVector();
        }
      }
      catch (std::exception &ee) {
        xp = ROL::nullPtr;
      }
    }
    return xp;
  }
}; // class FEMdata

#endif
