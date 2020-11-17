
#ifndef MISFIT_ROBJ_H
#define MISFIT_ROBJ_H

#include "ROL_Objective.hpp"
#include "ROL_SimController.hpp"
#include "femdata.hpp"

template <class Real>
class Misfit_Objective : public ROL::Objective<Real> {
private:
  // Advection Diffusion PDE
  const ROL::Ptr<FEMdata<Real>> fem_;

  // Vector Storage
  ROL::Ptr<ROL::SimController<Real>> stateStore_, adjointStore_;
  ROL::Ptr<ROL::Vector<Real>> state_, adjoint_, state_sens_, adjoint_sens_;
  ROL::Ptr<ROL::Vector<Real>> dualadjoint_;

  bool storage_;

  unsigned nstat_, nadjo_, nsens_, nsadj_;
  unsigned nupda_, nfval_, ngrad_, nhess_, nprec_;

public:
  Misfit_Objective(const ROL::Ptr<FEMdata<Real>> &fem,
                   ROL::ParameterList            &list)
    : fem_(fem), nstat_(0), nadjo_(0), nsens_(0), nsadj_(0),
      nupda_(0), nfval_(0), ngrad_(0), nhess_(0), nprec_(0) {
    stateStore_   = ROL::makePtr<ROL::SimController<Real>>();
    adjointStore_ = ROL::makePtr<ROL::SimController<Real>>();

    // Vector Storage
    state_        = fem_->createStateVector(list);
    adjoint_      = fem_->createStateVector(list);
    state_sens_   = fem_->createStateVector(list);
    adjoint_sens_ = fem_->createStateVector(list);
    dualadjoint_  = fem_->createResidualVector(list);

    storage_  = list.sublist("Problem").get("Use state storage", true);
  }

  const ROL::Ptr<Assembler<Real>> getAssembler(void) const {
    return fem_->getAssembler();
  }

  void update( const ROL::Vector<Real> &z, bool flag = true, int iter = -1 ) {
    nupda_++;
    stateStore_->objectiveUpdate(true);
    adjointStore_->objectiveUpdate(flag);
  }

  void solvePDE(ROL::Vector<Real> &u, const ROL::Vector<Real> &z) {
    update(z,false);
    solve_state_equation(u,z);
  }

  Real value( const ROL::Vector<Real> &z, Real &tol ) {
    nfval_++;
    solve_state_equation(*state_,z);
    fem_->applyObjectiveHessian(*dualadjoint_,*state_); // Hu
    dualadjoint_->scale(static_cast<Real>(0.5));        // 0.5 Hu
    fem_->addObjectiveGradient(*dualadjoint_);          // 0.5 Hu + g
    // 0.5 uHu + gu + c0
    return state_->dot(dualadjoint_->dual()) + fem_->getObjectiveConstant();
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &z, Real &tol ) {
    ngrad_++;
    solve_state_equation(*state_,z);
    solve_adjoint_equation(*adjoint_,*state_);
    fem_->applyControlJacobian(g,*adjoint_,true);
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    nhess_++;
    solve_state_sensitivity(*state_sens_,v);
    solve_adjoint_sensitivity(*adjoint_sens_,*state_sens_);
    fem_->applyControlJacobian(hv,*adjoint_sens_,true);
  }

  /** \brief Apply a reduced Hessian preconditioner.
  */
  virtual void precond( ROL::Vector<Real> &Pv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    nprec_++;
    Pv.set(v.dual());
  }

  void summarize(std::ostream &stream) const {
    stream << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << "  Misfit_Objective::summarize" << std::endl;
    stream << "    Number of calls to update:            " << nupda_ << std::endl;
    stream << "    Number of calls to value:             " << nfval_ << std::endl;
    stream << "    Number of calls to gradient:          " << ngrad_ << std::endl;
    stream << "    Number of calls to hessVec:           " << nhess_ << std::endl;
    stream << "    Number of calls to precond:           " << nprec_ << std::endl;
    stream << "    Number of state solves:               " << nstat_ << std::endl;
    stream << "    Number of adjoint solves:             " << nadjo_ << std::endl;
    stream << "    Number of state sensitivity solves:   " << nsens_ << std::endl;
    stream << "    Number of adjoint sensitivity solves: " << nsadj_ << std::endl;
    stream << std::string(114,'=') << std::endl;
    stream << std::endl;
    fem_->summarize(stream);
  }

private:

  void solve_state_equation(ROL::Vector<Real> &state, const ROL::Vector<Real> &control) {
    bool isComputed = false;
    if (storage_) {
      isComputed = stateStore_->get(state,ROL::Objective<Real>::getParameter());
    }
    if (!isComputed || !storage_) {
      fem_->applyControlJacobian(*dualadjoint_,control,false);  // Bz
      fem_->addPDErhs(*dualadjoint_);                           // Bz + f
      fem_->applyInversePDEJacobian(state,*dualadjoint_,false); // inv(A)(Bz + f)
      state.scale(static_cast<Real>(-1));                       // -inv(A)(Bz + f)
      nstat_++;
      if (storage_) {
        stateStore_->set(state,ROL::Objective<Real>::getParameter());
      }
    }
  }

  void solve_adjoint_equation(ROL::Vector<Real> &adjoint, const ROL::Vector<Real> &state) {
    bool isComputed = false;
    if (storage_) {
      isComputed = adjointStore_->get(adjoint,ROL::Objective<Real>::getParameter());
    }
    if (!isComputed || !storage_) {
      fem_->applyObjectiveHessian(*dualadjoint_,state);          // Hu
      fem_->addObjectiveGradient(*dualadjoint_);                 // Hu + g
      fem_->applyInversePDEJacobian(adjoint,*dualadjoint_,true); // inv(A')(Hu + g)
      adjoint.scale(static_cast<Real>(-1));                      // -inv(A')(Hu + g)
      nadjo_++;
      if (storage_) {
        adjointStore_->set(adjoint,ROL::Objective<Real>::getParameter());
      }
    }
  }

  void solve_state_sensitivity(ROL::Vector<Real> &state_sens, const ROL::Vector<Real> &v) {
    fem_->applyControlJacobian(*dualadjoint_,v,false);               // Bv
    fem_->applyInversePDEJacobian(state_sens,*dualadjoint_,false); // inv(A)Bv
    nsens_++;
  }

  void solve_adjoint_sensitivity(ROL::Vector<Real> &adjoint_sens, const ROL::Vector<Real> &state_sens) {
    fem_->applyObjectiveHessian(*dualadjoint_,state_sens);          // Hv
    fem_->applyInversePDEJacobian(adjoint_sens,*dualadjoint_,true); // inv(A')Hv
    nsadj_++;
  }
}; // class Misfit_Objective

#endif
