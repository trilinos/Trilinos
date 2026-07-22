// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_PRIMALDUALACTIVESETALGORITHM_HPP
#define ROL_TYPEB_PRIMALDUALACTIVESETALGORITHM_HPP

#include "ROL_TypeB_Algorithm.hpp"
#include "ROL_KrylovFactory.hpp"
#include "ROL_SecantFactory.hpp"

/** \class ROL::TypeB::PrimalDualActiveSetAlgorithm
    \brief Provides an interface to run the projected secant algorithm.
*/

namespace ROL {
namespace TypeB {

template<typename Real>
class PrimalDualActiveSetAlgorithm : public TypeB::Algorithm<Real> {
private:
  Ptr<Secant<Real>> secant_; ///< Secant object (used for quasi-Newton)
  ESecant esec_;             ///< Secant type
  std::string secantName_;   ///< Secant name

  Ptr<Krylov<Real>> krylov_; ///< Krylov solver object (used for inexact Newton)
  EKrylov ekv_;              ///< Krylov type
  std::string krylovName_;   ///< Krylov name

  int totalKrylov_; ///< Total number of Krylov iterations per PDAS iteration
  int iterKrylov_;  ///< Number of Krylov iterations (used for inexact Newton)
  int flagKrylov_;  ///< Termination flag for Krylov method (used for inexact Newton)
 
  bool useSecantHessVec_; ///< Whether or not to use to a secant approximation as the Hessian
  bool useSecantPrecond_; ///< Whether or not to use a secant approximation to precondition inexact Newton

  int maxit_;         ///< Maximum number of PDAS steps (default: 10)
  int iter_;          ///< PDAS iteration counter
  int flag_;          ///< PDAS termination flag
  Real stol_;         ///< PDAS minimum step size stopping tolerance (default: 1e-8)
  Real gtol_;         ///< PDAS gradient stopping tolerance (default: 1e-6)
  Real scale_;        ///< Scale for dual variables in the active set, \f$c\f$ (default: 1)
  Real neps_;         ///< \f$\epsilon\f$-active set parameter 
  Real itol_;         ///< Inexact Krylov tolerance
  Real atolKrylov_;   ///< Absolute tolerance for Krylov solve (default: 1e-4)
  Real rtolKrylov_;   ///< Relative tolerance for Krylov solve (default: 1e-2)
  int maxitKrylov_;   ///< Maximum number of Krylov iterations (default: 100)
  bool feasible_;     ///< Flag whether the current iterate is feasible or not

  int verbosity_;
  bool writeHeader_;
  bool hasPoly_;

  class HessianPDAS : public LinearOperator<Real> {
  private:
    const Ptr<Objective<Real>>       obj_;
    const Ptr<BoundConstraint<Real>> bnd_;
    const Ptr<const Vector<Real>>    x_;
    const Ptr<const Vector<Real>>    xlam_;
    const Real                       eps_;
    const Ptr<Secant<Real>>          secant_;
    const bool                       useSecant_;
    const Ptr<Vector<Real>>          pwa_;
  public:
    HessianPDAS(const Ptr<Objective<Real>>       &obj,
                const Ptr<BoundConstraint<Real>> &bnd,
                const Ptr<const Vector<Real>>    &x,
                const Ptr<const Vector<Real>>    &xlam,
                Real                              eps,
                const Ptr<Secant<Real>>          &secant,
                bool                              useSecant,
                const Ptr<Vector<Real>>          &pwa)
      : obj_(obj), bnd_(bnd), x_(x), xlam_(xlam), eps_(eps),
        secant_(secant), useSecant_(useSecant), pwa_(pwa) {}
    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      pwa_->set(v);
      bnd_->pruneActive(*pwa_,*xlam_,eps_);
      if (!useSecant_) obj_->hessVec(Hv,*pwa_,*x_,tol);
      else             secant_->applyB(Hv,*pwa_);
      bnd_->pruneActive(Hv,*xlam_,eps_);
    }
  };

  class PrecondPDAS : public LinearOperator<Real> {
  private:
    const Ptr<Objective<Real>>       obj_;
    const Ptr<BoundConstraint<Real>> bnd_;
    const Ptr<const Vector<Real>>    x_;
    const Ptr<const Vector<Real>>    xlam_;
    const Real                       eps_;
    const Ptr<Secant<Real>>          secant_;
    const bool                       useSecant_;
    const Ptr<Vector<Real>>          dwa_;
  public:
    PrecondPDAS(const Ptr<Objective<Real>>       &obj,
                const Ptr<BoundConstraint<Real>> &bnd,
                const Ptr<const Vector<Real>>    &x,
                const Ptr<const Vector<Real>>    &xlam,
                Real                              eps,
                const Ptr<Secant<Real>>          &secant,
                bool                              useSecant,
                const Ptr<Vector<Real>>          &dwa)
      : obj_(obj), bnd_(bnd), x_(x), xlam_(xlam), eps_(eps),
        secant_(secant), useSecant_(useSecant), dwa_(dwa) {}
    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      Hv.set(v.dual()); 
    }
    void applyInverse(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      dwa_->set(v);
      bnd_->pruneActive(*dwa_,*xlam_,eps_);
      if ( useSecant_ ) secant_->applyH(Hv,*dwa_);
      else              obj_->precond(Hv,*dwa_,*x_,tol);
      bnd_->pruneActive(Hv,*xlam_,eps_);
      dwa_->set(v);
      bnd_->pruneInactive(*dwa_,*xlam_,eps_);
      Hv.plus(dwa_->dual());
    }
  };

  class HessianPDAS_Poly : public LinearOperator<Real> {
  private:
    const Ptr<Objective<Real>>       obj_;
    const Ptr<BoundConstraint<Real>> bnd_;
    const Ptr<Constraint<Real>>      con_;
    const Ptr<const Vector<Real>>    x_;
    const Ptr<const Vector<Real>>    xlam_;
    const Real                       eps_;
    const Ptr<Secant<Real>>          secant_;
    const bool                       useSecant_;
    const Ptr<Vector<Real>>          pwa_, dwa_;
  public:
    HessianPDAS_Poly(const Ptr<Objective<Real>>       &obj,
                     const Ptr<BoundConstraint<Real>> &bnd,
                     const Ptr<Constraint<Real>>      &con,
                     const Ptr<const Vector<Real>>    &x,
                     const Ptr<const Vector<Real>>    &xlam,
                     Real                              eps,
                     const Ptr<Secant<Real>>          &secant,
                     bool                              useSecant,
                     const Ptr<Vector<Real>>          &pwa,
                     const Ptr<Vector<Real>>          &dwa)
      : obj_(obj), bnd_(bnd), con_(con), x_(x), xlam_(xlam), eps_(eps),
        secant_(secant), useSecant_(useSecant), pwa_(pwa), dwa_(dwa) {}
    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      PartitionedVector<Real> &Hvp = dynamic_cast<PartitionedVector<Real>&>(Hv);
      const PartitionedVector<Real> &vp = dynamic_cast<const PartitionedVector<Real>&>(v);
      pwa_->set(*vp.get(0));
      bnd_->pruneActive(*pwa_,*xlam_,eps_);
      if (!useSecant_) obj_->hessVec(*Hvp.get(0),*pwa_,*x_,tol);
      else             secant_->applyB(*Hvp.get(0),*pwa_);
      con_->applyAdjointJacobian(*dwa_,*vp.get(1),*x_,tol);
      Hvp.get(0)->plus(*dwa_);
      bnd_->pruneActive(*Hvp.get(0),*xlam_,eps_);
      con_->applyJacobian(*Hvp.get(1),*pwa_,*x_,tol);
    }
  };

  class PrecondPDAS_Poly : public LinearOperator<Real> {
  private:
    const Ptr<Objective<Real>>       obj_;
    const Ptr<BoundConstraint<Real>> bnd_;
    const Ptr<const Vector<Real>>    x_;
    const Ptr<const Vector<Real>>    xlam_;
    const Real                       eps_;
    const Ptr<Secant<Real>>          secant_;
    const bool                       useSecant_;
    const Ptr<Vector<Real>>          dwa_;
  public:
    PrecondPDAS_Poly(const Ptr<Objective<Real>>       &obj,
                     const Ptr<BoundConstraint<Real>> &bnd,
                     const Ptr<const Vector<Real>>    &x,
                     const Ptr<const Vector<Real>>    &xlam,
                     Real                              eps,
                     const Ptr<Secant<Real>>          &secant,
                     bool                              useSecant,
                     const Ptr<Vector<Real>>          &dwa)
      : obj_(obj), bnd_(bnd), x_(x), xlam_(xlam), eps_(eps),
        secant_(secant), useSecant_(useSecant), dwa_(dwa) {}
    void apply(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      Hv.set(v.dual()); 
    }
    void applyInverse(Vector<Real> &Hv, const Vector<Real> &v, Real &tol) const {
      PartitionedVector<Real> &Hvp = dynamic_cast<PartitionedVector<Real>&>(Hv);
      const PartitionedVector<Real> &vp = dynamic_cast<const PartitionedVector<Real>&>(v);
      dwa_->set(*vp.get(0));
      bnd_->pruneActive(*dwa_,*xlam_,eps_);
      if ( useSecant_ ) secant_->applyH(*Hvp.get(0),*dwa_);
      else              obj_->precond(*Hvp.get(0),*dwa_,*x_,tol);
      bnd_->pruneActive(*Hvp.get(0),*xlam_,eps_);
      dwa_->set(*vp.get(0));
      bnd_->pruneInactive(*dwa_,*xlam_,eps_);
      Hvp.get(0)->plus(dwa_->dual());
      Hvp.get(1)->set(vp.get(1)->dual());
    }
  };

  using TypeB::Algorithm<Real>::status_;
  using TypeB::Algorithm<Real>::state_;
  using TypeB::Algorithm<Real>::proj_;

  void initialize(Vector<Real>          &x,
                  const Vector<Real>    &g,
                  Objective<Real>       &obj,
                  BoundConstraint<Real> &bnd,
                  std::ostream &outStream = std::cout); 

public:

  PrimalDualActiveSetAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant = nullPtr);

  using TypeB::Algorithm<Real>::run;
  void run( Vector<Real>          &x,
            const Vector<Real>    &g, 
            Objective<Real>       &obj,
            BoundConstraint<Real> &bnd,
            std::ostream          &outStream = std::cout) override;

  void writeHeader( std::ostream& os ) const override;

  void writeName( std::ostream& os ) const override;

  void writeOutput( std::ostream& os, const bool write_header = false ) const override;

}; // class ROL::TypeB::PrimalDualActiveSetAlgorithm

} // namespace TypeB
} // namespace ROL

#include "ROL_TypeB_PrimalDualActiveSetAlgorithm_Def.hpp"

#endif
