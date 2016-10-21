#ifndef ROL_COLEMANLIMODEL_HPP
#define ROL_COLEMANLIMODEL_HPP

#include "ROL_TrustRegionModel.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Secant.hpp"

namespace ROL {

template<class Real>
class ColemanLiModel : public TrustRegionModel<Real> {
private:
  Teuchos::RCP<BoundConstraint<Real> > bnd_;
  Teuchos::RCP<Secant<Real> > sec_;
  Teuchos::RCP<Vector<Real> > prim_, dual_;
  Teuchos::RCP<Vector<Real> > di_;  // di_i=sqrt(abs(v_i))
  Teuchos::RCP<Vector<Real> > j_;   // "Jacobian" of v

  const bool useSecantPrecond_;
  const bool useSecantHessVec_;

  Elementwise::Multiply<Real> mult_;
  Elementwise::Divide<Real>   div_;

  void invert( Vector<Real> &x ) {
    const Real one(1);
    x.scale(-one);
    x.applyUnary(Elementwise::Shift<Real>(one));
  }

  void initialize(Objective<Real> &obj, BoundConstraint<Real> &bnd,
                  const Vector<Real> &x, const Vector<Real> &g) {
    const Real zero(0), one(1), INF(ROL_INF<Real>()), NINF(ROL_NINF<Real>());
    const int LESS_THAN    = 0;
    const int EQUAL_TO     = 1;
    const int GREATER_THAN = 2;
    
    bnd_ = Teuchos::rcpFromRef(bnd);

    prim_ = x.clone();
    dual_ = g.clone();
    di_   = x.clone();
    j_    = x.clone();

    Teuchos::RCP<Vector<Real> > l =  bnd.getLowerVectorRCP();
    Teuchos::RCP<Vector<Real> > u =  bnd.getUpperVectorRCP();

    Teuchos::RCP<Vector<Real> > m1 = x.clone();
    Teuchos::RCP<Vector<Real> > m2 = x.clone();

    j_->set(g);
    j_->applyUnary(Elementwise::Sign<Real>());

    di_->zero();

    // CASE (i)
    // Mask for negative gradient (m1 is 1 if g is negative and 0 otherwise)
    m1->applyBinary(Elementwise::ValueSet<Real>(zero, LESS_THAN),g);
    // Mask for finite upper bounds (m2 is 1 if u is finite and 0 otherwise)
    m2->applyBinary(Elementwise::ValueSet<Real>(INF, LESS_THAN),*u);
    // Zero out elements of Jacobian with u_i=inf
    j_->applyBinary(mult_,*m2);
    // Mask for g_i < 0 and u_i < inf
    m2->applyBinary(mult_,*m1);
    // prim_i = { u_i-x_i if g_i < 0 and u_i < inf
    //          { 0       otherwise
    prim_->set(*u); prim_->axpy(-one,x);
    prim_->applyBinary(mult_,*m2);
    // Add to D
    di_->plus(*prim_);

    // CASE (iii)
    // Mask for infinite upper bounds
    m2->applyBinary(Elementwise::ValueSet<Real>(INF, EQUAL_TO),*u);
    // Mask for g_i < 0 and u_i = inf
    m2->applyBinary(mult_,*m1);
    // prim_i = { -1 if g_i < 0 and u_i = inf
    //          { 0  otherwise
    prim_->applyUnary(Elementwise::Fill<Real>(-one)); 
    prim_->applyBinary(mult_,*m2);
    // Add to D
    di_->plus(*prim_);

    // CASE (ii)
    invert(*m1);
    // Mask for finite lower bounds
    m2->applyBinary(Elementwise::ValueSet<Real>(NINF, GREATER_THAN),*l);
    // Zero out elements of Jacobian with l_i=-inf
    j_->applyBinary(mult_,*m2);
    m2->applyBinary(mult_,*m1);  
    // prim_i = { x_i-l_i if g_i >= 0 and l_i > -inf
    //          { 0       otherwise
    prim_->set(x); prim_->axpy(-one,*l);
    prim_->applyBinary(mult_,*m2);
    // Add to D
    di_->plus(*prim_);

    // CASE (iv)
    // Mask for infinite lower bounds
    m2->applyBinary(Elementwise::ValueSet<Real>(NINF, EQUAL_TO),*l);
    // Mask for g_i>=0 and l_i=-inf
    m2->applyBinary(mult_,*m1);
    // prim_i = { 1 if g_i >= 0 and l_i = -inf
    //          { 0 otherwise
    prim_->applyUnary(Elementwise::Fill<Real>(one));
    prim_->applyBinary(mult_,*m2);
    // Add to D
    di_->plus(*prim_);
  
    // d_i = { u_i-x_i if g_i <  0, u_i<inf
    //       { -1      if g_i <  0, u_i=inf
    //       { x_i-l_i if g_i >= 0, l_i>-inf
    //       { 1       if g_i >= 0, l_i=-inf 
    di_->applyUnary(Elementwise::AbsoluteValue<Real>());
    di_->applyUnary(Elementwise::SquareRoot<Real>());
  }

 public:

  ColemanLiModel( Objective<Real> &obj, BoundConstraint<Real> &bnd,
                  const Vector<Real> &x, const Vector<Real> &g )
    : TrustRegionModel<Real>::TrustRegionModel(obj,x,g,false),
      sec_(Teuchos::null), useSecantPrecond_(false), useSecantHessVec_(false) {
    initialize(obj,bnd,x,g);
  }

  ColemanLiModel( Objective<Real> &obj, BoundConstraint<Real> &bnd,
                  const Vector<Real> &x, const Vector<Real> &g,
                  const Teuchos::RCP<Secant<Real> > &sec,
                  const bool useSecantPrecond, const bool useSecantHessVec)
    : TrustRegionModel<Real>::TrustRegionModel(obj,x,g,false),
      sec_(sec), useSecantPrecond_(useSecantPrecond), useSecantHessVec_(useSecantHessVec) {
    initialize(obj,bnd,x,g);
  }
 
  // Note that s is the \f$\hat{s}\f$ and \f$\psi\f$ is the $\hat\psi$ from the paper
  Real value( const Vector<Real> &s, Real &tol ) {
    const Teuchos::RCP<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    // Apply Hessian to s
    hessVec(*dual_, s, s, tol);
    dual_->scale(static_cast<Real>(0.5));
    // Form inv(D) * g
    prim_->set(gc->dual());
    prim_->applyBinary(mult_,*di_);
    // Add scaled gradient to Hessian in direction s
    dual_->plus(prim_->dual());
    return dual_->dot(s.dual());    
  }

  void gradient( Vector<Real> &g, const Vector<Real> &s, Real &tol ) {
    const Teuchos::RCP<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    hessVec(g, s, s, tol);
    dualTransform(*dual_,*gc);
    g.plus(*dual_);    
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &s, Real &tol ) {
    const Teuchos::RCP<const Vector<Real> > gc = TrustRegionModel<Real>::getGradient();
    // Build B = inv(D) * Hessian * inv(D)
    primalTransform(*prim_, v);
    if(useSecantHessVec_) {
      sec_->applyB(*dual_, *prim_);
    }
    else {
      const Teuchos::RCP<const Vector<Real> > xc = TrustRegionModel<Real>::getIterate();
      TrustRegionModel<Real>::getObjective()->hessVec(*dual_, *prim_, *xc, tol);   
    }
    dualTransform(hv, *dual_);
    // Build C = diag(g) J
    prim_->set(v);
    prim_->applyBinary(mult_, *j_);
    prim_->applyBinary(mult_, gc->dual());
    hv.plus(prim_->dual()); 
  }
  
  void dualTransform( Vector<Real> &tv, const Vector<Real> &v ) {
    tv.set(v);
    tv.applyBinary(mult_,*di_);
  }

  void primalTransform( Vector<Real> &tiv, const Vector<Real> &v ) { 
    tiv.set(v);
    tiv.applyBinary(mult_,*di_);
  }

  const Teuchos::RCP<BoundConstraint<Real> > getBoundConstraint(void) const {
    return bnd_;
  }

}; // class ColemanLiModel

}

#endif // ROL_COLEMANLIMODEL_HPP
