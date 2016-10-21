#ifndef ROL_COLEMANLIMODEL_HPP
#define ROL_COLEMANLIMODEL_HPP

#include "ROL_TrustRegionModel.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Secant.hpp"

namespace ROL {

template<class Real>
class ColemanLiModel : public TrustRegionModel<Real> {

  template <typename T> using RCP = Teuchos::RCP<T>;

  typedef Objective<Real>                  OBJ;
  typedef BoundConstraint<Real>            BND;
  typedef Secant<Real>                     SEC; 
  typedef Vector<Real>                     V;

  typedef Elementwise::AbsoluteValue<Real> ABS;
  typedef Elementwise::Divide<Real>        DIV;
  typedef Elementwise::Fill<Real>          FILL;
  typedef Elementwise::Multiply<Real>      MULT;
  typedef Elementwise::Power<Real>         POW;
  typedef Elementwise::Shift<Real>         SHIFT;
  typedef Elementwise::ValueSet<Real>      VALSET;

private:

  RCP<OBJ> obj_;
  RCP<BND> bnd_;
  RCP<SEC> sec_;

  RCP<V>   x_;   // Optimization vector
  RCP<V>   g_;   // Gradient
  RCP<V>   di_;  // di_i=sqrt(abs(v_i))
  RCP<V>   j_;   // "Jacobian" of v

  RCP<V>   z_;   // Scratch vector


  Real one_  = Real(1);
  Real zero_ = Real(0);
  const Real inf_;
  const Real ninf_;

  static const int LESS_THAN    = 0;
  static const int EQUAL_TO     = 1;
  static const int GREATER_THAN = 2;

  Real eps_;

  MULT mult_;
  DIV  div_;

  bool useSecant_;

  void applyB( V &Bv, const V &v ) {
    if(useSecant_) {
      sec_->applyB(Bv,v);
    }
    else {
      obj_->hessVec(Bv,v,*x_,eps_);   
    }
  } 

  void invert( V &x ) {
    x.scale(-one_);
    x.applyUnary(SHIFT(one_));
  }

 public:

  ColemanLiModel( OBJ &obj, BND &bnd, const V &x, const V &g ) :
    ColemanLiModel(obj,bnd,Teuchos::null,x,g) {
  }

  ColemanLiModel( OBJ &obj, BND &bnd, const RCP<SEC> &sec, const V &x, const V &g ) :
    sec_(sec),
    x_(x.clone()),
    g_(g.clone()),
    di_(x.clone()),
    j_(x.clone()),
    z_(x.clone()),
    inf_(ROL_INF<Real>()),
    ninf_(ROL_NINF<Real>()) {
    
    obj_ = Teuchos::rcpFromRef(obj);
    bnd_ = Teuchos::rcpFromRef(bnd);

    eps_ = std::sqrt(ROL_EPSILON<Real>());

    useSecant_ =  sec == Teuchos::null ? false : true;

    x_->set(x);

    RCP<V> l =  bnd.getLowerVectorRCP();
    RCP<V> u =  bnd.getUpperVectorRCP();

    RCP<V> m1 = x.clone();
    RCP<V> m2 = x.clone();

    j_->set(g);
    j_->applyUnary(Elementwise::Sign<Real>());

    di_->zero();

    // CASE (i)

    z_->set(*u);
    z_->axpy(-one_,x);

    //  m_i = { 1 if g_i <  0
    //        { 0 if g_i >= 0


    m1->applyBinary(VALSET(zero_, LESS_THAN),*g_);          // Mask for gradient
    m2->applyBinary(VALSET(inf_, LESS_THAN),*u); // Mask for finite upper bounds

    j_->applyBinary(mult_,*m2);  // Zero out elements of Jacobian with u_i=inf

    m2->applyBinary(mult_,*m1);  // Mask for g_i<0 and u_i<inf

    z_->applyBinary(mult_,*m2);
 
    // z_i = { u_i-x_i if g_i <  0 and u_i<inf
    //       { 0       otherwise
 
    di_->plus(*z_);

    // CASE (iii)

    z_->applyUnary(FILL(-one_)); 

    m2->applyBinary(VALSET(inf_, EQUAL_TO),*u); // Mask for infinite upper bounds
    m2->applyBinary(mult_,*m1);  // Mask for g_i<0 and u_i=inf
   
    z_->applyBinary(mult_,*m2);

    di_->plus(*z_);

    // CASE (ii)
  
    z_->set(x);
    z_->axpy(-one_,*l);

    invert(*m1);

    m2->applyBinary(VALSET(ninf_, GREATER_THAN),*l); // Mask for finite lower bounds

    j_->applyBinary(mult_,*m2); // Zero out elements of Jacobian with l_i=-inf

    m2->applyBinary(mult_,*m1);  


    z_->applyBinary(mult_,*m2);

    // z_i = { x_i-l_i if g_i >=  0 and l_i>-inf
    //       { 0       otherwise
 
    di_->plus(*z_);

    // CASE (iv)
    
    z_->applyUnary(FILL(one_));
    
    m2->applyBinary(VALSET(ninf_, EQUAL_TO),*l); // Mask for infinite lower bounds
    m2->applyBinary(mult_,*m1);  // Mask for g_i>=0 and l_i=-inf

    z_->applyBinary(mult_,*m2);

    di_->plus(*z_);

    // d_i = { u_i-x_i if g_i <  0, u_i<inf
    //       { -1      if g_i <  0, u_i=inf
    //       { x_i-l_i if g_i >= 0, l_i>-inf
    //       { 1       if g_i >= 0, l_i=-inf 
  
    di_->applyUnary(ABS());
    di_->applyUnary(POW(0.5));

    // Stored g is ghat
    forwardTransform(*g_,g);

  }
 
  // Note that s is the \f$\hat{s}\f$ and \f$\psi\f$ is the $\hat\psi$ from the paper
  Real value( const V& s, Real &tol ) {

    Real psi = g_->dot(s.dual()); // \f$\hat g^\top \hat s\f$

    forwardTransform(*z_,s);    // z = inv(D)*s
    applyB(*z_,*z_);            // y = B*inv(D)*s
    forwardTransform(*z_,*z_);  // z = inv(D)*B*inv(D)*s   
     
    psi += 0.5*z_->dot(s);          

    z_->set(s);                 // z = s
    z_->applyBinary(mult_,*j_);  // z = J*z
    z_->applyBinary(mult_,*g_);  // z = diag(g)*J*s
    
    psi += 0.5*s.dot(*z_); 

    return psi;    

  }


  void gradient( V& gs, const V& s, Real &tol ) {
    hessVec( gs, s, s, tol );
    gs.plus(*g_);    
  }

  void hessVec( V& hv, const V& v, const V& s, Real &tol ) {
    z_->set(v);
    forwardTransform(hv,*z_);
    applyB(*z_,hv);
    forwardTransform(hv,*z_);

    backwardTransform(*z_,v);
    z_->applyBinary(mult_,*j_);
    z_->applyBinary(mult_,*di_);
    hv.plus(*z_); 
    
  }
  
  void forwardTransform( Vector<Real> &tv, const Vector<Real> &v ) {
    tv.set(v);
    tv.applyBinary(mult_,*di_);
  }

  void backwardTransform( Vector<Real> &tiv, const Vector<Real> &v ) { 
    tiv.set(v);
    tiv.applyBinary(div_,*di_);
  }


}; // class ColemanLiModel


}


#endif // ROL_COLEMANLIMODEL_HPP
