// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* \file test_11.hpp
    \brief Verify that the Coleman-Li Model function produces the
           correct values for a test problem
 
    \f[ \min_x f(x) = \frac{1}{2}(x_1^2+2x_2^2),\quad x_1 \geq 1,\; x_2 <= -1 \f]

    The gradient is
    \f[  \nabla f(x) = (x_1,2 x_2 ) \f]
    
    and the Hessian is 
    \[f \nabla^2 f(x) = \begin{pmatrix} 1 & 0 \\ 0 & 2 \end{pmatrix} \f]

    The minimizer is \f$ x^\ast = (1,-1)\f$

    For feasible \f$x\f$, the Coleman-Li quantities of interest are 

    \f[ v_1 = x_1 - 1,\quad v_2 = -1 - x_2 \f]
    \f[ D^{-1} = \begin{pmatrix} \sqrt{|x_1-1|} & 0 \\ 0 & \sqrt{|x_2+1|} \end{pmatrix} \f]
    \f[ J=\begin{pmatrix} 1 & -1 \end{pmatrix} \f]
    \f[ \hat M_k = \begin{pmatrix} |x_1-1|^2+|x_1| & 0  \\
                                   0 & |x_2+1|^2+2|x_2| \end{pmatrix} \f]
    

*/

#include "ROL_Objective.hpp"
#include "ROL_StdVector.hpp"

template<class Real> 
class CLTestObjective : public ROL::Objective<Real> {


public:

  Real value( const ROL::Vector<Real> &x, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > xp = 
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
      return 0.5*((*xp)[0]*(*xp)[0] + 2*(*xp)[1]*(*xp)[1]);
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x, Real &tol ) {
    ROL::Ptr<std::vector<Real> > gp = 
      dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
    ROL::Ptr<const std::vector<Real> > xp = 
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    (*gp)[0] =   (*xp)[0];
    (*gp)[1] = 2*(*xp)[1];
  }

  void hessVec( ROL::Vector<Real> &hv, 
                const ROL::Vector<Real> &v, 
                const ROL::Vector<Real> &x,
                Real &tol ) {
    ROL::Ptr<std::vector<Real> > hvp = 
      dynamic_cast<ROL::StdVector<Real>&>(hv).getVector();
    ROL::Ptr<const std::vector<Real> > vp = 
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();
    (*hvp)[0] =   (*vp)[0];
    (*hvp)[1] = 2*(*vp)[1];
  }

}; // CLTestObjective

template<class Real>
class CLExactModel : public ROL::Objective<Real> {

ROL::Ptr<std::vector<Real> > x_;
const ROL::Ptr<const std::vector<Real> > l_;
const ROL::Ptr<const std::vector<Real> > u_;
ROL::Ptr<std::vector<Real> > g_;
ROL::Ptr<std::vector<Real> > di_;
ROL::Ptr<std::vector<Real> > j_;
ROL::Ptr<ROL::Objective<Real> > obj_;

public:

  CLExactModel( ROL::Ptr<std::vector<Real> > &xp,
                const ROL::Ptr<const std::vector<Real> > &lp,
                const ROL::Ptr<const std::vector<Real> > &up ) : 
                x_(xp), l_(lp), u_(up) { 
    g_  = ROL::makePtr<std::vector<double>(x_->size>());
    di_ = ROL::makePtr<std::vector<double>(x_->size>());
    j_  = ROL::makePtr<std::vector<double>(x_->size>());
 

    obj_ = ROL::makePtr<CLTestObjective<Real>>();

    ROL::StdVector<Real> g(g_);
    ROL::StdVector<Real> x(x_);
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    obj_->gradient(g,x,tol);

    std::vector<Real> v(2);
    
    for(int i=0; i<2;++i) {
      (*j_)[i] = 0;
      // Case (i)
      if( (*g_)[i]<0 && (*u_)[i] < ROL::ROL_INF<Real>() ) {
        v[i] = (*u_)[i]-(*x_)[i];
        (*j_)[i] = -1;
      }
      // Case (ii)
      else if( (*g_)[i]>=0 && (*l_)[i] > ROL::ROL_NINF<Real>() ) {
         v[i] = (*x_)[i] - (*l_)[i]; 
         (*j_)[i] = 1;
      }
      // Case (iii)
      else if( (*g_)[i]<=0 && (*u_)[i] == ROL::ROL_INF<Real>() ) {
        v[i] = -1;
      }
      // Case (iv)
      else {
        v[i] = 1;
      }
      (*di_)[i] = std::sqrt(std::abs(v[i]));
    }  

 
    std::cout << "x[0]  = " << (*x_)[0] << std::endl;
    std::cout << "x[1]  = " << (*x_)[1] << std::endl;
    std::cout << "g[0]  = " << (*g_)[0] << std::endl;
    std::cout << "g[0]  = " << (*g_)[1] << std::endl;
    std::cout << "di[0] = " << (*di_)[0] << std::endl;
    std::cout << "di[1] = " << (*di_)[1] << std::endl;

  } 

  void update( const ROL::Vector<Real> &x, bool flag = true, int iter=-1 ) {
    ROL::Ptr<const std::vector<Real> > xc = 
      dynamic_cast<const ROL::StdVector<Real>&>(x).getVector();
    (*x_)[0]  = (*xc)[0];
    (*x_)[1]  = (*xc)[1];

    std::vector<Real> v(2);
    ROL::StdVector<Real> g(g_);
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    obj_->gradient(g,x,tol);

    for(int i=0; i<2;++i) {
      (*j_)[i] = 0;
      // Case (i)
      if( (*g_)[i]<0 && (*u_)[i] < ROL::ROL_INF<Real>() ) {
        v[i] = (*u_)[i]-(*x_)[i];
        (*j_)[i] = -1;
      }
      // Case (ii)
      else if( (*g_)[i]>=0 && (*l_)[i] > ROL::ROL_NINF<Real>() ) {
         v[i] = (*x_)[i] - (*l_)[i]; 
         (*j_)[i] = 1;
      }
      // Case (iii)
      else if( (*g_)[i]<=0 && (*u_)[i] == ROL::ROL_INF<Real>() ) {
        v[i] = -1;
      }
      // Case (iv)
      else {
        v[i] = 1;
      }
      (*di_)[i] = std::sqrt(std::abs(v[i]));
    }  

    std::cout << "x[0]  = " << (*x_)[0] << std::endl;
    std::cout << "x[1]  = " << (*x_)[1] << std::endl;
    std::cout << "g[0]  = " << (*g_)[0] << std::endl;
    std::cout << "g[0]  = " << (*g_)[1] << std::endl;
    std::cout << "di[0] = " << (*di_)[0] << std::endl;
    std::cout << "di[1] = " << (*di_)[1] << std::endl;
  }

  Real value( const ROL::Vector<Real> &s, Real &tol ) {
    ROL::Ptr<const std::vector<Real> > sp = 
      dynamic_cast<const ROL::StdVector<Real>&>(s).getVector();

    ROL::Ptr<ROL::Vector<Real> > y = s.clone();
    hessVec(*y,s,s,tol);  
    Real result = 0.5*y->dot(s);
    result += (*di_)[0]*(*g_)[0]*(*sp)[0];
    result += (*di_)[1]*(*g_)[1]*(*sp)[1];
    return result; 
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &s, Real &tol ) {
    ROL::Ptr<std::vector<Real> > gp = 
      dynamic_cast<ROL::StdVector<Real>&>(g).getVector();
    hessVec(g,s,s,tol);
     
    (*gp)[0] += (*di_)[0]*(*g_)[0];
    (*gp)[1] += (*di_)[1]*(*g_)[1];
  }

  void hessVec( ROL::Vector<Real> &hv, 
                const ROL::Vector<Real> &v, 
                const ROL::Vector<Real> &s,
                Real &tol ) {

    ROL::Ptr<std::vector<Real> > hvp = 
      dynamic_cast<ROL::StdVector<Real>&>(hv).getVector();
    ROL::Ptr<const std::vector<Real> > vp = 
      dynamic_cast<const ROL::StdVector<Real>&>(v).getVector();

    obj_->hessVec(hv,v,s,tol);

    for(int i=0; i<2; ++i) {
      (*hvp)[i] *=  (*di_)[i]*(*di_)[i];
      (*hvp)[i] += (*g_)[i]*(*j_)[i]*(*vp)[i];
    }
    
  }


}; // CLExactModel 




