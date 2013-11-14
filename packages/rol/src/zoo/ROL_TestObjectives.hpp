// @HEADER
// ************************************************************************
//
// Questions? Contact Denis Ridzal (dridzal@sandia.gov), or
//                    Drew Kouri   (dpkouri@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions of test objective functions.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_TESTOBJECTIVES_HPP
#define ROL_TESTOBJECTIVES_HPP

#include "ROL_Types.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"

namespace ROL {

  /** \brief Rosenbrock's function.
   */
  template<class Real>
  class Objective_Rosenbrock : public ROL::Objective<Real> {
  private:
    Real alpha_;

  public:
    Objective_Rosenbrock(Real alpha = 100.0) : alpha_(alpha) {}

    Real value( const ROL::Vector<Real> &x ) {
      ROL::StdVector<Real> & ex =
        Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast <ROL::Vector<Real> &>(x));
      Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();

      int n = xp->size();
      Real val = 0;
      for( int i=0; i<n/2; i++ ) {
        val += alpha_ * pow(pow((*xp)[2*i],2) - (*xp)[2*i+1], 2);
        val += pow((*xp)[2*i] - 1.0, 2);
      }
 
      return val;
    }

    void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > gp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(g)).getVector());

      int n = xp->size();
      for( int i=0; i<n/2; i++ ) {
        (*gp)[2*i]   =  4.0*alpha_*(pow((*xp)[2*i],2) - (*xp)[2*i+1])*(*xp)[2*i] + 2.0*((*xp)[2*i]-1.0);
        (*gp)[2*i+1] = -2.0*alpha_*(pow((*xp)[2*i],2) - (*xp)[2*i+1]);
      }
    }

    void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());

      int n = xp->size();
      for( int i=0; i<n/2; i++ ) {
        Real h11 = 4.0*alpha_*(3.0*pow((*xp)[2*i],2)-(*xp)[2*i+1]) + 2.0;
        Real h12 = -4.0*alpha_*(*xp)[2*i];
        Real h22 = 2.0*alpha_;

        (*hvp)[2*i]   = h11*(*vp)[2*i] + h12*(*vp)[2*i+1];
        (*hvp)[2*i+1] = h12*(*vp)[2*i] + h22*(*vp)[2*i+1];
      }
    }

    void invHessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());

      int n = xp->size();
      for( int i=0; i<n/2; i++ ) {
        Real h11 = 4.0*alpha_*(3.0*pow((*xp)[2*i],2)-(*xp)[2*i+1]) + 2.0;
        Real h12 = -4.0*alpha_*(*xp)[2*i];
        Real h22 = 2.0*alpha_;
 
        (*hvp)[2*i]   = (1.0/(h11*h22-h12*h12))*( h22*(*vp)[2*i] - h12*(*vp)[2*i+1]);
        (*hvp)[2*i+1] = (1.0/(h11*h22-h12*h12))*(-h12*(*vp)[2*i] + h11*(*vp)[2*i+1]);
      }
    }
  };

  /** \brief Freudenstein and Roth's function.
   */
  template<class Real>
  class Objective_FreudensteinRoth : public ROL::Objective<Real> {
  public:
    Objective_FreudensteinRoth() {}

    Real value( const ROL::Vector<Real> &x ) {
      ROL::StdVector<Real> & ex =
        Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast <ROL::Vector<Real> &>(x));
      Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();

      Real f1 = -13.0 + (*xp)[0] + ((5.0-(*xp)[1])*(*xp)[1] - 2.0)*(*xp)[1];
      Real f2 = -29.0 + (*xp)[0] + (((*xp)[1]+1.0)*(*xp)[1] - 14.0)*(*xp)[1];
 
      return f1*f1+f2*f2;
    }

    void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > gp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(g)).getVector());
  
      Real f1 = -13.0 + (*xp)[0] + ((5.0-(*xp)[1])*(*xp)[1] - 2.0)*(*xp)[1];
      Real f2 = -29.0 + (*xp)[0] + (((*xp)[1]+1.0)*(*xp)[1] - 14.0)*(*xp)[1];
 
      Real f11 = 1.0;
      Real f12 = 10.0*(*xp)[1] - 3.0*(*xp)[1]*(*xp)[1] - 2.0;
      Real f21 = 1.0;
      Real f22 = 3.0*(*xp)[1]*(*xp)[1] + 2.0*(*xp)[1] - 14.0;  

      (*gp)[0] = 2.0*(f11*f1 + f21*f2);
      (*gp)[1] = 2.0*(f12*f1 + f22*f2);
    }

    void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());
  
      Real f1 = -13.0 + (*xp)[0] + ((5.0-(*xp)[1])*(*xp)[1] - 2.0)*(*xp)[1];
      Real f2 = -29.0 + (*xp)[0] + (((*xp)[1]+1.0)*(*xp)[1] - 14.0)*(*xp)[1];
 
      Real f11 = 1.0;
      Real f12 = 10.0*(*xp)[1] - 3.0*(*xp)[1]*(*xp)[1] - 2.0;
      Real f21 = 1.0;
      Real f22 = 3.0*(*xp)[1]*(*xp)[1] + 2.0*(*xp)[1] - 14.0;  

      Real f122 = 10.0 - 6.0*(*xp)[1];
      Real f222 = 6.0*(*xp)[1] + 2.0;

      Real h11 = 2.0*(f11*f11) + 2.0*(f21*f21);
      Real h12 = 2.0*(f12*f11) + 2.0*(f22*f21);
      Real h22 = 2.0*(f122*f1 + f12*f12) + 2.0*(f222*f2 + f22*f22);

      (*hvp)[0] = h11*(*vp)[0] + h12*(*vp)[1];
      (*hvp)[1] = h12*(*vp)[0] + h22*(*vp)[1];
    }

    void invHessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());
  
      Real f1 = -13.0 + (*xp)[0] + ((5.0-(*xp)[1])*(*xp)[1] - 2.0)*(*xp)[1];
      Real f2 = -29.0 + (*xp)[0] + (((*xp)[1]+1.0)*(*xp)[1] - 14.0)*(*xp)[1];
 
      Real f11 = 1.0;
      Real f12 = 10.0*(*xp)[1] - 3.0*(*xp)[1]*(*xp)[1] - 2.0;
      Real f21 = 1.0;
      Real f22 = 3.0*(*xp)[1]*(*xp)[1] + 2.0*(*xp)[1] - 14.0;  

      Real f122 = 10.0 - 6.0*(*xp)[1];
      Real f222 = 6.0*(*xp)[1] + 2.0;

      Real h11 = 2.0*(f11*f11) + 2.0*(f21*f21);
      Real h12 = 2.0*(f12*f11) + 2.0*(f22*f21);
      Real h22 = 2.0*(f122*f1 + f12*f12) + 2.0*(f222*f2 + f22*f22);

      (*hvp)[0] = (1.0/(h11*h22-h12*h12))*( h22*(*vp)[0] - h12*(*vp)[1]);
      (*hvp)[1] = (1.0/(h11*h22-h12*h12))*(-h12*(*vp)[0] + h11*(*vp)[1]);
    }
  };

  /** \brief Powell's badly scaled function.
   */
  template<class Real>
  class Objective_Powell : public ROL::Objective<Real> {
  public:
    Objective_Powell() {}

    Real value( const ROL::Vector<Real> &x ) {
      ROL::StdVector<Real> & ex =
        Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast <ROL::Vector<Real> &>(x));
      Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();

      Real f1 = 1.e4*(*xp)[0]*(*xp)[1] - 1.0;
      Real f2 = std::exp(-(*xp)[0]) + std::exp(-(*xp)[1]) - 1.0001;
 
      return f1*f1+f2*f2;
    }

    void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > gp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(g)).getVector());

      Real f1 = 1.e4*(*xp)[0]*(*xp)[1] - 1.0;
      Real f2 = std::exp(-(*xp)[0]) + std::exp(-(*xp)[1]) - 1.0001;
 
      Real f11 = 1.e4*(*xp)[1];
      Real f12 = 1.e4*(*xp)[0];
      Real f21 = -std::exp(-(*xp)[0]);
      Real f22 = -std::exp(-(*xp)[1]);  

      (*gp)[0] = 2.0*(f11*f1 + f21*f2);
      (*gp)[1] = 2.0*(f12*f1 + f22*f2);
    }

    void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());
  
      Real f1 = 1.e4*(*xp)[0]*(*xp)[1] - 1.0;
      Real f2 = std::exp(-(*xp)[0]) + std::exp(-(*xp)[1]) - 1.0001;
 
      Real f11 = 1.e4*(*xp)[1];
      Real f12 = 1.e4*(*xp)[0];
      Real f21 = -std::exp(-(*xp)[0]);
      Real f22 = -std::exp(-(*xp)[1]);  

      Real f111 = 0.0;
      Real f112 = 1.e4;
      Real f121 = 1.e4;
      Real f122 = 0.0;
      Real f211 = std::exp(-(*xp)[0]);
      Real f212 = 0.0;
      Real f221 = 0.0;
      Real f222 = std::exp(-(*xp)[1]);

      Real h11 = 2.0*(f111*f1 + f11*f11) + 2.0*(f211*f2 + f21*f21);
      Real h12 = 2.0*(f112*f1 + f11*f12) + 2.0*(f212*f2 + f21*f22);
      Real h21 = 2.0*(f121*f1 + f21*f11) + 2.0*(f221*f2 + f22*f21);
      Real h22 = 2.0*(f122*f1 + f12*f12) + 2.0*(f222*f2 + f22*f22);

      (*hvp)[0] = h11*(*vp)[0] + h12*(*vp)[1];
      (*hvp)[1] = h21*(*vp)[0] + h22*(*vp)[1];
    }

    void invHessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());
  
      Real f1 = 1.e4*(*xp)[0]*(*xp)[1] - 1.0;
      Real f2 = std::exp(-(*xp)[0]) + std::exp(-(*xp)[1]) - 1.0001;
 
      Real f11 = 1.e4*(*xp)[1];
      Real f12 = 1.e4*(*xp)[0];
      Real f21 = -std::exp(-(*xp)[0]);
      Real f22 = -std::exp(-(*xp)[1]);  

      Real f111 = 0.0;
      Real f112 = 1.e4;
      Real f121 = 1.e4;
      Real f122 = 0.0;
      Real f211 = std::exp(-(*xp)[0]);
      Real f212 = 0.0;
      Real f221 = 0.0;
      Real f222 = std::exp(-(*xp)[1]);

      Real h11 = 2.0*(f111*f1 + f11*f11) + 2.0*(f211*f2 + f21*f21);
      Real h12 = 2.0*(f112*f1 + f11*f12) + 2.0*(f212*f2 + f21*f22);
      Real h21 = 2.0*(f121*f1 + f21*f11) + 2.0*(f221*f2 + f22*f21);
      Real h22 = 2.0*(f122*f1 + f12*f12) + 2.0*(f222*f2 + f22*f22);

      (*hvp)[0] = (1.0/(h11*h22-h12*h21))*( h22*(*vp)[0] - h21*(*vp)[1]);
      (*hvp)[1] = (1.0/(h11*h22-h12*h21))*(-h12*(*vp)[0] + h11*(*vp)[1]);
    }
  };

  /** \brief Sum of squares function. 
   */
  template<class Real>
  class Objective_SumOfSquares : public ROL::Objective<Real> {
  public:
    Real value( const ROL::Vector<Real> &x ) {
      ROL::StdVector<Real> & ex =
        Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast <ROL::Vector<Real> &>(x));
      Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();

      int n = xp->size();
      Real val = 0;
      for (int i=0; i<n; i++) {
        val += pow((*xp)[i], 2);
      }

      return val;
    }

    void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > gp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(g)).getVector());

      int n = xp->size();
      for( int i=0; i<n; i++ ) {
        (*gp)[i] = 2.0*(*xp)[i];
      }
    }

    void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());

      int n = xp->size();
      for( int i=0; i<n; i++ ) {
        (*hvp)[i] = 2.0*(*vp)[i];
      }
    }

    void invHessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());

      int n = xp->size();
      for( int i=0; i<n; i++ ) {
        (*hvp)[i] = 0.5*(*vp)[i];
      }
    }
  };

  /** \brief Least squares function.
   */
  template<class Real>
  class Objective_LeastSquares : public ROL::Objective<Real> {
  public:
    Real value( const ROL::Vector<Real> &x ) {
      ROL::StdVector<Real> & ex =
        Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast <ROL::Vector<Real> &>(x));
      Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();

      int n    = xp->size();
      Real h   = 1.0/((Real)n+1.0);
      Real val = 0.0;
      Real res = 0.0;
      for (int i=0; i<n; i++) {
        if ( i == 0 ) {
          res = 2.0*h*(5.0/6.0) + 1.0/h*((*xp)[i+1]-2.0*(*xp)[i]);
        }
        else if ( i == n-1 ) {
          res = 2.0*h*(5.0/6.0) + 1.0/h*((*xp)[i-1]-2.0*(*xp)[i]);
        }
        else {
          res = 2.0*h + 1.0/h*((*xp)[i-1]-2.0*(*xp)[i]+(*xp)[i+1]);
        }
        val += 0.5*res*res;
      }

      return val;
   }

    void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > gp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(g)).getVector());

      int n  = xp->size();
      Real h = 1.0/((Real)n+1.0);
      std::vector<Real> res(n,0.0);
      for (int i=0; i<n; i++) {
        if ( i == 0 ) {
          res[i] = 2.0*h*(5.0/6.0) + 1.0/h*((*xp)[i+1]-2.0*(*xp)[i]);
        }
        else if ( i == n-1 ) {
          res[i] = 2.0*h*(5.0/6.0) + 1.0/h*((*xp)[i-1]-2.0*(*xp)[i]);
        }
        else {
          res[i] = 2.0*h + 1.0/h*((*xp)[i-1]-2.0*(*xp)[i]+(*xp)[i+1]);
        }
      }

      for (int i=0; i<n; i++) {
        if ( i == 0 ) {
          (*gp)[i] = 1.0/h*(res[i+1]-2.0*res[i]);
        }
        else if ( i == n-1 ) {
          (*gp)[i] = 1.0/h*(res[i-1]-2.0*res[i]);
        }
        else {
          (*gp)[i] = 1.0/h*(res[i-1]-2.0*res[i]+res[i+1]);
        }
      }
    }

    void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &x ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<ROL::StdVector<Real> >(const_cast<ROL::Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(hv)).getVector());

      int n  = xp->size();
      Real h = 1.0/((Real)n+1.0);
      std::vector<Real> res(n,0.0);
      for (int i=0; i<n; i++) {
        if ( i == 0 ) {
          res[i] = 1.0/h*((*vp)[i+1]-2.0*(*vp)[i]);
        }
        else if ( i == n-1 ) {
          res[i] = 1.0/h*((*vp)[i-1]-2.0*(*vp)[i]);
        }
        else {
          res[i] = 1.0/h*((*vp)[i-1]-2.0*(*vp)[i]+(*vp)[i+1]);
        }
      }

      for (int i=0; i<n; i++) {
        if ( i == 0 ) {
          (*hvp)[i] = 1.0/h*(res[i+1]-2.0*res[i]);
        }
        else if ( i == n-1 ) {
          (*hvp)[i] = 1.0/h*(res[i-1]-2.0*res[i]);
        }
        else {
          (*hvp)[i] = 1.0/h*(res[i-1]-2.0*res[i]+res[i+1]);
        }
      }
    }
  };

  template<class Real>
  void getTestObjectives( Teuchos::RCP<Objective<Real> > &obj, Vector<Real> &x0, Vector<Real> &x, 
                          const ETestObjectives test ) {
    // Cast Initial Guess and Solution Vectors
    Teuchos::RCP<std::vector<Real> > x0p =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(x0)).getVector());
    Teuchos::RCP<std::vector<Real> > xp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<ROL::StdVector<Real> >(x)).getVector());
    int n = xp->size();

    if ( test == TESTOBJECTIVES_ROSENBROCK ) {
      // Resize Vectors
      n = 100;
      x0p->resize(n);
      xp->resize(n);
      // Instantiate Objective Function
      obj = Teuchos::rcp( new Objective_Rosenbrock<Real> );
      // Get Initial Guess
      for (int i=0; i<n/2; i++) {
        (*x0p)[2*i]   = -1.2;
        (*x0p)[2*i+1] =  1.0;
      }
      // Get Solution
      for( int i=0; i<n; i++ ) {
        (*xp)[i] = 1.0;
      }
    }
    else if ( test == TESTOBJECTIVES_FREUDENSTEINANDROTH ) {
      // Resize Vectors
      n = 2;
      x0p->resize(n);
      xp->resize(n);
      // Instantiate Objective Function
      obj = Teuchos::rcp( new Objective_FreudensteinRoth<Real> );
      // Get Initial Guess
      (*x0p)[0] =  0.5;
      (*x0p)[1] = -2.0;
      // Get Solution
      (*xp)[0] = 0.0;
      (*xp)[1] = 0.0;
    }
    else if ( test == TESTOBJECTIVES_POWELL ) {
      // Resize Vectors
      n = 2;
      x0p->resize(n);
      xp->resize(n);
      // Instantiate Objective Function
      obj = Teuchos::rcp( new Objective_Powell<Real> );
      // Get Initial Guess
      (*x0p)[0] = 0.0;
      (*x0p)[1] = 1.0;
      // Get Solution
      (*xp)[0] = 1.098*1.e-5;
      (*xp)[1] = 9.106;
    }
    else if ( test == TESTOBJECTIVES_SUMOFSQUARES ) {
      // Resize Vectors
      n = 100;
      x0p->resize(n);
      xp->resize(n);
      // Instantiate Objective Function
      obj = Teuchos::rcp( new Objective_SumOfSquares<Real> );
      // Get Initial Guess
      for (int i=0; i<n; i++) {
        (*x0p)[i] = 1.0;
      }
      // Get Solution
      for( int i=0; i<n; i++ ) {
        (*xp)[i] = 0.0;
      }
    }
    else if ( test == TESTOBJECTIVES_LEASTSQUARES ) {
      // Resize Vectors
      n = 32;
      x0p->resize(n);
      xp->resize(n);
      // Instantiate Objective Function
      obj = Teuchos::rcp( new Objective_LeastSquares<Real> );
      // Get Initial Guess
      for (int i=0; i<n; i++) {
        (*x0p)[i] = 0.0;
      }
      // Get Solution
      Real h  = 1.0/((Real)n+1.0);
      Real pt = 0.0;
      for( int i=0; i<n; i++ ) {
        pt = (Real)(i+1)*h;
        (*xp)[i] = pt*(1.0-pt);
      }
    }
  }

} // namespace ROL

#endif
