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
  class Objective_Rosenbrock : public Objective<Real> {
  private:
    Real alpha_;

  public:
    Objective_Rosenbrock(Real alpha = 100.0) : alpha_(alpha) {}

    Real value( const Vector<Real> &x, Real &tol ) {
      StdVector<Real> & ex =
        Teuchos::dyn_cast<StdVector<Real> >(const_cast <Vector<Real> &>(x));
      Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();

      int n = xp->size();
      Real val = 0;
      for( int i=0; i<n/2; i++ ) {
        val += alpha_ * pow(pow((*xp)[2*i],2) - (*xp)[2*i+1], 2);
        val += pow((*xp)[2*i] - 1.0, 2);
      }
 
      return val;
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > gp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());

      int n = xp->size();
      for( int i=0; i<n/2; i++ ) {
        (*gp)[2*i]   =  4.0*alpha_*(pow((*xp)[2*i],2) - (*xp)[2*i+1])*(*xp)[2*i] + 2.0*((*xp)[2*i]-1.0);
        (*gp)[2*i+1] = -2.0*alpha_*(pow((*xp)[2*i],2) - (*xp)[2*i+1]);
      }
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

      int n = xp->size();
      for( int i=0; i<n/2; i++ ) {
        Real h11 = 4.0*alpha_*(3.0*pow((*xp)[2*i],2)-(*xp)[2*i+1]) + 2.0;
        Real h12 = -4.0*alpha_*(*xp)[2*i];
        Real h22 = 2.0*alpha_;

        (*hvp)[2*i]   = h11*(*vp)[2*i] + h12*(*vp)[2*i+1];
        (*hvp)[2*i+1] = h12*(*vp)[2*i] + h22*(*vp)[2*i+1];
      }
    }
#endif
    void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

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
  class Objective_FreudensteinRoth : public Objective<Real> {
  public:
    Objective_FreudensteinRoth() {}

    Real value( const Vector<Real> &x, Real &tol ) {
      StdVector<Real> & ex =
        Teuchos::dyn_cast<StdVector<Real> >(const_cast <Vector<Real> &>(x));
      Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();

      Real f1 = -13.0 + (*xp)[0] + ((5.0-(*xp)[1])*(*xp)[1] - 2.0)*(*xp)[1];
      Real f2 = -29.0 + (*xp)[0] + (((*xp)[1]+1.0)*(*xp)[1] - 14.0)*(*xp)[1];
 
      return f1*f1+f2*f2;
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > gp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());
  
      Real f1 = -13.0 + (*xp)[0] + ((5.0-(*xp)[1])*(*xp)[1] - 2.0)*(*xp)[1];
      Real f2 = -29.0 + (*xp)[0] + (((*xp)[1]+1.0)*(*xp)[1] - 14.0)*(*xp)[1];
 
      Real f11 = 1.0;
      Real f12 = 10.0*(*xp)[1] - 3.0*(*xp)[1]*(*xp)[1] - 2.0;
      Real f21 = 1.0;
      Real f22 = 3.0*(*xp)[1]*(*xp)[1] + 2.0*(*xp)[1] - 14.0;  

      (*gp)[0] = 2.0*(f11*f1 + f21*f2);
      (*gp)[1] = 2.0*(f12*f1 + f22*f2);
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());
  
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
#endif
    void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());
  
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

  /** \brief Beale's function.
   */
  template<class Real>
  class Objective_Beale : public Objective<Real> {
  private: 
    std::vector<Real> y_;

  public:
    Objective_Beale() {
      y_.clear();
      y_.push_back(1.5);
      y_.push_back(2.25);
      y_.push_back(2.625);
    }

    Real value( const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > ex =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();

      Real f1 = 1.5-(*ex)[0]*(1.0-(*ex)[1]);
      Real f2 = 2.25-(*ex)[0]*(1.0-pow((*ex)[1],2));
      Real f3 = 2.625-(*ex)[0]*(1.0-pow((*ex)[1],3));
  
      return pow(f1,2)+pow(f2,2)+pow(f3,2);
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > ex =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > eg =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());

      Real f1 = 1.5-(*ex)[0]*(1.0-(*ex)[1]);
      Real f2 = 2.25-(*ex)[0]*(1.0-pow((*ex)[1],2));
      Real f3 = 2.625-(*ex)[0]*(1.0-pow((*ex)[1],3));
      Real df1dx = -(1.0-(*ex)[1]);
      Real df1dy = (*ex)[0];
      Real df2dx = -(1.0-pow((*ex)[1],2));
      Real df2dy = 2.0*(*ex)[0]*(*ex)[1];
      Real df3dx = -(1.0-pow((*ex)[1],3));
      Real df3dy = 3.0*(*ex)[0]*pow((*ex)[1],2);

      (*eg)[0] = 2.0*df1dx*f1+2.0*df2dx*f2+2.0*df3dx*f3;
      (*eg)[1] = 2.0*df1dy*f1+2.0*df2dy*f2+2.0*df3dy*f3;
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > ex =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > ev =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > ehv =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

      Real f1 = 1.5-(*ex)[0]*(1.0-(*ex)[1]);
      Real f2 = 2.25-(*ex)[0]*(1.0-pow((*ex)[1],2));
      Real f3 = 2.625-(*ex)[0]*(1.0-pow((*ex)[1],3));
      Real df1dx = -(1.0-(*ex)[1]);
      Real df1dy = (*ex)[0];
      Real df2dx = -(1.0-pow((*ex)[1],2));
      Real df2dy = 2.0*(*ex)[0]*(*ex)[1];
      Real df3dx = -(1.0-pow((*ex)[1],3));
      Real df3dy = 3.0*(*ex)[0]*pow((*ex)[1],2);
      Real d2f1dx2 = 0.0;
      Real d2f1dy2 = 0.0;
      Real d2f1dxdy = 1.0;
      Real d2f2dx2 = 0.0;
      Real d2f2dy2 = 2.0*(*ex)[0];
      Real d2f2dxdy = 2.0*(*ex)[1];
      Real d2f3dx2 = 0.0;
      Real d2f3dy2 = 6.0*(*ex)[0]*(*ex)[1];
      Real d2f3dxdy = 3.0*pow((*ex)[1],2);

      Real H11 = 2.0*(d2f1dx2*f1+df1dx*df1dx)+2.0*(d2f2dx2*f2+df2dx*df2dx)
                  +2.0*(d2f3dx2*f3+df3dx*df3dx);
      Real H22 = 2.0*(d2f1dy2*f1+df1dy*df1dy)+2.0*(d2f2dy2*f2+df2dy*df2dy)
                  +2.0*(d2f3dy2*f3+df3dy*df3dy);
      Real H12 = 2.0*(d2f1dxdy*f1 + df1dx*df1dy)+2.0*(d2f2dxdy*f2 + df2dx*df2dy)
                  +2.0*(d2f3dxdy*f3 + df3dx*df3dy);

      (*ehv)[0] = H11*(*ev)[0]+H12*(*ev)[1];
      (*ehv)[1] = H12*(*ev)[0]+H22*(*ev)[1];
    }
#endif
    void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > ex =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > ev =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > ehv =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

      Real f1 = 1.5-(*ex)[0]*(1.0-(*ex)[1]);
      Real f2 = 2.25-(*ex)[0]*(1.0-pow((*ex)[1],2));
      Real f3 = 2.625-(*ex)[0]*(1.0-pow((*ex)[1],3));
      Real df1dx = -(1.0-(*ex)[1]);
      Real df1dy = (*ex)[0];
      Real df2dx = -(1.0-pow((*ex)[1],2));
      Real df2dy = 2.0*(*ex)[0]*(*ex)[1];
      Real df3dx = -(1.0-pow((*ex)[1],3));
      Real df3dy = 3.0*(*ex)[0]*pow((*ex)[1],2);
      Real d2f1dx2 = 0.0;
      Real d2f1dy2 = 0.0;
      Real d2f1dxdy = 1.0;
      Real d2f2dx2 = 0.0;
      Real d2f2dy2 = 2.0*(*ex)[0];
      Real d2f2dxdy = 2.0*(*ex)[1];
      Real d2f3dx2 = 0.0;
      Real d2f3dy2 = 6.0*(*ex)[0]*(*ex)[1];
      Real d2f3dxdy = 3.0*pow((*ex)[1],2);

      Real H11 = 2.0*(d2f1dx2*f1+df1dx*df1dx)+2.0*(d2f2dx2*f2+df2dx*df2dx)
                  +2.0*(d2f3dx2*f3+df3dx*df3dx);
      Real H22 = 2.0*(d2f1dy2*f1+df1dy*df1dy)+2.0*(d2f2dy2*f2+df2dy*df2dy)
                  +2.0*(d2f3dy2*f3+df3dy*df3dy);
      Real H12 = 2.0*(d2f1dxdy*f1 + df1dx*df1dy)+2.0*(d2f2dxdy*f2 + df2dx*df2dy)
                  +2.0*(d2f3dxdy*f3 + df3dx*df3dy);

      (*ehv)[0] = (1.0/(H11*H22-H12*H12))*( H22*(*ev)[0] - H12*(*ev)[1]);
      (*ehv)[1] = (1.0/(H11*H22-H12*H12))*(-H12*(*ev)[0] + H11*(*ev)[1]);
    }
  };

  /** \brief Powell's badly scaled function.
   */
  template<class Real>
  class Objective_Powell : public Objective<Real> {
  public:
    Objective_Powell() {}

    Real value( const Vector<Real> &x, Real &tol ) {
      StdVector<Real> & ex =
        Teuchos::dyn_cast<StdVector<Real> >(const_cast <Vector<Real> &>(x));
      Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();

      Real f1 = 1.e4*(*xp)[0]*(*xp)[1] - 1.0;
      Real f2 = std::exp(-(*xp)[0]) + std::exp(-(*xp)[1]) - 1.0001;
 
      return f1*f1+f2*f2;
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > gp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());

      Real f1 = 1.e4*(*xp)[0]*(*xp)[1] - 1.0;
      Real f2 = std::exp(-(*xp)[0]) + std::exp(-(*xp)[1]) - 1.0001;
 
      Real f11 = 1.e4*(*xp)[1];
      Real f12 = 1.e4*(*xp)[0];
      Real f21 = -std::exp(-(*xp)[0]);
      Real f22 = -std::exp(-(*xp)[1]);  

      (*gp)[0] = 2.0*(f11*f1 + f21*f2);
      (*gp)[1] = 2.0*(f12*f1 + f22*f2);
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());
  
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
#endif
    void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());
  
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
  class Objective_SumOfSquares : public Objective<Real> {
  public:
    Real value( const Vector<Real> &x, Real &tol ) {
      StdVector<Real> & ex =
        Teuchos::dyn_cast<StdVector<Real> >(const_cast <Vector<Real> &>(x));
      Teuchos::RCP<const std::vector<Real> > xp = ex.getVector();

      int n = xp->size();
      Real val = 0;
      for (int i=0; i<n; i++) {
        val += pow((*xp)[i], 2);
      }

      return val;
    }

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > gp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());

      int n = xp->size();
      for( int i=0; i<n; i++ ) {
        (*gp)[i] = 2.0*(*xp)[i];
      }
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

      int n = xp->size();
      for( int i=0; i<n; i++ ) {
        (*hvp)[i] = 2.0*(*vp)[i];
      }
    }
#endif
    void invHessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

      int n = xp->size();
      for( int i=0; i<n; i++ ) {
        (*hvp)[i] = 0.5*(*vp)[i];
      }
    }
  };

  /** \brief Least squares function.
   */
  template<class Real>
  class Objective_LeastSquares : public Objective<Real> {
  public:
    Real value( const Vector<Real> &x, Real &tol ) {
      StdVector<Real> & ex =
        Teuchos::dyn_cast<StdVector<Real> >(const_cast <Vector<Real> &>(x));
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

    void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<std::vector<Real> > gp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());

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
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > xp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

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
#endif
  };

  /** \brief Poisson distributed control.
   */
  template<class Real>
  class Objective_PoissonControl : public Objective<Real> {
  private: 
    Real alpha_;

  public:

    Objective_PoissonControl(Real alpha = 1.e-4) : alpha_(alpha) {}

    void apply_mass(Vector<Real> &Mz, const Vector<Real> &z ) {
      Teuchos::RCP<const std::vector<Real> > zp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
      Teuchos::RCP<std::vector<Real> > Mzp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(Mz)).getVector());

      int  n = zp->size();
      Real h = 1.0/((Real)n+1.0);
      for (int i=0; i<n; i++) {
        if ( i == 0 ) {
          (*Mzp)[i] = h/6.0*(4.0*(*zp)[i] + (*zp)[i+1]);
        }
        else if ( i == n-1 ) {
          (*Mzp)[i] = h/6.0*((*zp)[i-1] + 4.0*(*zp)[i]);
        }
        else {
          (*Mzp)[i] = h/6.0*((*zp)[i-1] + 4.0*(*zp)[i] + (*zp)[i+1]);
        }
      }
    }

    void solve_poisson(Vector<Real> & u, const Vector<Real> & z) {
      Teuchos::RCP<std::vector<Real> > up =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(u)).getVector());

      int  n = up->size();
      Real h = 1.0/((Real)n+1.0);
      StdVector<Real> b( Teuchos::rcp( new std::vector<Real>(n,0.0) ) );
      Teuchos::RCP<std::vector<Real> > bp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(b)).getVector());
      this->apply_mass(b,z);
  
      Real d   =  2.0/h;
      Real o   = -1.0/h;
      Real m   = 0.0;
      std::vector<Real> c(n,o);
      c[0]     = c[0]/d;
      (*up)[0] = (*bp)[0]/d;
      for ( int i = 1; i < n; i++ ) {
        m        = 1.0/(d - o*c[i-1]);
        c[i]     = c[i]*m;
        (*up)[i] = ( (*bp)[i] - o*(*up)[i-1] )*m;
      }
      for ( int i = n-1; i > 0; i-- ) {
        (*up)[i-1] = (*up)[i-1] - c[i-1]*(*up)[i];
      }
    }

    Real evaluate_target(Real x) {
      Real val = 1.0/3.0*std::pow(x,4.0) - 2.0/3.0*std::pow(x,3.0) + 1.0/3.0*x + 8.0*this->alpha_;
      return val;
    }

    Real value( const Vector<Real> &z, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > zp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
      int  n    = zp->size();
      Real h    = 1.0/((Real)n+1.0);
      // SOLVE STATE EQUATION
      StdVector<Real> u( Teuchos::rcp( new std::vector<Real>(n,0.0) ) );
      this->solve_poisson(u,z);      
      Teuchos::RCP<std::vector<Real> > up =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(u)).getVector());

      Real val  = 0.0;
      Real res  = 0.0;
      Real res1 = 0.0;
      Real res2 = 0.0;
      Real res3 = 0.0;
      for (int i=0; i<n; i++) {
	res = this->alpha_*(*zp)[i];
        if ( i == 0 ) {
          res *= h/6.0*(4.0*(*zp)[i] + (*zp)[i+1]);
          res1 = (*up)[i]-evaluate_target((Real)(i+1)*h);
          res2 = (*up)[i+1]-evaluate_target((Real)(i+2)*h);
          res += h/6.0*(4.0*res1 + res2)*res1;
        }
        else if ( i == n-1 ) {
          res *= h/6.0*((*zp)[i-1] + 4.0*(*zp)[i]);
          res1 = (*up)[i-1]-evaluate_target((Real)(i)*h);
          res2 = (*up)[i]-evaluate_target((Real)(i+1)*h);
          res += h/6.0*(res1 + 4.0*res2)*res2;
        }
        else {
          res *= h/6.0*((*zp)[i-1] + 4.0*(*zp)[i] + (*zp)[i+1]);
          res1 = (*up)[i-1]-evaluate_target((Real)(i)*h);
          res2 = (*up)[i]-evaluate_target((Real)(i+1)*h);
          res3 = (*up)[i+1]-evaluate_target((Real)(i+2)*h);
          res += h/6.0*(res1 + 4.0*res2 + res3)*res2;
        }
        val += 0.5*res;
      }
      return val;
   }

    void gradient( Vector<Real> &g, const Vector<Real> &z, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > zp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
      Teuchos::RCP<std::vector<Real> > gp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(g)).getVector());
      int  n = zp->size();
      Real h = 1.0/((Real)n+1.0);

      // SOLVE STATE EQUATION
      StdVector<Real> u( Teuchos::rcp( new std::vector<Real>(n,0.0) ) );
      this->solve_poisson(u,z);      
      Teuchos::RCP<std::vector<Real> > up =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(u)).getVector());

      // SOLVE ADJOINT EQUATION
      StdVector<Real> res( Teuchos::rcp( new std::vector<Real>(n,0.0) ) );
      Teuchos::RCP<std::vector<Real> > rp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(res)).getVector());
      for (int i=0; i<n; i++) {
        (*rp)[i] = -((*up)[i]-evaluate_target((Real)(i+1)*h));
      }
      StdVector<Real> p( Teuchos::rcp( new std::vector<Real>(n,0.0) ) );
      this->solve_poisson(p,res);      
      Teuchos::RCP<std::vector<Real> > pp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(p)).getVector());

      Real res1 = 0.0;
      Real res2 = 0.0;
      Real res3 = 0.0;
      for (int i=0; i<n; i++) {
        if ( i == 0 ) {
          res1 = this->alpha_*(*zp)[i] - (*pp)[i];
          res2 = this->alpha_*(*zp)[i+1] - (*pp)[i+1];
          (*gp)[i] = h/6.0*(4.0*res1 + res2);
        }
        else if ( i == n-1 ) {
          res1 = this->alpha_*(*zp)[i-1] - (*pp)[i-1];
          res2 = this->alpha_*(*zp)[i] - (*pp)[i];
          (*gp)[i] = h/6.0*(res1 + 4.0*res2);
        }
        else {
          res1 = this->alpha_*(*zp)[i-1] - (*pp)[i-1];
          res2 = this->alpha_*(*zp)[i] - (*pp)[i];
          res3 = this->alpha_*(*zp)[i+1] - (*pp)[i+1];
          (*gp)[i] = h/6.0*(res1 + 4.0*res2 + res3);
        }
      }
    }
#if USE_HESSVEC
    void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &z, Real &tol ) {
      Teuchos::RCP<const std::vector<Real> > zp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(z))).getVector();
      Teuchos::RCP<const std::vector<Real> > vp =
        (Teuchos::dyn_cast<StdVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      Teuchos::RCP<std::vector<Real> > hvp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(hv)).getVector());

      int  n = zp->size();
      Real h = 1.0/((Real)n+1.0);

      // SOLVE STATE EQUATION
      StdVector<Real> u( Teuchos::rcp( new std::vector<Real>(n,0.0) ) );
      this->solve_poisson(u,v);      
      Teuchos::RCP<std::vector<Real> > up =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(u)).getVector());

      // SOLVE ADJOINT EQUATION
      StdVector<Real> p( Teuchos::rcp( new std::vector<Real>(n,0.0) ) );
      this->solve_poisson(p,u);      
      Teuchos::RCP<std::vector<Real> > pp =
        Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(p)).getVector());

      Real res1 = 0.0;
      Real res2 = 0.0;
      Real res3 = 0.0;
      for (int i=0; i<n; i++) {
        if ( i == 0 ) {
          res1 = this->alpha_*(*vp)[i] + (*pp)[i];
          res2 = this->alpha_*(*vp)[i+1] + (*pp)[i+1];
          (*hvp)[i] = h/6.0*(4.0*res1 + res2);
        }
        else if ( i == n-1 ) {
          res1 = this->alpha_*(*vp)[i-1] + (*pp)[i-1];
          res2 = this->alpha_*(*vp)[i] + (*pp)[i];
          (*hvp)[i] = h/6.0*(res1 + 4.0*res2);
        }
        else {
          res1 = this->alpha_*(*vp)[i-1] + (*pp)[i-1];
          res2 = this->alpha_*(*vp)[i] + (*pp)[i];
          res3 = this->alpha_*(*vp)[i+1] + (*pp)[i+1];
          (*hvp)[i] = h/6.0*(res1 + 4.0*res2 + res3);
        }
      }
    }
#endif
  };

  template<class Real>
  void getTestObjectives( Teuchos::RCP<Objective<Real> > &obj, Vector<Real> &x0, Vector<Real> &x, 
                          const ETestObjectives test ) {
    // Cast Initial Guess and Solution Vectors
    Teuchos::RCP<std::vector<Real> > x0p =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x0)).getVector());
    Teuchos::RCP<std::vector<Real> > xp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x)).getVector());
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
    else if ( test == TESTOBJECTIVES_BEALE ) {
      // Resize Vectors
      n = 2;
      x0p->resize(n);
      xp->resize(n);
      // Instantiate Objective Function
      obj = Teuchos::rcp( new Objective_Beale<Real> );
      // Get Initial Guess
      (*x0p)[0] =  1.0;
      (*x0p)[1] =  1.0;
      // Get Solution
      (*xp)[0] = 3.0;
      (*xp)[1] = 0.5;
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
    else if ( test == TESTOBJECTIVES_POISSONCONTROL ) {
      // Resize Vectors
      n = 128;
      x0p->resize(n);
      xp->resize(n);
      // Instantiate Objective Function
      obj = Teuchos::rcp( new Objective_PoissonControl<Real> );
      // Get Initial Guess
      for (int i=0; i<n; i++) {
        (*x0p)[i] = 0.0;
      }
      // Get Solution
      Real h  = 1.0/((Real)n+1.0);
      Real pt = 0.0;
      for( int i=0; i<n; i++ ) {
        pt = (Real)(i+1)*h;
        (*xp)[i] = 4.0*pt*(1.0-pt);
      }
    }
  }

} // namespace ROL

#endif
