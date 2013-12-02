// @HEADER
// ************************************************************************
//
// Questions? Contact Denis Ridzal (dridzal@sandia.gov), or
//                    Drew Kouri   (dpkouri@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for Beale's function.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef USE_HESSVEC 
#define USE_HESSVEC 1
#endif

#ifndef ROL_BEALE_HPP
#define ROL_BEALE_HPP

#include "ROL_StdVector.hpp"
#include "ROL_Objective.hpp"

namespace ROL {

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

  template<class Real>
  void getBeale( Teuchos::RCP<Objective<Real> > &obj, Vector<Real> &x0, Vector<Real> &x ) {
    // Cast Initial Guess and Solution Vectors
    Teuchos::RCP<std::vector<Real> > x0p =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x0)).getVector());
    Teuchos::RCP<std::vector<Real> > xp =
      Teuchos::rcp_const_cast<std::vector<Real> >((Teuchos::dyn_cast<StdVector<Real> >(x)).getVector());
    int n = xp->size();
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


}// End ROL Namespace

#endif
