/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef PMMMsqMatrix_hpp
#define PMMMsqMatrix_hpp


#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

#include <MsqMatrix.hpp>

namespace stk {
  namespace percept {


    
    inline double my_sqr_Frobenius( const Mesquite::MsqMatrix<3,3>& m)
    {
      double sum=0.0;
#define R(i,j)   sum += m(i,j)*m(i,j)
      R(0,0);
      R(0,1);
      R(0,2);

      R(1,0);
      R(1,1);
      R(1,2);

      R(2,0);
      R(2,1);
      R(2,2);
#undef R
      return sum;
    }

    inline void inverse( const Mesquite::MsqMatrix<3,3>& m, Mesquite::MsqMatrix<3,3>& result )
    {
      //return adj(m) * (1.0 / det(m));
      const double detInv = 1.0/det(m);
      result(0,0) = (m(1,1)*m(2,2) - m(1,2)*m(2,1))*detInv;
      result(0,1) = (m(0,2)*m(2,1) - m(0,1)*m(2,2))*detInv;
      result(0,2) = (m(0,1)*m(1,2) - m(0,2)*m(1,1))*detInv;
  
      result(1,0) = (m(1,2)*m(2,0) - m(1,0)*m(2,2))*detInv;
      result(1,1) = (m(0,0)*m(2,2) - m(0,2)*m(2,0))*detInv;
      result(1,2) = (m(0,2)*m(1,0) - m(0,0)*m(1,2))*detInv;
  
      result(2,0) = (m(1,0)*m(2,1) - m(1,1)*m(2,0))*detInv;
      result(2,1) = (m(0,1)*m(2,0) - m(0,0)*m(2,1))*detInv;
      result(2,2) = (m(0,0)*m(1,1) - m(0,1)*m(1,0))*detInv;
    }

    // z = x * y
    inline void product( const Mesquite::MsqMatrix<3,3>& x, const Mesquite::MsqMatrix<3,3>& y, Mesquite::MsqMatrix<3,3>& z )
    {
#define R(i,j)   z(i,j) = x(i,0)*y(0,j) + x(i,1)*y(1,j) + x(i,2)*y(2,j)
      R(0,0);
      R(0,1);
      R(0,2);

      R(1,0);
      R(1,1);
      R(1,2);

      R(2,0);
      R(2,1);
      R(2,2);
#undef R

    }

    // z = x + y
    inline void sum( const Mesquite::MsqMatrix<3,3>& x, const Mesquite::MsqMatrix<3,3>& y, Mesquite::MsqMatrix<3,3>& z )
    {
#define R(i,j)   z(i,j) = x(i,j)+y(i,j)
      R(0,0);
      R(0,1);
      R(0,2);

      R(1,0);
      R(1,1);
      R(1,2);

      R(2,0);
      R(2,1);
      R(2,2);
#undef R

    }

    // z = x - y
    inline void difference( const Mesquite::MsqMatrix<3,3>& x, const Mesquite::MsqMatrix<3,3>& y, Mesquite::MsqMatrix<3,3>& z )
    {
#define R(i,j)   z(i,j) = x(i,j)-y(i,j)
      R(0,0);
      R(0,1);
      R(0,2);

      R(1,0);
      R(1,1);
      R(1,2);

      R(2,0);
      R(2,1);
      R(2,2);
#undef R

    }

    inline void identity( Mesquite::MsqMatrix<3,3>& I )
    {
      I(0,0) = 1.0;
      I(0,1) = 0.0;
      I(0,2) = 0.0;

      I(1,0) = 0.0;
      I(1,1) = 1.0;
      I(1,2) = 0.0;

      I(2,0) = 0.0;
      I(2,1) = 0.0;
      I(2,2) = 1.0;
    }

  }
}

#endif
#endif
