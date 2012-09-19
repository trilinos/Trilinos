#ifndef STK_MMS_stk_mms_h
#define STK_MMS_stk_mms_h

#include <Sacado.hpp>

namespace stk
{

namespace mms
{

  typedef Sacado::Fad::DFad<double> FAD_Type;
  typedef Sacado::Fad::DFad<Sacado::Fad::DFad<double> > FAD2_Type;
  
  typedef FAD_Type (scalar_FAD_function)(const FAD_Type & x,
					 const FAD_Type & y,
					 const FAD_Type & z,
					 const FAD_Type & t
					 );

  typedef FAD2_Type (scalar_FAD2_function)(const FAD2_Type & x,
					   const FAD2_Type & y,
					   const FAD2_Type & z,
					   const FAD2_Type & t
					   );

  enum PARTIAL_DERIVATIVES {_X=0, _Y, _Z, _T};
  
  // first partial derivatives for first order FADs
  double DX(const FAD_Type &f);
  double DY(const FAD_Type &f);
  double DZ(const FAD_Type &f);
  double DT(const FAD_Type &f);
  
  double   D(const  FAD_Type &f, const int i);

  FAD_Type D(const FAD2_Type &f, const int i);

  // first order partial derivatives for second order FADs
  FAD_Type DX2(const FAD2_Type &f);
  FAD_Type DY2(const FAD2_Type &f);
  FAD_Type DZ2(const FAD2_Type &f);
  FAD_Type DT2(const FAD2_Type &f);

  FAD_Type D2(const FAD2_Type &f, const int i);

  // NOTE I would like to overload using DX like this:
  //   FAD_Type DX(const FAD2_Type &f) { return f.dx(0);}
  // but that code does not compile

/*   FAD_Type DX2(const FAD2_Type &f) { return f.dx(0);} */
/*   FAD_Type DY2(const FAD2_Type &f) { return f.dx(1);} */
/*   FAD_Type DZ2(const FAD2_Type &f) { return f.dx(2);} */
/*   FAD_Type DT2(const FAD2_Type &f) { return f.dx(3);} */
  
/*   // second order partial derivatives for second order FADs */
/*   double DXX(const FAD2_Type &f) { return DX(DX2(f));} */
/*   double DXY(const FAD2_Type &f) { return DX(DY2(f));} */
/*   double DXZ(const FAD2_Type &f) { return DX(DZ2(f));} */

/*   double DYX(const FAD2_Type &f) { return DY(DX2(f));} */
/*   double DYY(const FAD2_Type &f) { return DY(DY2(f));} */
/*   double DYZ(const FAD2_Type &f) { return DY(DZ2(f));} */

/*   double DZX(const FAD2_Type &f) { return DZ(DX2(f));} */
/*   double DZY(const FAD2_Type &f) { return DZ(DY2(f));} */
/*   double DZZ(const FAD2_Type &f) { return DZ(DZ2(f));} */
  
/*   // examples of scalar differential operators */
/*   double LAPLACE2D(const FAD2_Type &f) { return DXX(f)+DYY(f);} */
/*   double LAPLACE3D(const FAD2_Type &f) { return DXX(f)+DYY(f)+DZZ(f);} */
  
  // API for scalar source terms  
  class ScalarMMSBase
  {
  public:
    
    ScalarMMSBase() {}
    
    virtual ~ScalarMMSBase() {}
    
    // scalar source term interface
    virtual double source(const FAD_Type & x,
			  const FAD_Type & y,
			  const FAD_Type & z,
			  const FAD_Type & time) = 0;
  };
  
}
}

#endif // STK_MMS_stk_mms_h
