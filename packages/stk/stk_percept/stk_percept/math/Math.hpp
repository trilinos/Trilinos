#ifndef stk_percept_math_TransformationMatrix_hpp
#define stk_percept_math_TransformationMatrix_hpp

#include <stdexcept>
#include <sstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <typeinfo>

#include <math.h>

#include <Teuchos_ScalarTraits.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace stk {
  namespace percept {

    namespace ublas =  boost::numeric::ublas;

    class Math
    {
    public:

      typedef ublas::c_matrix<double,3,3> Matrix;

      typedef ublas::c_vector<double,3> Vector;

      typedef ublas::c_vector<double,3> ubvec;

      class MyVector : public ubvec
      {
      public:

        MyVector(double x=0.0) : ubvec()
        {
          (*this)(0) = x;
          (*this)(1) = x;
          (*this)(2) = x;
        }

        MyVector(double *x) : ubvec()
        {
          (*this)(0) = x[0];
          (*this)(1) = x[1];
          (*this)(2) = x[2];
        }
        //Vector(const ubvec& v) : ubvec(v) {}

        MyVector& operator=(const ubvec& v) 
        { 
          //ubvec& v0 = ubvec::operator=(v); 
          (*this)(0) = v(0);
          (*this)(1) = v(1);
          (*this)(2) = v(2);
          return *this;
        }

        //v = ublas::prod(m_rotMat, v);

      };


      static double random01()
      {
        double rnd = Teuchos::ScalarTraits<double>::random();
        return (rnd+1.0)/2.0;
      }

      static Matrix rotationMatrix(int axis, double angle_degrees) 
      {
        Matrix rm;
        rm.clear();
        double theta = M_PI * angle_degrees / 180.0;
        double cost = std::cos(theta);
        double sint = std::sin(theta);
        if (axis == 2)
          {
            rm(0,0) = cost; rm(0,1) = -sint;
            rm(1,0) = sint; rm(1,1) = cost;
            rm(2,2) = 1.0;
          }
        else if (axis == 1)
          {
            rm(0,0) = cost; rm(0,2) = -sint;
            rm(2,0) = sint; rm(2,2) = cost;
            rm(1,1) = 1.0;
          }
        else if (axis == 0)
          {
            rm(1,1) = cost; rm(1,2) = -sint;
            rm(2,1) = sint; rm(2,2) = cost;
            rm(0,0) = 1.0;
          }
        return rm;
      }

      static Matrix scalingMatrix(int axis, double scale)
      {
        Matrix sm;
        sm.clear();
        sm(0,0)=1.0;
        sm(1,1)=1.0;
        sm(2,2)=1.0;
        sm(axis,axis)=scale;
        return sm;
      }

      static Matrix scalingMatrix( double scale)
      {
        Matrix sm;
        sm.clear();
        sm(0,0)=scale;
        sm(1,1)=scale;
        sm(2,2)=scale;
        return sm;
      }
      
    };

    inline Math::Vector operator*(Math::Matrix& mat, Math::Vector& vec) { return ublas::prod(mat, vec); }
    inline Math::Matrix operator*(Math::Matrix& mat, Math::Matrix& mat2) { return ublas::prod(mat, mat2); }

  }
}
#endif
