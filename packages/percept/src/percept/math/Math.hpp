// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_math_TransformationMatrix_hpp
#define percept_math_TransformationMatrix_hpp

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

#include <percept/Stacktrace.hpp>
#include <percept/Util.hpp>

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

    static double my_abs_hi(double x, double eps=1.e-6) { return std::sqrt(x*x + eps*eps); }
    static double my_min_hi(double x, double y, double eps=1.e-6) { return 0.5*(x+y - my_abs_hi(x-y,eps)); }
    static double my_max_hi(double x, double y, double eps=1.e-6) { return 0.5*(x+y + my_abs_hi(x-y,eps)); }
    // heavyside
    static double heavy_smooth(double x, double x0, double eps=1.e-6) { return 1./(1.+std::exp(-2.0*(x-x0)/eps)); }

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

    static double norm_3d(const double * vec)
    {
      double norm = std::sqrt(vec[0]*vec[0]+
                              vec[1]*vec[1]+
                              vec[2]*vec[2]);
      return norm;
    }

    static void normalize_3d(double * vec)
    {
      double norm = norm_3d(vec);
      if (norm > 0.0)
        {
          vec[0] /= norm;
          vec[1] /= norm;
          vec[2] /= norm;
        }
      else
        {
          std::cout << "norm= " << vec[0] << ", " << vec[1] << ", " << vec[2] << "\n"
                    << Stacktrace::demangled_stacktrace(30) << std::endl;
          VERIFY_MSG("bad norm");
        }
    }
    static void cross_3d(const double * a, const double * b, double * axb)
    {
      axb[0] = (a[1]*b[2]-a[2]*b[1]);
      axb[1] = -(a[0]*b[2]-a[2]*b[0]);
      axb[2] = (a[0]*b[1]-a[1]*b[0]);
    }
    static double dot_3d(const double *a, const double *b)
    {
      return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }
    static double distance_squared_3d(const double *a, const double *b)
    {
      return (a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2]);
    }
    static double distance_3d(const double *a, const double *b)
    {
      return std::sqrt(distance_squared_3d(a,b));
    }
    static void copy_3d(double *a, const double *b)
    {
      a[0] = b[0];
      a[1] = b[1];
      a[2] = b[2];
    }
    static void copy_2d(double *a, const double *b)
    {
      a[0] = b[0];
      a[1] = b[1];
    }
    static std::string print_3d(const double *a, int prec=6)
    {
      std::ostringstream ostr;
      ostr << std::setprecision(prec);
      ostr << " " << a[0] << ", " << a[1] << ", " << a[2];
      return ostr.str();
    }
    static std::string print_2d(const double *a, int prec=6)
    {
      std::ostringstream ostr;
      ostr << std::setprecision(prec);
      ostr << " " << a[0] << ", " << a[1];
      return ostr.str();
    }
    static void subtract_3d(double *a, const double *b)
    {
      a[0] -= b[0];
      a[1] -= b[1];
      a[2] -= b[2];
    }

    static double project_to_line_3d(const double *a, const double *b, double *xyz, double& u)
    {
      double xma[3] = {0,0,0}, bma[3] = {0,0,0}, bman[3] = {0,0,0}, xyz_in[3] = {0,0,0};
      copy_3d(xyz_in, xyz);
      copy_3d(bma, b);
      subtract_3d(bma, a);
      copy_3d(xma, xyz);
      subtract_3d(xma, a);
      copy_3d(bman, bma);
      normalize_3d(bman);
      double dd = dot_3d(xma, bman);
      double ba_len = distance_3d(a, b);
      if (dd < 0.0)
        {
          u = 0.0;
          copy_3d(xyz, a);
        }
      else if (dd > ba_len)
        {
          u = 1.0;
          copy_3d(xyz, b);
        }
      else
        {
          u = dd/ba_len;
          xyz[0] = a[0] + u*(b[0] - a[0]);
          xyz[1] = a[1] + u*(b[1] - a[1]);
          xyz[2] = a[2] + u*(b[2] - a[2]);
        }
      return distance_3d(xyz, xyz_in);
    }
  };

  inline Math::Vector operator*(Math::Matrix& mat, Math::Vector& vec) { return ublas::prod(mat, vec); }
  inline Math::Matrix operator*(Math::Matrix& mat, Math::Matrix& mat2) { return ublas::prod(mat, mat2); }

}

#endif
