// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_MasterElementBasis_h
#define Akri_MasterElementBasis_h

#define ATTR_RESTRICT __restrict__

namespace krino {

class Basis
{
public:
  Basis(unsigned deg) : my_degree(deg) {}

  virtual ~Basis() {}
  virtual void nodal_parametric_coordinates(double* ATTR_RESTRICT p_coords) const = 0;
  virtual double parametric_volume() const = 0;
  virtual void shape_fcn(const int nint, const double* p_coords, double* result) const = 0;
  virtual void shape_fcn_deriv(const int nint, const double* p_coords, double* result ) const = 0;
  unsigned degree() const { return my_degree; }

private:
  unsigned my_degree;
};

class Basis_LINE_2 : public Basis
{
public:
  Basis_LINE_2() : Basis(1) {}

  double parametric_volume() const override { return 2.; }

  void nodal_parametric_coordinates(double* ATTR_RESTRICT p_coords) const override
  {
    p_coords[0] = -1.0;
    p_coords[1] =  1.0;
  }
  void shape_fcn(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for (int ip(0); ip < nint; ++ip)
    {
      const double x = p_coords[ip];
      result[ip*2 + 0] = ( 1.0 - x ) * 0.5;
      result[ip*2 + 1] = ( 1.0 + x ) * 0.5;
    }
  }
  void shape_fcn_deriv(const int nint, const double* ATTR_RESTRICT /*p_coords*/, double* ATTR_RESTRICT result) const override
  {
    for ( int ip(0); ip < nint; ++ip )
    {
      result[ip*2 + 0] = -0.5;
      result[ip*2 + 1] =  0.5;
    }
  }
};

class Basis_LINE_3 : public Basis
{
public:
  Basis_LINE_3() : Basis(2) {}

  double parametric_volume() const override { return 2.; }

  void nodal_parametric_coordinates(double* ATTR_RESTRICT p_coords) const override
  {
    p_coords[0] = -1.0;
    p_coords[1] =  1.0;
    p_coords[2] =  0.0;
  }
  void shape_fcn(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for (int ip(0); ip < nint; ++ip)
    {
      const double x = p_coords[ip];
      result[ip*3 + 0] = -x * ( 1.0 - x ) * 0.5;
      result[ip*3 + 1] =  x * ( 1.0 + x ) * 0.5;
      result[ip*3 + 2] = ( 1.0 - x ) * ( 1.0 + x );
    }
  }
  void shape_fcn_deriv(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for ( int ip(0); ip < nint; ++ip )
    {
      const double x = p_coords[ip];
      result[ip*3 + 0] = x;
      result[ip*3 + 1] = x;
      result[ip*3 + 2] = -2.0*x;
    }
  }
};

class Basis_TRI_3 : public Basis
{
public:
  Basis_TRI_3() : Basis(1) {}

  double parametric_volume() const override { return 0.5; }

  void nodal_parametric_coordinates(double* ATTR_RESTRICT p_coords) const override
  {
    p_coords[0] = 0.0; p_coords[1] = 0.0;
    p_coords[2] = 1.0; p_coords[3] = 0.0;
    p_coords[4] = 0.0; p_coords[5] = 1.0;
  }
  void shape_fcn(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for (int ip(0); ip < nint; ++ip)
    {
      const double x = p_coords[2*ip + 0];
      const double y = p_coords[2*ip + 1];
      result[ip*3 + 0] = 1.0 - x - y;
      result[ip*3 + 1] = x;
      result[ip*3 + 2] = y;
    }
  }
  void shape_fcn_deriv(const int nint, const double* ATTR_RESTRICT /*p_coords*/, double* ATTR_RESTRICT result) const override
  {
    for ( int ip(0); ip < nint; ++ip )
    {
      result[ip*6 + 0] = -1.0;
      result[ip*6 + 1] = -1.0;
      result[ip*6 + 2] =  1.0;
      result[ip*6 + 3] =  0.0;
      result[ip*6 + 4] =  0.0;
      result[ip*6 + 5] =  1.0;
    }
  }
};

class Basis_TRI_6 : public Basis
{
public:
  Basis_TRI_6() : Basis(2) {}

  double parametric_volume() const override { return 0.5; }

  void nodal_parametric_coordinates(double* ATTR_RESTRICT p_coords) const override
  {
    p_coords[ 0] = 0.0; p_coords[ 1] = 0.0;
    p_coords[ 2] = 1.0; p_coords[ 3] = 0.0;
    p_coords[ 4] = 0.0; p_coords[ 5] = 1.0;
    p_coords[ 6] = 0.5; p_coords[ 7] = 0.0;
    p_coords[ 8] = 0.5; p_coords[ 9] = 0.5;
    p_coords[10] = 0.0; p_coords[11] = 0.5;
  }
  void shape_fcn(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for (int ip(0); ip < nint; ++ip)
    {
      const double x = p_coords[2*ip + 0];
      const double y = p_coords[2*ip + 1];
      result[ip*6 + 0] = (x + y - 1.0)*(2.0*x + 2.0*y - 1.0);
      result[ip*6 + 1] = x*(2.0*x - 1.0);
      result[ip*6 + 2] = y*(2.0*y - 1.0);
      result[ip*6 + 3] = -4.0*x*(x + y - 1.0);
      result[ip*6 + 4] = 4.0*x*y;
      result[ip*6 + 5] = -4.0*y*(x + y - 1.0);
    }
  }
  void shape_fcn_deriv(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for ( int ip(0); ip < nint; ++ip )
    {
      const double x = p_coords[2*ip + 0];
      const double y = p_coords[2*ip + 1];
      result[ip*12 +  0] =  4.0*x + 4.0*y - 3.0;
      result[ip*12 +  1] =  4.0*x + 4.0*y - 3.0;
      result[ip*12 +  2] =  4.0*x - 1.0;
      result[ip*12 +  3] =  0.0;
      result[ip*12 +  4] =  0.0;
      result[ip*12 +  5] =  4.0*y - 1.0;
      result[ip*12 +  6] = -4.0*(2.0*x + y - 1.0);
      result[ip*12 +  7] = -4.0*x;
      result[ip*12 +  8] =  4.0*y;
      result[ip*12 +  9] =  4.0*x;
      result[ip*12 + 10] = -4.0*y;
      result[ip*12 + 11] = -4.0*(x + 2.0*y - 1.0);
    }
  }
};

class Basis_QUAD_4 : public Basis
{
public:
  Basis_QUAD_4() : Basis(1) {}

  double parametric_volume() const override { return 4.0; }

  void nodal_parametric_coordinates(double* ATTR_RESTRICT p_coords) const override
  {
    p_coords[ 0] = -1.0; p_coords[ 1] = -1.0;
    p_coords[ 2] =  1.0; p_coords[ 3] = -1.0;
    p_coords[ 4] =  1.0; p_coords[ 5] =  1.0;
    p_coords[ 6] = -1.0; p_coords[ 7] =  1.0;
  }
  void shape_fcn(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for (int ip(0); ip < nint; ++ip)
    {
      const double x = p_coords[2*ip + 0];
      const double y = p_coords[2*ip + 1];
      result[ip*4 + 0] = 0.25*(1.0 - x)*(1.0 - y);
      result[ip*4 + 1] = 0.25*(1.0 + x)*(1.0 - y);
      result[ip*4 + 2] = 0.25*(1.0 + x)*(1.0 + y);
      result[ip*4 + 3] = 0.25*(1.0 - x)*(1.0 + y);
    }
  }
  void shape_fcn_deriv(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for ( int ip(0); ip < nint; ++ip )
    {
      const double x = p_coords[2*ip + 0];
      const double y = p_coords[2*ip + 1];
      result[ip*8 + 0] = -0.25*(1.0 - y);
      result[ip*8 + 1] = -0.25*(1.0 - x);
      result[ip*8 + 2] =  0.25*(1.0 - y);
      result[ip*8 + 3] = -0.25*(1.0 + x);
      result[ip*8 + 4] =  0.25*(1.0 + y);
      result[ip*8 + 5] =  0.25*(1.0 + x);
      result[ip*8 + 6] = -0.25*(1.0 + y);
      result[ip*8 + 7] =  0.25*(1.0 - x);
    }
  }
};

class Basis_QUAD_9 : public Basis
{
public:
  Basis_QUAD_9() : Basis(2) {}

  double parametric_volume() const override { return 4.0; }

  void nodal_parametric_coordinates(double* ATTR_RESTRICT p_coords) const override
  {
    p_coords[ 0] = -1.0; p_coords[ 1] = -1.0;
    p_coords[ 2] =  1.0; p_coords[ 3] = -1.0;
    p_coords[ 4] =  1.0; p_coords[ 5] =  1.0;
    p_coords[ 6] = -1.0; p_coords[ 7] =  1.0;
    p_coords[ 8] =  0.0; p_coords[ 9] = -1.0;
    p_coords[10] =  1.0; p_coords[11] =  0.0;
    p_coords[12] =  0.0; p_coords[13] =  1.0;
    p_coords[14] = -1.0; p_coords[15] =  0.0;
    p_coords[16] =  0.0; p_coords[17] =  0.0;
  }
  void shape_fcn(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for (int ip(0); ip < nint; ++ip)
    {
      const double x = p_coords[2*ip + 0];
      const double y = p_coords[2*ip + 1];
      result[ip*9 + 0] = 0.25*x*(x - 1.0)*y*(y - 1.0);
      result[ip*9 + 1] = 0.25*x*(x + 1.0)*y*(y - 1.0);
      result[ip*9 + 2] = 0.25*x*(x + 1.0)*y*(y + 1.0);
      result[ip*9 + 3] = 0.25*x*(x - 1.0)*y*(y + 1.0);
      result[ip*9 + 4] = 0.5*(1.0 - x)*(1.0 + x)*y*(y - 1.0);
      result[ip*9 + 5] = 0.5*x*(x + 1.0)*(1.0 - y)*(1.0 + y);
      result[ip*9 + 6] = 0.5*(1.0 - x)*(1.0 + x)*y*(y + 1.0);
      result[ip*9 + 7] = 0.5*x*(x - 1.0)*(1.0 - y)*(1.0 + y);
      result[ip*9 + 8] = (1.0 - x)*(1.0 + x)*(1.0 - y)*(1.0 + y);
    }
  }
  void shape_fcn_deriv(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for ( int ip(0); ip < nint; ++ip )
    {
      const double x = p_coords[2*ip + 0];
      const double y = p_coords[2*ip + 1];
      result[ip*18 +  0] = (-0.25 + 0.5*x)*(-1. + y)*y;
      result[ip*18 +  1] = (-1.0 + x)*x*(-0.25 + 0.5*y);
      result[ip*18 +  2] = (0.25 + 0.5*x)*(-1. + y)*y;
      result[ip*18 +  3] = x*(1. + x)*(-0.25 + 0.5*y);
      result[ip*18 +  4] = (0.25 + 0.5*x)*y*(1. + y);
      result[ip*18 +  5] = x*(1. + x)*(0.25 + 0.5*y);
      result[ip*18 +  6] = (-0.25 + 0.5*x)*y*(1. + y);
      result[ip*18 +  7] = (-1. + x)*x*(0.25 + 0.5*y);
      result[ip*18 +  8] = x*(1.0 - y)*y;
      result[ip*18 +  9] = 0.5*(1.0 - x)*(1.0 + x)*(-1.0 + 2.0*y);
      result[ip*18 + 10] = 0.5*(1.0 - y)*(1.0 + y)*(1.0 + 2.0*x);
      result[ip*18 + 11] = -x*(1.0 + x)*y;
      result[ip*18 + 12] = -y*(1.0 + y)*x;
      result[ip*18 + 13] = 0.5*(1.0 - x)*(1.0 + x)*(1.0 + 2.0*y);
      result[ip*18 + 14] = 0.5*(1.0 - y)*(1.0+ y)*(-1.0 + 2.0*x);
      result[ip*18 + 15] = (1.0 - x)*x*y;
      result[ip*18 + 16] = -2.0*(1.0 - y)*(1.0 + y)*x;
      result[ip*18 + 17] = -2.0*(1.0 - x)*(1.0 + x)*y;
    }
  }
};

class Basis_TET_4 : public Basis
{
public:
  Basis_TET_4() : Basis(1) {}

  double parametric_volume() const override { return 1./6.; }

  void nodal_parametric_coordinates(double* ATTR_RESTRICT p_coords) const override
  {
    p_coords[ 0] = 0.0; p_coords[ 1] = 0.0; p_coords[ 2] = 0.0;
    p_coords[ 3] = 1.0; p_coords[ 4] = 0.0; p_coords[ 5] = 0.0;
    p_coords[ 6] = 0.0; p_coords[ 7] = 1.0; p_coords[ 8] = 0.0;
    p_coords[ 9] = 0.0; p_coords[10] = 0.0; p_coords[11] = 1.0;
  }
  void shape_fcn(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for (int ip(0); ip < nint; ++ip)
    {
      const double x = p_coords[3*ip + 0];
      const double y = p_coords[3*ip + 1];
      const double z = p_coords[3*ip + 2];
      result[ip*4 + 0] = 1.0 - x - y - z;
      result[ip*4 + 1] = x;
      result[ip*4 + 2] = y;
      result[ip*4 + 3] = z;
    }
  }
  void shape_fcn_deriv(const int nint, const double* ATTR_RESTRICT /*p_coords*/, double* ATTR_RESTRICT result) const override
  {
    for ( int ip(0); ip < nint; ++ip )
    {
      result[ip*12 +  0] = -1.0;
      result[ip*12 +  1] = -1.0;
      result[ip*12 +  2] = -1.0;
      result[ip*12 +  3] =  1.0;
      result[ip*12 +  4] =  0.0;
      result[ip*12 +  5] =  0.0;
      result[ip*12 +  6] =  0.0;
      result[ip*12 +  7] =  1.0;
      result[ip*12 +  8] =  0.0;
      result[ip*12 +  9] =  0.0;
      result[ip*12 + 10] =  0.0;
      result[ip*12 + 11] =  1.0;
    }
  }
};

class Basis_TET_10 : public Basis
{
public:
  Basis_TET_10() : Basis(2) {}

  double parametric_volume() const override { return 1./6.; }

  void nodal_parametric_coordinates(double* ATTR_RESTRICT p_coords) const override
  {
    p_coords[ 0] = 0.0; p_coords[ 1] = 0.0; p_coords[ 2] = 0.0;
    p_coords[ 3] = 1.0; p_coords[ 4] = 0.0; p_coords[ 5] = 0.0;
    p_coords[ 6] = 0.0; p_coords[ 7] = 1.0; p_coords[ 8] = 0.0;
    p_coords[ 9] = 0.0; p_coords[10] = 0.0; p_coords[11] = 1.0;
    p_coords[12] = 0.5; p_coords[13] = 0.0; p_coords[14] = 0.0;
    p_coords[15] = 0.5; p_coords[16] = 0.5; p_coords[17] = 0.0;
    p_coords[18] = 0.0; p_coords[19] = 0.5; p_coords[20] = 0.0;
    p_coords[21] = 0.0; p_coords[22] = 0.0; p_coords[23] = 0.5;
    p_coords[24] = 0.5; p_coords[25] = 0.0; p_coords[26] = 0.5;
    p_coords[27] = 0.0; p_coords[28] = 0.5; p_coords[29] = 0.5;
  }
  void shape_fcn(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for (int ip(0); ip < nint; ++ip)
    {
      const double x = p_coords[3*ip + 0];
      const double y = p_coords[3*ip + 1];
      const double z = p_coords[3*ip + 2];
      result[ip*10 + 0] = (-1. + x + y + z)*(-1. + 2.*x + 2.*y + 2.*z);
      result[ip*10 + 1] = x*(-1. + 2.*x);
      result[ip*10 + 2] = y*(-1. + 2.*y);
      result[ip*10 + 3] = z*(-1. + 2.*z);
      result[ip*10 + 4] = -4.*x*(-1. + x + y + z);
      result[ip*10 + 5] = 4.*x*y;
      result[ip*10 + 6] = -4.*y*(-1. + x + y + z);
      result[ip*10 + 7] = -4.*z*(-1. + x + y + z);
      result[ip*10 + 8] = 4.*x*z;
      result[ip*10 + 9] = 4.*y*z;
    }
  }
  void shape_fcn_deriv(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for ( int ip(0); ip < nint; ++ip )
    {
      const double x = p_coords[3*ip + 0];
      const double y = p_coords[3*ip + 1];
      const double z = p_coords[3*ip + 2];
      result[ip*30 +  0] = -3.+ 4.*x + 4.*y + 4.*z;
      result[ip*30 +  1] = -3.+ 4.*x + 4.*y + 4.*z;
      result[ip*30 +  2] = -3.+ 4.*x + 4.*y + 4.*z;
      result[ip*30 +  3] = -1.+ 4.*x;
      result[ip*30 +  4] =  0.;
      result[ip*30 +  5] =  0.;
      result[ip*30 +  6] =  0.;
      result[ip*30 +  7] = -1.+ 4.*y;
      result[ip*30 +  8] =  0.;
      result[ip*30 +  9] =  0.;
      result[ip*30 + 10] =  0.;
      result[ip*30 + 11] = -1.+ 4.*z;
      result[ip*30 + 12] = -4.*(-1.+ 2*x + y + z);
      result[ip*30 + 13] = -4.*x;
      result[ip*30 + 14] = -4.*x;
      result[ip*30 + 15] =  4.*y;
      result[ip*30 + 16] =  4.*x;
      result[ip*30 + 17] =  0.;
      result[ip*30 + 18] = -4.*y;
      result[ip*30 + 19] = -4.*(-1.+ x + 2*y + z);
      result[ip*30 + 20] = -4.*y;
      result[ip*30 + 21] = -4.*z;
      result[ip*30 + 22] = -4.*z;
      result[ip*30 + 23] = -4.*(-1.+ x + y + 2*z);
      result[ip*30 + 24] =  4.*z;
      result[ip*30 + 25] =  0.;
      result[ip*30 + 26] =  4.*x;
      result[ip*30 + 27] =  0.;
      result[ip*30 + 28] =  4.*z;
      result[ip*30 + 29] =  4.*y;
    }
  }
};

class Basis_HEX_8 : public Basis
{
public:
  Basis_HEX_8() : Basis(1) {}

  double parametric_volume() const override { return 8.; }

  void nodal_parametric_coordinates(double* ATTR_RESTRICT p_coords) const override
  {
    p_coords[ 0] = -1.0; p_coords[ 1] = -1.0; p_coords[ 2] = -1.0;
    p_coords[ 3] =  1.0; p_coords[ 4] = -1.0; p_coords[ 5] = -1.0;
    p_coords[ 6] =  1.0; p_coords[ 7] =  1.0; p_coords[ 8] = -1.0;
    p_coords[ 9] = -1.0; p_coords[10] =  1.0; p_coords[11] = -1.0;
    p_coords[12] = -1.0; p_coords[13] = -1.0; p_coords[14] =  1.0;
    p_coords[15] =  1.0; p_coords[16] = -1.0; p_coords[17] =  1.0;
    p_coords[18] =  1.0; p_coords[19] =  1.0; p_coords[20] =  1.0;
    p_coords[21] = -1.0; p_coords[22] =  1.0; p_coords[23] =  1.0;
  }
  void shape_fcn(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for (int ip(0); ip < nint; ++ip)
    {
      const double x = p_coords[3*ip + 0];
      const double y = p_coords[3*ip + 1];
      const double z = p_coords[3*ip + 2];
      result[ip*8 + 0] = 0.125*(1.0 - x)*(1.0 - y)*(1.0 - z);
      result[ip*8 + 1] = 0.125*(1.0 + x)*(1.0 - y)*(1.0 - z);
      result[ip*8 + 2] = 0.125*(1.0 + x)*(1.0 + y)*(1.0 - z);
      result[ip*8 + 3] = 0.125*(1.0 - x)*(1.0 + y)*(1.0 - z);
      result[ip*8 + 4] = 0.125*(1.0 - x)*(1.0 - y)*(1.0 + z);
      result[ip*8 + 5] = 0.125*(1.0 + x)*(1.0 - y)*(1.0 + z);
      result[ip*8 + 6] = 0.125*(1.0 + x)*(1.0 + y)*(1.0 + z);
      result[ip*8 + 7] = 0.125*(1.0 - x)*(1.0 + y)*(1.0 + z);
    }
  }
  void shape_fcn_deriv(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for ( int ip(0); ip < nint; ++ip )
    {
      const double x = p_coords[3*ip + 0];
      const double y = p_coords[3*ip + 1];
      const double z = p_coords[3*ip + 2];
      result[ip*24 +  0] = -(1.0 - y)*(1.0 - z)*0.125;
      result[ip*24 +  1] = -(1.0 - x)*(1.0 - z)*0.125;
      result[ip*24 +  2] = -(1.0 - x)*(1.0 - y)*0.125;
      result[ip*24 +  3] =  (1.0 - y)*(1.0 - z)*0.125;
      result[ip*24 +  4] = -(1.0 + x)*(1.0 - z)*0.125;
      result[ip*24 +  5] = -(1.0 + x)*(1.0 - y)*0.125;
      result[ip*24 +  6] =  (1.0 + y)*(1.0 - z)*0.125;
      result[ip*24 +  7] =  (1.0 + x)*(1.0 - z)*0.125;
      result[ip*24 +  8] = -(1.0 + x)*(1.0 + y)*0.125;
      result[ip*24 +  9] = -(1.0 + y)*(1.0 - z)*0.125;
      result[ip*24 + 10] =  (1.0 - x)*(1.0 - z)*0.125;
      result[ip*24 + 11] = -(1.0 - x)*(1.0 + y)*0.125;
      result[ip*24 + 12] = -(1.0 - y)*(1.0 + z)*0.125;
      result[ip*24 + 13] = -(1.0 - x)*(1.0 + z)*0.125;
      result[ip*24 + 14] =  (1.0 - x)*(1.0 - y)*0.125;
      result[ip*24 + 15] =  (1.0 - y)*(1.0 + z)*0.125;
      result[ip*24 + 16] = -(1.0 + x)*(1.0 + z)*0.125;
      result[ip*24 + 17] =  (1.0 + x)*(1.0 - y)*0.125;
      result[ip*24 + 18] =  (1.0 + y)*(1.0 + z)*0.125;
      result[ip*24 + 19] =  (1.0 + x)*(1.0 + z)*0.125;
      result[ip*24 + 20] =  (1.0 + x)*(1.0 + y)*0.125;
      result[ip*24 + 21] = -(1.0 + y)*(1.0 + z)*0.125;
      result[ip*24 + 22] =  (1.0 - x)*(1.0 + z)*0.125;
      result[ip*24 + 23] =  (1.0 - x)*(1.0 + y)*0.125;
    }
  }
};

class Basis_HEX_27 : public Basis
{
public:
  Basis_HEX_27() : Basis(1) {}

  double parametric_volume() const override { return 8.; }

  void nodal_parametric_coordinates(double* ATTR_RESTRICT p_coords) const override
  {
    p_coords[ 0] = -1.0; p_coords[ 1] = -1.0; p_coords[ 2] = -1.0;
    p_coords[ 3] =  1.0; p_coords[ 4] = -1.0; p_coords[ 5] = -1.0;
    p_coords[ 6] =  1.0; p_coords[ 7] =  1.0; p_coords[ 8] = -1.0;
    p_coords[ 9] = -1.0; p_coords[10] =  1.0; p_coords[11] = -1.0;
    p_coords[12] = -1.0; p_coords[13] = -1.0; p_coords[14] =  1.0;
    p_coords[15] =  1.0; p_coords[16] = -1.0; p_coords[17] =  1.0;
    p_coords[18] =  1.0; p_coords[19] =  1.0; p_coords[20] =  1.0;
    p_coords[21] = -1.0; p_coords[22] =  1.0; p_coords[23] =  1.0;
    p_coords[24] =  0.0; p_coords[25] = -1.0; p_coords[26] = -1.0;
    p_coords[27] =  1.0; p_coords[28] =  0.0; p_coords[29] = -1.0;
    p_coords[30] =  0.0; p_coords[31] =  1.0; p_coords[32] = -1.0;
    p_coords[33] = -1.0; p_coords[34] =  0.0; p_coords[35] = -1.0;
    p_coords[36] = -1.0; p_coords[37] = -1.0; p_coords[38] =  0.0;
    p_coords[39] =  1.0; p_coords[40] = -1.0; p_coords[41] =  0.0;
    p_coords[42] =  1.0; p_coords[43] =  1.0; p_coords[44] =  0.0;
    p_coords[45] = -1.0; p_coords[46] =  1.0; p_coords[47] =  0.0;
    p_coords[48] =  0.0; p_coords[49] = -1.0; p_coords[50] =  1.0;
    p_coords[51] =  1.0; p_coords[52] =  0.0; p_coords[53] =  1.0;
    p_coords[54] =  0.0; p_coords[55] =  1.0; p_coords[56] =  1.0;
    p_coords[57] = -1.0; p_coords[58] =  0.0; p_coords[59] =  1.0;
    p_coords[60] =  0.0; p_coords[61] =  0.0; p_coords[62] =  0.0;
    p_coords[63] =  0.0; p_coords[64] =  0.0; p_coords[65] = -1.0;
    p_coords[66] =  0.0; p_coords[67] =  0.0; p_coords[68] =  1.0;
    p_coords[69] = -1.0; p_coords[70] =  0.0; p_coords[71] =  0.0;
    p_coords[72] =  1.0; p_coords[73] =  0.0; p_coords[74] =  0.0;
    p_coords[75] =  0.0; p_coords[76] = -1.0; p_coords[77] =  0.0;
    p_coords[78] =  0.0; p_coords[79] =  1.0; p_coords[80] =  0.0;
  }
  void shape_fcn(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for (int ip(0); ip < nint; ++ip)
    {
      const double x = p_coords[3*ip + 0];
      const double y = p_coords[3*ip + 1];
      const double z = p_coords[3*ip + 2];
      result[ip*27 +  0] = 0.125*(-1. + x)*x*(-1. + y)*y*(-1. + z)*z;
      result[ip*27 +  1] = 0.125*x*(1.+ x)*(-1. + y)*y*(-1. + z)*z;
      result[ip*27 +  2] = 0.125*x*(1.+ x)*y*(1.+ y)*(-1. + z)*z;
      result[ip*27 +  3] = 0.125*(-1. + x)*x*y*(1.+ y)*(-1. + z)*z;
      result[ip*27 +  4] = 0.125*(-1. + x)*x*(-1. + y)*y*z*(1.+ z);
      result[ip*27 +  5] = 0.125*x*(1.+ x)*(-1. + y)*y*z*(1.+ z);
      result[ip*27 +  6] = 0.125*x*(1.+ x)*y*(1.+ y)*z*(1.+ z);
      result[ip*27 +  7] = 0.125*(-1. + x)*x*y*(1.+ y)*z*(1.+ z);
      result[ip*27 +  8] = 0.25*(1. - x)*(1. + x)*(-1. + y)*y*(-1. + z)*z;
      result[ip*27 +  9] = 0.25*x*(1.+ x)*(1. - y)*(1. + y)*(-1. + z)*z;
      result[ip*27 + 10] = 0.25*(1. - x)*(1. + x)*y*(1.+ y)*(-1. + z)*z;
      result[ip*27 + 11] = 0.25*(-1. + x)*x*(1. - y)*(1. + y)*(-1. + z)*z;
      result[ip*27 + 12] = 0.25*(-1. + x)*x*(-1. + y)*y*(1. - z)*(1. + z);
      result[ip*27 + 13] = 0.25*x*(1.+ x)*(-1. + y)*y*(1. - z)*(1. + z);
      result[ip*27 + 14] = 0.25*x*(1.+ x)*y*(1.+ y)*(1. - z)*(1. + z);
      result[ip*27 + 15] = 0.25*(-1. + x)*x*y*(1.+ y)*(1. - z)*(1. + z);
      result[ip*27 + 16] = 0.25*(1. - x)*(1. + x)*(-1. + y)*y*z*(1.+ z);
      result[ip*27 + 17] = 0.25*x*(1.+ x)*(1. - y)*(1. + y)*z*(1.+ z);
      result[ip*27 + 18] = 0.25*(1. - x)*(1. + x)*y*(1.+ y)*z*(1.+ z);
      result[ip*27 + 19] = 0.25*(-1. + x)*x*(1. - y)*(1. + y)*z*(1.+ z);
      result[ip*27 + 20] = (1. - x)*(1. + x)*(1. - y)*(1. + y)*(1. - z)*(1. + z);
      result[ip*27 + 21] = 0.5*(1. - x)*(1. + x)*(1. - y)*(1. + y)*(-1. + z)*z;
      result[ip*27 + 22] = 0.5*(1. - x)*(1. + x)*(1. - y)*(1. + y)*z*(1.+ z);
      result[ip*27 + 23] = 0.5*(-1. + x)*x*(1. - y)*(1. + y)*(1. - z)*(1. + z);
      result[ip*27 + 24] = 0.5*x*(1.+ x)*(1. - y)*(1. + y)*(1. - z)*(1. + z);
      result[ip*27 + 25] = 0.5*(1. - x)*(1. + x)*(-1. + y)*y*(1. - z)*(1. + z);
      result[ip*27 + 26] = 0.5*(1. - x)*(1. + x)*y*(1.+ y)*(1. - z)*(1. + z);
    }
  }
  void shape_fcn_deriv(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for ( int ip(0); ip < nint; ++ip )
    {
      const double x = p_coords[3*ip + 0];
      const double y = p_coords[3*ip + 1];
      const double z = p_coords[3*ip + 2];
      result[ip*81 +  0] = (-0.125 + 0.25*x)*(-1. + y)*y*(-1. + z)*z;
      result[ip*81 +  1] = (-1. + x)*x*(-0.125 + 0.25*y)*(-1. + z)*z;
      result[ip*81 +  2] = (-1. + x)*x*(-1. + y)*y*(-0.125 + 0.25*z);
      result[ip*81 +  3] = (0.125 + 0.25*x)*(-1. + y)*y*(-1. + z)*z;
      result[ip*81 +  4] = x*(1. + x)*(-0.125 + 0.25*y)*(-1. + z)*z;
      result[ip*81 +  5] = x*(1. + x)*(-1. + y)*y*(-0.125 + 0.25*z);
      result[ip*81 +  6] = (0.125 + 0.25*x)*y*(1. + y)*(-1. + z)*z;
      result[ip*81 +  7] = x*(1. + x)*(0.125 + 0.25*y)*(-1. + z)*z;
      result[ip*81 +  8] = x*(1. + x)*y*(1. + y)*(-0.125 + 0.25*z);
      result[ip*81 +  9] = (-0.125 + 0.25*x)*y*(1. + y)*(-1. + z)*z;
      result[ip*81 + 10] = (-1. + x)*x*(0.125 + 0.25*y)*(-1. + z)*z;
      result[ip*81 + 11] = (-1. + x)*x*y*(1. + y)*(-0.125 + 0.25*z);
      result[ip*81 + 12] = (-0.125 + 0.25*x)*(-1. + y)*y*z*(1. + z);
      result[ip*81 + 13] = (-1. + x)*x*(-0.125 + 0.25*y)*z*(1. + z);
      result[ip*81 + 14] = (-1. + x)*x*(-1. + y)*y*(0.125 + 0.25*z);
      result[ip*81 + 15] = (0.125 + 0.25*x)*(-1. + y)*y*z*(1. + z);
      result[ip*81 + 16] = x*(1. + x)*(-0.125 + 0.25*y)*z*(1. + z);
      result[ip*81 + 17] = x*(1. + x)*(-1. + y)*y*(0.125 + 0.25*z);
      result[ip*81 + 18] = (0.125 + 0.25*x)*y*(1. + y)*z*(1. + z);
      result[ip*81 + 19] = x*(1. + x)*(0.125 + 0.25*y)*z*(1. + z);
      result[ip*81 + 20] = x*(1. + x)*y*(1. + y)*(0.125 + 0.25*z);
      result[ip*81 + 21] = (-0.125 + 0.25*x)*y*(1. + y)*z*(1. + z);
      result[ip*81 + 22] = (-1. + x)*x*(0.125 + 0.25*y)*z*(1. + z);
      result[ip*81 + 23] = (-1. + x)*x*y*(1. + y)*(0.125 + 0.25*z);
      result[ip*81 + 24] = -0.5*x*(-1. + y)*y*(-1. + z)*z;
      result[ip*81 + 25] = (1. - x)*(1. + x)*(-0.25 + 0.5*y)*(-1. + z)*z;
      result[ip*81 + 26] = (1. - x)*(1. + x)*(-1. + y)*y*(-0.25 + 0.5*z);
      result[ip*81 + 27] = (0.25 + 0.5*x)*(1. - y)*(1. + y)*(-1. + z)*z;
      result[ip*81 + 28] = x*(1. + x)*(-0.5*y)*(-1. + z)*z;
      result[ip*81 + 29] = x*(1. + x)*(1. - y)*(1. + y)*(-0.25 + 0.5*z);
      result[ip*81 + 30] = -0.5*x*y*(1. + y)*(-1. + z)*z;
      result[ip*81 + 31] = (1. - x)*(1. + x)*(0.25 + 0.5*y)*(-1. + z)*z;
      result[ip*81 + 32] = (1. - x)*(1. + x)*y*(1. + y)*(-0.25 + 0.5*z);
      result[ip*81 + 33] = (-0.25 + 0.5*x)*(1. - y)*(1. + y)*(-1. + z)*z;
      result[ip*81 + 34] = (-1. + x)*x*(-0.5*y)*(-1. + z)*z;
      result[ip*81 + 35] = (-1. + x)*x*(1. - y)*(1. + y)*(-0.25 + 0.5*z);
      result[ip*81 + 36] = (-0.25 + 0.5*x)*(-1. + y)*y*(1. - z)*(1. + z);
      result[ip*81 + 37] = (-1. + x)*x*(-0.25 + 0.5*y)*(1. - z)*(1. + z);
      result[ip*81 + 38] = (-1. + x)*x*(-1. + y)*y*(-0.5*z);
      result[ip*81 + 39] = (0.25 + 0.5*x)*(-1. + y)*y*(1. - z)*(1. + z);
      result[ip*81 + 40] = x*(1. + x)*(-0.25 + 0.5*y)*(1. - z)*(1. + z);
      result[ip*81 + 41] = x*(1. + x)*(-1. + y)*y*(-0.5*z);
      result[ip*81 + 42] = (0.25 + 0.5*x)*y*(1. + y)*(1. - z)*(1. + z);
      result[ip*81 + 43] = x*(1. + x)*(0.25 + 0.5*y)*(1. - z)*(1. + z);
      result[ip*81 + 44] = x*(1. + x)*y*(1. + y)*(-0.5*z);
      result[ip*81 + 45] = (-0.25 + 0.5*x)*y*(1. + y)*(1. - z)*(1. + z);
      result[ip*81 + 46] = (-1. + x)*x*(0.25 + 0.5*y)*(1. - z)*(1. + z);
      result[ip*81 + 47] = (-1. + x)*x*y*(1. + y)*(-0.5*z);
      result[ip*81 + 48] = -0.5*x*(-1. + y)*y*z*(1. + z);
      result[ip*81 + 49] = (1. - x)*(1. + x)*(-0.25 + 0.5*y)*z*(1. + z);
      result[ip*81 + 50] = (1. - x)*(1. + x)*(-1. + y)*y*(0.25 + 0.5*z);
      result[ip*81 + 51] = (0.25 + 0.5*x)*(1. - y)*(1. + y)*z*(1. + z);
      result[ip*81 + 52] = x*(1. + x)*(-0.5*y)*z*(1. + z);
      result[ip*81 + 53] = x*(1. + x)*(1. - y)*(1. + y)*(0.25 + 0.5*z);
      result[ip*81 + 54] = -0.5*x*y*(1. + y)*z*(1. + z);
      result[ip*81 + 55] = (1. - x)*(1. + x)*(0.25 + 0.5*y)*z*(1. + z);
      result[ip*81 + 56] = (1. - x)*(1. + x)*y*(1. + y)*(0.25 + 0.5*z);
      result[ip*81 + 57] = (-0.25 + 0.5*x)*(1. - y)*(1. + y)*z*(1. + z);
      result[ip*81 + 58] = (-1. + x)*x*(-0.5*y)*z*(1. + z);
      result[ip*81 + 59] = (-1. + x)*x*(1. - y)*(1. + y)*(0.25 + 0.5*z);
      result[ip*81 + 60] = -2.*x*(1. - y)*(1. + y)*(1. - z)*(1. + z);
      result[ip*81 + 61] = (1. - x)*(1. + x)*(-2.*y)*(1. - z)*(1. + z);
      result[ip*81 + 62] = (1. - x)*(1. + x)*(1. - y)*(1. + y)*(-2.*z);
      result[ip*81 + 63] = -x*(1. - y)*(1. + y)*(-1. + z)*z;
      result[ip*81 + 64] = (1. - x)*(1. + x)*(-y)*(-1. + z)*z;
      result[ip*81 + 65] = (1. - x)*(1. + x)*(1. - y)*(1. + y)*(-0.5 + z);
      result[ip*81 + 66] = -x*(1. - y)*(1. + y)*z*(1. + z);
      result[ip*81 + 67] = (1. - x)*(1. + x)*(-y)*z*(1. + z);
      result[ip*81 + 68] = (1. - x)*(1. + x)*(1. - y)*(1. + y)*(0.5 + z);
      result[ip*81 + 69] = (-0.5 + x)*(1. - y)*(1. + y)*(1. - z)*(1. + z);
      result[ip*81 + 70] = (-1. + x)*x*(-y)*(1. - z)*(1. + z);
      result[ip*81 + 71] = (-1. + x)*x*(1. - y)*(1. + y)*(-z);
      result[ip*81 + 72] = (0.5 + x)*(1. - y)*(1. + y)*(1. - z)*(1. + z);
      result[ip*81 + 73] = x*(1. + x)*(-y)*(1. - z)*(1. + z);
      result[ip*81 + 74] = x*(1. + x)*(1. - y)*(1. + y)*(-z);
      result[ip*81 + 75] = -x*(-1. + y)*y*(1. - z)*(1. + z);
      result[ip*81 + 76] = (1. - x)*(1. + x)*(-0.5 + y)*(1. - z)*(1. + z);
      result[ip*81 + 77] = (1. - x)*(1. + x)*(-1. + y)*y*(-z);
      result[ip*81 + 78] = -x*y*(1. + y)*(1. - z)*(1. + z);
      result[ip*81 + 79] = (1. - x)*(1. + x)*(0.5 + y)*(1. - z)*(1. + z);
      result[ip*81 + 80] = (1. - x)*(1. + x)*y*(1. + y)*(-z);
    }
  }
};

class Basis_WEDGE_6 : public Basis
{
public:
  Basis_WEDGE_6() : Basis(1) {}

  double parametric_volume() const override { return 1.; }

  void nodal_parametric_coordinates(double* ATTR_RESTRICT p_coords) const override
  {
    p_coords[ 0] =  0.0; p_coords[ 1] =  0.0; p_coords[ 2] = -1.0;
    p_coords[ 3] =  1.0; p_coords[ 4] =  0.0; p_coords[ 5] = -1.0;
    p_coords[ 6] =  0.0; p_coords[ 7] =  1.0; p_coords[ 8] = -1.0;
    p_coords[ 9] =  0.0; p_coords[10] =  0.0; p_coords[11] =  1.0;
    p_coords[12] =  1.0; p_coords[13] =  0.0; p_coords[14] =  1.0;
    p_coords[15] =  0.0; p_coords[16] =  1.0; p_coords[17] =  1.0;
  }
  void shape_fcn(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for (int ip(0); ip < nint; ++ip)
    {
      const double x = p_coords[3*ip + 0];
      const double y = p_coords[3*ip + 1];
      const double z = p_coords[3*ip + 2];
      const double t = 1. - x - y;

      result[ip*6 + 0] = 0.5 * t * (1.0 - z);
      result[ip*6 + 1] = 0.5 * x * (1.0 - z);
      result[ip*6 + 2] = 0.5 * y * (1.0 - z);
      result[ip*6 + 3] = 0.5 * t * (1.0 + z);
      result[ip*6 + 4] = 0.5 * x * (1.0 + z);
      result[ip*6 + 5] = 0.5 * y * (1.0 + z);
    }
  }
  void shape_fcn_deriv(const int nint, const double* ATTR_RESTRICT p_coords, double* ATTR_RESTRICT result) const override
  {
    for ( int ip(0); ip < nint; ++ip )
    {
      const double x = p_coords[3*ip + 0];
      const double y = p_coords[3*ip + 1];
      const double z = p_coords[3*ip + 2];
      const double t = 1. - x - y;

      result[ip*18 +  0] = -0.5 * (1.0 - z);
      result[ip*18 +  1] = -0.5 * (1.0 - z);
      result[ip*18 +  2] = -0.5 * t;
      result[ip*18 +  3] =  0.5 * (1.0 - z);
      result[ip*18 +  4] =  0.;
      result[ip*18 +  5] = -0.5 * x;
      result[ip*18 +  6] =  0.;
      result[ip*18 +  7] =  0.5 * (1.0 - z);
      result[ip*18 +  8] = -0.5 * y;
      result[ip*18 +  9] = -0.5 * (1.0 + z);
      result[ip*18 + 10] = -0.5 * (1.0 + z);
      result[ip*18 + 11] =  0.5 * t;
      result[ip*18 + 12] =  0.5 * (1.0 + z);
      result[ip*18 + 13] =  0.;
      result[ip*18 + 14] =  0.5 * x;
      result[ip*18 + 15] =  0.;
      result[ip*18 + 16] =  0.5 * (1.0 + z);
      result[ip*18 + 17] =  0.5 * y;
    }
  }
};

}  // end namespace krino

#endif // Akri_MasterElementBasis_h
