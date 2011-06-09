/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Glen Hansen (Glen.Hansen@inl.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file mrtr_functions.H
 *
 * \brief A couple of shape functions for 1D and 2D elements
 *
 * \date Last update do Doxygen: 20-March-06
 *
 */
#ifndef MOERTEL_FUNCTIONS_H
#define MOERTEL_FUNCTIONS_H

#include <ctime>
#include <iostream>
#include <cstdlib>

#include "mrtr_function.H"

/*!
\brief MOERTEL: namespace of the Moertel package

The Moertel package depends on \ref Epetra, \ref EpetraExt, \ref Teuchos,
\ref Amesos, \ref ML and \ref AztecOO:<br>
Use at least the following lines in the configure of Trilinos:<br>
\code
--enable-moertel 
--enable-epetra 
--enable-epetraext
--enable-teuchos 
--enable-ml
--enable-aztecoo --enable-aztecoo-teuchos 
--enable-amesos
\endcode

*/
namespace MOERTEL
{

/*!
\class Function_Constant1D

\brief <b> A 1D constant shape function for a 2-noded 1D Segment</b>

phi_1 = 1 <br>
phi_2 = 1 <br>
phi_1,xi = 0 <br>
phi_2,xi = 0 <br>

\author Glen Hansen (Glen.Hansen@inl.gov)

*/
class Function_Constant1D : public MOERTEL::Function 
{
public:
  
  /*!
  \brief Constructor
  
  \param out : Level of output information written to stdout ( 0 - 10 )
  */
  Function_Constant1D(int out) :
  MOERTEL::Function(out,MOERTEL::Function::func_Constant1D)
  {}
  
  /*!
  \brief Copy-Constructor
  */
  Function_Constant1D(const MOERTEL::Function_Constant1D& old) :
  MOERTEL::Function(old.outputlevel_,old.type_)
  {}

  /*!
  \brief Clone method
  
  Makes a deep copy of this instance and returns a pointer to it.
  
  */
  MOERTEL::Function* Clone() const;

  /*!
  \brief Destructor
  */
  virtual ~Function_Constant1D() {}

  /*!
  \brief Evaluate the function and derivatives at a given local coordinate
  
  \param xi (in) : Local coordinate where to evaluate the function
  \param val (out) : Function values at xi, dimension 2, if NULL on input no evaluation
  \param valdim (in) : Dimension of val should be 2
  \param deriv (out) : Derivative of functions at xi, dimension 4, if NULL on input no evaluation
  */
  bool EvaluateFunction(const MOERTEL::Segment& seg, const double* xi, 
                        double* val, const int valdim, double* deriv);

};

/*!
\class Function_Linear1D

\brief <b> A 1D linear shape function </b>

phi_1 = 0.5*(1-xi) <br>
phi_2 = 0.5*(1+xi) <br>
phi_1,xi = -0.5 <br>
phi_2,xi = 0.5 <br>

\author Glen Hansen (Glen.Hansen@inl.gov)

*/
class Function_Linear1D : public MOERTEL::Function 
{
public:
  
  /*!
  \brief Constructor
  
  \param out : Level of output information written to stdout ( 0 - 10 )
  */
  Function_Linear1D(int out) :
  MOERTEL::Function(out,MOERTEL::Function::func_Linear1D)
  {}
  
  /*!
  \brief Copy-Constructor
  */
  Function_Linear1D(const MOERTEL::Function_Linear1D& old) :
  MOERTEL::Function(old.outputlevel_,old.type_)
  {}

  /*!
  \brief Clone method
  
  Makes a deep copy of this instance and returns a pointer to it.
  
  */
  MOERTEL::Function* Clone() const;

  /*!
  \brief Destructor
  */
  virtual ~Function_Linear1D() {}

  /*!
  \brief Evaluate the function and derivatives at a given local coordinate
  
  \param xi (in) : Local coordinate where to evaluate the function
  \param val (out) : Function values at xi, dimension 2, if NULL on input no evaluation
  \param valdim (in) : Dimension of val should be 2
  \param deriv (out) : Derivative of functions at xi, dimension 4, if NULL on input no evaluation
  */
  bool EvaluateFunction(const MOERTEL::Segment& seg, const double* xi, 
                        double* val, const int valdim, double* deriv);

};

/*!
\class Function_DualLinear1D

\brief <b> A 1D dual linear shape function </b>

phi_1 = -1.5*xi+0.5 <br>
phi_2 = 1.5*xi+0.5 <br>
phi_1,xi = -1.5 <br>
phi_2,xi = 1.5 <br>

\author Glen Hansen (Glen.Hansen@inl.gov)

*/
class Function_DualLinear1D : public MOERTEL::Function 
{
public:
  
  /*!
  \brief Constructor
  
  \param out : Level of output information written to stdout ( 0 - 10 )
  */
  Function_DualLinear1D(int out) :
  MOERTEL::Function(out,MOERTEL::Function::func_DualLinear1D)
  {}
  
  /*!
  \brief Copy-Constructor
  */
  Function_DualLinear1D(const MOERTEL::Function_DualLinear1D& old) :
  MOERTEL::Function(old.outputlevel_,old.type_)
  {}

  /*!
  \brief Clone method
  
  Makes a deep copy of this instance and returns a pointer to it.
  
  */
  virtual MOERTEL::Function* Clone() const;

  /*!
  \brief Destructor
  */
  virtual ~Function_DualLinear1D() {}

  /*!
  \brief Evaluate the function and derivatives at a given local coordinate
  
  \param xi (in) : Local coordinate where to evaluate the function
  \param val (out) : Function values at xi, dimension 2, if NULL on input no evaluation
  \param valdim (in) : Dimension of val should be 2
  \param deriv (out) : Derivative of functions at xi, dimension 4, if NULL on input no evaluation
  */
  bool EvaluateFunction(const MOERTEL::Segment& seg, const double* xi, 
                        double* val, const int valdim, double* deriv);

};

/*!
\class Function_LinearTri

\brief <b> A 2D linear shape function for a 3-noded triangle</b>

phi_1 = 1 - xi_1 - xi_2 <br>
phi_2 = xi_1 <br>
phi_3 = xi_2 <br>
phi_1,xi_1 = -1 <br>
phi_1,xi_2 = -1 <br>
phi_2,xi_1 = 1 <br>
phi_2,xi_2 = 0 <br>
phi_3,xi_1 = 0 <br>
phi_3,xi_2 = 1 <br>

\author Glen Hansen (Glen.Hansen@inl.gov)

*/
class Function_LinearTri : public MOERTEL::Function 
{
public:
  
  /*!
  \brief Constructor
  
  \param out : Level of output information written to stdout ( 0 - 10 )
  */
  Function_LinearTri(int out) :
  MOERTEL::Function(out,MOERTEL::Function::func_LinearTri)
  {}
  
  /*!
  \brief Copy-Constructor
  */
  Function_LinearTri(const MOERTEL::Function_LinearTri& old) :
  MOERTEL::Function(old.outputlevel_,old.type_)
  {}

  /*!
  \brief Clone method
  
  Makes a deep copy of this instance and returns a pointer to it.
  
  */
  virtual MOERTEL::Function* Clone() const;

  /*!
  \brief Destructor
  */
  virtual ~Function_LinearTri() {}

  /*!
  \brief Evaluate the function and derivatives at a given local coordinate
  
  \param xi (in) : Local coordinate where to evaluate the function (dimension 2)
  \param val (out) : Function values at xi, dimension 3, if NULL on input no evaluation
  \param valdim (in) : Dimension of val should be 3
  \param deriv (out) : Derivative of functions at xi, dimension 6, if NULL on input no evaluation
  */
  bool EvaluateFunction(const MOERTEL::Segment& seg, const double* xi, 
                        double* val, const int valdim, double* deriv);
};

/*!
\class Function_DualLinearTri

\brief <b> A 2D dual linear shape function for a 3-noded triangle</b>

phi_1 = 3 - 2 * xi_1 - 2 * xi_2 <br>
phi_2 = 4 * xi_1 - 1<br>
phi_3 = 4 * xi_2 - 1 <br>
phi_1,xi_1 = -2 <br>
phi_1,xi_2 = -2 <br>
phi_2,xi_1 = 4 <br>
phi_2,xi_2 = 0 <br>
phi_3,xi_1 = 0 <br>
phi_3,xi_2 = 4 <br>

\author Glen Hansen (Glen.Hansen@inl.gov)

*/
class Function_DualLinearTri : public MOERTEL::Function 
{
public:
  
  /*!
  \brief Constructor
  
  \param out : Level of output information written to stdout ( 0 - 10 )
  */
  Function_DualLinearTri(int out) :
  MOERTEL::Function(out,MOERTEL::Function::func_DualLinearTri)
  {}
  
  /*!
  \brief Copy-Constructor
  */
  Function_DualLinearTri(const MOERTEL::Function_DualLinearTri& old) :
  MOERTEL::Function(old.outputlevel_,old.type_)
  {}

  /*!
  \brief Clone method
  
  Makes a deep copy of this instance and returns a pointer to it.
  
  */
  virtual MOERTEL::Function* Clone() const;

  /*!
  \brief Destructor
  */
  virtual ~Function_DualLinearTri() {}

  /*!
  \brief Evaluate the function and derivatives at a given local coordinate
  
  \param xi (in) : Local coordinate where to evaluate the function (dimension 2)
  \param val (out) : Function values at xi, dimension 3, if NULL on input no evaluation
  \param valdim (in) : Dimension of val should be 3
  \param deriv (out) : Derivative of functions at xi, dimension 6, if NULL on input no evaluation
  */
  bool EvaluateFunction(const MOERTEL::Segment& seg, const double* xi, 
                        double* val, const int valdim, double* deriv);

};

/*!
\class Function_ConstantTri

\brief <b> A 2D constant shape function for a 3-noded triangle</b>

phi_1 = 1 <br>
phi_2 = 1<br>
phi_3 = 1 <br>
phi_1,xi_1 = 0 <br>
phi_1,xi_2 = 0 <br>
phi_2,xi_1 = 0 <br>
phi_2,xi_2 = 0 <br>
phi_3,xi_1 = 0 <br>
phi_3,xi_2 = 0 <br>

\author Glen Hansen (Glen.Hansen@inl.gov)

*/
class Function_ConstantTri : public MOERTEL::Function 
{
public:
  
  /*!
  \brief Constructor
  
  \param out : Level of output information written to stdout ( 0 - 10 )
  */
  Function_ConstantTri(int out) :
  MOERTEL::Function(out,MOERTEL::Function::func_ConstantTri)
  {}
  
  /*!
  \brief Copy-Constructor
  */
  Function_ConstantTri(const MOERTEL::Function_ConstantTri& old) :
  MOERTEL::Function(old.outputlevel_,old.type_)
  {}

  /*!
  \brief Clone method
  
  Makes a deep copy of this instance and returns a pointer to it.
  
  */
  virtual MOERTEL::Function* Clone() const;

  /*!
  \brief Destructor
  */
  virtual ~Function_ConstantTri() {}

  /*!
  \brief Evaluate the function and derivatives at a given local coordinate
  
  \param xi (in) : Local coordinate where to evaluate the function (dimension 2)
  \param val (out) : Function values at xi, dimension 3, if NULL on input no evaluation
  \param valdim (in) : Dimension of val should be 3
  \param deriv (out) : Derivative of functions at xi, dimension 6, if NULL on input no evaluation
  */
  bool EvaluateFunction(const MOERTEL::Segment& seg, const double* xi, 
                        double* val, const int valdim, double* deriv);

};


class Function_BiLinearQuad : public MOERTEL::Function 
{
public:
  
  /*!
  \brief Constructor
  
  \param out : Level of output information written to stdout ( 0 - 10 )
  */
  Function_BiLinearQuad(int out) :
  MOERTEL::Function(out,MOERTEL::Function::func_BiLinearQuad)
  {}
  
  /*!
  \brief Copy-Constructor
  */
  Function_BiLinearQuad(const MOERTEL::Function_BiLinearQuad& old) :
  MOERTEL::Function(old.outputlevel_,old.type_)
  {}

  /*!
  \brief Clone method
  
  Makes a deep copy of this instance and returns a pointer to it.
  
  */
  virtual MOERTEL::Function* Clone() const;

  /*!
  \brief Destructor
  */
  virtual ~Function_BiLinearQuad() {}

  /*!
  \brief Evaluate the function and derivatives at a given local coordinate
  
  \param xi (in) : Local coordinate where to evaluate the function (dimension 2)
  \param val (out) : Function values at xi, dimension 4, if NULL on input no evaluation
  \param valdim (in) : Dimension of val should be 4
  \param deriv (out) : Derivative of functions at xi, dimension 8, if NULL on input no evaluation
  */
  bool EvaluateFunction(const MOERTEL::Segment& seg, const double* xi, 
                        double* val, const int valdim, double* deriv);
};



class Function_DualBiLinearQuad : public MOERTEL::Function 
{
public:
  
  /*!
  \brief Constructor
  
  \param out : Level of output information written to stdout ( 0 - 10 )
  */
  Function_DualBiLinearQuad(int out) :
  MOERTEL::Function(out,MOERTEL::Function::func_DualBiLinearQuad)
  {}
  
  /*!
  \brief Copy-Constructor
  */
  Function_DualBiLinearQuad(const MOERTEL::Function_DualBiLinearQuad& old) :
  MOERTEL::Function(old.outputlevel_,old.type_)
  {}

  /*!
  \brief Clone method
  
  Makes a deep copy of this instance and returns a pointer to it.
  
  */
  virtual MOERTEL::Function* Clone() const;

  /*!
  \brief Destructor
  */
  virtual ~Function_DualBiLinearQuad() {}

  /*!
  \brief Evaluate the function and derivatives at a given local coordinate
  
  \param xi (in) : Local coordinate where to evaluate the function (dimension 2)
  \param val (out) : Function values at xi, dimension 4, if NULL on input no evaluation
  \param valdim (in) : Dimension of val should be 4
  \param deriv (out) : Derivative of functions at xi, dimension 8, if NULL on input no evaluation
  */
  bool EvaluateFunction(const MOERTEL::Segment& seg, const double* xi, 
                        double* val, const int valdim, double* deriv);
};




} // namespace MOERTEL


#endif // MOERTEL_FUNCTION_LINEAR1D_H
