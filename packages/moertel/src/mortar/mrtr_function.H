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
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file mrtr_function.H
 *
 * \class MOERTEL::Function
 *
 * \brief A virtual class to form a shape function of some type
 *
 * \date Last update do Doxygen: 20-March-06
 *
 */
#ifndef MOERTEL_FUNCTION_H
#define MOERTEL_FUNCTION_H

#include <ctime>
#include <iostream>

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

class Segment;

/*!
\class Function

\brief <b> A virtual class to form a shape function of some type </b>



\author Glen Hansen (gahanse@sandia.gov)

*/
class Function 
{
public:

  /*!
  \brief Type of function
  
  \param func_none : default value
  \param func_Constant1D : Constant function on a 1D Segment
  \param func_Linear1D : Linear function on a 1D Segment
  \param func_DualLinear1D : Dual linear function on a 1D Segment
  \param func_LinearTri : Linear function on a 3-noded triangle
  \param func_DualLinearTri : Dual linear function on a 3-noded triangle
  \param func_ConstantTri : Constant function on a 3-noded triangle
  */
  enum FunctionType
  {
    func_none,
    func_Constant1D,
    func_Linear1D,
    func_DualLinear1D,
    func_LinearTri,
    func_DualLinearTri,
    func_ConstantTri,
    func_BiLinearQuad,
    func_DualBiLinearQuad
  };
  
  // @{ \name Constructors and destructors

  /*!
  \brief Constructor
  
  Constructs an instance of this base class.<br>
  
  \param outlevel : Level of output information written to stdout ( 0 - 10 )
  \param type : Function type
  */
  Function(int outlevel,MOERTEL::Function::FunctionType type) : 
  type_(type),
  outputlevel_(outlevel) 
  {};
  
  /*!
  \brief Copy-Constructor
  
  Makes a deep copy<br>
  
  \param old : Instance to be copied
  */
  Function(const MOERTEL::Function& old) :
  type_(old.type_),
  outputlevel_(old.outputlevel_)
  {}
  
  /*!
  \brief Destructor
  */
  virtual ~Function() {}

  /*!
  \brief Clone method
  
  The derived class implementing this virtual method is supposed to make a
  deep copy of itself and return a pointer to the copy.
  
  */
  virtual MOERTEL::Function* Clone() const = 0;


  //@}
  // @{ \name Public members

  /*!
  \brief Return the level of output written to stdout ( 0 - 10 )
  
  */
  int OutLevel() { return outputlevel_; }
  
  /*!
  \brief Return the function type
  
  */
  MOERTEL::Function::FunctionType Type() const { return type_;}
  
  /*!
  \brief Evaluate the function and derivatives at a given local coordinate
  
  \param xi (in) : Local coordinate where to evaluate the function
  \param val (out) : Function values at xi
  \param valdim (in) : Dimension of val
  \param deriv (out) : Derivative of functions at xi
  */
  virtual bool EvaluateFunction(const MOERTEL::Segment& seg, const double* xi, 
                                double* val, const int valdim,double* deriv) = 0;

  //@}

protected:

  MOERTEL::Function::FunctionType type_;
  int                             outputlevel_;

  // Note that this base class does not hold any data.
  // If Your derived class needs to hold data, make sure it's carefully
  // taken care of in the copy-ctor!

};


} // namespace MOERTEL


#endif // MOERTEL_FUNCTION_H
