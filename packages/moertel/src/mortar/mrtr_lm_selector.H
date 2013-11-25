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
 * \file mrtr_lm_selector.H
 *
 * \class MOERTEL::Lmselector
 *
 * \brief A virtual class that supports user activation of LMs
 *
 * \date Last update do Doxygen: 6-Aug-10
 *
 */
#ifndef MOERTEL_LMSELECTOR_H
#define MOERTEL_LMSELECTOR_H

#include "Teuchos_RefCountPtr.hpp"
#include "mrtr_node.H"

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
\class Lmselector

\brief <b> A virtual class that supports user activation of LMs </b>



\author Glen Hansen (gahanse@sandia.gov)

*/
class Lmselector
{
public:

  // @{ \name Constructors and destructors

  /*!
  \brief Constructor
  
  Constructs an instance of this base class.<br>
  
  */
  Lmselector(){}
  
  /*!
  \brief Copy-Constructor
  
  Makes a deep copy<br>
  
  \param old : Instance to be copied
  */
  Lmselector(const MOERTEL::Lmselector& old){}
  
  /*!
  \brief Destructor
  */
  virtual ~Lmselector() {}

  //@}
  // @{ \name Public members

  /*!
  \brief Function that determines if the LMs should be active for this node and dof
  
  \param node (in) : the node hosting the LM
  \param LMDof (in) : the DOF of the LM in question

  \return bool : true if the M and D values of this node should be accumulated into the rhs, false if not.
  */

  virtual bool EvaluateLM(Teuchos::RCP<MOERTEL::Node> node, const int LMDof) = 0;

  /*!
  \brief Functions that accumulate values into the RHS depending on the state of the contact problem
  
  \param node (in) : the slave node of interest
  \param rhs (out) : the JFNK residual vector
  \param soln (in) : the current state vector

  Note that this class can (and should) be used to accumulate particular values of the LM into rhs using soln.
  This function is called for the slave interface nodes (that host the LMs).
  */

  virtual void AccumulateRHS(Teuchos::RCP<MOERTEL::Node> node) = 0;

  /*!
  \brief Functions that accumulate values into the RHS depending on the state of the contact problem
  
  \param row (in) : Row number of the D (or M) matrix to assemble
  \param col (in) : Column number of the D (or M) matrix to assemble
  \param val (in) : value of D (or M) to assemble

  This function is used to accumulate values of D and M into the JFNK rhs using soln. 
  */

  virtual void AssembleNodeVal(int row, int col, double val) = 0;


  //@}

protected:

  // Note that this base class does not hold any data.
  // If Your derived class needs to hold data, make sure it's carefully
  // taken care of in the copy-ctor!

};


} // namespace MOERTEL


#endif // MOERTEL_LMSELECTOR_H
