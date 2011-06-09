/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2010) Sandia Corporation
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



\author Glen Hansen (Glen.Hansen@inl.gov)

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
