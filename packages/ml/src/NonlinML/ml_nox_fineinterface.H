/*
#@HEADER
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)

#ifndef ML_NOX_FINEINTERFACE_H
#define ML_NOX_FINEINTERFACE_H

#include <iostream>
// ML headers
#include "ml_common.h"
#include "ml_include.h"
// Epetra headers
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
// NOX interfaces
#include "NOX_EpetraNew_Interface_Required.H"
#include "NOX_EpetraNew_Interface_Jacobian.H"
#include "NOX_EpetraNew_Interface_Preconditioner.H"

namespace ML_NOX
{
class  Ml_Nox_Fineinterface : public NOX::EpetraNew::Interface::Required,
                              public NOX::EpetraNew::Interface::Jacobian,
                              public NOX::EpetraNew::Interface::Preconditioner
{
public:

  Ml_Nox_Fineinterface() 
  { isnewJacobian_=false; numJacobian_=0; printlevel_=10; return;}

  virtual ~Ml_Nox_Fineinterface() {};

  //! Compute and return F (derived from NOX::EpetraNew::Interface::Required)
  virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag) = 0;

  // Compute Jacobian (derived from NOX::EpetraNew::Interface::Jacobian)
  virtual bool computeJacobian(const Epetra_Vector& x) = 0;
  
  // Compute Preocnditioner (derived from NOX::EpetraNew::Interface::Preconditioner)
  virtual bool computePreconditioner(const Epetra_Vector& x,
				     NOX::Parameter::List* precParams = 0) = 0;

  // return the Jacobian from the problem
  virtual Epetra_CrsMatrix* getJacobian()  = 0;
/*
  {Epetra_CrsMatrix& Jac = problem.getJacobian(); 
   return(&Jac);}
*/  
  // return the graph from the problem
  virtual const Epetra_CrsGraph* getGraph() = 0;
/*
  {const Epetra_CrsGraph& Graph = problem.getGraph();
   return(&Graph);}
*/                                     
  // return the current solution from the problem
  virtual const Epetra_Vector* getSolution() = 0;
  //{const Epetra_Vector& vec = problem.getSolution(); return &vec; }                                     
                               
  // return the Map from the problem
  virtual const Epetra_Map& getMap() = 0;
  //{ return problem.getMap(); }                                     
                               
  // return the Jacobian flag (Jacobian valid or not)
  virtual bool isnewJacobian() = 0;
  //{ return  (isnewJacobian_); }     
  
  // return the number of Jacobians evaluated so far
  virtual int getnumJacobian() = 0;
  //{ return (numJacobian_); }                         

  // get nullspace
  virtual double* Get_Nullspace(const int nummyrows, const int numpde, 
                                const int dim_nullsp) = 0;
  
  // get block information
  virtual bool getBlockInfo(int *nblocks, int** blocks, int **block_pde) = 0;

private:  
  // don't want copy-ctor and = operator
  Ml_Nox_Fineinterface(ML_NOX::Ml_Nox_Fineinterface& old);
  Ml_Nox_Fineinterface operator = (const ML_NOX::Ml_Nox_Fineinterface& old);
  
  //! Application: Class of the the user's evaluation routines.
  /*! This is used to point to the actual routines and to store 
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is type of 
   *  passdown class design by the application code.
   *  The user should invlude a reference to this in his ML_NOX::Ml_Nox_Fineinterface
   *  derived class
   */ 
  //Ml_Nox_Nlnproblem& problem;
  
public:

  // flag to indicate, that the Jacobian has been evaluated
  // the L_NOX::Ml_Nox_Fineinterface should set this to true, whenever
  // the evaluated Jacobian matches the F in the computeF.
  // It should be set to false if there is either no Jacobian evaluated or
  // it does not match the F from the computeF function
  bool isnewJacobian_;
  
  // count number of evaluated Jacobians
  // The user's ML_NOX::Ml_Nox_Fineinterface derived class should use this to
  // count the number of Jacobian evaluations for statistical and debug-purpose
  int numJacobian_;
  
  // printlevel
  // The user's ML_NOX::Ml_Nox_Fineinterface derived class should set a
  // printlevel value between 0 and 10
  int printlevel_;
  
};
} // namespace ML_NOX


#endif // ML_NOX_FINEINTERFACE_H
#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
