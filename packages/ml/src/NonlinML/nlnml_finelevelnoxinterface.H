/*
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
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file nlnml_finelevelnoxinterface.H
 *
 * \class NLNML_finelevelnoxinterface
 *
 * \brief virtual class NLNML_finelevelnoxinterface used by the nonlinear preconditioner
 *
 * \date Last update do Doxygen: 31-Mar-05
 *
 */
// ML headers
#include "ml_common.h"
#include "ml_include.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)
#ifndef NLNML_FINELEVELNOXINTERFACE_H
#define NLNML_FINELEVELNOXINTERFACE_H

#include <iostream>
// Epetra headers
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
// NOX interface
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Required.H"

namespace NLNML
{

class  NLNML_FineLevelNoxInterface : public virtual NOX::Epetra::Interface::Required,
                                     public virtual NOX::Epetra::Interface::Jacobian
{
public:

  NLNML_FineLevelNoxInterface()
  { isnewJacobian_=false; numJacobian_=0; printlevel_=0; t_ = 0.; ncalls_computeF_=0; return;}

  virtual ~NLNML_FineLevelNoxInterface() {}
  
  virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F, 
			const FillType fillFlag) = 0;
  
  virtual bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac) = 0;
  
  virtual Epetra_CrsMatrix* getJacobian()  = 0;
  
  virtual const Epetra_CrsGraph* getGraph() = 0;
  virtual const Epetra_CrsGraph* getModifiedGraph() = 0;
  
  virtual const Epetra_Vector* getSolution() = 0;
  
  virtual const Epetra_Map& getMap() = 0;

  virtual double* Get_Nullspace(const int nummyrows, const int numpde, 
                                const int dim_nullsp) = 0;

  virtual bool getBlockInfo(int *nblocks, vector<int>& blocks, vector<int>& block_pde) = 0;

  virtual void ApplyAllConstraints(Epetra_Vector& gradient, int level) 
  { return;}

  virtual bool isnewJacobian() { return isnewJacobian_; }

  virtual int getnumJacobian() { return numJacobian_; }

  double getsumtime() {return t_;}

  void resetsumtime() { t_ = 0.; return; }

  int getnumcallscomputeF() { return ncalls_computeF_; }

  bool setnumcallscomputeF(int ncalls) { ncalls_computeF_=ncalls; return true; }

  inline int OutLevel() { return printlevel_; }

private:  
  
  NLNML_FineLevelNoxInterface(NLNML_FineLevelNoxInterface& old);

protected:  
  
  bool isnewJacobian_;
  int numJacobian_;
  int printlevel_;
  double t_;
  int ncalls_computeF_;
  
}; // class NLNML_FineLevelNoxInterface
} // namespace NLNML

#endif
#endif
