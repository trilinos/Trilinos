//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//                                                                                 
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA                                                                                
// Questions? Contact Michael Gee (mwgee@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER
                                                                                
//-----------------------------------------------------------------------------
#ifndef Problem_Interface_H
#define Problem_Interface_H


// ML headers
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)
#include "ml_include.h"
#include "nlnml_finelevelnoxinterface.H"
// ---------- Standard Includes ----------
#include <iostream>
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include "FiniteElementProblem.H"
// ---------- Forward Declarations ----------
class FiniteElementProblem;

class  Problem_Interface : public NLNML::NLNML_FineLevelNoxInterface
{
public:
  explicit Problem_Interface(FiniteElementProblem& Problem);
  
  virtual ~Problem_Interface();

  //! Compute and return F
  bool computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag);

  //! Compute an explicit Jacobian
  bool computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac);

  // return the Jacobian from the problem
  Epetra_CrsMatrix* getJacobian() 
  {
    Epetra_CrsMatrix& Jac = problem.getJacobian(); 
    return(&Jac);
  }
  
  // return the graph from the problem
  const Epetra_CrsGraph* getGraph() {const Epetra_CrsGraph& Graph = problem.getGraph();
                                     return(&Graph);}
  const Epetra_CrsGraph* getModifiedGraph() {const Epetra_CrsGraph& Graph = problem.getGraph();
                                     return(&Graph);}
                                     
  // return the current solution from the problem
  const Epetra_Vector* getSolution() {const Epetra_Vector& vec = problem.getSolution(); return &vec; }                                     
                               
  // return the Map from the problem
  const Epetra_Map& getMap() { return problem.getMap(); }                                     

  // get nullspace
  double* Get_Nullspace(const int nummyrows, const int numpde, const int dim_nullsp)
  { 
    return NULL; 
  }
  
  // get block information
  bool getBlockInfo(int *nblocks, std::vector<int>& blocks, std::vector<int>& block_pde)
  { *nblocks = 0; blocks.clear(); block_pde.clear(); return false; }

  // derived from NOX::Parameter::PrePostOperator
  const std::string& getType() const { return type_; }

  //! Application Operator: Object that points to the user's evaluation routines.
  /*! This is used to point to the actual routines and to store 
   *  auxiliary data required by the user's application for function/Jacobian
   *  evaluations that NOX does not need to know about.  This is type of 
   *  passdown class design by the application code.
   */ 

  FiniteElementProblem& problem;
  std::string type_;
};

#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO)
#endif
