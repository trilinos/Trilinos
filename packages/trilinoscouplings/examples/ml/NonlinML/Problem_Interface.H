// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
                                                                                
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
