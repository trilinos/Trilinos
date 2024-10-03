// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Transfer Operator for use with 1D linear elements

#ifndef _NOX_MULTIPHYSICS_XFER_OPERATOR
#define _NOX_MULTIPHYSICS_XFER_OPERATOR

//#include "GenericEpetraProblem.H"
class GenericEpetraProblem;
class Epetra_Vector;

#include <map>

class XferOp {

 public:

  // Constructor
  XferOp(GenericEpetraProblem& probA, const GenericEpetraProblem& probB);

  // Destructor
  ~XferOp();

  // Transfers data from one problem to another using a previously constructed
  // mapping.
  void transferField(Epetra_Vector& vecTo, Epetra_Vector& vecFrom);

  // Return a const reference to the dependent nodes map
  multimap<int,int>& getDependentNodesMap()
    { return dependentNodes; };

private:

  std::multimap<int,int> dependentNodes;
  std::multimap<int,double> dependentWeights;

};
#endif
