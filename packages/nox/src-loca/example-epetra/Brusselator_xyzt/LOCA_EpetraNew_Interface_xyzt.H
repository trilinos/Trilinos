//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_EPETRANEW_INTERFACE_XYZT_H
#define LOCA_EPETRANEW_INTERFACE_XYZT_H

#include "NOX_Common.H"
#include "LOCA_EpetraNew_Interface_Required.H"
#include "NOX_EpetraNew_Interface_Jacobian.H"
#include "LOCA_EpetraNew_Interface_MassMatrix.H"
  
#ifdef HAVE_MPI
#ifdef HAVE_NOX_EPETRAEXT

#include <iostream>
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_DistObject.h"
#include "Epetra_Import.h"
#include "Epetra_MpiComm.h"

#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockCrsMatrix.h"

namespace LOCA {
namespace EpetraNew {

  /*! 
    \brief Provides a set of interfaces for users to provide information about the nonlinear problem to NOX.  

    Contains interfaces for the user to supply (1) the evaluation of the nonlinear equations, (2) the Jacobian, and (3) any preconditioning if required.
  */
namespace Interface {

  /*!
    \brief Interface for space-time (xyzt) problems.
 
    Code that takes a standard NOX/LOCA problem interface for spatially-discretized
    problems, and creates a NOX/LOCA interface for space-time problems.  This
    interface can be called by multiple replicas of the spatially discretized
    problem for parallelism of the time domain. Each replica can also own on 
    one or more time steps. This interface assembles a single space-time system
    including the solution and residual BlockVectors and the Jacobian BlockCrsMatrix.
    This interface relies heavily on the EpetraExt_Block* classes.
  */
class xyzt : public LOCA::EpetraNew::Interface::Required,
             public NOX::EpetraNew::Interface::Jacobian
{

public:
 
  //! Constructor
  xyzt(LOCA::EpetraNew::Interface::Required &iReq, NOX::EpetraNew::Interface::Jacobian &iJac,
       LOCA::EpetraNew::Interface::MassMatrix &iMass, Epetra_Vector &splitVec,
       Epetra_RowMatrix &splitJac, Epetra_RowMatrix &splitMass, 
       Epetra_MpiComm &globalComm, int replica, int timeStepsPerProc);

  //! Destructor
  virtual ~xyzt();

  //! NOX Interface function for computing the function, F, given the vector x.
  //Returns true if computation was successful.
  virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F, 
			const FillType fillFlag);
  
  //! NOX Interface function for computing the Jacobian matrix
  //Returns true if computation was successful.
  virtual bool computeJacobian(const Epetra_Vector& x);

  //! LOCA Interface function for setting the continuation parameter
  virtual void setParameters(const LOCA::ParameterVector& param);

  //! LOCA Interface function for printing the solution
  virtual void printSolution(const Epetra_Vector& x_, double conParam);

  //! Accessor for the global solution BlockVector
  virtual EpetraExt::BlockVector& getSolution();
  //! Accessor for the global BlockCrsMatrix
  virtual EpetraExt::BlockCrsMatrix& getJacobian();

private:
  
  // Info from split (3D) problem, passed in to the constructor.
  LOCA::EpetraNew::Interface::Required &iReq;
  NOX::EpetraNew::Interface::Jacobian &iJac;
  LOCA::EpetraNew::Interface::MassMatrix &iMass;
  Epetra_Vector& splitVec;
  Epetra_RowMatrix& splitJac;
  Epetra_RowMatrix& splitMass; // can be same object as splitJac
  Epetra_MpiComm& globalComm;
  int replica;
  int timeStepsPerProc;

  // xyzt objects created in constructor
  Epetra_Vector splitRes; //tmp space
  Epetra_Vector splitVecOld; //tmp space
  EpetraExt::BlockCrsMatrix* jacobian;
  EpetraExt::BlockVector* solution;
  EpetraExt::BlockVector* solutionOverlap;
  Epetra_Import* overlapImporter;
  int numReplicas;
  std::vector< std::vector<int> >* rowStencil;
  std::vector<int>* rowIndex;

};
} // namespace Interface
} // namespace EpetraNew
} // namespace LOCA

#endif
#endif
#endif
