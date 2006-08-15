//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_EPETRA_INTERFACE_XYZT_H
#define LOCA_EPETRA_INTERFACE_XYZT_H

#include "LOCA.H"
#include "NOX_Common.H"
#include "LOCA_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "LOCA_Epetra_Interface_MassMatrix.H"
#include "LOCA_Epetra_xyztPrec.H"

#ifdef HAVE_NOX_EPETRAEXT

#include <iostream>
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_DistObject.h"
#include "Epetra_Import.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "EpetraExt_MultiMpiComm.h"
#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockCrsMatrix.h"

namespace LOCA {
namespace Epetra {

namespace Interface {

  /*!
    \brief Interface for space-time (xyzt) problems.
 
    Code that takes a standard NOX/LOCA problem interfaces for spatially-discretized
    problems, and creates a NOX/LOCA interface for space-time problems.  This
    interface is called by multiple replicas of the spatially discretized
    problem for parallelism of the time domain. Each replica can also own on 
    one or more time steps. This interface assembles a single space-time system
    including the solution and residual BlockVectors and the Jacobian BlockCrsMatrix.
    This interface relies heavily on the EpetraExt_Block* classes. For usage,
    see the example problem:  nox/src-loca/example-epetra/Brusselator_xyzt
  */
class xyzt : public LOCA::Epetra::Interface::Required,
             public NOX::Epetra::Interface::Jacobian
{

public:
 
  //! Constructor
  /*! The interface, vector, and matrix arguments are all for the spatially
   *  decomposed problem. The MultiMpiCOmm communicator gives this
   *  class the information needed to construct the space-time system. 
   */

  xyzt(const Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required> &iReq,
       const Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> &iJac,
       const Teuchos::RefCountPtr<LOCA::Epetra::Interface::MassMatrix> &iMass,
       const Epetra_MultiVector &splitMultiVec,
       const Teuchos::RefCountPtr<Epetra_RowMatrix> &splitJac,
       const Teuchos::RefCountPtr<Epetra_RowMatrix> &splitMass,
       const Teuchos::RefCountPtr<EpetraExt::MultiMpiComm> &globalComm,
       Teuchos::ParameterList *precPrintParams = 0, 
       Teuchos::ParameterList *precLSParams = 0);

  //! Destructor
  virtual ~xyzt();

  //! NOX Interface function for computing the function, F, given the vector x.
  //Returns true if computation was successful.
  virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F, 
                        const FillType fillFlag);
  
  //! NOX Interface function for computing the Jacobian matrix
  //Returns true if computation was successful.
  virtual bool computeJacobian(const Epetra_Vector& x,
			       Epetra_Operator& Jac);

  //! LOCA Interface function for setting the continuation parameter
  virtual void setParameters(const LOCA::ParameterVector& param);

  //! LOCA Interface function for printing the solution
  virtual void printSolution(const Epetra_Vector& x_, double conParam);

  virtual void throwError(const string& functionName, 
			  const string& errorMsg) const;

  //! Accessor for the global solution BlockVector
  virtual EpetraExt::BlockVector& getSolution();
  //! Accessor for the global BlockCrsMatrix
  virtual EpetraExt::BlockCrsMatrix& getJacobian();
  //! Accessor for the xyzt Preconditioner
  virtual LOCA::Epetra::xyztPrec& getPreconditioner();

private:
  
  // Info from split (spatial) problem, passed in to the constructor.
  Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required> iReq;
  Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac;
  Teuchos::RefCountPtr<LOCA::Epetra::Interface::MassMatrix> iMass;
  Teuchos::RefCountPtr<Epetra_RowMatrix> splitJac;
  Teuchos::RefCountPtr<Epetra_RowMatrix> splitMass; // can be same object as splitJac
  Teuchos::RefCountPtr<EpetraExt::MultiMpiComm> globalComm;

  // xyzt objects created in constructor
  Epetra_Vector splitVec; //tmp space
  Epetra_Vector splitRes; //tmp space
  Epetra_Vector splitVecOld; //tmp space
  EpetraExt::BlockCrsMatrix* jacobian;
  EpetraExt::BlockVector* solution;
  EpetraExt::BlockVector* solutionOverlap;
  Epetra_Import* overlapImporter;
  int timeStepsOnTimeDomain; 
  int numTimeDomains;
  int timeDomain;
  int conStep;
  std::vector< std::vector<int> >* rowStencil;
  std::vector<int>* rowIndex;
  // xyztPrec object created in constructor
  Teuchos::ParameterList* precPrintParams;
  Teuchos::ParameterList* precLSParams;
  LOCA::Epetra::xyztPrec* preconditioner;
  Epetra_CrsMatrix* splitJacCrs;
  Epetra_CrsMatrix* splitMassCrs;
  bool  isPeriodic;

};
} // namespace Interface
} // namespace Epetra
} // namespace LOCA

#endif
#endif
