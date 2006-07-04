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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#ifndef LOCA_EPETRA_XYZTPREC_H_
#define LOCA_EPETRA_XYZTPREC_H_

#include "NOX_Common.H"
#include "Teuchos_ParameterList.hpp"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_Epetra_Vector.H"
#include "EpetraExt_MultiMpiComm.h"
#include "EpetraExt_BlockVector.h"
#include "EpetraExt_BlockCrsMatrix.h"
#include "Epetra_Operator.h"
#include "Epetra_Import.h"
#include "NOX_Epetra_LinearSystem_AztecOO.H"

// Forward Declarations
namespace LOCA {
  namespace Epetra {

  /*! 
    \brief Abstract base class for preconditioners to be used with LOCA::Epetra::xyzt.
  */
class xyztPrec : public Epetra_Operator, 
		 public NOX::Epetra::Interface::Required,
		 public NOX::Epetra::Interface::Jacobian,
		 public NOX::Epetra::Interface::Preconditioner
{

 public:

  //! Constructor
  xyztPrec(EpetraExt::BlockCrsMatrix& jacobian, 
	   Epetra_CrsMatrix &splitJac,
	   Epetra_CrsMatrix &splitMass,
	   EpetraExt::BlockVector& solution,
	   EpetraExt::BlockVector& solutionOverlap,
	   Epetra_Import &overlapImporter,
	   Teuchos::ParameterList& precPrintParams, 
	   Teuchos::ParameterList& precLSParams, 
	   const Teuchos::RefCountPtr<EpetraExt::MultiMpiComm> globalComm_);

  //! Destructor
  virtual ~xyztPrec();

  // Inherited from Epetra_Operator
  virtual int SetUseTranspose(bool UseTranspose);

  virtual int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  virtual int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  virtual double NormInf() const;

  virtual const char* Label() const;

  virtual bool UseTranspose() const;

  virtual bool HasNormInf() const;

  virtual const Epetra_Comm& Comm() const;

  virtual const Epetra_Map& OperatorDomainMap() const;

  virtual const Epetra_Map& OperatorRangeMap() const;

  // Inherited from NOX::Epetra::Interface::Required
  virtual bool computeF(const Epetra_Vector&, 
			Epetra_Vector&, 
			NOX::Epetra::Interface::Required::FillType);

  // Inherited from NOX::Epetra::Interface::Jacobian
  virtual bool computeJacobian(const Epetra_Vector&, Epetra_Operator&);

  // Inherited from NOX::Epetra::Interface::Preconditioner
  virtual bool computePreconditioner(const Epetra_Vector& x,
				     Epetra_Operator& Prec,
				     Teuchos::ParameterList* p = 0);


  // New in LOCA::Epetra::xyztPrec
  virtual void throwError(const string& functionName, 
			  const string& errorMsg) const;


private:

  // Info passed in to construtor
  EpetraExt::BlockCrsMatrix& jacobian;
  Epetra_CrsMatrix &splitJac;
  Epetra_CrsMatrix &splitMass;
  EpetraExt::BlockVector& solution;
  EpetraExt::BlockVector& solutionOverlap;
  Epetra_Import& overlapImporter;
  Teuchos::ParameterList& printParams;
  Teuchos::ParameterList& lsParams;
  const Teuchos::RefCountPtr<EpetraExt::MultiMpiComm> globalComm;

  // Objects created in constructor
  //Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linSys;
  std::vector<NOX::Epetra::LinearSystemAztecOO *> linSys;
  std::vector<Teuchos::RefCountPtr<Epetra_CrsMatrix> > jacobianBlock;
  std::vector<Teuchos::RefCountPtr<Epetra_CrsMatrix> > massBlock;
  std::vector<Teuchos::RefCountPtr<Epetra_Vector> > diagBlockSubdiag;
  EpetraExt::BlockVector *residual;
  Epetra_Vector *splitVec; //tmp space
  Epetra_Vector *splitRes; //tmp space
  Epetra_Vector *splitVecOld; //tmp space
  NOX::Epetra::Vector *splitRes_NEV; //view of splitRes
  NOX::Epetra::Vector *splitVec_NEV; //view of splitVec
  bool  isPeriodic;

  //! Label for this operator
  string label;
  
};

} // namespace Epetra
} // namespace LOCA

#endif /* LOCA_EPETRA_XYZTPREC_H_ */
