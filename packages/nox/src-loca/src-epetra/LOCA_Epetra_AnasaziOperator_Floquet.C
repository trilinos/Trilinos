//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include "LOCA_Epetra_AnasaziOperator_Floquet.H"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::Epetra::AnasaziOperator::Floquet::Floquet(
	const Teuchos::RCP<LOCA::GlobalData>& global_data,
	const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	const Teuchos::RCP<Teuchos::ParameterList>& eigenParams_,
	const Teuchos::RCP<Teuchos::ParameterList>& solverParams_,
	const Teuchos::RCP<NOX::Abstract::Group>& grp_)
  : globalData(global_data),
    myLabel("Floquet Transformation"),
    eigenParams(eigenParams_),
    solverParams(solverParams_),
    grp(grp_),
    xyztInterface()
{
  std::string callingFunction =
    "LOCA::Epetra::AnasaziOperator::Floquet::Floquet()";

  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  Teuchos::RCP<NOX::Epetra::Group> NEGrp
      = Teuchos::rcp_dynamic_cast<NOX::Epetra::Group>(grp);
  if (NEGrp == Teuchos::null) std::cout << callingFunction << "  NEGrp cast failed " << std::endl;
  else std::cout << callingFunction << "  NEGrp cast succeeded." << std::endl;

  xyztInterface = Teuchos::rcp_dynamic_cast<LOCA::Epetra::Interface::xyzt>(NEGrp->getRequiredInterface());
  if (xyztInterface == Teuchos::null) std::cout << callingFunction << "  xyztInterface cast failed " << std::endl;
  else std::cout << callingFunction << "  xyztInterface cast succeeded." << std::endl;

  // make sure Jacobian is up-to-date
  
  xyztInterface->setFloquetFillFlag(true);  //Thhis setting is undone in destructor
  status = grp->computeJacobian();

  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                                                           callingFunction);
}

LOCA::Epetra::AnasaziOperator::Floquet::~Floquet()
{
  xyztInterface->setFloquetFillFlag(false);
}

const std::string&
LOCA::Epetra::AnasaziOperator::Floquet::label() const
{
  return myLabel;
}

void
LOCA::Epetra::AnasaziOperator::Floquet::apply(const NOX::Abstract::MultiVector& input, 
				     NOX::Abstract::MultiVector& output) const
{

  // Apply first part of monodromy operator on input vector
  Teuchos::RCP<NOX::Abstract::MultiVector> tmpVec =  input.clone();
  for (int i=0; i < input.numVectors(); i++) {
    NOX::Abstract::Vector& nAV = tmpVec->operator[](i);
    NOX::Epetra::Vector& nEV = dynamic_cast<NOX::Epetra::Vector&>(nAV);
    Epetra_Vector& eV = nEV.getEpetraVector();

    xyztInterface->beginFloquetOperatorApplication(eV);
  }
    
  // Now apply the main part of the monodromy matrix
  NOX::Abstract::Group::ReturnType status =
    grp->applyJacobianInverseMultiVector(*solverParams, *(tmpVec.get()), output);
  globalData->locaErrorCheck->checkReturnType(status,
                       "LOCA::Epetra::AnasaziOperator::Floquet::apply()");


  for (int i=0; i < input.numVectors(); i++) {
    NOX::Abstract::Vector& nAV = output.operator[](i);
    NOX::Epetra::Vector& nEV = dynamic_cast<NOX::Epetra::Vector&>(nAV);
    Epetra_Vector& eV = nEV.getEpetraVector();

    xyztInterface->finishFloquetOperatorApplication(eV);
  }

// Was this needed?

// TESTING:  Doubling the call to this routine resulted in the
//  squaring of the Floquet multipliers, as they should.
//  Replacing the apply function so the operator is diagonal
//  with entries 1/(i+2) led to the Floquet multipliers.

/*
  std::cout << " Fixing apply so Floquets at 1/2 1/3 1/4 ... " << std::endl;

  Teuchos::RCP<NOX::Abstract::MultiVector> tmpVec =  input.clone();
  for (int i=0; i < input.numVectors(); i++) {
    NOX::Abstract::Vector& nAV = output.operator[](i);
    NOX::Epetra::Vector& nEV = dynamic_cast<NOX::Epetra::Vector&>(nAV);
    Epetra_Vector& oV = nEV.getEpetraVector();

    NOX::Abstract::Vector& nAV2 = tmpVec->operator[](i);
    NOX::Epetra::Vector& nEV2 = dynamic_cast<NOX::Epetra::Vector&>(nAV2);
    Epetra_Vector& iV = nEV2.getEpetraVector();

    for (int j=0; j < iV.MyLength(); j++) {
     oV[j] = iV[j] / (j + 2.0);
    }
  }
*/
}

void
LOCA::Epetra::AnasaziOperator::Floquet::transformEigenvalue(double& ev_r, 
						   double& ev_i) const
{
  // Floquet multipliers need no transformations
}

NOX::Abstract::Group::ReturnType 
LOCA::Epetra::AnasaziOperator::Floquet::rayleighQuotient(
				         NOX::Abstract::Vector& evec_r,
					 NOX::Abstract::Vector& evec_i,
					 double& rq_r, double& rq_i) const
{
  std::string callingFunction = 
    "LOCA::Epetra::AnasaziOperator::Floquet::rayleighQuotient()";

  // create two-column  multivector of two eigenvectors
  Teuchos::RCP<NOX::Abstract::MultiVector> z = evec_r.createMultiVector(2, NOX::DeepCopy);
  z->operator[](1) = evec_i;
  Teuchos::RCP<NOX::Abstract::MultiVector> Az = z->clone(NOX::ShapeCopy);

  apply(*(z.get()), *(Az.get()));
  NOX::Abstract::Vector&  Ax = Az->operator[](0);
  NOX::Abstract::Vector&  Ay = Az->operator[](1);

  double mag = evec_r.innerProduct(evec_r) +  evec_i.innerProduct(evec_i) ;
  rq_r = (Ax.innerProduct(evec_r) +  Ay.innerProduct(evec_i)) / mag ;
  rq_i = (Ay.innerProduct(evec_r) -  Ax.innerProduct(evec_i)) / mag ;

  return NOX::Abstract::Group::Ok; 
}
