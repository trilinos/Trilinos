// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "NOX_Abstract_Vector.H"

// Included multivector declarations
#include "NOX_MultiVector.H"

void NOX::Abstract::Vector::print(std::ostream& /* stream */) const
{
  return;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::Abstract::Vector::createMultiVector(
                    const NOX::Abstract::Vector* const* vecs,
                    int numVecs, NOX::CopyType type) const
{
  if (numVecs < 0) {
    std::cerr << "NOX::Abstract::Vector::createMultiVector:  Error!  Multivector"
     << " must have positive number of columns!" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  const NOX::Abstract::Vector** tmp =
    new const NOX::Abstract::Vector*[numVecs+1];

  tmp[0] = this;
  for (int i=0; i<numVecs; i++)
    tmp[i+1] = vecs[i];

  Teuchos::RCP<NOX::MultiVector> mv =
    Teuchos::rcp(new NOX::MultiVector(tmp, numVecs+1, type));

  delete [] tmp;

  return mv;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::Abstract::Vector::createMultiVector(int numVecs, NOX::CopyType type) const
{
  if (numVecs <= 0) {
    std::cerr << "NOX::Abstract::Vector::createMultiVector:  Error!  Multivector"
     << " must have positive number of columns!" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  return Teuchos::rcp(new NOX::MultiVector(*this, numVecs, type));
}
