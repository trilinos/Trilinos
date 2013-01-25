// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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

#include "NOX_Epetra_Vector.H"
#include "NOX_Epetra_MultiVector.H"
#include "Epetra_Vector.h"
#include "NOX_Epetra_VectorSpace_L2.H"

NOX::Epetra::Vector::
Vector(const Teuchos::RCP<Epetra_Vector>& source,
       NOX::Epetra::Vector::MemoryType memoryType, 
       NOX::CopyType type,
       Teuchos::RCP<NOX::Epetra::VectorSpace> vs)
{
  if (Teuchos::is_null(vs))
    vectorSpace = Teuchos::rcp(new NOX::Epetra::VectorSpaceL2);
  else
    vectorSpace = vs;

  if (memoryType == NOX::Epetra::Vector::CreateView)
    epetraVec = source;
  else {

    switch (type) {
      
    case DeepCopy:		// default behavior
      
      epetraVec = Teuchos::rcp(new Epetra_Vector(*source));
      break;
      
    case ShapeCopy:
      
      epetraVec = Teuchos::rcp(new Epetra_Vector(source->Map())); 
      break;  
    }

  }
}

NOX::Epetra::Vector::Vector(const Epetra_Vector& source, NOX::CopyType type,
			    Teuchos::RCP<NOX::Epetra::VectorSpace> vs)
{
  if (Teuchos::is_null(vs))
    vectorSpace = Teuchos::rcp(new NOX::Epetra::VectorSpaceL2);
  else
    vectorSpace = vs;

  switch (type) {

  case DeepCopy:		// default behavior

    epetraVec = Teuchos::rcp(new Epetra_Vector(source)); 
    break;

  case ShapeCopy:

    epetraVec = Teuchos::rcp(new Epetra_Vector(source.Map())); 
    break;  

  }
}

NOX::Epetra::Vector::Vector(const NOX::Epetra::Vector& source, 
			    NOX::CopyType type)
{
  vectorSpace = source.vectorSpace;

  switch (type) {

  case DeepCopy:		// default behavior

    epetraVec = Teuchos::rcp(new Epetra_Vector(source.getEpetraVector())); 
    break;

  case ShapeCopy:

    epetraVec = 
      Teuchos::rcp(new Epetra_Vector(source.getEpetraVector().Map())); 
    break;  

  }
}

NOX::Epetra::Vector::~Vector()
{

}

NOX::Abstract::Vector& NOX::Epetra::Vector::operator=(const Epetra_Vector& source)
{
  epetraVec->Scale(1.0, source);
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::operator=(const NOX::Abstract::Vector& source)
{
  return operator=(dynamic_cast<const NOX::Epetra::Vector&>(source));
}

NOX::Abstract::Vector& NOX::Epetra::Vector::operator=(const NOX::Epetra::Vector& source)
{
  epetraVec->Scale(1.0, source.getEpetraVector());
  return *this;
}

Epetra_Vector& NOX::Epetra::Vector::getEpetraVector()
{
  return *epetraVec;
}

const Epetra_Vector& NOX::Epetra::Vector::getEpetraVector() const
{
  return *epetraVec;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::init(double value)
{
  epetraVec->PutScalar(value);
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::random(bool useSeed, int seed)
{
  if (useSeed)
    epetraVec->SetSeed(seed);
  epetraVec->Random();
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::abs(const NOX::Abstract::Vector& base)
{
  return abs(dynamic_cast<const NOX::Epetra::Vector&>(base));
}

NOX::Abstract::Vector& NOX::Epetra::Vector::abs(const NOX::Epetra::Vector& base)
{
  epetraVec->Abs(base.getEpetraVector());
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::reciprocal(const NOX::Abstract::Vector& base)
{
  return reciprocal(dynamic_cast<const NOX::Epetra::Vector&>(base));
}

NOX::Abstract::Vector& NOX::Epetra::Vector::reciprocal(const NOX::Epetra::Vector& base)
{
  epetraVec->Reciprocal(base.getEpetraVector());
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::scale(double alpha)
{
  epetraVec->Scale(alpha);
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::
update(double alpha, const NOX::Abstract::Vector& a, double gamma)
{
  return update(alpha, dynamic_cast<const NOX::Epetra::Vector&>(a), gamma);
}

NOX::Abstract::Vector& NOX::Epetra::Vector::
update(double alpha, const NOX::Epetra::Vector& a, double gamma)
{
  epetraVec->Update(alpha, a.getEpetraVector(), gamma);
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::
update(double alpha, const NOX::Abstract::Vector& a, 
       double beta, const NOX::Abstract::Vector& b,
       double gamma)
{
  return update(alpha, dynamic_cast<const NOX::Epetra::Vector&>(a), 
		beta, dynamic_cast<const NOX::Epetra::Vector&>(b), gamma);
}

NOX::Abstract::Vector& NOX::Epetra::Vector::
update(double alpha, const NOX::Epetra::Vector& a, 
       double beta, const NOX::Epetra::Vector& b,
       double gamma)
{
  epetraVec->Update(alpha, a.getEpetraVector(), beta, b.getEpetraVector(), gamma);
  return *this;
}

NOX::Abstract::Vector& NOX::Epetra::Vector::scale(const NOX::Abstract::Vector& a)
{  
  return scale(dynamic_cast<const Epetra::Vector&>(a));
}

NOX::Abstract::Vector& NOX::Epetra::Vector::scale(const NOX::Epetra::Vector& a)
{  
  epetraVec->Multiply(1.0, *epetraVec, a.getEpetraVector(), 0.0);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Vector> NOX::Epetra::Vector::
clone(CopyType type) const
{
  Teuchos::RCP<NOX::Abstract::Vector> newVec = 
    Teuchos::rcp(new NOX::Epetra::Vector(*epetraVec, type, vectorSpace));
  return newVec;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::Epetra::Vector::createMultiVector(
				    const NOX::Abstract::Vector* const* vecs,
				    int numVecs, NOX::CopyType type) const
{
  if (numVecs < 0) {
    std::cerr << "NOX::Epetra::Vector::createMultiVector:  Error!  Multivector" 
	 << " must have postive number of columns!" << std::endl;
    throw "NOX Error";
  }

  double** v = new double*[numVecs+1];
  const Epetra_BlockMap& map = epetraVec->Map();
  const NOX::Epetra::Vector* noxEpetraVecPtr;
  const Epetra_Vector* vecPtr;
 
  epetraVec->ExtractView(&(v[0]));
  for (int i=0; i<numVecs; i++) {
    noxEpetraVecPtr = dynamic_cast<const NOX::Epetra::Vector*>(vecs[i]);
    vecPtr = &(noxEpetraVecPtr->getEpetraVector());
    vecPtr->ExtractView(&(v[i+1]));
  }

  Epetra_MultiVector epetra_mv(View, map, v, numVecs+1);

  Teuchos::RCP<NOX::Epetra::MultiVector> mv = 
    Teuchos::rcp(new NOX::Epetra::MultiVector(epetra_mv, type));

  delete [] v;

  return mv;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::Epetra::Vector::createMultiVector(int numVecs, NOX::CopyType type) const
{
  if (numVecs <= 0) {
    std::cerr << "NOX::Epetra::Vector::createMultiVector:  Error!  Multivector" 
	 << " must have postive number of columns!" << std::endl;
    throw "NOX Error";
  }

  const Epetra_BlockMap& map = epetraVec->Map();
  Epetra_MultiVector *epetra_mv;
  
  if (type == NOX::ShapeCopy)
    epetra_mv = new Epetra_MultiVector(map, numVecs, true);
  else {
    epetra_mv = new Epetra_MultiVector(map, numVecs, false);
    Epetra_Vector* v;
    for (int i=0; i<numVecs; i++) {
      v = (*epetra_mv)(i);
      *v = *epetraVec;
    }
  }

  Teuchos::RCP<NOX::Epetra::MultiVector> mv = 
    Teuchos::rcp(new NOX::Epetra::MultiVector(*epetra_mv, type));

  delete epetra_mv;

  return mv;
}

double NOX::Epetra::Vector::norm(NOX::Abstract::Vector::NormType type) const
{
  return vectorSpace->norm(*epetraVec, type);
}

double NOX::Epetra::Vector::norm(const NOX::Abstract::Vector& weights) const
{
  return norm(dynamic_cast<const NOX::Epetra::Vector&>(weights));
}

double NOX::Epetra::Vector::norm(const NOX::Epetra::Vector& weights) const
{
    std::cerr << "NOX::Epetra::Vector - Weighted norm not supported" << std::endl;
    throw "NOX-Epetra Error";
}

double NOX::Epetra::Vector::innerProduct(const NOX::Abstract::Vector& y) const
{
  return innerProduct(dynamic_cast<const NOX::Epetra::Vector&>(y));
}

double NOX::Epetra::Vector::innerProduct(const NOX::Epetra::Vector& y) const
{
  return vectorSpace->innerProduct(*epetraVec, y.getEpetraVector());
}

int NOX::Epetra::Vector::length() const
{
  return epetraVec->GlobalLength();
}

void NOX::Epetra::Vector::print(std::ostream& stream) const
{
  epetraVec->Print(stream);
  return;
}

Teuchos::RCP<NOX::Epetra::VectorSpace> 
NOX::Epetra::Vector::getVectorSpace() const
{
  return vectorSpace;
}
