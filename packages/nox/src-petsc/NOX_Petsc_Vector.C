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

#include "petscvec.h" // Petsc Vec header 

#include "NOX_Common.H"
#include "NOX_Petsc_Vector.H"

using namespace NOX;
using namespace NOX::Petsc;

Vector::Vector(const Vec& source, CopyType type)
{
  allocate(source, type);
}

Vector::Vector(const Vector& source, CopyType type)
{
  allocate(source.getPetscVector(), type);
}

Vector::Vector(const Vec& source, std::string Name, CopyType type) :
  name(Name)
{
  allocate(source, type);
}

Vector::~Vector()
{
  if(isAlloc)
  {
    //cout << "\n\t\tVecDestroy called ....  " << name << "  " << this << "   " 
    //     << &petscVec << std::endl;
    //VecDestroy(petscVec);
    isAlloc = false;
  }
}

Abstract::Vector& 
Vector::operator=(const Vec& source)
{
  VecCopy(source, petscVec);
  name = "Unnamed";
  isAlloc = true;
  return *this;
}

Abstract::Vector& 
Vector::operator=(const Abstract::Vector& source)
{
  return operator=(dynamic_cast<const Vector&>(source));
}

Abstract::Vector& 
Vector::operator=(const Vector& source)
{
  VecCopy(source.getPetscVector(), petscVec);
  isAlloc = source.isAlloc;
  name = source.name;
  return *this;
}

int 
Vector::allocate(const Vec& source, CopyType type)
{
 
  int ierr = VecDuplicate(source, &petscVec);
  isAlloc = true;

  switch (type) {

    case DeepCopy:                // default behavior

    ierr += VecCopy(source, petscVec);
    break;

  case ShapeCopy:

    break;
  }

  if(ierr)
    std::cout << "ERROR: value " << ierr << " returned during "
         << "NOX::Petsc::Vector allocation !!" << std::endl;

  return ierr;
}

Vec& 
Vector::getPetscVector()
{
  return petscVec;
}

const Vec& 
Vector::getPetscVector() const
{
  return petscVec;
}

Abstract::Vector& 
Vector::init(double value)
{
  VecSet(petscVec, value);
  return *this;
}

Abstract::Vector& 
Vector::abs(const Abstract::Vector& base)
{
  return abs(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& 
Vector::abs(const Vector& base)
{
  VecCopy(base.getPetscVector(), petscVec);
  VecAbs(petscVec);
  return *this;
}

Abstract::Vector& 
Vector::reciprocal(const Abstract::Vector& base)
{
  return reciprocal(dynamic_cast<const Vector&>(base));
}

Abstract::Vector& 
Vector::reciprocal(const Vector& base)
{
  VecCopy(base.getPetscVector(), petscVec);
  VecReciprocal(petscVec);
  return *this;
}

Abstract::Vector& 
Vector::scale(double alpha)
{
  VecScale(petscVec, alpha);
  return *this;
}

Abstract::Vector& 
Vector::scale(const Abstract::Vector& a)
{
  return scale(dynamic_cast<const Petsc::Vector&>(a));
}
  
Abstract::Vector& 
Vector::scale(const Vector& a)
{
  VecPointwiseMult(petscVec, a.getPetscVector(), petscVec);
  return *this;
} 
  
Abstract::Vector& 
Vector::update(double alpha, const Abstract::Vector& a, 
				 double gammaval)
{
  return update(alpha, dynamic_cast<const Vector&>(a), gammaval);
}

Abstract::Vector& 
Vector::update(double alpha, const Vector& a, 
				 double gammaval)
{
  VecAXPBY(petscVec, alpha, gammaval, a.getPetscVector());
  return *this;
}

Abstract::Vector& 
Vector::update(double alpha, const Abstract::Vector& a, 
				 double beta, const Abstract::Vector& b,
				 double gammaval)
{
  return update(alpha, dynamic_cast<const Vector&>(a), 
		beta, dynamic_cast<const Vector&>(b), gammaval);
}

Abstract::Vector& 
Vector::update(double alpha, const Vector& a, 
				 double beta, const Vector& b,
				 double gammaval)
{
  VecAXPBY(petscVec, alpha, gammaval, a.getPetscVector());
  VecAXPY(petscVec, beta, b.getPetscVector());
  return *this;
}


Teuchos::RCP<NOX::Abstract::Vector> 
Vector::clone(CopyType type) const
{
  Teuchos::RCP<NOX::Abstract::Vector> newVec = 
    Teuchos::rcp(new NOX::Petsc::Vector(petscVec, type));
  return newVec;
}

double 
Vector::norm(Abstract::Vector::NormType type) const
{
  double n;
  switch (type) {
  case MaxNorm:
    VecNorm(petscVec, NORM_INFINITY, &n);
    break;
  case OneNorm:
    VecNorm(petscVec, NORM_1, &n);
    break;
  case TwoNorm:
  default:
   VecNorm(petscVec, NORM_2, &n);
   break;
  }
  return n;
}

double 
Vector::norm(const Abstract::Vector& weights) const
{
  return norm(dynamic_cast<const Vector&>(weights));
}

double 
Vector::norm(const Vector& weights) const
{
  double n = 0.0;
  std::cerr << "Norm type not supported for weighted norm" << std::endl;
  throw;
  return n;
}

double 
Vector::innerProduct(const Abstract::Vector& y) const
{
  return innerProduct(dynamic_cast<const Vector&>(y));
}

double 
Vector::innerProduct(const Vector& y) const
{
  double dotprod;
  VecDot(y.getPetscVector(), petscVec, &dotprod);
  return dotprod;
}

int 
Vector::length() const
{
  int size;
  VecGetSize(petscVec, &size);
  return size;
}
