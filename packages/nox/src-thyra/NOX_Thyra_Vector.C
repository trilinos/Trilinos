// $Id$ 
// $Source$ 

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

#include "NOX_Thyra_Vector.H"
#include "NOX_Thyra_MultiVector.H"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"

NOX::Thyra::Vector::
Vector(const Teuchos::RCP< ::Thyra::VectorBase<double> >& source) :
  thyraVec(source)
{
}

NOX::Thyra::Vector::
Vector(const ::Thyra::VectorBase<double>& source) :
  thyraVec(source.clone_v())
{
}

NOX::Thyra::Vector::
Vector(const NOX::Thyra::Vector& source, NOX::CopyType type) :
  thyraVec(source.thyraVec->clone_v())
{
}

NOX::Thyra::Vector::
~Vector()
{
}

NOX::Abstract::Vector& 
NOX::Thyra::Vector::
operator=(const NOX::Abstract::Vector& src)
{
  const NOX::Thyra::Vector& source = 
    dynamic_cast<const NOX::Thyra::Vector&>(src);
  this->operator=(source);
  return *this;
}

NOX::Abstract::Vector& 
NOX::Thyra::Vector::
operator=(const NOX::Thyra::Vector& src)
{
  ::Thyra::copy(*src.thyraVec, thyraVec.get());
  return *this;
}

::Thyra::VectorBase<double>& 
NOX::Thyra::Vector::
getThyraVector()
{
  return *thyraVec;
}

const ::Thyra::VectorBase<double>& 
NOX::Thyra::Vector::
getThyraVector() const
{
  return *thyraVec;
}

Teuchos::RCP< ::Thyra::VectorBase<double> > 
NOX::Thyra::Vector::
getThyraRCPVector()
{
  return thyraVec;
}

Teuchos::RCP<const ::Thyra::VectorBase<double> > 
NOX::Thyra::Vector::
getThyraRCPVector() const
{
  return thyraVec;
}

NOX::Abstract::Vector& 
NOX::Thyra::Vector::
init(double value)
{
  ::Thyra::put_scalar(value, thyraVec.get());
  return *this;
}

NOX::Abstract::Vector& 
NOX::Thyra::Vector::
random(bool useSeed, int seed)
{
  if (useSeed)
    ::Thyra::seed_randomize<double>(seed);
  ::Thyra::randomize(-1.0, 1.0, thyraVec.get());
  return *this;
}

NOX::Abstract::Vector& 
NOX::Thyra::Vector::
abs(const NOX::Abstract::Vector& src)
{
  const NOX::Thyra::Vector& source = 
    dynamic_cast<const NOX::Thyra::Vector&>(src);
  ::Thyra::abs(thyraVec.get(), *source.thyraVec);
  return *this;
}

NOX::Abstract::Vector& 
NOX::Thyra::Vector::
reciprocal(const NOX::Abstract::Vector& src)
{
  const NOX::Thyra::Vector& source = 
    dynamic_cast<const NOX::Thyra::Vector&>(src);
  ::Thyra::reciprocal(thyraVec.get(), *source.thyraVec);
  return *this;
}

NOX::Abstract::Vector& 
NOX::Thyra::Vector::
scale(double alpha)
{
  ::Thyra::scale(alpha, thyraVec.get());
  return *this;
}

NOX::Abstract::Vector& 
NOX::Thyra::Vector::
scale(const NOX::Abstract::Vector& src)
{  
  const NOX::Thyra::Vector& source = 
    dynamic_cast<const NOX::Thyra::Vector&>(src);
  ::Thyra::Vt_StV(thyraVec.get(), 1.0, *source.thyraVec);
  return *this;
}

NOX::Abstract::Vector& 
NOX::Thyra::Vector::
update(double alpha, const NOX::Abstract::Vector& a, double gamma)
{
  const NOX::Thyra::Vector& aa = 
    dynamic_cast<const NOX::Thyra::Vector&>(a);
  const ::Thyra::VectorBase<double> *tv = aa.thyraVec.get();
  ::Thyra::linear_combination(1, &alpha, &tv, gamma, thyraVec.get());
  return *this;
}

NOX::Abstract::Vector& NOX::Thyra::Vector::
update(double alpha, const NOX::Abstract::Vector& x, 
       double beta, const NOX::Abstract::Vector& y,
       double gamma)
{
  const NOX::Thyra::Vector& xx = 
    dynamic_cast<const NOX::Thyra::Vector&>(x);
  const NOX::Thyra::Vector& yy = 
    dynamic_cast<const NOX::Thyra::Vector&>(y);
  double a[] = {alpha, beta};
  const ::Thyra::VectorBase<double>* z[] = {xx.thyraVec.get(), 
					    yy.thyraVec.get()};
  ::Thyra::linear_combination(2, a, z, gamma, thyraVec.get());
  return *this;
}

Teuchos::RCP<NOX::Abstract::Vector> 
NOX::Thyra::Vector::
clone(CopyType type) const
{
  return Teuchos::rcp(new NOX::Thyra::Vector(*this, type));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::Thyra::Vector::
createMultiVector(const NOX::Abstract::Vector* const* vecs,
		  int numVecs, NOX::CopyType type) const
{
  // Get vector space
  Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > space = 
    thyraVec->space();

  // Create Thyra multivector
  Teuchos::RCP< ::Thyra::MultiVectorBase<double> > mv = 
    ::Thyra::createMembers(*space, numVecs+1);

  // Copy vectors
  if (type == NOX::DeepCopy) {
    Teuchos::RCP< ::Thyra::VectorBase<double> > v = mv->col(0);
    ::Thyra::copy(*thyraVec, v.get());
    for (int i=0; i<numVecs; i++) {
      const NOX::Thyra::Vector* tv = 
	dynamic_cast<const NOX::Thyra::Vector*>(vecs[i]);
      v = mv->col(i+1);
      ::Thyra::copy(*(tv->thyraVec), v.get());
    }
  }
  
  // Create multi-vector
  Teuchos::RCP<NOX::Thyra::MultiVector> nmv = 
    Teuchos::rcp(new NOX::Thyra::MultiVector(mv));

  return nmv;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::Thyra::Vector::
createMultiVector(int numVecs, NOX::CopyType type) const
{
  // Get vector space
  Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > space = 
    thyraVec->space();

  // Create Thyra multivector
  Teuchos::RCP< ::Thyra::MultiVectorBase<double> > mv = 
    createMembers(*space, numVecs);

  // Copy vectors
  if (type == NOX::DeepCopy) {
    for (int i=0; i<numVecs; i++) {
      Teuchos::RCP< ::Thyra::VectorBase<double> > v = mv->col(i);
      ::Thyra::copy(*thyraVec, v.get());
    }
  }
  
  // Create multi-vector
  Teuchos::RCP<NOX::Thyra::MultiVector> nmv = 
    Teuchos::rcp(new NOX::Thyra::MultiVector(mv));

  return nmv;
}

double 
NOX::Thyra::Vector::
norm(NOX::Abstract::Vector::NormType type) const
{
  if (type == NOX::Abstract::Vector::TwoNorm)
    return ::Thyra::norm_2(*thyraVec);
  else if (type == NOX::Abstract::Vector::OneNorm)
    return ::Thyra::norm_1(*thyraVec);
  else
    return ::Thyra::norm_inf(*thyraVec);
}

double 
NOX::Thyra::Vector::
norm(const NOX::Abstract::Vector& weights) const
{
  const NOX::Thyra::Vector& w = 
    dynamic_cast<const NOX::Thyra::Vector&>(weights);
  return ::Thyra::norm_2(*w.thyraVec, *thyraVec);
}

double 
NOX::Thyra::Vector::
innerProduct(const NOX::Abstract::Vector& y) const
{
  const NOX::Thyra::Vector& yy = 
    dynamic_cast<const NOX::Thyra::Vector&>(y);
  return thyraVec->space()->scalarProd(*thyraVec, *yy.thyraVec);
}

int 
NOX::Thyra::Vector::
length() const
{
  return thyraVec->space()->dim();
}

void 
NOX::Thyra::Vector::
print(std::ostream& stream) const
{
  //stream << *thyraVec;
  thyraVec->describe(
    *Teuchos::getFancyOStream(Teuchos::rcpFromRef(stream)),
    Teuchos::VERB_HIGH);
  return;
}

