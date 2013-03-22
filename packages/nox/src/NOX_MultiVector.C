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

#include "NOX_MultiVector.H"

NOX::MultiVector::MultiVector(int numVecs) :
  vecs(numVecs)
{
  if (numVecs <= 0) {
    std::cerr << "NOX::MultiVector:  Error!  Multivector" 
	 << " must have postive number of columns!" << std::endl;
    throw "NOX Error";
  }
}

NOX::MultiVector::MultiVector(const NOX::Abstract::Vector& v, int numVecs,
			      NOX::CopyType type) :
  vecs(numVecs)
{
  if (numVecs <= 0) {
    std::cerr << "NOX::MultiVector:  Error!  Multivector" 
	 << " must have postive number of columns!" << std::endl;
    throw "NOX Error";
  }

  for (int i=0; i<numVecs; i++) {
    vecs[i] = v.clone(type);
  }
}

NOX::MultiVector::MultiVector(const NOX::Abstract::Vector* const* vs,
			      int numVecs,
			      NOX::CopyType type) :
  vecs(numVecs)
{
  if (numVecs <= 0) {
    std::cerr << "NOX::MultiVector:  Error!  Multivector" 
	 << " must have postive number of columns!" << std::endl;
    throw "NOX Error";
  }

  for (int i=0; i<numVecs; i++) {
    vecs[i] = vs[i]->clone(type);
  }
}

NOX::MultiVector::MultiVector(const NOX::MultiVector& source,
			      NOX::CopyType type) :
  vecs(source.vecs.size())
{
  for (unsigned int i=0; i<source.vecs.size(); i++) {
    vecs[i] = source.vecs[i]->clone(type);
  }
}

NOX::MultiVector::~MultiVector()
{

}

NOX::Abstract::MultiVector& 
NOX::MultiVector::init(double gamma)
{
  for (unsigned int i=0; i<vecs.size(); i++)
    vecs[i]->init(gamma);
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::MultiVector::random(bool useSeed, int seed)
{
  if (vecs.size() > 0)
    vecs[0]->random(useSeed, seed);
  for (unsigned int i=1; i<vecs.size(); i++)
    vecs[i]->random();
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::MultiVector::operator=(const NOX::Abstract::MultiVector& source)
{
  return operator=(dynamic_cast<const NOX::MultiVector&>(source));
}

NOX::Abstract::MultiVector& 
NOX::MultiVector::operator=(const NOX::MultiVector& source)
{
  if (this != &source) {
    checkSize(source.vecs.size());
    for (unsigned int i=0; i<vecs.size(); i++)
      *(vecs[i]) = *(source.vecs[i]);
  }

  return *this;
}

NOX::Abstract::MultiVector&
NOX::MultiVector::setBlock(const NOX::Abstract::MultiVector& source, 
			   const std::vector<int>& index)
{
  return setBlock(dynamic_cast<const NOX::MultiVector&>(source), index);
}

NOX::Abstract::MultiVector&
NOX::MultiVector::setBlock(const NOX::MultiVector& source, 
			   const std::vector<int>& index)
{
  int ind;

  source.checkIndex(index.size()-1);
  for (unsigned int i=0; i<index.size(); i++) {
    ind = index[i];
    checkIndex(ind);
    *(vecs[ind]) = *(source.vecs[i]);
  }

  return *this;
}

NOX::Abstract::MultiVector&
NOX::MultiVector::augment(const NOX::Abstract::MultiVector& source)
{
  return augment(dynamic_cast<const NOX::MultiVector&>(source));
}

NOX::Abstract::MultiVector&
NOX::MultiVector::augment(const NOX::MultiVector& source)
{
  int sz = vecs.size();
  int newsize = sz + source.vecs.size();
  vecs.resize(newsize);

  for (unsigned int i=0; i<source.vecs.size(); i++) {
    vecs[sz+i] = source.vecs[i]->clone(NOX::DeepCopy);
  }

  return *this;
}

NOX::Abstract::Vector&
NOX::MultiVector::operator [] (int i)
{
  checkIndex(i);
  return *(vecs[i]);
}

const NOX::Abstract::Vector&
NOX::MultiVector::operator [] (int i) const
{
  checkIndex(i);
  return *(vecs[i]);
}

NOX::Abstract::MultiVector& 
NOX::MultiVector::scale(double gamma)
{
  for (unsigned int i=0; i<vecs.size(); i++)
    vecs[i]->scale(gamma);
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::MultiVector::update(double alpha, const NOX::Abstract::MultiVector& a, 
			 double gamma)
{
  return update(alpha, dynamic_cast<const NOX::MultiVector&>(a), gamma);
}

NOX::Abstract::MultiVector& 
NOX::MultiVector::update(double alpha, const NOX::MultiVector& a, 
			 double gamma)
{
  checkSize(a.vecs.size());
  for (unsigned int i=0; i<vecs.size(); i++) 
    vecs[i]->update(alpha, *(a.vecs[i]), gamma);
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::MultiVector::update(double alpha, const NOX::Abstract::MultiVector& a, 
			 double beta, const NOX::Abstract::MultiVector& b,
			 double gamma)
{
  return update(alpha, dynamic_cast<const NOX::MultiVector&>(a),
		beta, dynamic_cast<const NOX::MultiVector&>(b), gamma);
}

NOX::Abstract::MultiVector& 
NOX::MultiVector::update(double alpha, const NOX::MultiVector& a, 
			 double beta, const NOX::MultiVector& b,
			 double gamma)
{
  checkSize(a.vecs.size());
  checkSize(b.vecs.size());
  for (unsigned int i=0; i<vecs.size(); i++) 
    vecs[i]->update(alpha, *(a.vecs[i]), beta, *(b.vecs[i]), gamma);
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::MultiVector::update(Teuchos::ETransp transb, double alpha, 
			 const NOX::Abstract::MultiVector& a, 
			 const NOX::Abstract::MultiVector::DenseMatrix& b, 
			 double gamma)
{
  return update(transb, alpha, dynamic_cast<const NOX::MultiVector&>(a), b, 
		gamma);
}

NOX::Abstract::MultiVector& 
NOX::MultiVector::update(Teuchos::ETransp transb, double alpha, 
			 const NOX::MultiVector& a, 
			 const NOX::Abstract::MultiVector::DenseMatrix& b, 
			 double gamma)
{
  if (transb == Teuchos::NO_TRANS) {
    a.checkSize(b.numRows());
    checkSize(b.numCols());
  }
  else {
    a.checkSize(b.numCols());
    checkSize(b.numRows());
  }
  
  int sz_a = a.vecs.size();
  int p = sz_a / 2;
  int q = sz_a - 2*p;

  if (transb == Teuchos::NO_TRANS) {
    for (unsigned int i=0; i<vecs.size(); i++) {

      if (p == 0)
	vecs[i]->update(alpha*b(0,i), *(a.vecs[0]), gamma);
      else {
	vecs[i]->update(alpha*b(0,i), *(a.vecs[0]), 
			alpha*b(1,i), *(a.vecs[1]), gamma);
	
	for (int j=1; j<p; j++) 
	  vecs[i]->update(alpha*b(2*j,i), *(a.vecs[2*j]), 
			  alpha*b(2*j+1,i), *(a.vecs[2*j+1]), 1.0);
	
	if (q > 0)
	  vecs[i]->update(alpha*b(sz_a-1,i), *(a.vecs[sz_a-1]), 1.0);
      }
      
    }

  }
  else {
    
    for (unsigned int i=0; i<vecs.size(); i++) {

      if (p == 0)
	vecs[i]->update(alpha*b(i,0), *(a.vecs[0]), gamma);
      else {
	vecs[i]->update(alpha*b(i,0), *(a.vecs[0]), 
			alpha*b(i,1), *(a.vecs[1]), gamma);
	
	for (int j=1; j<p; j++) 
	  vecs[i]->update(alpha*b(i,2*j), *(a.vecs[2*j]), 
			  alpha*b(i,2*j+1), *(a.vecs[2*j+1]), 1.0);
	
	if (q > 0)
	  vecs[i]->update(alpha*b(i,sz_a-1), *(a.vecs[sz_a-1]), 1.0);
      }
      
    }

  }
  return *this;
}

Teuchos::RCP<NOX::Abstract::MultiVector> 
NOX::MultiVector::clone(NOX::CopyType type) const
{
  Teuchos::RCP<NOX::Abstract::MultiVector> tmp = 
    Teuchos::rcp(new NOX::MultiVector(*this, type));
  return tmp;
}

Teuchos::RCP<NOX::Abstract::MultiVector> 
NOX::MultiVector::clone(int numvecs) const
{
  Teuchos::RCP<NOX::MultiVector> tmp = 
    Teuchos::rcp(new NOX::MultiVector(numvecs));
  for (int i=0; i<numvecs; i++) {
    tmp->vecs[i] = vecs[0]->clone(NOX::ShapeCopy);
  }
  return tmp;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::MultiVector::subCopy(const std::vector<int>& index) const
{
  Teuchos::RCP<NOX::MultiVector> tmp = 
    Teuchos::rcp(new NOX::MultiVector(index.size()));
  int ind;
  for (unsigned int i=0; i<index.size(); i++) {
    ind = index[i];
    checkIndex(ind);
    tmp->vecs[i] = vecs[ind]->clone(NOX::DeepCopy);
  }
  return tmp;
}

Teuchos::RCP<NOX::Abstract::MultiVector> 
NOX::MultiVector::subView(const std::vector<int>& index) const
{
  Teuchos::RCP<NOX::MultiVector> tmp = 
    Teuchos::rcp(new NOX::MultiVector(index.size()));
  int ind;
  for (unsigned int i=0; i<index.size(); i++) {
    ind = index[i];
    checkIndex(ind);
    tmp->vecs[i] = vecs[ind];
  }
  return tmp;
}

void
NOX::MultiVector::norm(std::vector<double>& result,
		       NOX::Abstract::Vector::NormType type) const
{
  if (result.size() != vecs.size())
    result.resize(vecs.size());

  for (unsigned int i=0; i<vecs.size(); i++)
    result[i] = vecs[i]->norm(type);
}

void 
NOX::MultiVector::multiply(double alpha, const NOX::Abstract::MultiVector& y,
			   NOX::Abstract::MultiVector::DenseMatrix& b) const
{
  multiply(alpha, dynamic_cast<const NOX::MultiVector&>(y), b);
}

void 
NOX::MultiVector::multiply(double alpha, const NOX::MultiVector& y,
			   NOX::Abstract::MultiVector::DenseMatrix& b) const
{
  for (unsigned int i=0; i<y.vecs.size(); i++) {
    for (unsigned int j=0; j<vecs.size(); j++) {
      b(i,j) = alpha*(y.vecs[i]->innerProduct(*(vecs[j])));
    }
  }
}

int 
NOX::MultiVector::length() const
{
  return vecs[0]->length();
}

int 
NOX::MultiVector::numVectors() const
{
  return vecs.size();
}

void 
NOX::MultiVector::print(std::ostream& stream) const
{
  for (unsigned int i=0; i<vecs.size(); i++)
    vecs[i]->print(stream);
}

void 
NOX::MultiVector::checkIndex(int idx) const 
{
  if ( idx < 0 || idx >= static_cast<int>(vecs.size()) ) {
    std::cerr << "NOX::MultiVector:  Error!  Invalid index " << idx << std::endl;
    throw "NOX Error";
  }
}

void
NOX::MultiVector::checkSize(int sz) const
{
  if (static_cast<int>(vecs.size()) != sz) {
    std::cerr << "NOX::MultiVector:  Error!  Size of supplied multivector is"
	 << " incompatible with this multivector" << std::endl;
    throw "NOX Error";
  }
}
