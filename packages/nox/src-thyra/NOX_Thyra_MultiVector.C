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

#include "NOX_Thyra_MultiVector.H"
#include "NOX_Thyra_Vector.H"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"

NOX::Thyra::MultiVector::
MultiVector(
       const Teuchos::RCP< ::Thyra::MultiVectorBase<double> >& source)
  : thyraMultiVec(source),
    noxThyraVectors(thyraMultiVec->domain()->dim())
{
}

NOX::Thyra::MultiVector::
MultiVector(const ::Thyra::MultiVectorBase<double>& source) 
  : thyraMultiVec(source.clone_mv()),
    noxThyraVectors(thyraMultiVec->domain()->dim())
{
}

NOX::Thyra::MultiVector::
MultiVector(const NOX::Thyra::MultiVector& source, NOX::CopyType type) 
  : thyraMultiVec(source.thyraMultiVec->clone_mv()),
    noxThyraVectors(thyraMultiVec->domain()->dim())
{
}

NOX::Thyra::MultiVector::
~MultiVector()
{
}

NOX::Abstract::MultiVector& 
NOX::Thyra::MultiVector::
operator=(const NOX::Abstract::MultiVector& src)
{
  const NOX::Thyra::MultiVector& source = 
    dynamic_cast<const NOX::Thyra::MultiVector&>(src);
  ::Thyra::assign(thyraMultiVec.get(), *source.thyraMultiVec);
  return *this;
}

Teuchos::RCP< ::Thyra::MultiVectorBase<double> > 
NOX::Thyra::MultiVector::
getThyraMultiVector()
{
  return thyraMultiVec;
}

Teuchos::RCP< const ::Thyra::MultiVectorBase<double> >
NOX::Thyra::MultiVector::
getThyraMultiVector() const
{
  return thyraMultiVec;
}

NOX::Abstract::MultiVector& 
NOX::Thyra::MultiVector::
init(double value)
{
  ::Thyra::assign(thyraMultiVec.get(), value);
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Thyra::MultiVector::
random(bool useSeed, int seed)
{
  if (useSeed)
    ::Thyra::seed_randomize<double>(seed);
  ::Thyra::randomize(-1.0, 1.0, thyraMultiVec.get());
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Thyra::MultiVector::
setBlock(const NOX::Abstract::MultiVector& src, const vector<int>& index) 
{
  const NOX::Thyra::MultiVector& source = 
    dynamic_cast<const NOX::Thyra::MultiVector&>(src);

  // Create view
  Teuchos::RCP< ::Thyra::MultiVectorBase<double> > v;
  if (isContiguous(index))
    v = thyraMultiVec->subView(::Thyra::Range1D(index[0],
						index[index.size()-1]));
  else
    v = thyraMultiVec->subView(index.size(), &index[0]);

  // Assign
  ::Thyra::assign(v.get(), *source.thyraMultiVec);

  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Thyra::MultiVector::
augment(const NOX::Abstract::MultiVector& src) 
{
  const NOX::Thyra::MultiVector& source = 
    dynamic_cast<const NOX::Thyra::MultiVector&>(src);

  // Total number of columns
  int num_cols = thyraMultiVec->domain()->dim() + 
    source.thyraMultiVec->domain()->dim();

  // Create new multivector with num_cols columns
  Teuchos::RCP< ::Thyra::MultiVectorBase<double> > new_mv = 
    ::Thyra::createMembers(thyraMultiVec->range(), num_cols);

  // Create views
  Teuchos::RCP< ::Thyra::MultiVectorBase<double> > v1 = 
    new_mv->subView(::Thyra::Range1D(0,thyraMultiVec->domain()->dim()-1));
  Teuchos::RCP< ::Thyra::MultiVectorBase<double> > v2 = 
    new_mv->subView(::Thyra::Range1D(thyraMultiVec->domain()->dim(),
				     num_cols-1));

  // Assign
  ::Thyra::assign(v1.get(), *thyraMultiVec);
  ::Thyra::assign(v2.get(), *source.thyraMultiVec);

  // Set new multivector
  thyraMultiVec = new_mv;
					     
  return *this;
}

NOX::Abstract::Vector&
NOX::Thyra::MultiVector::
operator [] (int i)
{
  if (noxThyraVectors[i] == Teuchos::null)
    noxThyraVectors[i] = 
      Teuchos::rcp(new NOX::Thyra::Vector(thyraMultiVec->col(i)));
  return *noxThyraVectors[i];
}

const NOX::Abstract::Vector&
NOX::Thyra::MultiVector::
operator [] (int i) const
{
  if (noxThyraVectors[i] == Teuchos::null)
    noxThyraVectors[i] = 
      Teuchos::rcp(new NOX::Thyra::Vector(thyraMultiVec->col(i)));
  return *noxThyraVectors[i];
}

NOX::Abstract::MultiVector& 
NOX::Thyra::MultiVector::
scale(double gamma)
{
  ::Thyra::scale(gamma, thyraMultiVec.get());
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Thyra::MultiVector::
update(double alpha, const NOX::Abstract::MultiVector& a, double gamma)
{
  const NOX::Thyra::MultiVector& aa = 
    dynamic_cast<const NOX::Thyra::MultiVector&>(a);
  const ::Thyra::MultiVectorBase<double> *tv = aa.thyraMultiVec.get();
  ::Thyra::linear_combination(1, &alpha, &tv, gamma, thyraMultiVec.get());
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Thyra::MultiVector::
update(double alpha, 
       const NOX::Abstract::MultiVector& a, 
       double beta, 
       const NOX::Abstract::MultiVector& b,
       double gamma)
{
  const NOX::Thyra::MultiVector& aa = 
    dynamic_cast<const NOX::Thyra::MultiVector&>(a);
  const NOX::Thyra::MultiVector& bb = 
    dynamic_cast<const NOX::Thyra::MultiVector&>(b);
  double c[] = {alpha, beta};
  const ::Thyra::MultiVectorBase<double>* z[] = {aa.thyraMultiVec.get(), 
						 bb.thyraMultiVec.get()};
  ::Thyra::linear_combination(2, c, z, gamma, thyraMultiVec.get());
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Thyra::MultiVector::
update(Teuchos::ETransp transb, double alpha, 
       const NOX::Abstract::MultiVector& a, 
       const NOX::Abstract::MultiVector::DenseMatrix& b, double gamma)
{
  const NOX::Thyra::MultiVector& aa = 
    dynamic_cast<const NOX::Thyra::MultiVector&>(a);
  
  int m = b.numRows();
  int n = b.numCols();
  Teuchos::RCP<const ::Thyra::MultiVectorBase<double> > bb;
  if (transb == Teuchos::NO_TRANS) {
    bb = ::Thyra::createMembersView(
	aa.thyraMultiVec->domain(), 
	RTOpPack::ConstSubMultiVectorView<double>(0,m,0,n,&b(0,0),b.stride()));
  }
  else {
    // Make a copy of the transpose of b
    NOX::Abstract::MultiVector::DenseMatrix btrans(n, m);
    for (int i=0; i<m; i++)
      for (int j=0; j<n; j++)
	btrans(j,i) = b(i,j);
    bb = ::Thyra::createMembersView(
	aa.thyraMultiVec->domain(), 
	RTOpPack::ConstSubMultiVectorView<double>(0,n,0,n,&btrans(0,0),
						  btrans.stride()));
  }
  aa.thyraMultiVec->apply(::Thyra::NONCONJ_ELE, *bb, thyraMultiVec.get(),
			  alpha, gamma);

  return *this;
}

Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::Thyra::MultiVector::
clone(CopyType type) const
{
  return Teuchos::rcp(new NOX::Thyra::MultiVector(*this, type));
}

Teuchos::RCP<NOX::Abstract::MultiVector> 
NOX::Thyra::MultiVector::
clone(int numvecs) const
{
  Teuchos::RCP< ::Thyra::MultiVectorBase<double> > mv =
    ::Thyra::createMembers(thyraMultiVec->range(), numvecs);
  return Teuchos::rcp(new NOX::Thyra::MultiVector(mv));
}

Teuchos::RCP<NOX::Abstract::MultiVector> 
NOX::Thyra::MultiVector::
subCopy(const vector<int>& index) const
{
  Teuchos::RCP< ::Thyra::MultiVectorBase<double> > mv;
  if (isContiguous(index))
    mv = thyraMultiVec->subView(::Thyra::Range1D(index[0],
						 index[index.size()-1]));
  else
    mv = thyraMultiVec->subView(index.size(), &index[0]);

  return Teuchos::rcp(new NOX::Thyra::MultiVector(*mv));
}

Teuchos::RCP<NOX::Abstract::MultiVector>
NOX::Thyra::MultiVector::
subView(const vector<int>& index) const
{
  Teuchos::RCP< ::Thyra::MultiVectorBase<double> > mv;
  if (isContiguous(index))
    mv = thyraMultiVec->subView(::Thyra::Range1D(index[0],
						 index[index.size()-1]));
  else
     mv = thyraMultiVec->subView(index.size(), &index[0]);
  
  return Teuchos::rcp(new NOX::Thyra::MultiVector(mv));
}

void 
NOX::Thyra::MultiVector::
norm(vector<double>& result, NOX::Abstract::Vector::NormType type) const
{
  if (type == NOX::Abstract::Vector::TwoNorm)
    ::Thyra::norms_2(*thyraMultiVec, &result[0]);
  else if (type == NOX::Abstract::Vector::OneNorm)
    ::Thyra::norms_1(*thyraMultiVec, &result[0]);
  else
    ::Thyra::norms_inf(*thyraMultiVec, &result[0]);
}

void 
NOX::Thyra::MultiVector::
multiply(double alpha, 
	 const NOX::Abstract::MultiVector& y,
	 NOX::Abstract::MultiVector::DenseMatrix& b) const
{
  const NOX::Thyra::MultiVector& yy = 
    dynamic_cast<const NOX::Thyra::MultiVector&>(y);

  int m = b.numRows();
  int n = b.numCols();
  Teuchos::RCP< ::Thyra::MultiVectorBase<double> > bb =
    ::Thyra::createMembersView(
	yy.thyraMultiVec->domain(), 
	RTOpPack::SubMultiVectorView<double>(0,m,0,n,&b(0,0),b.stride()));
    
  yy.thyraMultiVec->applyTranspose(::Thyra::NONCONJ_ELE, *thyraMultiVec, 
				   bb.get(), alpha, 0.0);
}

int 
NOX::Thyra::MultiVector::
length() const
{
  return thyraMultiVec->range()->dim();
}

int 
NOX::Thyra::MultiVector::
numVectors() const
{
  return thyraMultiVec->domain()->dim();
}

void 
NOX::Thyra::MultiVector::
print(std::ostream& stream) const
{
  //stream << *thyraMultiVec;
  thyraMultiVec->describe(
    *Teuchos::getFancyOStream(Teuchos::rcpFromRef(stream)),
    Teuchos::VERB_EXTREME);
  return;
}

bool
NOX::Thyra::MultiVector::
isContiguous(const vector<int>& index) const
{
  if (index.size()==0)
    return true;
  int i0 = index[0];
  for (std::size_t i=1; i<index.size(); i++)
    if (index[i] != static_cast<int>(i0+i))
      return false;
  return true;
}
