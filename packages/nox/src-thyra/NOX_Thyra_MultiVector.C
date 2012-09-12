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

#include "NOX_Thyra_MultiVector.H"
#include "NOX_Thyra_Vector.H"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
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
  ::Thyra::assign(thyraMultiVec.ptr(), *source.thyraMultiVec);
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
  ::Thyra::assign(thyraMultiVec.ptr(), value);
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Thyra::MultiVector::
random(bool useSeed, int seed)
{
  if (useSeed)
    ::Thyra::seed_randomize<double>(seed);
  ::Thyra::randomize(-1.0, 1.0, thyraMultiVec.ptr());
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
    v = thyraMultiVec->subView(Teuchos::arrayViewFromVector(index));

  // Assign
  ::Thyra::assign(v.ptr(), *source.thyraMultiVec);

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
  ::Thyra::assign(v1.ptr(), *thyraMultiVec);
  ::Thyra::assign(v2.ptr(), *source.thyraMultiVec);

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
  ::Thyra::scale(gamma, thyraMultiVec.ptr());
  return *this;
}

NOX::Abstract::MultiVector& 
NOX::Thyra::MultiVector::
update(double alpha, const NOX::Abstract::MultiVector& a, double gamma)
{
  using Teuchos::tuple;
  const NOX::Thyra::MultiVector& aa = 
    dynamic_cast<const NOX::Thyra::MultiVector&>(a);
  ::Thyra::linear_combination<double>(tuple(alpha)().getConst(),
    tuple(aa.thyraMultiVec.ptr().getConst())(), gamma,
    thyraMultiVec.ptr());
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
  using Teuchos::tuple;
  const NOX::Thyra::MultiVector& aa = 
    dynamic_cast<const NOX::Thyra::MultiVector&>(a);
  const NOX::Thyra::MultiVector& bb = 
    dynamic_cast<const NOX::Thyra::MultiVector&>(b);
  ::Thyra::linear_combination<double>(tuple<double>(alpha, beta)().getConst(),
    tuple(aa.thyraMultiVec.ptr().getConst(), bb.thyraMultiVec.ptr().getConst())(),
    gamma, thyraMultiVec.ptr());
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
      RTOpPack::ConstSubMultiVectorView<double>(
        0, m, 0, n,
        Teuchos::arcp(b.values(), 0, b.stride()*b.numCols(), false),
        b.stride())
      );
    // rabartl: ToDo: This conversion from a NOX::Abstract::DenseMatrix to an
    // RTOpPack::ConstSubMultiVectorView objects needs a utility function very
    // badly!
  }
  else {
    // Make a copy of the transpose of b
    NOX::Abstract::MultiVector::DenseMatrix btrans(n, m);
    for (int i=0; i<m; i++) {
      for (int j=0; j<n; j++) {
        btrans(j,i) = b(i,j);
      }
    }
    bb = ::Thyra::createMembersView(
      aa.thyraMultiVec->domain(), 
      RTOpPack::ConstSubMultiVectorView<double>(0, n, 0, n,
        Teuchos::arcp(btrans.values(), 0, btrans.stride() * btrans.numCols(), false),
        btrans.stride())
      );
  }
  ::Thyra::apply(*aa.thyraMultiVec, ::Thyra::NOTRANS, *bb, thyraMultiVec.ptr(),
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
    mv = thyraMultiVec->subView(Teuchos::arrayViewFromVector(index));

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
    mv = thyraMultiVec->subView(Teuchos::arrayViewFromVector(index));
  
  return Teuchos::rcp(new NOX::Thyra::MultiVector(mv));
}

void 
NOX::Thyra::MultiVector::
norm(vector<double>& result, NOX::Abstract::Vector::NormType type) const
{
  using Teuchos::arrayViewFromVector;
  if (type == NOX::Abstract::Vector::TwoNorm)
    ::Thyra::norms_2<double>(*thyraMultiVec, arrayViewFromVector(result));
  else if (type == NOX::Abstract::Vector::OneNorm)
    ::Thyra::norms_1<double>(*thyraMultiVec, arrayViewFromVector(result));
  else
    ::Thyra::norms_inf<double>(*thyraMultiVec, arrayViewFromVector(result));
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
      RTOpPack::SubMultiVectorView<double>(0, m, 0, n,
        Teuchos::arcp(b.values(), 0, b.stride()*b.numCols(), false),
        b.stride()
        )
      );
  
  ::Thyra::apply(*yy.thyraMultiVec, ::Thyra::CONJTRANS, *thyraMultiVec, 
				   bb.ptr(), alpha, 0.0);
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
