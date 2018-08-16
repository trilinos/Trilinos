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

#include "NOX_Thyra_Vector.H"
#include "NOX_Thyra_MultiVector.H"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Teuchos_Tuple.hpp"

NOX::Thyra::Vector::
Vector(const Teuchos::RCP< ::Thyra::VectorBase<double> >& source) :
  thyraVec(source),
  do_implicit_weighting_(false)
{
}

NOX::Thyra::Vector::
Vector(const ::Thyra::VectorBase<double>& source) :
  thyraVec(source.clone_v()),
  do_implicit_weighting_(false)
{
}

NOX::Thyra::Vector::
Vector(const NOX::Thyra::Vector& source, NOX::CopyType /* type */) :
  thyraVec(source.thyraVec->clone_v()),
  do_implicit_weighting_(false)
{
  if (nonnull(source.weightVec_)) {
    weightVec_ = source.weightVec_;
    tmpVec_ = source.tmpVec_->clone_v();
    do_implicit_weighting_ = source.do_implicit_weighting_;
  }
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
  // NOTE that we only copy over values, we never copy over weighting.
  // Weighting property only transfers via clone or copy ctor!
  ::Thyra::copy(*src.thyraVec, thyraVec.ptr());
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
  ::Thyra::put_scalar(value, outArg(*thyraVec));
  return *this;
}

NOX::Abstract::Vector&
NOX::Thyra::Vector::
random(bool useSeed, int seed)
{
  if (useSeed)
    ::Thyra::seed_randomize<double>(seed);
  ::Thyra::randomize(-1.0, 1.0, outArg(*thyraVec));
  return *this;
}

NOX::Abstract::Vector&
NOX::Thyra::Vector::
abs(const NOX::Abstract::Vector& src)
{
  const NOX::Thyra::Vector& source =
    dynamic_cast<const NOX::Thyra::Vector&>(src);
  ::Thyra::abs(*source.thyraVec, outArg(*thyraVec));
  return *this;
}

NOX::Abstract::Vector&
NOX::Thyra::Vector::
reciprocal(const NOX::Abstract::Vector& src)
{
  const NOX::Thyra::Vector& source =
    dynamic_cast<const NOX::Thyra::Vector&>(src);
  ::Thyra::reciprocal(*source.thyraVec, outArg(*thyraVec));
  return *this;
}

NOX::Abstract::Vector&
NOX::Thyra::Vector::
scale(double alpha)
{
  ::Thyra::scale(alpha, outArg(*thyraVec));
  return *this;
}

NOX::Abstract::Vector&
NOX::Thyra::Vector::
scale(const NOX::Abstract::Vector& src)
{
  const NOX::Thyra::Vector& source =
    dynamic_cast<const NOX::Thyra::Vector&>(src);
  ::Thyra::ele_wise_scale(*source.thyraVec, outArg(*thyraVec));
  return *this;
}

NOX::Abstract::Vector&
NOX::Thyra::Vector::
update(double alpha, const NOX::Abstract::Vector& a, double gamma)
{
  const NOX::Thyra::Vector& aa =
    dynamic_cast<const NOX::Thyra::Vector&>(a);

  ::Thyra::linear_combination<double>(
      Teuchos::tuple<double>(alpha)(),
      Teuchos::tuple<Teuchos::Ptr<const ::Thyra::VectorBase<double> > >(aa.thyraVec.ptr())(),
      gamma,
      outArg(*thyraVec)
      );

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

  ::Thyra::linear_combination<double>(
      Teuchos::tuple<double>(alpha,beta)(),
      Teuchos::tuple<Teuchos::Ptr<const ::Thyra::VectorBase<double> > >(xx.thyraVec.ptr(),yy.thyraVec.ptr())(),
      gamma,
      outArg(*thyraVec)
      );

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
  TEUCHOS_TEST_FOR_EXCEPTION(nonnull(weightVec_), std::logic_error,
                 "Can NOT create NOX::Thyra::MultiVector from a NOX::Thyra::Vector that contains a weighting vector! ");

  // Get vector space
  Teuchos::RCP<const ::Thyra::VectorSpaceBase<double> > space =
    thyraVec->space();

  // Create Thyra multivector
  Teuchos::RCP< ::Thyra::MultiVectorBase<double> > mv =
    ::Thyra::createMembers(*space, numVecs+1);

  // Copy vectors
  if (type == NOX::DeepCopy) {
    Teuchos::RCP< ::Thyra::VectorBase<double> > v = mv->col(0);
    ::Thyra::copy(*thyraVec, v.ptr());
    for (int i=0; i<numVecs; i++) {
      const NOX::Thyra::Vector & tv =
        dynamic_cast<const NOX::Thyra::Vector&>(*vecs[i]);
      v = mv->col(i+1);
      ::Thyra::copy(*(tv.thyraVec), v.ptr());
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
      ::Thyra::copy(*thyraVec, v.ptr());
    }
  }

  // Create multi-vector
  Teuchos::RCP<NOX::Thyra::MultiVector> nmv =
    Teuchos::rcp(new NOX::Thyra::MultiVector(mv));

  nmv->setImplicitWeighting(do_implicit_weighting_);
  nmv->setWeightVector(weightVec_);

  return nmv;
}

double
NOX::Thyra::Vector::
norm(NOX::Abstract::Vector::NormType type) const
{
  if (nonnull(weightVec_) && do_implicit_weighting_) {
    ::Thyra::copy(*thyraVec, outArg(*tmpVec_));
    ::Thyra::ele_wise_scale(*weightVec_, outArg(*tmpVec_));

    if (type == NOX::Abstract::Vector::TwoNorm)
      return ::Thyra::norm_2(*tmpVec_);
    else if (type == NOX::Abstract::Vector::OneNorm)
      return ::Thyra::norm_1(*tmpVec_);
    else
      return ::Thyra::norm_inf(*tmpVec_);
  }
  else {
    if (type == NOX::Abstract::Vector::TwoNorm)
      return ::Thyra::norm_2(*thyraVec);
    else if (type == NOX::Abstract::Vector::OneNorm)
      return ::Thyra::norm_1(*thyraVec);
    else
      return ::Thyra::norm_inf(*thyraVec);
  }
}

double
NOX::Thyra::Vector::
norm(const NOX::Abstract::Vector& weights) const
{
  const NOX::Thyra::Vector& w =
    dynamic_cast<const NOX::Thyra::Vector&>(weights);

  if (nonnull(weightVec_) && do_implicit_weighting_) {
    ::Thyra::copy(*thyraVec, outArg(*tmpVec_));
    ::Thyra::ele_wise_scale(*weightVec_, outArg(*tmpVec_));
    return ::Thyra::norm_2(*w.thyraVec, *tmpVec_);
  }

  return ::Thyra::norm_2(*w.thyraVec, *thyraVec);
}

double
NOX::Thyra::Vector::
innerProduct(const NOX::Abstract::Vector& y) const
{
  const NOX::Thyra::Vector& yy =
    dynamic_cast<const NOX::Thyra::Vector&>(y);

  if (nonnull(weightVec_) && do_implicit_weighting_) {
    ::Thyra::copy(*thyraVec, tmpVec_.ptr());
    // double the scaling one for each vector in product
    ::Thyra::ele_wise_scale(*weightVec_, outArg(*tmpVec_));
    ::Thyra::ele_wise_scale(*weightVec_, outArg(*tmpVec_));
    return thyraVec->space()->scalarProd(*tmpVec_, *yy.thyraVec);
  }

  return thyraVec->space()->scalarProd(*thyraVec, *yy.thyraVec);
}

NOX::size_type
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


void
NOX::Thyra::Vector::
setWeightVector(const Teuchos::RCP<const ::Thyra::VectorBase<double> >& weightVec)
{
  weightVec_ = weightVec;
  tmpVec_ = weightVec_->clone_v();
  do_implicit_weighting_ = true;
}

bool NOX::Thyra::Vector::hasWeightVector() const
{
  return Teuchos::nonnull(weightVec_);
}

Teuchos::RCP<const ::Thyra::VectorBase<double> >
NOX::Thyra::Vector::getWeightVector() const
{
  return weightVec_;
}

bool NOX::Thyra::Vector::getImplicitWeighting() const
{
  return do_implicit_weighting_;

}

void NOX::Thyra::Vector::setImplicitWeighting(bool do_implicit_weighting)
{
  do_implicit_weighting_ = do_implicit_weighting;
}
