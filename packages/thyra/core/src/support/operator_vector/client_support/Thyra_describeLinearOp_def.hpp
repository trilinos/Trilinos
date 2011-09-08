// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DESCRIBE_LINEAR_OP_HPP
#define THYRA_DESCRIBE_LINEAR_OP_HPP

#include "Thyra_describeLinearOp_def.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_AssertOp.hpp"


template<class Scalar>
void Thyra::describeLinearOp(
  const LinearOpBase<Scalar> &A,
  Teuchos::FancyOStream &out_arg,
  const Teuchos::EVerbosityLevel verbLevel
  )
{
  using Teuchos::RCP;
  using Teuchos::FancyOStream;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> DST;

  RCP<FancyOStream> out = rcp(&out_arg,false);
  OSTab tab(out);
  *out << A.description() << "\n";

  const Teuchos::RCP<const VectorSpaceBase<Scalar> >
    range = A.range();
  const Teuchos::RCP<const VectorSpaceBase<Scalar> >
    domain = A.domain();

  if(!range.get()) {
    return;
  }

  const Ordinal dimDomain = domain->dim(), dimRange = range->dim();
  if ( dimDomain > 0 && dimRange > 0 && verbLevel >= Teuchos::VERB_EXTREME ) {
    // Copy into dense matrix by column
    Teuchos::RCP<VectorBase<Scalar> > e_j = createMember(domain);
    Teuchos::RCP<VectorBase<Scalar> > t = createMember(range); // temp column
    RTOpPack::ConstSubVectorView<Scalar> sv;
    Array<Scalar>  Md(dimRange*dimDomain); // Column major
    const Ordinal
      cs = 1,         // stride for columns or rows 
      rs = dimRange;  // stride for rows or columns
    Ordinal i, j;
    OSTab tab2(out);
    for( j = 0; j < dimDomain; ++j ) {
      Thyra::assign( e_j.ptr(), DST::zero() );
      Thyra::set_ele( j, DST::one(), e_j.ptr() );
      Thyra::apply<Scalar>(A, NOTRANS, *e_j, t.ptr());  // extract the ith column or row
      t->acquireDetachedView(Range1D(),&sv);
      for( i = 0; i < dimRange; ++i ) Md[ i*cs + j*rs ] = sv(i);
      t->releaseDetachedView(&sv);
    }
    // Print the matrix
    for( i = 0; i < dimRange; ++i ) {
      for( j = 0; j < dimDomain; ++j )
        *out << " " << i << ":" << j << ":" << Md[ i + j*dimRange ];
      *out << std::endl;
    }
  }

}


//
// Explicit instant macro
//

#define THYRA_DESCRIBE_LINEAR_INSTANT(SCALAR) \
    \
  template void describeLinearOp(  \
    const LinearOpBase<SCALAR > &A,  \
    Teuchos::FancyOStream &out_arg,  \
    const Teuchos::EVerbosityLevel verbLevel  \
    );  \



#endif // THYRA_DESCRIBE_LINEAR_OP_HPP
