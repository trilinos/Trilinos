// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DESCRIBE_LINEAR_OP_HPP
#define THYRA_DESCRIBE_LINEAR_OP_HPP

#include "Thyra_describeLinearOpDecl.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_AssertOp.hpp"
#include "Thyra_AssertOp.hpp"

template<class RangeScalar, class DomainScalar>
std::ostream& Thyra::describeLinearOp(
  const LinearOpBase<RangeScalar,DomainScalar>   &A
  ,std::ostream                                  &out
  ,const Teuchos::EVerbosityLevel                verbLevel
  ,const std::string                             leadingIndent
  ,const std::string                             indentSpacer
  )
{
  typedef Teuchos::ScalarTraits<DomainScalar> DST;
  const Index dimDomain = A.domain()->dim(), dimRange = A.range()->dim();
  out << leadingIndent << indentSpacer << "type = \'" << A.description()
      << "\', rangeDim = " << dimRange
      << ", domainDim = " << dimDomain << "\n";
  if(verbLevel >= Teuchos::VERB_EXTREME) {
    // Copy into dense matrix by column
    Teuchos::RefCountPtr<VectorBase<DomainScalar> > e_j = createMember(A.domain());
    Teuchos::RefCountPtr<VectorBase<RangeScalar> >  t   = createMember(A.range()); // temp column
    RTOpPack::SubVectorT<RangeScalar> sv;
    std::vector<RangeScalar>  Md( dimRange * dimDomain ); // Column major
    const Index
      cs = 1,         // stride for columns or rows 
      rs = dimRange;  // stride for rows or columns
    Index i, j;
    for( j = 1; j <= dimDomain; ++j ) {
      Thyra::assign( e_j.get(), DST::zero() );
      Thyra::set_ele( j, DST::one(), e_j.get() );
      apply(A,NONCONJ_ELE,*e_j,t.get());  // extract the ith column or row
      t->getSubVector(Range1D(),&sv);
      for( i = 1; i <= dimRange; ++i ) Md[ (i-1)*cs + (j-1)*rs ] = sv(i);
      t->freeSubVector(&sv);
    }
    // Print the matrix
    for( i = 1; i <= dimRange; ++i ) {
      out << leadingIndent << indentSpacer << indentSpacer;
      for( j = 1; j <= dimDomain; ++j )
        out << " " << i << ":" << j << ":" << Md[ (i-1) + (j-1)*dimRange ];
      out << std::endl;
    }
  }
  return out;
}

#endif // THYRA_DESCRIBE_LINEAR_OP_HPP
