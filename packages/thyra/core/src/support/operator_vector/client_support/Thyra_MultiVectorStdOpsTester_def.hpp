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

#ifndef THYRA_MULTI_VECTOR_STD_OPS_TESTER_HPP
#define THYRA_MULTI_VECTOR_STD_OPS_TESTER_HPP

#include "Thyra_MultiVectorStdOpsTester_decl.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"

namespace Thyra {

// MultiVectorStdOpsTester

template <class Scalar>
MultiVectorStdOpsTester<Scalar>::MultiVectorStdOpsTester(
  const ScalarMag    &warning_tol_in
  ,const ScalarMag   &error_tol_in
  ,const int         num_mv_cols_in
  )
  :warning_tol_(warning_tol_in)
  ,error_tol_(error_tol_in)
  ,num_mv_cols_(num_mv_cols_in)
{}

template <class Scalar>
bool MultiVectorStdOpsTester<Scalar>::checkStdOps(
  const VectorSpaceBase<Scalar>    &vecSpc
  ,std::ostream                    *out
  ,const bool                      &dumpAll
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;

  if(out)
    *out << "\n*** Entering MultiVectorStdOpsTester<"<<ST::name()<<">::checkStdOps(...) ...\n"
         << "using a \'" << vecSpc.description() << "\' object ...\n";

  bool success = true;
  if(out) *out << "\nvecSpc.dim() = " << vecSpc.dim() << std::endl;

  if(out) *out << "\nCreating MultiVectorBase objects V1, V2, V3 and Z ...\n";
  Teuchos::RCP<MultiVectorBase<Scalar> >
    V1 = createMembers(vecSpc,num_mv_cols()),
    V2 = createMembers(vecSpc,num_mv_cols()),
    V3 = createMembers(vecSpc,num_mv_cols()),
    Z  = createMembers(vecSpc,num_mv_cols());

  if(out) *out << "\nassign(&*V1,-2.0);\n";
  assign(&*V1,Scalar(-2.0));
  if(out) *out << "\nassign(&*V2,-3.0);\n";
  assign(&*V2,Scalar(-3.0));
  if(out) *out << "\nassign(&*V3,-4.0);\n";
  assign(&*V3,Scalar(-4.0));

  // ToDo: Fill in the tests!

  if(out) *out
    << "\n*** Leaving MultiVectorStdOpsTester<"<<ST::name()<<">::checkStdOps(...) ...\n";

  return success;

}

} // namespace Thyra

#endif // THYRA_MULTI_VECTOR_STD_OPS_TESTER_HPP
