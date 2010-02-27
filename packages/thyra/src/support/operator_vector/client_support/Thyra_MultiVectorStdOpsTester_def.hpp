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
