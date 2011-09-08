/*
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
*/

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "EpetraThyraAdaptersTestHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


//
// Unit Tests
//


TEUCHOS_UNIT_TEST( get_Epetra_MultiVector, singleBlockProductVector )
{
  using Teuchos::Comm;
  typedef Teuchos_Ordinal Ordinal;
  using Thyra::VectorSpaceBase;
  using Thyra::MultiVectorBase;

  const RCP<const Epetra_Comm> epetra_comm = getEpetraComm();
  const RCP<const Comm<Ordinal> > comm = Thyra::create_Comm(epetra_comm);
  
  const RCP<const Epetra_Map> epetra_map = rcp(new Epetra_Map(g_localDim, 0, *epetra_comm));
  const RCP<const VectorSpaceBase<double> > vs =  Thyra::create_VectorSpace(epetra_map);

  const RCP<const VectorSpaceBase<double> > pvs = Thyra::productVectorSpace(vs, 1);

  const RCP<MultiVectorBase<double> > pmv = Thyra::createMembers(pvs, 1);

  const double alpha = 3.5;
  Thyra::assign<double>( pmv.ptr(), alpha );

  const RCP<Epetra_MultiVector> epetra_mv =
    Thyra::get_Epetra_MultiVector(*epetra_map, pmv);

  const RCP<MultiVectorBase<double> > mv2 =
    Thyra::create_MultiVector(epetra_mv, pvs);

  Thyra::testRelNormDiffErr<double>(
    "*pmv->col(0)", *pmv->col(0),
    "*mv2->col(0)", *mv2->col(0),
    "max-error", 0.0,
    "max-warning", 0.0,
    &out
    );
   
}


} // namespace
