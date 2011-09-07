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

//#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_TpetraThyraWrappers_decl.hpp"


#include "Teuchos_DefaultSerialComm.hpp"
#ifdef HAVE_MPI
#  include "Teuchos_DefaultMpiComm.hpp"
#endif


Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> >
Thyra::convertTpetraToThyraComm(const RCP<const Teuchos::Comm<int> > &tpetraComm)
{

  using Teuchos::rcp_dynamic_cast;

#ifdef HAVE_MPI
  const RCP<const Teuchos::MpiComm<int> > tpetraMpiComm = 
    rcp_dynamic_cast<const Teuchos::MpiComm<int> >(tpetraComm);
  if (nonnull(tpetraMpiComm)) {
    return Teuchos::createMpiComm<Ordinal>(tpetraMpiComm->getRawMpiComm());
  }
#endif // HAVE_MPI

  // Assert conversion to Teuchos::SerialComm as a last resort (or throw)
  rcp_dynamic_cast<const Teuchos::SerialComm<int> >(tpetraComm, true);
   return Teuchos::createSerialComm<Ordinal>();

  // NOTE: Above will throw if the type is not Teuchos::SerialComm.  In this
  // case, the type could not be converted.  We need to either get rid of the
  // Ordinal templating on Comm or we need to use the same ordinal type for
  // Tpetra and Thyra so this conversion function goes away!

}
