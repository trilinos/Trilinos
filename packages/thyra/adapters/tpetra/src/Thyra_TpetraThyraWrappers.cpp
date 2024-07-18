// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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
    return Teuchos::createMpiComm<Ordinal>(tpetraMpiComm->getRawMpiComm(),tpetraMpiComm->getTag());
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
