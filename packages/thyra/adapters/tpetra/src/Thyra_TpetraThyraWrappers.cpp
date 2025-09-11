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

template<typename InputOrdinalType, typename OutputOrdinalType>
Teuchos::RCP<const Teuchos::Comm<OutputOrdinalType> >
convertCommunicatorType(const Teuchos::RCP<const Teuchos::Comm<InputOrdinalType> > &inputComm)
{
  using Teuchos::rcp_dynamic_cast;

#ifdef HAVE_MPI
  const Teuchos::RCP<const Teuchos::MpiComm<InputOrdinalType> > inputMpiComm = 
    rcp_dynamic_cast<const Teuchos::MpiComm<InputOrdinalType> >(inputComm);
  if (Teuchos::nonnull(inputMpiComm)) {
    return Teuchos::createMpiComm<OutputOrdinalType>(inputMpiComm->getRawMpiComm(),inputMpiComm->getTag());
  }
#endif // HAVE_MPI

  // Assert conversion to Teuchos::SerialComm as a last resort (or throw)
  rcp_dynamic_cast<const Teuchos::SerialComm<InputOrdinalType> >(inputComm, true);
  return Teuchos::createSerialComm<OutputOrdinalType>();

  // NOTE: Above will throw if the type is not Teuchos::SerialComm.  In this
  // case, the type could not be converted.  We need to either get rid of the
  // Ordinal templating on Comm or we need to use the same ordinal type for
  // Tpetra and Thyra so this conversion function goes away!
}

Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> >
Thyra::convertTpetraToThyraComm(const RCP<const Teuchos::Comm<int> > &tpetraComm)
{
  return convertCommunicatorType<int, Thyra::Ordinal>(tpetraComm);
}

Teuchos::RCP<const Teuchos::Comm<int> >
Thyra::convertThyraToTpetraComm(const RCP<const Teuchos::Comm<Thyra::Ordinal> > &thyraComm)
{
  return convertCommunicatorType<Thyra::Ordinal, int>(thyraComm);
}
