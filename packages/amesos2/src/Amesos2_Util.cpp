// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Amesos2_Util.hpp"

#include <Teuchos_RCP.hpp>
#include <Teuchos_DefaultSerialComm.hpp>

#ifdef HAVE_MPI
#  include <Teuchos_DefaultMpiComm.hpp>
#  ifdef HAVE_AMESOS2_EPETRA
#    include <Teuchos_OpaqueWrapper.hpp>
#  endif // HAVE_AMESOS2_EPETRA
#endif // HAVE_MPI

#ifdef HAVE_AMESOS2_EPETRA
#  include <Epetra_Map.h>
#  ifdef HAVE_MPI
#    include <Epetra_MpiComm.h>
#  endif // HAVE_MPI
#  include <Epetra_SerialComm.h>
#endif // HAVE_AMESOS2_EPETRA

#ifdef HAVE_AMESOS2_EPETRA
const Teuchos::RCP<const Teuchos::Comm<int> >
Amesos2::Util::to_teuchos_comm(Teuchos::RCP<const Epetra_Comm> c)
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::set_extra_data;

#ifdef HAVE_MPI
  Teuchos::RCP<const Epetra_MpiComm>
    mpiEpetraComm = rcp_dynamic_cast<const Epetra_MpiComm>(c);
  if( mpiEpetraComm.get() ) {
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >
      rawMpiComm = Teuchos::opaqueWrapper(mpiEpetraComm->Comm());
    return Teuchos::createMpiComm<int>(rawMpiComm);
  } else
#endif // HAVE_MPI
    if( rcp_dynamic_cast<const Epetra_SerialComm>(c) != Teuchos::null )
      return Teuchos::createSerialComm<int>();
    else
      return(Teuchos::null);
}

const Teuchos::RCP<const Epetra_Comm>
Amesos2::Util::to_epetra_comm(Teuchos::RCP<const Teuchos::Comm<int> > c)
{
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::set_extra_data;

#ifdef HAVE_MPI
  Teuchos::RCP<const Teuchos::MpiComm<int> >
    mpiTeuchosComm = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(c);
  if( mpiTeuchosComm.get() ){
    Teuchos::RCP<const Teuchos::OpaqueWrapper<MPI_Comm> >
      rawMpiComm = mpiTeuchosComm->getRawMpiComm();
    Teuchos::RCP<const Epetra_MpiComm>
      mpiComm = rcp(new Epetra_MpiComm(*rawMpiComm()));
    return mpiComm;
  }
#else // NOT HAVE_MPI
  Teuchos::RCP<const Teuchos::SerialComm<int> >
    serialTeuchosComm = rcp_dynamic_cast<const Teuchos::SerialComm<int> >(c);

  if( serialTeuchosComm.get() ){
    Teuchos::RCP<const Epetra_SerialComm> serialComm = rcp(new Epetra_SerialComm());
    return serialComm;
  }
#endif // HAVE_MPI

  return Teuchos::null;
}
#endif // HAVE_AMESOS2_EPETRA

/// Prints a line of 80 "-"s on out.
void Amesos2::Util::printLine( Teuchos::FancyOStream& out )
{
  out << "----------------------------------------"
      << "----------------------------------------"
      << std::endl;
}


