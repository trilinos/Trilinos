// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "Xpetra_ConfigDefs.hpp"

#ifdef HAVE_XPETRA_EPETRA

#include "Xpetra_EpetraUtils.hpp"

// header files for comm objects conversion
#ifdef HAVE_MPI
#include <mpi.h>
#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#endif
#include <Teuchos_DefaultSerialComm.hpp>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#endif
#include <Epetra_SerialComm.h>

#include "Xpetra_Exceptions.hpp"


namespace Xpetra {

  using Teuchos::RCP;

  const RCP<const Epetra_Comm> toEpetra(const RCP<const Teuchos::Comm<int> > & comm) {
#ifdef HAVE_MPI
    const RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
    if (mpiComm != Teuchos::null) {
      return Teuchos::rcp(new Epetra_MpiComm(*mpiComm->getRawMpiComm()));
    }  else
#endif
      if ((Teuchos::rcp_dynamic_cast<const Teuchos::SerialComm<int> >(comm) != Teuchos::null))
        return Teuchos::rcp(new Epetra_SerialComm());
      else
        TEUCHOS_TEST_FOR_EXCEPTION(1,Xpetra::Exceptions::BadCast,"Cannot convert a Teuchos::Comm to an Epetra_Comm: The exact type of the Teuchos::Comm object is unknown");
  }

  const RCP<const Teuchos::Comm<int> > toXpetra(const Epetra_Comm & comm) {
#ifdef HAVE_MPI
    try {
      const Epetra_MpiComm& mpiComm = dynamic_cast<const Epetra_MpiComm&>(comm);
      return Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(mpiComm.Comm())));
    } catch (std::bad_cast & b) {}
#endif
    try {
      const Epetra_SerialComm& serialComm = dynamic_cast<const Epetra_SerialComm&>(comm);
      serialComm.NumProc(); // avoid compilation warning
      return Teuchos::rcp(new Teuchos::SerialComm<int>());
    } catch (std::bad_cast & b) {
      TEUCHOS_TEST_FOR_EXCEPTION(1,Xpetra::Exceptions::BadCast,"Cannot convert an Epetra_Comm to a Teuchos::Comm: The exact type of the Epetra_Comm object is unknown");
    }
  }

  bool toEpetra(Teuchos::ETransp trans) {
    if (trans == Teuchos::NO_TRANS)
      return false;
    else if (trans == Teuchos::TRANS)
      return true;
    else {
      TEUCHOS_TEST_FOR_EXCEPTION((trans != Teuchos::NO_TRANS) && (trans == Teuchos::TRANS), Xpetra::Exceptions::NotImplemented, "Cannot convert Teuchos::ETransp to a boolean.");
    }

    return false; // to skip a compilation warning msg.
  }

}

#endif
