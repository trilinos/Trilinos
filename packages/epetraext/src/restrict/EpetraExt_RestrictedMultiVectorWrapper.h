//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef EPETRAEXT_RESTRICTEDMULTIVECTORWRAPPER_H
#define EPETRAEXT_RESTRICTEDMULTIVECTORWRAPPER_H

#if defined(EpetraExt_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The EpetraExt package is deprecated"
#endif
#endif
#include <EpetraExt_ConfigDefs.h>
#include <Teuchos_RCP.hpp>
#include <mpi.h>

#ifdef HAVE_MPI
class Epetra_Comm;
class Epetra_MpiComm;
class Epetra_BlockMap;
class Epetra_MultiVector;

namespace EpetraExt {
class RestrictedMultiVectorWrapper{
  
public:
  RestrictedMultiVectorWrapper();
  
  ~RestrictedMultiVectorWrapper();
    
  /// \brief Set the subcommunicator over which to restrict multivectors.
  ///
  /// This sets the subcommunicator only if it has not already been
  /// set.  If it has been set, this routine does nothing.
  ///
  /// \param MPI_Comm [in] MPI communicator representing the
  ///   subcommunicator over which the multivector is to be
  ///   restricted.
  ///
  /// \return 0 if the subcommunicator has not already been set, else
  ///   -1 if it has already been set.
  int SetMPISubComm(MPI_Comm MPI_SubComm);

  //! Get the MPI communicator corresponding to the restricted multivector.
  MPI_Comm GetMPISubComm(){return MPI_SubComm_;}

  //! Get the Epetra communicator corresponding to the restricted multivector.
  const Epetra_MpiComm & RestrictedComm(){return *RestrictedComm_;}    

  //! Whether the calling process is active for the restricted multivector.
  bool RestrictedProcIsActive(){return proc_is_active;}

  //! The input multivector.
  Teuchos::RCP<Epetra_MultiVector> InputMultiVector(){return input_mv_;}

  //! The input multivector, restricted to the active processes.
  Teuchos::RCP<Epetra_MultiVector> RestrictedMultiVector(){return restricted_mv_;}

  /// \brief Restrict the given multivector.
  ///
  /// If the MPISubComm is not set, this method sets it to include
  /// only the active processes.  Then this method restricts the
  /// multivector.  If MPISubComm is set in a fashion compatible with
  /// the existing multivector (that is, all processes with active
  /// rows must be in the subcommunicator), it creates a multivector
  /// using that MPISubComm.  This method also sets the input
  /// multivector (what InputMultiVector() returns) to the given
  /// multivector.
  ///
  /// \return 0 if sucessful, -1 if the input_multivector is deficient
  ///   (e.g., its Epetra_Comm is not an Epetra_MpiComm), or -2 if the
  ///   MPI_Comm object set with SetMPISubComm is not consistent with
  ///   the input multivector or if splitting the communicator failed.
  int restrict_comm(Teuchos::RCP<Epetra_MultiVector> input_mv);

private:
  //! Whether the calling process is part of the restricted multivector.
  bool proc_is_active;
  //! Whether the subcommunicator has been set (using SetMPISubComm()).
  bool subcomm_is_set;
  
  MPI_Comm MPI_SubComm_;
  Epetra_MpiComm *RestrictedComm_;
  Epetra_BlockMap *ResMap_;
  
  Teuchos::RCP<Epetra_MultiVector> input_mv_;
  Teuchos::RCP<Epetra_MultiVector> restricted_mv_;
};


} // namespace EpetraExt


#endif
#endif /* EPETRAEXT_RESTRICTEDMULTIVECTORWRAPPER_H */
