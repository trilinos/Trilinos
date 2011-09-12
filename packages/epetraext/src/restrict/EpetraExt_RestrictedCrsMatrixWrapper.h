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

#ifndef EPETRAEXT_RESTRICTEDCRSMATRIXWRAPPER_H
#define EPETRAEXT_RESTRICTEDCRSMATRIXWRAPPER_H
#include <EpetraExt_ConfigDefs.h>
#include <Teuchos_RCP.hpp>
#include <mpi.h>

#ifdef HAVE_MPI
class Epetra_Comm;
class Epetra_MpiComm;
class Epetra_Map;
class Epetra_CrsMatrix;

namespace EpetraExt {
class RestrictedCrsMatrixWrapper{
  
public:
  RestrictedCrsMatrixWrapper();
  
  ~RestrictedCrsMatrixWrapper();
    
  //! Set the MPI communicator corresponding to the restricted matrix.
  /*! Sets the MPISubComm if it hasn't been set already.
    \param MPI_Comm (In) An MPI_Comm object representing the subcommunicator.

    \return Returns 0 if the SubComm hasn't ready been set , -1 if it has */  
  int SetMPISubComm(MPI_Comm MPI_SubComm);

  //! Get the MPI communicator corresponding to the restricted matrix.
  MPI_Comm GetMPISubComm(){return MPI_SubComm_;}

  //! Get the Epetra communicator corresponding to the restricted matrix.
  const Epetra_MpiComm & RestrictedComm(){return *RestrictedComm_;}    

  //! Notes whether or not the proc is active for the restricted matrix.
  bool RestrictedProcIsActive(){return proc_is_active;}

  //! Returns the input matrix.
  Teuchos::RCP<Epetra_CrsMatrix> InputMatrix(){return input_matrix_;}

  //! Returns the input matrix, restricted to the active procs.
  Teuchos::RCP<Epetra_CrsMatrix> RestrictedMatrix(){return restricted_matrix_;}

  //! Restricts the matrix.
  /*! Restricts the matrix.  If the MPISubComm is not set, restrict sets it to
     include only the active processors and then restricts the matrix.  If
     MPISubComm is set in a fashion compatible with the existing matrix (aka all
     procs with active rows must be in the subcomm), it
     creates a matrix using that MPISubComm

     \return 0 if sucessful, -1 if the input_matrix is deficient, and -2 if the
  MPI_Comm object set with SetMPISubComm is not consistent with the input matrix.
  */
  int restrict_comm(Teuchos::RCP<Epetra_CrsMatrix> input_matrix);

private:
  bool proc_is_active;
  bool subcomm_is_set;
  
  MPI_Comm MPI_SubComm_;
  Epetra_MpiComm *RestrictedComm_;       
  Epetra_Map *ResRowMap_;
  Epetra_Map *ResColMap_;

  Teuchos::RCP<Epetra_CrsMatrix> input_matrix_;
  Teuchos::RCP<Epetra_CrsMatrix> restricted_matrix_;  
};


} // namespace EpetraExt


#endif
#endif /* EPETRAEXT_RESTRICTEDCRSMATRIXWRAPPER_H */
