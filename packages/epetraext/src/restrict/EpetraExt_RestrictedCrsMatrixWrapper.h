//@HEADER
// ***********************************************************************
// 
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2001) Sandia Corporation
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
