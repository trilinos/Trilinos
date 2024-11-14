#ifndef MLAPI_WORKSPACE_H
#define MLAPI_WORKSPACE_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_Workspace.h

\brief Collection of utilities for workspace.

\author Marzio Sala, D-INFK/ETHZ.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_include.h"
#include "Epetra_Comm.h"
#ifdef HAVE_MPI
#include <mpi.h>
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "MLAPI_Error.h"

//! MLAPI: Default namespace for all MLAPI objects and functions.
namespace MLAPI {

/*!
\file MLAPI_Workspace

\brief Basic functions to initialize, use and finalize the MLAPI workspace.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

//! Returns a pointer to the ML_Comm object defined on MPI_COMM_WORLD.
ML_Comm* GetML_Comm();

//! Returns a reference to the Epetra_Comm object defined on MPI_COMM_WORLD.
Epetra_Comm& GetEpetra_Comm();

//! Calls Mpi_Barrier() if MPI is enabled.
void Barrier();

//! Returns the ID of the calling process.
int GetMyPID();

//! Returns the total number of processes in the computation.
int GetNumProcs();

//! Retutns the level of output (always 0 if MyPID() != 0).
int GetPrintLevel();

//! Sets the level of output (from 0 to 10, 0 being verbose).
void SetPrintLevel(int Level);

//! Initialize the MLAPI workspace.
void Init(USR_COMM comm = USR_COMM_WORLD);

//! Destroys the MLAPI workspace.
void Finalize();

std::string GetString(const int& x);

std::string GetString(const double& x);

int GetMatrixType();

} // namespace MLAPI

#endif
