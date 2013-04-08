/*
//@HEADER
// ************************************************************************
// 
//               Epetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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
// ************************************************************************
//@HEADER
*/

#ifndef EPETRA_DISTRIBUTOR_H
#define EPETRA_DISTRIBUTOR_H

//! Epetra_Distributor:  The Epetra Gather/Scatter Setup Base Class.
/*! The Epetra_Distributor class is an interface that encapsulates the general
  information and services needed for other Epetra classes to perform gather/scatter
  operations on a parallel computer.
  An Epetra_Distributor object is actually produced by calling a method in the Epetra_Comm class.
  
  Epetra_Distributor has default implementations, via Epetra_SerialDistributor and 
  Epetra_MpiDistributor, for both serial execution and MPI
  distributed memory execution.  It is meant to insulate the user from
  the specifics of communication that are not required for normal
  manipulation of linear algebra objects..
*/

#include "Epetra_Object.h"
class Epetra_Distributor {
    
  public:
    //! @name Constructor and Destructor
  //@{ 
  //! Epetra_Distributor clone constructor.
  virtual Epetra_Distributor * Clone() = 0;
  //! Epetra_Distributor Destructor.
  virtual ~Epetra_Distributor(){};

  //! Create and extract the reverse version of the distributor
  /*! This is not a const method since a reverse distributor might need to be created.
   */
  virtual Epetra_Distributor * GetReverseDistributor() = 0;


  //@}

  
  //! @name Gather/Scatter Constructors
  //@{ 
  //! Create Distributor object using list of process IDs to which we export
  /*! Take a list of Process IDs and construct a plan for efficiently scattering to these processes.
      Return the number of IDs being sent to me.
    \param NumExportIDs In
           Number of IDs that need to be sent from this processor.
    \param ExportPIDs In
           List of processors that will get the exported IDs.
    \param Deterministic In
           No op.
    \param NumRemoteIDs Out
           Number of IDs this processor will be receiving.
  */
  virtual int CreateFromSends( const int & NumExportIDs,
                               const int * ExportPIDs,
			       bool Deterministic,
                               int & NumRemoteIDs ) = 0;

  //! Create Distributor object using list of Remote global IDs and corresponding PIDs
  /*! Take a list of global IDs and construct a plan for efficiently scattering to these processes.
      Return the number and list of IDs being sent by me.
    \param NumRemoteIDs In
           Number of IDs this processor will be receiving.
    \param RemoteGIDs In
           List of IDs that this processor wants.
    \param RemotePIDs In
           List of processors that will send the remote IDs.
    \param Deterministic In
           No op.
    \param NumExportIDs Out
           Number of IDs that need to be sent from this processor.
    \param ExportPIDs Out
           List of processors that will get the exported IDs.
  */
  virtual int CreateFromRecvs( const int & NumRemoteIDs,
                               const int * RemoteGIDs,
			       const int * RemotePIDs,
			       bool Deterministic,
			       int & NumExportIDs,
			       int *& ExportGIDs,
			       int *& ExportPIDs) = 0;

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  virtual int CreateFromRecvs( const int & NumRemoteIDs,
                               const long long * RemoteGIDs,
			       const int * RemotePIDs,
			       bool Deterministic,
			       int & NumExportIDs,
			       long long *& ExportGIDs,
			       int *& ExportPIDs) = 0;
#endif
  //@}

  //! @name Execute Gather/Scatter Operations (Constant size objects)
  //@{ 

  //! Execute plan on buffer of export objects in a single step
  virtual int Do( char * export_objs,
                  int obj_size,
                  int & len_import_objs,
                  char *& import_objs) = 0;

  //! Execute reverse of plan on buffer of export objects in a single step
  virtual int DoReverse( char * export_objs,
                         int obj_size,
                         int & len_import_objs,
                         char *& import_objs ) = 0;

  //! Post buffer of export objects (can do other local work before executing Waits)
  virtual int DoPosts( char * export_objs, 
                       int obj_size,
                       int & len_import_objs,
                       char *& import_objs ) = 0;

  //! Wait on a set of posts
  virtual int DoWaits() = 0;

  //! Do reverse post of buffer of export objects (can do other local work before executing Waits)
  virtual int DoReversePosts( char * export_objs,
                              int obj_size,
                              int & len_import_objs,
                              char *& import_objs) = 0;

  //! Wait on a reverse set of posts
  virtual int DoReverseWaits() = 0;
  //@}

  //! @name Execute Gather/Scatter Operations (Non-constant size objects)
  //@{ 

  //! Execute plan on buffer of export objects in a single step (object size may vary)
  virtual int Do( char * export_objs,
                  int obj_size,
                  int *& sizes,
                  int & len_import_objs,
                  char *& import_objs) = 0;

  //! Execute reverse of plan on buffer of export objects in a single step (object size may vary)
  virtual int DoReverse( char * export_objs,
                         int obj_size,
                         int *& sizes,
                         int & len_import_objs,
                         char *& import_objs) = 0;

  //! Post buffer of export objects (can do other local work before executing Waits)
  virtual int DoPosts( char * export_objs,
                       int obj_size,
                       int *& sizes,
                       int & len_import_objs,
                       char *& import_objs) = 0;

  //! Do reverse post of buffer of export objects (can do other local work before executing Waits)
  virtual int DoReversePosts( char * export_objs,
                              int obj_size,
                              int *& sizes,
                              int & len_import_objs,
                              char *& import_objs) = 0;
  //@}

  //! @name Print object to an output stream
  //@{ 
  virtual void Print(ostream & os) const = 0;
  //@}
};
#endif /* EPETRA_DISTRIBUTOR_H */
