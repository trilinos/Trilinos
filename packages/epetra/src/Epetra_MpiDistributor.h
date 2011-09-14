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

#ifndef EPETRA_MPIDISTRIBUTOR_H
#define EPETRA_MPIDISTRIBUTOR_H
#include "Epetra_Object.h"
#include "Epetra_Distributor.h"
class Epetra_MpiComm;
#include <mpi.h>

//! Epetra_MpiDistributor:  The Epetra MPI implementation of the Epetra_Distributor Gather/Scatter Setup Class.
/*! The Epetra_MpiDistributor class is an MPI implement of Epetra_Distributor that encapsulates the general
  information and services needed for other Epetra classes to perform gather/scatter
  operations on a parallel computer.
  An Epetra_MpiDistributor object is actually produced by calling a method in the Epetra_MpiComm class.
  
*/

class Epetra_MpiDistributor: public Epetra_Object, public virtual Epetra_Distributor {
    
  public:

    //! @name Constructors/Destructor
  //@{ 

  //! Epetra_Comm Default Constructor.
  Epetra_MpiDistributor(const Epetra_MpiComm & Comm);

  //! Epetra_Comm Copy Constructor.
  Epetra_MpiDistributor(const Epetra_MpiDistributor & Distributor);

  //! Clone method
  Epetra_Distributor * Clone(){return(dynamic_cast<Epetra_Distributor *>(new Epetra_MpiDistributor(*this)));};

  //! Epetra_Comm Destructor.
  virtual ~Epetra_MpiDistributor();
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
           No Op.
    \param NumRemoteIDs Out
           Number of IDs this processor will be receiving.
  */
  int CreateFromSends( const int & NumExportIDs,
                       const int * ExportPIDs,
		       bool Deterministic,
                       int & NumRemoteIDs );

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
           No Op.
    \param NumExportIDs Out
           Number of IDs that need to be sent from this processor.
    \param ExportGIDs Out
           List of processors that will get the exported IDs.
    \param ExportPIDs Out
           List of processors that will get the exported IDs.
  */
  int CreateFromRecvs( const int & NumRemoteIDs,
                       const int * RemoteGIDs,
                       const int * RemotePIDs,
		       bool Deterministic,
                       int & NumExportIDs,
                       int *& ExportGIDs,
                       int *& ExportPIDs);
  //@}

  //! @name Execute Gather/Scatter Operations
  //@{ 

  //! Execute plan on buffer of export objects in a single step
  int Do( char * export_objs,
          int obj_size,
          int & len_import_objs,
          char *& import_objs );

  //! Execute reverse of plan on buffer of export objects in a single step
  int DoReverse( char * export_objs,
                 int obj_size,
                 int & len_import_objs,
                 char *& import_objs );

  //! Post buffer of export objects (can do other local work before executing Waits)
  int DoPosts( char * export_objs,
               int obj_size,
               int & len_import_objs,
               char *& import_objs );
  //! Wait on a set of posts
  int DoWaits();

  //! Do reverse post of buffer of export objects (can do other local work before executing Waits)
  int DoReversePosts( char * export_objs,
                      int obj_size,
                      int & len_import_objs,
                      char *& import_objs );

  //! Wait on a reverse set of posts
  int DoReverseWaits();
  //@}

  //! @name Execute Gather/Scatter Operations (Non-constant size objects)
  //@{ 

  //! Execute plan on buffer of export objects in a single step (object size may vary)
  int Do( char * export_objs,
          int obj_size,
          int *& sizes,
          int & len_import_objs,
          char *& import_objs );
  
  //! Execute reverse of plan on buffer of export objects in a single step (object size may vary)
  int DoReverse( char * export_objs,
                 int obj_size,
                 int *& sizes,
                 int & len_import_objs,
                 char *& import_objs );
  
  //! Post buffer of export objects (can do other local work before executing Waits)
  int DoPosts( char * export_objs,
               int obj_size,
               int *& sizes,
               int & len_import_objs,
               char *& import_objs);
  
  //! Do reverse post of buffer of export objects (can do other local work before executing Waits)
  int DoReversePosts( char * export_objs,
                      int obj_size,
                      int *& sizes,
                      int & len_import_objs,
                      char *& import_objs );
  //@}
  
  //! @name Print object to an output stream
  //@{ 
  void Print(ostream & os) const;
  //@}
  private:

    int ComputeRecvs_( int my_proc,
	               int nprocs );

    int ComputeSends_( int num_imports,
		       const int *& import_ids,
		       const int *& import_procs,
		       int & num_exports,
		       int *& export_ids,
		       int *& export_procs,
		       int my_proc );


    int Resize_(int *sizes);

    int Sort_ints_( int *vals, int *other, int nvals );

  private:
    Epetra_MpiDistributor& operator=(const Epetra_MpiDistributor& src);

    int * lengths_to_;
    int * procs_to_;
    int * indices_to_;
    int   size_indices_to_;

    int * lengths_from_;
    int * procs_from_;
    int * indices_from_;
    int   size_indices_from_;

    bool  resized_;
    int * sizes_;

    int * sizes_to_;
    int * starts_to_;
    int * starts_to_ptr_;
    int * indices_to_ptr_;

    int * sizes_from_;
    int * starts_from_;
    int * starts_from_ptr_;
    int * indices_from_ptr_;

    int   nrecvs_;
    int   nsends_;
    int   nexports_;

    int   self_msg_;

    int   max_send_length_;
    int   total_recv_length_;

    int   tag_;

    const Epetra_MpiComm * epComm_;
    const MPI_Comm      comm_;

    MPI_Request * request_;
    MPI_Status  * status_;

    bool no_delete_;

    char * send_array_;
    int send_array_size_;

    Epetra_MpiDistributor * comm_plan_reverse_;

};
#endif /* EPETRA_MPIDISTRIBUTOR_H */
