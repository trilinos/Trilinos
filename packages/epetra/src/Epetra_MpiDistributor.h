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

/// \class Epetra_MpiDistributor
/// \brief MPI implementation of Epetra_Distributor.
///
/// This class is an MPI implementation of \c Epetra_Distributor.  It
/// encapsulates the general information and services needed for other
/// Epetra classes to perform gather/scatter operations on a parallel
/// computer.  An Epetra_MpiDistributor instance is actually produced
/// by calling a method in the Epetra_MpiComm class.
class Epetra_MpiDistributor: public Epetra_Object, public virtual Epetra_Distributor {
    
  public:

  //! @name Constructors/Destructor
  //@{ 

  //! Default constructor.
  Epetra_MpiDistributor(const Epetra_MpiComm & Comm);

  //! Copy constructor.
  Epetra_MpiDistributor(const Epetra_MpiDistributor & Distributor);

  //! Clone method
  Epetra_Distributor * Clone(){return(dynamic_cast<Epetra_Distributor *>(new Epetra_MpiDistributor(*this)));};

  //! Destructor (declared virtual for memory safety).
  virtual ~Epetra_MpiDistributor();
  //@}

  
  //! @name Gather/Scatter Constructors
  //@{ 

  /// \brief Create a communication plan from send list.
  ///
  /// Given a list of process IDs to which to send the given number of
  /// data IDs, construct a communication plan for efficiently
  /// scattering data to these processes.
  ///
  /// \return The number of data IDs being sent to me.
  ///
  /// \param NumExportIDs [in] Number of data IDs that need to be sent
  ///   from the calling process.
  /// \param ExportPIDs [in] List of process IDs that will get the
  ///   exported data IDs.
  /// \param Deterministic [in] Currently has no effect.
  /// \param NumRemoteIDs [out] Number of data IDs the calling process
  ///   will be receiving.
  int CreateFromSends( const int & NumExportIDs,
                       const int * ExportPIDs,
		       bool Deterministic,
                       int & NumRemoteIDs );

  /// \brief Create a communication plan from receive list.
  ///
  /// Given a list of remote data IDs and corresponding process IDs
  /// from which to receive data, construct a communication plan for
  /// efficiently scattering data to these processes.
  ///
  /// \return The number and list of data IDs being sent by me.
  ///
  /// \param NumRemoteIDs [in] Number of data IDs the calling process
  ///   will be receiving.
  /// \param RemoteGIDs [in] List of data IDs that the calling process
  ///   wants to receive.
  /// \param RemotePIDs [in] List of IDs of the processes that will
  ///   send the remote data IDs to the calling process.
  /// \param Deterministic [in] Currently has no effect.
  /// \param NumExportIDs [out] Number of data IDs that need to be
  ///   sent from the calling process.
  /// \param ExportGIDs [out] List of data IDs that the calling
  ///   process will send out.
  /// \param ExportPIDs [out] List of IDs of the processes that will
  ///   receive the data IDs sent by the calling process.  
  ///
  /// \note This method allocates the output arrays using \c new.  The
  ///   caller is responsible for deallocating them after use.
  int CreateFromRecvs( const int & NumRemoteIDs,
                       const int * RemoteGIDs,
                       const int * RemotePIDs,
		       bool Deterministic,
                       int & NumExportIDs,
                       int *& ExportGIDs,
                       int *& ExportPIDs);

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  int CreateFromRecvs( const int & NumRemoteIDs,
                       const long long * RemoteGIDs,
                       const int * RemotePIDs,
		       bool Deterministic,
                       int & NumExportIDs,
                       long long *& ExportGIDs,
                       int *& ExportPIDs);
#endif
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
  

  //! @name Attribute Accessor Methods
  //@{ 
  //! The number of procs from which we will receive data
  int NumReceives() const {return nrecvs_;}

  //! The number of procs to which we will send data
  int NumSends() const {return nsends_;}

  //! Maximum number of values that this proc is sending to another single proc.
  int MaxSendLength() const {return max_send_length_;}

  //! Total number of values that this proc is receiving from other procs.
  int TotalReceiveLength() const { return total_recv_length_;}

  //! A list of procs sending values to this proc.
  const int * ProcsFrom() const {return procs_from_;}
  
  //! A list of procs to which this proc is sending values. 
  const int * ProcsTo() const {return procs_to_;}

  //! Number of values we're receiving from each proc. 
  /*! We will receive <tt>LengthsFrom[i]</tt> values from proc <tt>ProcsFrom[i]</tt>. */
  const int * LengthsFrom() const {return lengths_from_;}

  //! Number of values we're sending to each proc. 
  /*! We will send <tt>LengthsTo[i]</tt> values to procs <tt>ProcsTo[i]</tt>. */
  const int * LengthsTo() const {return lengths_to_;}

  //@}

  //! @name Print object to an output stream
  //@{ 
  void Print(ostream & os) const;
  //@}
  private:

    int ComputeRecvs_( int my_proc,
	               int nprocs );

	template<typename id_type>
    int ComputeSends_( int num_imports,
		       const id_type *& import_ids,
		       const int *& import_procs,
		       int & num_exports,
		       id_type *& export_ids,
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
