
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _EPETRA_MPIDISTRIBUTOR_H_
#define _EPETRA_MPIDISTRIBUTOR_H_
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

  //@{ \name Constructors/Destructor

  //! Epetra_Comm Default Constructor.
  Epetra_MpiDistributor(const Epetra_MpiComm & Comm);

  //! Epetra_Comm Copy Constructor.
  Epetra_MpiDistributor(const Epetra_MpiDistributor & Distributor);

  //! Clone method
  Epetra_Distributor * Clone(){return(dynamic_cast<Epetra_Distributor *>(new Epetra_MpiDistributor(*this)));};

  //! Epetra_Comm Destructor.
  virtual ~Epetra_MpiDistributor();
  //@}

  
  //@{ \name Gather/Scatter Constructors
  //! Create Distributor object using list of process IDs to which we export
  /*! Take a list of Process IDs and construct a plan for efficiently scattering to these processes.
      Return the number of IDs being sent to me.
    \param NumExportIDs In
           Number of IDs that need to be sent from this processor.
    \param ExportPIDs In
           List of processors that will get the exported IDs.
    \param Deterministic In
           If set to true, communication will be deterministic (repeatable) from call to call.
    \param NumRemoteIDs Out
           Number of IDs this processor will be receiving.
  */
  int CreateFromSends( const int & NumExportIDs,const int * ExportPIDs,
			const bool & Deterministic, int & NumRemoteIDs );

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
           If set to true, communication will be deterministic (repeatable) from call to call.
    \param NumExportIDs Out
           Number of IDs that need to be sent from this processor.
    \param ExportPIDs Out
           List of processors that will get the exported IDs.
  */
  int CreateFromRecvs( const int & NumRemoteIDs, const int * RemoteGIDs, const int * RemotePIDs,
			const bool & Deterministic,int & NumExportIDs,
			int *& ExportGIDs, int *& ExportPIDs);
  //@}

  //@{ \name Execute Gather/Scatter Operations

  //! Execute plan on buffer of export objects in a single step
  int Do       (char * export_objs,const int & obj_size, char * import_objs);

  //! Execute reverse of plan on buffer of export objects in a single step
  int DoReverse(char * export_objs,const int & obj_size, char * import_objs);

  //! Post buffer of export objects (can do other local work before executing Waits)
  int DoPosts(char * export_objs,const int & obj_size, char * import_objs);
  //! Wait on a set of posts
  int DoWaits(char * export_objs,const int & obj_size, char * import_objs);

  //! Do reverse post of buffer of export objects (can do other local work before executing Waits)
  int DoReversePosts(char * export_objs,const int & obj_size, char * import_objs);

  //! Wait on a reverse set of posts
  int DoReverseWaits(char * export_objs,const int & obj_size, char * import_objs);
  //@}

  //@{ \name Execute Gather/Scatter Operations (Non-constant size objects: NOT IMPLEMENTED)

  //! Execute plan on buffer of export objects in a single step (object size may vary)
  int Do       (char * export_objs, const int * & obj_size, char * import_objs);
  
  //! Execute reverse of plan on buffer of export objects in a single step (object size may vary)
  int DoReverse(char * export_objs, const int * & obj_size, char * import_objs);
  
  //! Post buffer of export objects (can do other local work before executing Waits)
  int DoPosts(char * export_objs, const int * & obj_size, char * import_objs);
  //! Wait on a set of posts
  int DoWaits(char * export_objs, const int * & obj_size, char * import_objs);
  
  //! Do reverse post of buffer of export objects (can do other local work before executing Waits)
  int DoReversePosts(char * export_objs, const int * & obj_size, char * import_objs);
  
  //! Wait on a reverse set of posts
  int DoReverseWaits(char * export_objs, const int * & obj_size, char * import_objs);
  //@}
  
  //@{ \name Print object to an output stream
  void Print(ostream & os) const;
  //@}
  private:

    int ComputeRecvs( const int & my_proc,
	               const int & nprocs,
	               const bool & Deterministic );

    int ComputeSends( const int & num_imports,
		       const int * import_ids,
		       const int * import_procs,
		       int & num_exports,
		       int *& export_ids,
		       int *& export_procs,
		       const int & my_proc );

  private:

    int * lengths_to_;
    int * procs_to_;
    int * indices_to_;
    int  size_indices_to_;
    int * lengths_from_;
    int * procs_from_;
    int * indices_from_;
    int  size_indices_from_;
    int   nrecvs_;
    int   nsends_;
    int   self_msg_;
    int   max_send_length_;
    int   total_recv_length_;
    int   tag_;
    const Epetra_MpiComm * epComm_;
    const MPI_Comm      comm_;
    MPI_Request * request_;
    MPI_Status  * status_;

    bool no_delete_;

    char * recv_array_;
    char * send_array_;

    Epetra_MpiDistributor * comm_plan_reverse_;

};
#endif /* _EPETRA_MPIDISTRIBUTOR_H_ */
