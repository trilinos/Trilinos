// @HEADER
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
// @HEADER

#ifndef _EPETRAEXT_ZOLTANMPIDISTRIBUTOR_H_
#define _EPETRAEXT_ZOLTANMPIDISTRIBUTOR_H_

#include "Epetra_Object.h"
#include "Epetra_Distributor.h"

#include "zoltan_comm.h"

namespace EpetraExt {

class ZoltanMpiComm;

//! EpetraExt::ZoltanMpiDistributor:  The Zoltanized Epetra MPI implementation of the Epetra_Distributor Gather/Scatter Setup Class.
/*! The Epetra_MpiDistributor class is an MPI implement of Epetra_Distributor
  that encapsulates the general information and services needed for other Epetra
  classes to perform gather/scatter operations on a parallel computer.  An
  Epetra_MpiDistributor object is actually produced by calling a method in the
  Epetra_MpiComm class.
*/

class ZoltanMpiDistributor: public Epetra_Object, public virtual Epetra_Distributor {
    
  public:

  //@{ \name Constructors/Destructor

  //! Default Constructor.
  ZoltanMpiDistributor(const ZoltanMpiComm &Comm);

  //! Copy Constructor.
  ZoltanMpiDistributor(const ZoltanMpiDistributor &Distributor);

  //! Clone method
  Epetra_Distributor *Clone() {
     return dynamic_cast<Epetra_Distributor*>(new ZoltanMpiDistributor(*this));
  };

  //! Epetra_Comm Destructor.
  virtual ~ZoltanMpiDistributor();
  //@}

  
  //@{ \name Gather/Scatter Constructors
  //! Create Distributor object using list of process IDs to which we export
  /*! Take a list of Process IDs and construct a plan for efficiently scattering
      to these processes.  Return the number of IDs being sent to me.
    \param NumExportIDs In
           Number of IDs that need to be sent from this processor.
    \param ExportPIDs In
           List of processors that will get the exported IDs.
    \param Deterministic In
           If set to true, communication will be deterministic (repeatable) from
           call to call.
    \param NumRemoteIDs Out
           Number of IDs this processor will be receiving.
  */
  int CreateFromSends( const int &NumExportIDs,
                       const int *ExportPIDs,
                       bool Deterministic,
                       int &NumRemoteIDs );

  //! Create Distributor object using list of Remote global IDs and corresponding PIDs
  /*! Take a list of global IDs and construct a plan for efficiently scattering
      to these processes.  Return the number and list of IDs being sent by me.
    \param NumRemoteIDs In
           Number of IDs this processor will be receiving.
    \param RemoteGIDs In
           List of IDs that this processor wants.
    \param RemotePIDs In
           List of processors that will send the remote IDs.
    \param Deterministic In
           If set to true, communication will be deterministic (repeatable) from
           call to call.
    \param NumExportIDs Out
           Number of IDs that need to be sent from this processor.
    \param ExportPIDs Out
           List of processors that will get the exported IDs.
  */
  int CreateFromRecvs( const int &NumRemoteIDs,
                       const int *RemoteGIDs,
                       const int *RemotePIDs,
                       bool Deterministic,
		       int &NumExportIDs,
                       int *&ExportGIDs,
		       int *&ExportPIDs);
  //@}

  //@{ \name Execute Gather/Scatter Operations

  //! Execute plan on buffer of export objects in a single step
  int Do( char *export_objs,
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

  //@{ \name Execute Gather/Scatter Operations (Non-constant size objects)

  //! Execute plan on buffer of export objects in a single step (object size may vary)
  int Do( char *export_objs,
          int obj_size,
          int *& sizes,
          int & len_import_objs,
          char *& import_objs );
  
  //! Execute reverse of plan on buffer of export objects in a single step (object size may vary)
  int DoReverse( char *export_objs,
                 int obj_size,
                 int *& sizes,
                 int & len_import_objs,
                 char *& import_objs );
  
  //! Post buffer of export objects (can do other local work before executing Waits)
  int DoPosts( char *export_objs,
               int obj_size,
               int *& sizes,
               int & len_import_objs,
               char *& import_objs );
  
  //! Do reverse post of buffer of export objects (can do other local work before executing Waits)
  int DoReversePosts( char *export_objs,
                      int obj_size,
                      int *& sizes,
                      int & len_import_objs,
                      char *& import_objs );
  
  //@}
  
  //@{ \name Print object to an output stream
  void Print(ostream &os) const;
  //@}
 
  private:

    //! Resize method allow for variable message lengths
    int Resize (int *sizes, int *sum_recv_sizes);

    Zoltan_Comm_Obj *plan_;

    int   tag_;

    const ZoltanMpiComm *epComm_;
    const MPI_Comm comm_;

};

} //namespace EpetraExt

#endif /* _EPETRAEXT_ZOLTANMPIDISTRIBUTOR_H_ */
