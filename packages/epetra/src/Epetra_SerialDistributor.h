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

#ifndef EPETRA_SERIALDISTRIBUTOR_H
#define EPETRA_SERIALDISTRIBUTOR_H

#include "Epetra_Object.h"
#include "Epetra_Distributor.h"
class Epetra_SerialComm;

//! Epetra_SerialDistributor:  The Epetra Serial implementation of the Epetra_Distributor Gather/Scatter Setup Class.
/*! The Epetra_SerialDistributor class is an Serial implement of Epetra_Distributor that is essentially a trivial class
    since a serial machine is a trivial parallel machine.
  An Epetra_SerialDistributor object is actually produced by calling a method in the Epetra_SerialComm class.
  
*/

class Epetra_SerialDistributor: public Epetra_Object, public virtual Epetra_Distributor {
    
  public:

    //! @name Constructor/Destructor
  //@{ 

  //! Constructor.
  Epetra_SerialDistributor(const Epetra_SerialComm & Comm);

  //! Epetra_SerialDistributor Copy Constructor.
  Epetra_SerialDistributor(const Epetra_SerialDistributor & Plan);

  //! Clone method
  Epetra_Distributor * Clone(){return(dynamic_cast<Epetra_Distributor *>(new Epetra_SerialDistributor(*this)));};

  //! Epetra_Comm Destructor.
  virtual ~Epetra_SerialDistributor();

  //! Create and extract the reverse version of the distributor
  /*! \warning This is not implemented for Epetra_SerialDistributor.   
   */
  Epetra_Distributor * GetReverseDistributor() {return 0;}

  //@}

  
  int CreateFromSends( const int & NumExportIDs,
                       const int * ExportPIDs,
		       bool Deterministic,
                       int & NumRemoteIDs );

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

  int Do       (char * export_objs,
                int obj_size,
                int & len_import_objs,
                char *& import_objs);
  int DoReverse(char * export_objs,
                int obj_size,
                int & len_import_objs,
                char *& import_objs);

  int DoPosts(char * export_objs,
              int obj_size,
              int & len_import_objs,
              char *& import_objs);
  int DoWaits();

  int DoReversePosts(char * export_objs,
                     int obj_size,
                     int & len_import_objs,
                     char *& import_objs);
  int DoReverseWaits();


  int Do       (char * export_objs,
                int obj_size,
                int *& sizes,
                int & len_import_objs,
                char *& import_objs);
  int DoReverse(char * export_objs,
                int obj_size,
                int *& sizes,
                int & len_import_objs,
                char *& import_objs);

  int DoPosts(char * export_objs,
              int obj_size,
              int *& sizes,
              int & len_import_objs,
              char *& import_objs);

  int DoReversePosts(char * export_objs,
                     int obj_size,
                     int *& sizes,
                     int & len_import_objs,
                     char *& import_objs);

  virtual void Print(ostream & os) const;

 private:
  int nrecvs_;
  int nsends_;
};
#endif /* EPETRA_SERIALDISTRIBUTOR_H */
