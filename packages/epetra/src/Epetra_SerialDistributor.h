
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

#ifndef _EPETRA_SERIALDISTRIBUTOR_H_
#define _EPETRA_SERIALDISTRIBUTOR_H_

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

  //@{ \name Constructor/Destructor

  //! Epetra_Comm Default Constructor.
  Epetra_SerialDistributor(const Epetra_SerialComm & Comm);

  //! Epetra_Comm Copy Constructor.
  Epetra_SerialDistributor(const Epetra_SerialDistributor & Plan);

  //! Clone method
  Epetra_Distributor * Clone(){return(dynamic_cast<Epetra_Distributor *>(new Epetra_SerialDistributor(*this)));};

  //! Epetra_Comm Destructor.
  virtual ~Epetra_SerialDistributor();
  //@}

  
  int CreateFromSends( const int & NumExportIDs,const int * ExportPIDs,
			const bool & Deterministic, int & NumRemoteIDs );
  int CreateFromRecvs( const int & NumRemoteIDs, const int * RemoteGIDs, const int * RemotePIDs,
			const bool & Deterministic,int & NumExportIDs,
			int *& ExportGIDs, int *& ExportPIDs);


  int Do       (char * export_objs, const int & obj_size, char * import_objs);
  int DoReverse(char * export_objs, const int & obj_size, char * import_objs);

  int DoPosts(char * export_objs, const int & obj_size, char * import_objs);
  int DoWaits(char * export_objs, const int & obj_size, char * import_objs);

  int DoReversePosts(char * export_objs, const int & obj_size, char * import_objs);
  int DoReverseWaits(char * export_objs, const int & obj_size, char * import_objs);


  int Do       (char * export_objs, const int * & obj_size, char * import_objs);
  int DoReverse(char * export_objs, const int * & obj_size, char * import_objs);

  int DoPosts(char * export_objs, const int * & obj_size, char * import_objs);
  int DoWaits(char * export_objs, const int * & obj_size, char * import_objs);

  int DoReversePosts(char * export_objs, const int * & obj_size, char * import_objs);
  int DoReverseWaits(char * export_objs, const int * & obj_size, char * import_objs);

  virtual void Print(ostream & os) const;
};
#endif /* _EPETRA_SERIALDISTRIBUTOR_H_ */
