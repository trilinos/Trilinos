
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

#ifndef _EPETRA_DISTOBJECT_H_
#define _EPETRA_DISTOBJECT_H_
#include "Epetra_Object.h"
#include "Epetra_SrcDistObject.h"
class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_Import;
class Epetra_Export;
class Epetra_Distributor;


//! Epetra_DistObject: A class for constructing and using dense multi-vectors, vectors and matrices in parallel.

/*! The Epetra_DistObject is a base class for all Epetra distributed global objects.  It provides the basic
    mechanisms and interface specifications for importing and exporting operations using Epetra_Import and
    Epetra_Export objects.

<b> Distributed Global vs. Replicated Local.</b>

<ul>
  <li> Distributed Global objects - In most instances, a distributed object will be partitioned
       across multiple memory images associated with multiple processors.  In this case, there is 
       a unique copy of each element and elements are spread across all processors specified by 
       the Epetra_Comm communicator.
  <li> Replicated Local Objects - Some algorithms use objects that are too small to
       be distributed across all processors, the Hessenberg matrix in a GMRES
       computation.  In other cases, such as with block iterative methods,  block dot product 
       functions produce small
       dense matrices that are required by all processors.  Replicated local objectss handle
       these types of situation.
</ul>

*/

//==========================================================================
class Epetra_DistObject: public Epetra_Object, public virtual Epetra_SrcDistObject {

  public:
  //@{ \name Constructors/Destructor.
  //! Basic Epetra_DistObject constuctor.
  /*! Creates a Epetra_DistObject object.  

    \param In 
           Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

	   \warning Note that, because Epetra_LocalMap
	   derives from Epetra_Map and Epetra_Map derives from Epetra_BlockMap, this constructor works
	   for all three types of Epetra map classes.

    \return Pointer to a Epetra_DistObject.

  */
  Epetra_DistObject(const Epetra_BlockMap& Map);

  /*! Creates a Epetra_DistObject object.  

    \param In 
           Map - A Epetra_LocalMap, Epetra_Map or Epetra_BlockMap.

	   \warning Note that, because Epetra_LocalMap
	   derives from Epetra_Map and Epetra_Map derives from Epetra_BlockMap, this constructor works
	   for all three types of Epetra map classes.
    \param In 
           Label - An identifier for this object.  By default, set to the name of the object class.

    \return Pointer to a Epetra_DistObject.

  */
  Epetra_DistObject(const Epetra_BlockMap& Map, const char * const Label);

  //! Epetra_DistObject copy constructor.
  
  Epetra_DistObject(const Epetra_DistObject& Source);
  
  
  //! Epetra_DistObject destructor.  
  virtual ~Epetra_DistObject();
  //@}

  //@{ \name Import/Export Methods.

  //! Imports an Epetra_DistObject using the Epetra_Import object.
  /*!
    \param In
           Source - Distributed object that will be imported into the "\e this" object.
    \param In
           Importer - A Epetra_Import object specifying the communication required.

    \param In
           CombineMode - A Epetra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
  int Import(const Epetra_SrcDistObject& A, const Epetra_Import & Importer, Epetra_CombineMode CombineMode);

  //! Imports an Epetra_DistObject using the Epetra_Export object.
  /*!
    \param In
           Source - Distributed object that will be imported into the "\e this" object.
    \param In
           Exporter - A Epetra_Export object specifying the communication required.

    \param In
           CombineMode - A Epetra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
  int Import(const Epetra_SrcDistObject& A, const Epetra_Export & Exporter, Epetra_CombineMode CombineMode);

  //! Exports an Epetra_DistObject using the Epetra_Import object.
  /*!
    \param In
           Source - Distributed object that will be exported to the "\e this" object.
    \param In
           Importer - A Epetra_Import object specifying the communication required.

    \param In
           CombineMode - A Epetra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
  int Export(const Epetra_SrcDistObject& A, const Epetra_Import & Importer, Epetra_CombineMode CombineMode);

  //! Exports an Epetra_DistObject using the Epetra_Export object.
  /*!
    \param In
           Source - Distributed object that will be exported to the "\e this" multivector.
    \param In
           Exporter - A Epetra_Export object specifying the communication required.

    \param In
           CombineMode - A Epetra_CombineMode enumerated type specifying how results should be combined on the 
	   receiving processor.

    \return Integer error code, set to 0 if successful.
  */
  int Export(const Epetra_SrcDistObject& A, const Epetra_Export & Exporter, Epetra_CombineMode CombineMode);
  //@}
  
  //@{ \name Attribute accessor methods.
  //! Returns the address of the Epetra_BlockMap for this multi-vector.
  const Epetra_BlockMap & Map() const {return(*Map_);};

  //! Returns the address of the Epetra_BlockMap for this multi-vector.
  const Epetra_Comm & Comm() const {return(*Comm_);};

  //! Returns true if this multi-vector is distributed global, i.e., not local replicated.
  bool DistributedGlobal() const {return(DistributedGlobal_);};
  //@}

  //@{ \name Miscellaneous
  //! Print method
  virtual void Print(ostream & os) const;
  //@}

 protected:


  //@{ \name Internal utilities  
  //! Perform actual transfer (redistribution) of data across memory images, using Epetra_Distributor object.
  virtual int DoTransfer(const Epetra_SrcDistObject& A, Epetra_CombineMode CombineMode,
			 int NumSameIDs, int NumPermuteIDs, int NumRemoteIDs, int NumExportIDs, 
			 int *PermuteToLIDs, int *PermuteFromLIDs, int *RemoteLIDs, int * ExportLIDs,
			 int Nsend, int Nrecv,
			 int & LenExports, char * & Exports,
			 int & LenImports, char * & Imports,
			 Epetra_Distributor & Distor, 
			 bool DoReverse);
  //@}

  // These methods must be implemented by derived class

  //@{ \name Virtual methods to be implemented by derived class.
  //! Allows the source and target (\e this) objects to be compared for compatibility, return nonzero if not.
  virtual int CheckSizes(const Epetra_SrcDistObject& Source) = 0;
  //! Perform ID copies and permutations that are on processor.
  virtual int CopyAndPermute(const Epetra_SrcDistObject& Source, int NumSameIDs, 
			     int NumPermuteIDs, int * PermuteToLIDs, int * PermuteFromLIDs) = 0;

  //! Perform any packing or preparation required for call to DoTransfer().
  virtual int PackAndPrepare(const Epetra_SrcDistObject& Source, int NumExportIDs, int * ExportLIDs,
			     int Nsend, int Nrecv,
			     int & LenExports, char * & Exports, int & LenImports, 
			     char * & Imports, 
			     int & SizeOfPacket, Epetra_Distributor & Distor) = 0;
  
  //! Perform any unpacking and combining after call to DoTransfer().
  virtual int UnpackAndCombine(const Epetra_SrcDistObject & Source, 
			       int NumImportIDs, int * ImportLIDs, 
			       char * Imports, int & SizeOfPacket, 
			       Epetra_Distributor & Distor, Epetra_CombineMode CombineMode ) = 0;

  //@}
  const Epetra_BlockMap* Map_;
  bool DistributedGlobal_;
  char * Exports_;
  char * Imports_;
  int LenExports_;
  int LenImports_;
  const Epetra_Comm * Comm_;

};

#endif /* _EPETRA_DISTOBJECT_H_ */
