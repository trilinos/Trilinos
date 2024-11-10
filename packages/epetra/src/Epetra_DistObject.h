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

#ifndef EPETRA_DISTOBJECT_H
#define EPETRA_DISTOBJECT_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif


#include "Epetra_ConfigDefs.h"
#include "Epetra_Object.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_BlockMap.h"
class Epetra_Comm;
class Epetra_Import;
class Epetra_Export;
class Epetra_Distributor;
class Epetra_OffsetIndex;

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
class EPETRA_LIB_DLL_EXPORT Epetra_DistObject: public Epetra_Object, public virtual Epetra_SrcDistObject {

  public:
    //! @name Constructors/Destructor
  //@{
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
  Epetra_DistObject(const Epetra_BlockMap& Map, const char* const Label);

  //! Epetra_DistObject copy constructor.

  Epetra_DistObject(const Epetra_DistObject& Source);


  //! Epetra_DistObject destructor.
  virtual ~Epetra_DistObject();
  //@}

  //! @name Import/Export Methods
  //@{

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
  int Import(const Epetra_SrcDistObject& A, const Epetra_Import& Importer, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0);

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
  int Import(const Epetra_SrcDistObject& A, const Epetra_Export& Exporter, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0);

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
  int Export(const Epetra_SrcDistObject& A, const Epetra_Import & Importer, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0);

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
  int Export(const Epetra_SrcDistObject& A, const Epetra_Export& Exporter, Epetra_CombineMode CombineMode, const Epetra_OffsetIndex * Indexor = 0);
  //@}

  //! @name Attribute accessor methods
  //@{
  //! Returns the address of the Epetra_BlockMap for this multi-vector.
  const Epetra_BlockMap& Map() const {return(Map_);};

  //! Returns the address of the Epetra_Comm for this multi-vector.
  const Epetra_Comm& Comm() const {return(*Comm_);};

  //! Returns true if this multi-vector is distributed global, i.e., not local replicated.
  bool DistributedGlobal() const {return(Map_.DistributedGlobal());};
  //@}

  //! @name Miscellaneous
  //@{
  //! Print method
  virtual void Print(std::ostream& os) const;
  //@}

 protected:


  //! @name Internal utilities
  //@{
  //! Perform actual transfer (redistribution) of data across memory images, using Epetra_Distributor object.
  virtual int DoTransfer(const Epetra_SrcDistObject& A,
                         Epetra_CombineMode CombineMode,
                         int NumSameIDs,
                         int NumPermuteIDs,
                         int NumRemoteIDs,
                         int NumExportIDs,
                         int* PermuteToLIDs,
                         int* PermuteFromLIDs,
                         int* RemoteLIDs,
                         int* ExportLIDs,
                         int& LenExports,
                         char*& Exports,
                         int& LenImports,
                         char*& Imports,
                         Epetra_Distributor& Distor,
                         bool DoReverse,
                         const Epetra_OffsetIndex * Indexor );
  //@}

  // These methods must be implemented by derived class

  //! @name Virtual methods to be implemented by derived class
  //@{
  //! Allows the source and target (\e this) objects to be compared for compatibility, return nonzero if not.
  virtual int CheckSizes(const Epetra_SrcDistObject& Source) = 0;
  //! Perform ID copies and permutations that are on processor.
  virtual int CopyAndPermute(const Epetra_SrcDistObject& Source,
                             int NumSameIDs,
                             int NumPermuteIDs,
                             int * PermuteToLIDs,
                             int * PermuteFromLIDs,
                             const Epetra_OffsetIndex * Indexor,
                             Epetra_CombineMode CombineMode = Zero) = 0;

  //! Perform any packing or preparation required for call to DoTransfer().
  virtual int PackAndPrepare(const Epetra_SrcDistObject& Source,
                             int NumExportIDs,
                             int* ExportLIDs,
                             int& LenExports,
                             char*& Exports,
                             int& SizeOfPacket,
                             int* Sizes,
                             bool & VarSizes,
                             Epetra_Distributor& Distor) = 0;

  //! Perform any unpacking and combining after call to DoTransfer().
  virtual int UnpackAndCombine(const Epetra_SrcDistObject& Source,
                               int NumImportIDs,
                               int* ImportLIDs,
                               int LenImports,
                               char* Imports,
                               int& SizeOfPacket,
                               Epetra_Distributor& Distor,
                               Epetra_CombineMode CombineMode,
                               const Epetra_OffsetIndex * Indexor) = 0;

  //@}
  Epetra_BlockMap Map_;
  const Epetra_Comm* Comm_;
  char* Exports_;
  char* Imports_;
  int LenExports_;
  int LenImports_;
  int *Sizes_;

 private:
  Epetra_DistObject& operator=(const Epetra_DistObject& src);

};

#endif /* EPETRA_DISTOBJECT_H */
