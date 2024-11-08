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

#ifndef EPETRA_IMPORT_H
#define EPETRA_IMPORT_H

#if defined(Epetra_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Epetra package is deprecated"
#endif
#endif



#include "Epetra_Object.h"
#include "Epetra_BlockMap.h"
class Epetra_Distributor;
class Epetra_Export;

//! Epetra_Import: This class builds an import object for efficient importing of off-processor elements.

/*! Epetra_Import is used to construct a communication plan that can be called repeatedly by computational
    classes such the Epetra matrix, vector and multivector classes to efficiently obtain off-processor
    elements.

    This class currently has one constructor, taking two Epetra_Map or Epetra_BlockMap objects.
    The first map specifies the global IDs of elements that we want to import later. The
    second map specifies the global IDs that are owned by the calling processor.
*/

class EPETRA_LIB_DLL_EXPORT Epetra_Import: public Epetra_Object {
  friend class Epetra_Export;
  public:

  //! Constructs a Epetra_Import object from the source and target maps.
  /*! This constructor builds an Epetra_Import object by comparing the GID lists of the source and
      target maps.
      \param TargetMap (In) Map containing the GIDs from which data should be imported to each processor from
             the source map whenever an import operation is performed using this importer.
      \param  SourceMap (In) Map containing the GIDs that should be used for importing data.

      \warning Note that the SourceMap \e must have GIDs uniquely owned, each GID of the source map can occur only once.


  Builds an import object that will transfer objects built with SourceMap to objects built with TargetMap.

    A Epetra_Import object categorizes the elements of the target map into three sets as follows:
    <ol>
    <li> All elements in the target map that have the same GID as the corresponding element of the source map,
         starting with the first
         element in the target map, going up to the first element that is different from the source map.  The number of
	 these IDs is returned by NumSameIDs().
    <li> All elements that are local to the processor, but are not part of the first set of elements.  These elements
         have GIDs that are owned by the calling processor, but at least the first element of this list is permuted.
	 Even if subsequent elements are not permuted, they are included in this list.  The number of permuted elements
	 is returned by NumPermutedIDs().  The list of elements (local IDs) in the source map that are permuted can be
	 found in the list PermuteFromLIDs().  The list of elements (local IDs) in the target map that are the new locations
	 of the source elements can be found in the list PermuteToLIDs().
    <li> All remaining elements of the target map correspond to global IDs that are owned by remote processors.  The number
         of these elements is returned by NumRemoteIDs() and the list of these is returned by RemoteLIDs().
    </ol>

Given the above information, the Epetra_Import constructor builds a list of elements that must be communicated to other
processors as a result of import requests.  The number of exported elements (where multiple sends of the same element
to different processors is counted) is returned by NumExportIDs().  The local IDs to be sent are returned by the list
ExportLIDs().  The processors to which each of the elements will be sent in returned in a list of the same length by
ExportPIDs().

The total number of elements that will be sent by the calling processor is returned by NumSend().  The total number of
elements that will be received is returned by NumRecv().


The following example illustrates the basic concepts.

Assume we have 3 processors and 9 global elements with each processor owning 3 elements as follows
\verbatim
 PE 0 Elements |  PE 1 Elements  |  PE 2 Elements
    0  1  2          3  4  5           6  7  8
\endverbatim

The above layout essentially defines the source map argument of the import object.

This could correspond to a 9 by 9 matrix with the first three rows on PE 0, and so on.  Suppose that this matrix
is periodic tridiagonal having the following sparsity pattern:

\verbatim

PE 0 Rows:

  X  X  0  0  0  0  0  0  X
  X  X  X  0  0  0  0  0  0
  0  X  X  X  0  0  0  0  0

PE 1 Rows:

  0  0  X  X  X  0  0  0  0
  0  0  0  X  X  X  0  0  0
  0  0  0  0  X  X  X  0  0

PE 2 Rows:

  0  0  0  0  0  X  X  X  0
  0  0  0  0  0  0  X  X  X
  X  0  0  0  0  0  0  X  X

\endverbatim

To perform a matrix vector multiplication operation y = A*x (assuming that x has the same distribution as the
rows of the matrix A) each processor will need to import elements of x that
are not local.  To do this, we build a target map on each processor as follows:
\verbatim
    PE 0 Elements    |  PE 1 Elements    |  PE 2 Elements
    0  1  2  3  8       2  3  4  5  6       0  5  6  7  8
\endverbatim

The above list is the elements that will be needed to perform the matrix vector multiplication locally on each processor.
Note that the ordering of the elements on each processor is not unique, but has been chosen for illustration.

With these two maps passed into the Epetra_Import constructor, we get the following attribute definitions:

On PE 0:

\verbatim
NumSameIDs      = 3

NumPermuteIDs   = 0
PermuteToLIDs   = 0
PermuteFromLIDs = 0

NumRemoteIDs    = 2
RemoteLIDs      = [3, 4]

NumExportIDs    = 2
ExportLIDs      = [0, 2]
ExportPIDs      = [1, 2]

NumSend         = 2
NumRecv         = 2

\endverbatim

On PE 1:

\verbatim
NumSameIDs      = 0

NumPermuteIDs   = 3
PermuteFromLIDs = [0, 1, 2]
PermuteToLIDs   = [1, 2, 3]

NumRemoteIDs    = 2
RemoteLIDs      = [0, 4]

NumExportIDs    = 2
ExportLIDs      = [0, 2]
ExportPIDs      = [0, 2]

NumSend         = 2
NumRecv         = 2

\endverbatim

On PE 2:

\verbatim
NumSameIDs      = 0

NumPermuteIDs   = 3
PermuteFromLIDs = [0, 1, 2]
PermuteToLIDs   = [2, 3, 4]

NumRemoteIDs    = 2
RemoteLIDs      = [0, 1]

NumExportIDs    = 2
ExportLIDs      = [0, 2]
ExportPIDs      = [0, 1]

NumSend         = 2
NumRecv         = 2

\endverbatim


<b> Using Epetra_Import Objects </b>

Once a Epetra_Import object has been constructed, it can be used by any of the Epetra classes that support distributed global
objects, namely Epetra_Vector, Epetra_MultiVector, Epetra_CrsGraph, Epetra_CrsMatrix and Epetra_VbrMatrix.
All of these classes have Import and Export methods that will fill new objects whose distribution is described by
the target map, taking elements from the source object whose distribution is described by the source map.  Details of usage
for each class is given in the appropriate class documentation.

Note that the reverse operation, an export, using this importer is also possible and appropriate in some instances.
For example, if we compute y = A^Tx, the transpose matrix-multiplication operation, then we can use the importer we constructed
in the above example to do an export operation to y, adding the contributions that come from multiple processors.

  */

  Epetra_Import( const Epetra_BlockMap & TargetMap, const Epetra_BlockMap & SourceMap );

  //! Expert-only import constructor.
  /*! The additional RemotePIDs argument should be filled with the owning PIDs (from the SourceMap) of the remote GIDs
    in the TargetMap.  The normal Import constructor computes this for you with a call to RemoteIDList.  However in
    some cases (MakeImportExport) we already have this information.

    WARNING: THIS METHOD IS FOR INTERNAL USE ONLY.  USERS SHOULD NOT CALL THIS CONSTRUCTOR */
  Epetra_Import( const Epetra_BlockMap & TargetMap, const Epetra_BlockMap & SourceMap, int NumRemotePIDs,const int * RemotePIDs);

  //! Expert-only import constructor
  /*! The RemotePIDs argument should be filled with the owning PIDs (from the SourceMap) of the remote GIDs
    in the TargetMap.  The normal Import constructor computes this for you with a call to RemoteIDList.  However in
    some cases (MakeImportExport) we already have this information.  We also require information on the Export PIDs/GIDs so
    we can use the Distributor's CreateFromSendsAndReceives method.
    WARNING: THIS METHOD IS FOR INTERNAL USE ONLY.  USERS SHOULD NOT CALL THIS CONSTRUCTOR */
  Epetra_Import( const Epetra_BlockMap & TargetMap, const Epetra_BlockMap & SourceMap, int NumRemotePIDs,const int * RemotePIDs,		
		 const int & NumExportIDs, const int * ExportLIDs,  const int * ExportPIDs);


  //! Epetra_Import copy constructor.
  Epetra_Import(const Epetra_Import& Importer);

  //! Epetra_Import pseudo-copy constructor.  Creates an Epetra_Import in the reverse direction of the Epetra_Export argument.
  Epetra_Import(const Epetra_Export& Exporter);


  //! Epetra_Import destructor.
  virtual ~Epetra_Import(void);

  //! Returns the number of elements that are identical between the source and target maps, up to the first different ID
  int NumSameIDs() const {return(NumSameIDs_);};

  //! Returns the number of elements that are local to the calling processor, but not part of the first NumSameIDs() elements.
  int NumPermuteIDs() const {return(NumPermuteIDs_);};

  //! List of elements in the source map that are permuted.
  int * PermuteFromLIDs () const {return(PermuteFromLIDs_);};
  //! List of elements in the target map that are permuted.
  int * PermuteToLIDs () const {return(PermuteToLIDs_);};

  //! Returns the number of elements that are not on the calling processor.
  int NumRemoteIDs() const {return(NumRemoteIDs_);};

  //! List of elements in the target map that are coming from other processors.
  int * RemoteLIDs() const {return(RemoteLIDs_);};

  //! Returns the number of elements that must be sent by the calling processor to other processors.
  int  NumExportIDs () const {return(NumExportIDs_);};

  //! List of elements that will be sent to other processors.
  int * ExportLIDs () const {return(ExportLIDs_);};

  //! List of processors to which elements will be sent, ExportLIDs() [i] will be sent to processor ExportPIDs() [i].
  int * ExportPIDs () const {return(ExportPIDs_);};

  //! Total number of elements to be sent.
  int NumSend() const {return(NumSend_);};

  //! Total number of elements to be received.
  int NumRecv() const {return(NumRecv_);};

  //! Returns the SourceMap used to construct this importer
  const Epetra_BlockMap & SourceMap() const {return(SourceMap_);};

  //! Returns the TargetMap used to construct this importer
  const Epetra_BlockMap & TargetMap() const {return(TargetMap_);};

  Epetra_Distributor & Distributor() const {return(*Distor_);};

  const Epetra_Distributor * DistributorPtr() const {return(Distor_);}

  //! @name Print object to an output stream
  //@{
  virtual void Print(std::ostream & os) const;
  //@}
 protected:

 friend class Epetra_BlockMap;

 private:
 Epetra_Import& operator=(const Epetra_Import& src)
   {
     (void)src;
     //not currently supported
     bool throw_error = true;
     if (throw_error) {
       throw ReportError("Epetra_Import::operator= not supported.",-1);
     }
     return(*this);
   }

 Epetra_BlockMap TargetMap_;
 Epetra_BlockMap SourceMap_;

  int  NumSameIDs_;
  int  NumPermuteIDs_;
  int * PermuteToLIDs_;
  int * PermuteFromLIDs_;
  int  NumRemoteIDs_;
  int * RemoteLIDs_;

  int  NumExportIDs_;
  int * ExportLIDs_;
  int * ExportPIDs_;

  int NumSend_;
  int NumRecv_;

  Epetra_Distributor * Distor_;


  template<typename int_type>
  void Construct( const Epetra_BlockMap &  targetMap, const Epetra_BlockMap & sourceMap, int NumRemotePIDs=-1, const int * UserRemotePIDs=0);

  template<typename int_type>
    void Construct_Expert( const Epetra_BlockMap & TargetMap, const Epetra_BlockMap & SourceMap, int NumRemotePIDs,const int * RemotePIDs, const int & NumExportIDs, const int * ExportLIDs,  const int * ExportPIDs);

};

#endif /* EPETRA_IMPORT_H */
