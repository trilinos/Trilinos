#ifndef _PETRA_IMPORT_H_
#define _PETRA_IMPORT_H_

//! Petra_Import: This class builds an import object for efficient importing of off-processor elements.

/*! Petra_Import is used to construct a communication plan that can be called repeatedly by computational
    classes such the Petra matrix, vector and multivector classes to efficiently obtain off-processor
    elements.

    This class currently has one constructor, taking two Petra_Map or Petra_BlockMap objects.  
    The first map specifies the global IDs of elements that we want to import later. The
    secpnd map specifies the global IDs that are owned by the calling processor.  
*/

#include "Petra_Petra.h"
#include "Petra_BlockMap.h"
#include "Petra_Map.h"

#ifdef PETRA_MPI
#include "GSComm_Plan.h"
#include "GSComm_Comm.h"
#endif
class Petra_Import{
    
  public:

  //! Petra_Import constructor
  /*! Builds an import object that will transfer object built with SourceMap to objects built with TargetMap.

    A Petra_Import object categorizes the elements of the target map into three sets as follows:
    <ol>
    <li> All elements in the target map that have the same GID as the corresponding element of the source map, 
         starting with the first 
         element in the target map, going up to the first element that is different from the source map.  The number of
	 these IDs is returned by NumSameIDs().
    <li> All elements that are local to the processor, but are not part of the first set of elements.  These elements
         have GIDs that are owned by the calling processor, but at least the first element of this list is permuted.
	 Even if subsequent elements are not permuted, they are included in this list.  The number of permuted elements
	 is returned by NumPermutedIDs().  The list of elements (local IDs) in the source map that are permuted can be
	 found in the list PermuteToLIDs().  The list of elements (local IDs) in the target map that are the new locations
	 of the source elements can be found in the list PermuteFromLIDs().
    <li> All remaining elements of the target map correspond to global IDs that are owned by remote processors.  The number 
         of these elements is returned by NumRemoteIDs() and the list of these is returned by RemoteLIDs().
    </ol>

Given the above information, the Petra_Import constructor builds a list of elements that must be communicated to other
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
    0  1  2  3  8        2  3  4  5         0  5  6  7  8
\endverbatim

The above list is the elements that will be needed to perform the matrix vector multiplication locally on each processor.
Note that the ordering of the elements on each processor is not unique, but has been chosen for illustration.

With these two maps passed into the Petra_Import constructor, we get the following attribute definitions:

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
PermuteToLIDs   = [0, 1, 2]
PermuteFromLIDs = [1, 2, 3]

NumRemoteIDs    = 1
RemoteLIDs      = [0]

NumExportIDs    = 2
ExportLIDs      = [0, 2]
ExportPIDs      = [0, 2]

NumSend         = 2
NumRecv         = 1

\endverbatim

On PE 2:

\verbatim
NumSameIDs      = 0

NumPermuteIDs   = 3
PermuteToLIDs   = [0, 1, 2]
PermuteFromLIDs = [2, 3, 4]

NumRemoteIDs    = 2
RemoteLIDs      = [0, 1]

NumExportIDs    = 2
ExportLIDs      = [0, 2]
ExportPIDs      = [0, 1]

NumSend         = 2
NumRecv         = 2

\endverbatim


<b> Using Petra_Import Objects </b>

Once a Petra_Import object has been constructed, it can be used by any of the Petra classes that support distributed global
objects, namely Petra_RDP_Vector, Petra_RDP_MultiVector, Petra_CRS_Graph, Petra_RDP_CRS_Matrix and Petra_RDP_VBR_Matrix.  
All of these classes have Import and Export methods that will fill new objects whose distribution is described by
the target map, taking elements from the source object whose distribution is described by the source map.  Details of usage
for each class is given in the appropriate class documentation.

  */ 

  //! Constructs a Petra_Import object from the source and target maps.
  Petra_Import( const Petra_BlockMap & TargetMap, const Petra_BlockMap & SourceMap );
  
  //! Petra_Import copy constructor. 
  Petra_Import(const Petra_Import& Importer);
  
  //! Petra_Import destructor.
  
  virtual ~Petra_Import(void);
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
  const Petra_BlockMap & SourceMap() const {return(SourceMap_);};

  //! Returns the TargetMap used to construct this importer
  const Petra_BlockMap & TargetMap() const {return(TargetMap_);};

#ifdef PETRA_MPI
  GSComm_Plan & GSPlan() const {return(*GSPlan_);};
#endif

 protected:

 friend class Petra_BlockMap;

 friend ostream & operator<<( ostream & os, const Petra_Import & Importer );

 private:

  const Petra_BlockMap & TargetMap_;
  const Petra_BlockMap & SourceMap_;

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

#ifdef PETRA_MPI
  GSComm_Plan * GSPlan_;
#endif
  

};

#endif /* _PETRA_IMPORT_H_ */
