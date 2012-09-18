/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#ifndef _fei_VectorSpace_hpp_
#define _fei_VectorSpace_hpp_

#include <fei_macros.hpp>
#include <fei_constants.hpp>
#include <fei_fwd.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_Logger.hpp>
#include <fei_utils.hpp>
#include <fei_CommUtils.hpp>
#include <fei_FieldDofMap.hpp>
#include <fei_ctg_set.hpp>
#include <snl_fei_RaggedTable.hpp>

namespace fei {
  class FieldMask;
  class Lookup_Impl;
  class Pattern;
  template<typename GlobalIDType> class Record;
  template<typename GlobalIDType> class Record_Operator;
  template<typename GlobalIDType> class SharedIDs;

  /** Class containing the methods for defining a solution-space (a set of
      degrees-of-freedom) and mapping that space to a globally unique set of
      indices. The set of indices can then be used in defining vectors or in
      defining either the row-space or the column-space of matrices.

      Example: define a displacement field over a set of node-identifiers,
      and map that set to a set of equation-numbers.

      There are multiple ways to use an instance of this interface:

      For generating vectors:
      <ol>
      <li>Define fields and identifier-types
      <li>Initialize active fields over sets of identifiers
      <li>Obtain index offset and range information via the methods
      'getGlobalIndexOffsets()', 'getIndices_Owned()',
      'getIndices_SharedAndOwned()', etc., to use in constructing or
      initializing a vector object.
      </ol>

      For generating matrices:
      <ol>
      <li>Define fields and identifier-types
      <li>Construct an instance of fei::MatrixGraph (using this
      VectorSpace as a contructor argument) and then
      initialize connectivities and other structural attributes on
      the MatrixGraph object.
      <li>Obtain matrix-graph information from the fei::MatrixGraph
      object to use in constructing or initializing a matrix object.
      </ol>
  */
  class VectorSpace : private fei::Logger {
  public:
    /** VectorSpace Factory interface */
    class Factory {
    public:
      /** Destructor */
      virtual ~Factory(){}

     /** Produce an instance of a VectorSpace. name may be NULL. */
     virtual fei::SharedPtr<VectorSpace> createVectorSpace(MPI_Comm,
                                                           const char* name);
    };

    /** Constructor.
        @param comm MPI_Communicator
        @param name String to be used in the name of a debug-log file,
        if any. This is an optional argument, defaults to NULL.
    */
    VectorSpace(MPI_Comm comm, const char* name = NULL);

    /** Destructor */
    virtual ~VectorSpace();

    //@{ \name Setup/initialization

    /** Set parameter values from a fei::ParameterSet object.
       (Need to add documentation here regarding what parameters and values
       are valid and meaningful to this object...)
    */
    void setParameters(const fei::ParameterSet& paramset);

    /** Define fields that will occur in this solution space. <br>
        Example: a temperature field might be defined as fieldID 0, size 1.<br>
        Example: a velocity field might be defined as fieldID 5, size 3.<br>

        @param numFields Input. Length of the fieldIDs and fieldSizes lists.
        @param fieldIDs Input. List of user-supplied field-identifiers.
        Convention: Active solution-space fields should generally be denoted
        by non-negative field-identifiers, while "other" fields (such as
        geometric coordinates) should be denoted by negative field-identifiers.
        @param fieldSizes Input. List of user-specified field-sizes. A
        field-size is the number of scalar components that make up a field.
    */
    void defineFields(int numFields,
                      const int* fieldIDs,
                      const int* fieldSizes,
                      const int* fieldTypes = NULL);

    /** Define identifier-types in this solution space. <br>
        For example, define node-identifiers to be type 0, edge-identifiers to
        be type 1, lagrange-multiplier identifiers to be type 2, etc.<br>
        identifier-types need not be zero-based or contiguous.

        @param numIDTypes Number of distinct identifier-types
        @param idTypes User-supplied list of identifier-types
        @return error-code 0 if successful
    */
    void defineIDTypes(int numIDTypes,
                       const int* idTypes);

    void setIDMap(int idType,
                  const int* localIDs_begin, const int* localIDs_end,
                  const int* globalIDs_begin, const int* globalIDs_end);

    /** Add a set of identifiers to this solution-space. These solution-space
        entries consist of fields residing at identifiers.<br>
        Example: temperature field at a set of finite-element nodes.

        @param fieldID Input. The field-identifier to be added. Must be one of
        the fieldIDs defined previously via 'defineFields()'.
        @param idType Input. The identifier-type over which the active field is
        being initialized. Must be one of the idTypes defined previously via
        'defineIDTypes()'.
        @param numIDs Input. Length of the IDs list.
        @param IDs Input List of identifiers over which 'fieldID' is active.
        @return error-code 0 if successful
     */
    int addDOFs(int fieldID,
                int idType,
                int numIDs,
                const int* IDs);

    /** Add a set of identifiers to the solution-space. These solution-space
        entries consist of identifiers that don't have associated fields.<br>
        Example: Lagrange-multiplier constraint identifiers.<br>
        This method may also be used for initializing a finite-element
        solution-space where the user knows that the entire problem contains only
        one scalar field (e.g., temperature) and so it is sufficient to define
        a solution space on identifiers without associating fields with those
        identifiers. (This will achieve a performance gain for the structure-
        definition, graph-generation and matrix/vector assembly.)

        @param idType Input. The identifier-type over which the solution-
        space is being defined. Must be one of the idTypes defined previously via
        'defineIDTypes()'.
        @param numIDs Input. Number of identifiers being added to the solution-
        space.
        @param IDs Input. List of length numIDs. Identifiers being added to the
        solution-space.
        @return error-code 0 if successful
    */
    int addDOFs(int idType,
                int numIDs,
                const int* IDs);

    /** Specify a set of identifiers that are shared with other processors.
        The shared ids must be identified in a globally symmetric way. i.e., if
        the local processor identifies id x as being shared with processor p,
        then processor p MUST identify id x as being shared with the local
        processor.

        @param numShared Input. Length of the lists sharedIDs and 
                                numSharingProcsPerID.
        @param idType Input. The identifier-type of the ids that are being 
                                identified as shared.
        @param sharedIDs Input. List of shared identifiers.
        @param numSharingProcsPerID Input. List of length numShared, and the
          i-th entry gives the number of processors being identified as sharing
          the i-th sharedID.
        @param sharingProcs Input.
        Packed list of length sum(numSharingProcsPerID), containing the sharing
        processor ranks.
        @return error-code 0 if successful
    */
    int initSharedIDs(int numShared,
                      int idType,
                      const int* sharedIDs,
                      const int* numSharingProcsPerID,
                      const int* sharingProcs);

    /** Specify a set of identifiers that are shared with other processors.
        The shared ids must be identified in a globally symmetric way. i.e., if
        the local processor identifies id x as being shared with processor p,
        then processor p MUST identify id x as being shared with the local
        processor.

        @param numShared Input. Length of the lists sharedIDs and 
                                numSharingProcsPerID.
        @param idType Input. The identifier-type of the ids that are being 
                                identified as shared.
        @param sharedIDs Input. List of shared identifiers.
        @param numSharingProcsPerID Input. List of length numShared, and the
          i-th entry gives the number of processors being identified as sharing
          the i-th sharedID.
        @param sharingProcs Input.
        Table with 'numShared' rows, and each row is of length
        numSharingProcsPerID. This table contains the sharing processor ranks.
        @return error-code 0 if successful
    */
    int initSharedIDs(int numShared,
                      int idType,
                      const int* sharedIDs,
                      const int* numSharingProcsPerID,
                      const int* const* sharingProcs);

    /** Add the contents of another VectorSpace object to this one.
     */
    int addVectorSpace(fei::VectorSpace* inputSpace);

    /** Indicate that initialization is complete. This is a collective function,
        must be called on all processors. At this time ownership of shared IDs
        will be assigned, and the global index space calculated.

        @return error-code 0 if successful
    */
    int initComplete();

    /** Return true if initComplete() has already been called.
    */
    bool initCompleteAlreadyCalled() const { return initCompleteAlreadyCalled_; }

    //@}

    //@{ \name Attribute query methods

    /** Return the MPI_Comm held by this vector-space. When built/run in
        serial mode, MPI_Comm is #defined to be int.
    */
    MPI_Comm getCommunicator() const;

    /** Given a particular degree-of-freedom, request the corresponding global
        index. A particular degree-of-freedom is specified as a component of a
        particular field, residing at a particular location (ID).

        @param idType Input. Identifier-type of the location at which the 
        specified degree-of-freedom resides. Must be one of the identifier-types
        previously defined via a call to 'defineIDTypes()'.

        @param ID Input. Identifier for the location being specified, such as a
        node-identifier, etc.

        @param fieldID Input. Identifier for the field being specified.

        @param fieldOffset Input. In case there is more than one field with the
        specified fieldID residing at the specified ID, this provides an offset
        into those fields. If only one field with specified fieldID, then this
        parameter is 0.

        @param whichComponentOfField Input. Specifies a scalar component within
        the field. If the field only has 1 scalar component, then this parameter
        is 0.

        @param globalIndex Output. This is the global index of the specified
        degree-of-freedom. Not referenced if the specified degree-of-freedom is
        not found.

        @return error-code 0 if successful. If the specified degree-of-freedom
        is not found, -1 is returned.
     */
    int getGlobalIndex(int idType,
                       int ID,
                       int fieldID,
                       int fieldOffset,
                       int whichComponentOfField,
                       int& globalIndex);

    /** Given a particular degree-of-freedom, request the corresponding global
        index. A particular degree-of-freedom is specified as a component of a
        particular field, residing at a particular location (ID).

        @param idType Input. Identifier-type of the location at which the 
        specified degree-of-freedom resides. Must be one of the identifier-types
        previously defined via a call to 'defineIDTypes()'.

        @param ID Input. Identifier for the location being specified, such as a
        node-identifier, etc.

        @param fieldID Input. Identifier for the field being specified.

        @param globalIndex Output. This is the global index of the first component
        of the specified field.
         Not referenced if the specified field is not found on the given ID.

        @return error-code 0 if successful. If the specified degree-of-freedom
        is not found, -1 is returned.
     */
    int getGlobalIndex(int idType,
                       int ID,
                       int fieldID,
                       int& globalIndex);

    /** Given a particular identifier, request the corresponding global
        block-index.

        @param idType Input. Identifier-type of the identifier being queried.

        @param ID Input. Identifier for which a block-index is being requested.

        @param globalBlkIndex Output. This is the global block-index of the
        specified identifier.

        @return error-code 0 if successful. If the specified degree-of-freedom
        is not found, -1 is returned.
     */
    int getGlobalBlkIndex(int idType,
                          int ID,
                          int& globalBlkIndex);

    /** Given a list of IDs, fill an output-list of the global-indices that
        correspond to the first instance of the specified field at each ID.

        @param numIDs Input. Length of the IDs.
        @param IDs Input. User-provided list of identifiers.
        @param idType Input. Type of the IDs for which indices are being
        requested.
        @param fieldID Input. Specified field
        @param globalIndices Output. User-allocated list which, on exit, will
        contain the requested indices. Note that the length of this list is
        assumed to be numIDs*getFieldSize(fieldID).

        @return error-code 0 if successful Note that for any IDs that are not
        found, or IDs which don't have the specified field, the corresponding
        global-index will be -1.
    */
    int getGlobalIndices(int numIDs,
                         const int* IDs,
                         int idType,
                         int fieldID,
                         int* globalIndices);

    /** Given a list of localIDs, fill an output-list of the global-indices that
        correspond to the first instance of the specified field at each localID.

        @param numIDs Input. Length of the IDs.
        @param IDs Input. User-provided list of identifiers.
        @param idType Input. Type of the IDs for which indices are being
        requested.
        @param fieldID Input. Specified field
        @param globalIndices Output. User-allocated list which, on exit, will
        contain the requested indices. Note that the length of this list is
        assumed to be numIDs*getFieldSize(fieldID).

        @return error-code 0 if successful Note that for any IDs that are not
        found, or IDs which don't have the specified field, the corresponding
        global-index will be -1.
    */
    int getGlobalIndicesLocalIDs(int numIDs,
                         const int* localIDs,
                         int idType,
                         int fieldID,
                         int* globalIndices);

    /** Given a list of IDs, fill an output-list of the global-block-indices
        that correspond to each ID.

        @param numIDs Input. Length of the IDs list and of the globalBlkIndices
        list.
        @param IDs Input. User-provided list of identifiers.
        @param idType Input. Type of the IDs for which block-indices are being
        requested.
        @param globalBlkIndices Output. User-allocated list which, on exit,
        will contain the requested indices. Note that the length of this list
        is assumed to be numIDs.

        @return error-code 0 if successful Note that for any IDs that are not
        found, the corresponding global-index will be -1.
    */
    int getGlobalBlkIndices(int numIDs,
                         const int* IDs,
                         int idType,
                         int* globalBlkIndices);

    /** Given a list of IDs, fill an output-list of the global-indices that
        correspond to the first instance of the specified field at each ID.
        Somewhat more general version of the getGlobalIndices() method above.

        @param numIDs Input. Length of the IDs list.
        @param IDs Input. User-provided list of identifiers.
        @param idTypes Input. List of length numIDs, specifying the types of
        the IDs for which indices are being requested.
        @param fieldIDs Input. List of length numIDs, specifying a field at
        each ID.
        @param globalIndices Output. User-allocated list which, on exit, will
        contain the requested indices. Note that the length of this list is
        assumed to be numIDs*getFieldSize(fieldID).

        @return error-code 0 if successful Note that for any IDs that are not
        found, or IDs which don't have the specified field, the corresponding
        global-index will be -1.
    */
    int getGlobalIndices(int numIDs,
                         const int* IDs,
                         const int* idTypes,
                         const int* fieldIDs,                         
                         int* globalIndices);

    /** Given a particular degree-of-freedom, request the corresponding global
        index. In this case, the degree-of-freedom is specified simply by an
        identifier and identifier-type, without specifying a field. This is
        intended to be used for requesting global indices for constraint-
        identifiers or other identifiers which don't have associated fields.
        If the specified identifier actually does have associated fields, then
        the output globalIndex will be the global-index corresponding to the
        first component of the first associated field.

        @param idType Input. Identifier-type of the location at which the 
        specified degree-of-freedom resides. Must be one of the identifier-types
        previously defined via a call to 'defineIDTypes()'.

        @param ID Input. Identifier for the location being specified, such as a
        node-identifier, etc.

        @param globalIndex Output. This is the global index of the specified
        degree-of-freedom. Not referenced if the specified degree-of-freedom is
        not found.

        @return error-code 0 if successful. If the specified degree-of-freedom
        is not found, -1 is returned.
     */
    int getGlobalIndex(int idType,
                       int ID,
                       int& globalIndex);

    /** Given a particular identifier, request the number of scalar degrees-of-
        freedom that are associated with that identifier.
    */
    int getNumDegreesOfFreedom(int idType,
                               int ID);

    /** Query the number of fields defined for this vector-space.
     */
    int getNumFields();

    /** Fill a std::vector with fieldIDs defined for this vector-space.

        @param fieldIDs On exit, contains fieldIDs.
     */
    void getFields(std::vector<int>& fieldIDs);

    /** Given a particular identifier, request the number of fields that are
        associated with that identifier.
    */
    int getNumFields(int idType, int ID);

    /** Given a particular identifier, request the list of fields that are
        associated with that identifier.

        @param idType Identifier-type
        @param ID Specified identifier
        @param lenFieldIDs Input. Length of user-allocated 'fieldIDs' list.
        @param fieldIDs Input. User-allocated list, length must be at least 
        as large as the value produced by getNumFields() for this ID.
        @param numFields Output. Number of fields. If numFields > lenFieldIDs,
        then fieldIDs will contain the first 'lenFieldIDs' field identifiers.
    */
    void getFields(int idType, int ID, std::vector<int>& fieldIDs);

    /** Query for the number of identifier-types defined for this vector-space.
     */
    size_t getNumIDTypes();

    /** Query for the list of identifier-types defined for this vector-space.

        @param idTypes Output, on exit contents will contain id-types that are
         defined for this vector-space.
     */
    void getIDTypes(std::vector<int>& idTypes) const;

    /** Request the global index offsets. Indices are zero-based.

        @param globalOffsets Output. On exit, contains global-offsets.<br>
        globalOffsets[i] is first global offset on processor i,
        for i in 0 .. numPartitions - 1<br>
        globalOffsets[i+1] - globalOffsets[i] is the number of indices on the
        i-th processor
    */
    void getGlobalIndexOffsets(std::vector<int>& globalOffsets) const;

    /** Request the global block-index offsets. Indices are zero-based.

        @param globalBlkOffsets Output. On exit, contains global-block-offsets.<br>
        globalBlkOffsets[i] is first global block-offset on processor i,
        for i in 0 .. numPartitions - 1<br>
        globalBlkOffsets[i+1] - globalBlkOffsets[i] is the number of
        block-indices on the i-th processor
    */
    void getGlobalBlkIndexOffsets(std::vector<int>& globalBlkOffsets) const;

    /** Given a global index in the point-equation space, return the
        owning processor. If the global index is not in the equation space,
        return -1.
    */
    int getOwnerProcPtIndex(int globalIndex);

    /** Given a global index in the block-equation space, return the
        owning processor. If the global index is not in the equation space,
        return -1.
    */
    int getOwnerProcBlkIndex(int globalIndex);

    /** Given an identifier (with identifier-type), return true if it resides
        on the local processor, false if not. This will return true if the 
        identifier is either owned or shared by the local processor.
    */
    bool isLocal(int idType, int ID);

    /** Given an identifier (with identifier-type), return true if it resides
        on the local processor and is locally owned, false if not.
    */
    bool isLocallyOwned(int idType, int ID);

    /** Request the field-size for a specified field-identifier. If the specified
        field-identifier is not found, std::runtime_error is thrown.
        @param fieldID Input. Specified field-identifier
    */
    unsigned getFieldSize(int fieldID);

    /** Query the number of locally owned-or-shared identifiers. */
    int getNumOwnedAndSharedIDs(int idType);

    /** Query the number of locally-owned identifiers. */
    int getNumOwnedIDs(int idType);

    /** Obtain a list of the local identifiers. Note that this includes
     identifiers that are locally shared but not owned. */
    int getOwnedAndSharedIDs(int idtype,
                    int lenList,
                    int* IDs,
                    int& numOwnedAndSharedIDs);

    /** Obtain a list of the locally owned identifiers.
     */
    int getOwnedIDs(int idtype,
                           int lenList,
                           int* IDs,
                           int& numLocalIDs);

    /** Query number of indices on local processor, including ones that are
        locally owned as well as shared-but-not-owned. Only available after
        initComplete has been called. (returns 0 before that)
    */
    int getNumIndices_SharedAndOwned() const;

    /** Obtain list of global indices on local processor, including ones that
        are locally owned as well as shared-but-not-owned. Only available
        after initComplete has been called.

        @param globalIndices On output, will contain all
        indices owned or shared by local processor.
    */
    int getIndices_SharedAndOwned(std::vector<int>& globalIndices) const;

    /** Query number of block indices on local processor, including ones that
        are locally owned as well as shared-but-not-owned. Only available after
        initComplete has been called.
    */
    int getNumBlkIndices_SharedAndOwned(int& numBlkIndices) const;

    /** Obtain list of global block indices on local processor, including ones
        that are locally owned as well as shared-but-not-owned. Only available
        after initComplete has been called.

        @param lenBlkIndices Input. Length of user-allocated 'globalBlkIndices'
        list.
        @param globalBlkIndices User-allocated list. On output, will contain all
        indices owned or shared by local processor.
        @param blkSizes User-allocated list. On output, will contain the number
         of scalars (point-indices) associated with each corresponding
         block-index.
        @param numBlkIndices Output. Number of indices. If 'numBlkIndices' is
        different than 'lenBlkIndices', then globalBlkIndices will contain
        'min(lenBlkIndices, numBlkIndices)' of the local processor's indices.
    */
    int getBlkIndices_SharedAndOwned(int lenBlkIndices,
                                     int* globalBlkIndices,
                                     int* blkSizes,
                                     int& numBlkIndices);

    /** Query number of indices owned by local processor.
    */
    int getNumIndices_Owned() const;

    /** Obtain list of global indices owned by local processor. Only
        available after initComplete has been called.

        @param lenIndices Input. Length of user-allocated 'globalIndices' list.
        @param globalIndices User-allocated list. On output, will contain all
        indices owned by local processor.
        @param numIndices Output. Number of indices. If 'numIndices' is different
        than 'lenIndices', then globalIndices will contain
        'min(lenIndices, numIndices)' of the local processor's indices.
    */
    int getIndices_Owned(std::vector<int>& globalIndices) const;
    int getIndices_Owned(int lenIndices, int* globalIndices, int& numIndices) const;

    /** Query number of block indices owned by local processor.
    */
    int getNumBlkIndices_Owned() const;

    /** Obtain list of global block indices owned by local processor. Only
        available after        initComplete has been called.

        @param lenBlkIndices Input. Length of user-allocated 'globalBlkIndices'
        list.
        @param globalBlkIndices User-allocated list. On output, will contain all
        indices owned by local processor.
        @param blkSizes User-allocated list. On output, will contain the number of
        scalars (point-indices) associated with each corresponding block-index.
        @param numBlkIndices Output. Number of indices. If 'numBlkIndices' is
        different than 'lenBlkIndices', then globalBlkIndices will contain
        'min(lenBlkIndices, numBlkIndices)' of the local processor's indices.
    */
    int getBlkIndices_Owned(int lenBlkIndices,
                            int* globalBlkIndices,
                            int* blkSizes,
                            int& numBlkIndices);

    /** Query the number of shared identifiers of a given id-type. */
    int getNumSharedIDs(int idType, int& numShared);

    /** Query the number of eqn indices across all processors.
     */
    int getGlobalNumIndices() const;

    /** Query the number of block-eqn indices across all processors.
     */
    int getGlobalNumBlkIndices() const;

    /** Intended to be used by other fei classes.
    */
    int getRecordCollection(int idType, snl_fei::RecordCollection*& records);

    /** Intended to be used by other fei classes.
    */
    int getRecordCollection(int idType, const snl_fei::RecordCollection*& records) const;

    /** Intended to be used only by other fei classes.
    */
    std::vector<int>& getEqnNumbers();

    /** Intended to be used only by other fei classes.
    */
    const std::vector<int>& getEqnNumbers() const;

    /** Intended to be used by other implementation classes.
    */
    snl_fei::PointBlockMap* getPointBlockMap();
    const snl_fei::PointBlockMap* getPointBlockMap() const;

    fei::FieldDofMap<int>& getFieldDofMap();

    void getGlobalIndices(const fei::Pattern* pattern,
                          const fei::Record<int>*const* records,
                          std::vector<int>& indices);

    void getGlobalIndicesL(const fei::Pattern* pattern,
                          const int* records,
                          std::vector<int>& indices);

    void getGlobalBlkIndices(const fei::Pattern* pattern,
                             const fei::Record<int>*const* records,
                             std::vector<int>& indices);

    void getGlobalIndices(int numRecords,
                          const fei::Record<int>*const* records,
                          int fieldID,
                          int fieldSize,
                          int indicesAllocLen,
                          int* indices,
                          int& numIndices);

    void getGlobalIndicesL(int numRecords,
                          const snl_fei::RecordCollection*const* recordCollections,
                          const int* records,
                          int fieldID,
                          int fieldSize,
                          int indicesAllocLen,
                          int* indices,
                          int& numIndices);

    void getGlobalIndices(int numRecords,
                          const fei::Record<int>*const* records,
                          const int* numFieldsPerID,
                          const int* fieldIDs,
                          const int* fieldSizes,
                          int indicesAllocLen,
                          int* indices,
                          int& numIndices);

    void getGlobalIndicesL(int numRecords,
                          const snl_fei::RecordCollection*const* recordCollections,
                          const int* records,
                          const int* numFieldsPerID,
                          const int* fieldIDs,
                          const int* fieldSizes,
                          int indicesAllocLen,
                          int* indices,
                          int& numIndices);

    void getGlobalBlkIndices(int numRecords,
                             const fei::Record<int>*const* records,
                             int indicesAllocLen,
                             int* indices,
                             int& numIndices);

    void getGlobalBlkIndicesL(int numRecords,
                             const snl_fei::RecordCollection*const* recordCollections,
                             const int* records,
                             int indicesAllocLen,
                             int* indices,
                             int& numIndices);

    int addDOFs(int fieldID,
                            int idType,
                            int numIDs,
                            const int* IDs,
                            int* records);

    int addDOFs(int idType,
                            int numIDs,
                            const int* IDs,
                            int* records);

    std::vector<fei::FieldMask*> fieldMasks_;

    void getSendProcs(std::vector<int>& sendProcs) const;

    fei::SharedIDs<int>& getSharedIDs(int idType);

  private:
    friend class fei::Lookup_Impl;

  private:
    VectorSpace(const VectorSpace& src);
    VectorSpace& operator=(const VectorSpace& src);

    void compute_shared_ids();

    inline void check_version() { fei::utils::version(); }

    void setOwners_lowestSharing();

    int calculateGlobalIndices();

    void runRecords(fei::Record_Operator<int>& record_op);

    int synchronizeSharedRecords();

    int setLocalEqnNumbers();

    int exchangeGlobalIndices();

    int exchangeFieldInfo(fei::comm_map* ownerPattern,
                          fei::comm_map* sharerPattern,
                          snl_fei::RecordCollection* recordCollection,
                          std::vector<fei::FieldMask*>& fieldMasks);

    void setName(const char* name);

  private:
    MPI_Comm comm_;

    std::vector<int> idTypes_;
    std::map<int,unsigned> fieldDatabase_;
    fei::FieldDofMap<int> fieldDofMap_;
    int maxFieldSize_;
    std::vector<snl_fei::RecordCollection*> recordCollections_;

    std::map<int, fei::SharedIDs<int> > sharedIDTables_;
    std::map<int, fei::comm_map*> ownerPatterns_;
    std::map<int, fei::comm_map*> sharerPatterns_;

    bool sharedRecordsSynchronized_;

    snl_fei::PointBlockMap* ptBlkMap_;

    std::vector<int> globalOffsets_;
    std::vector<int> globalIDOffsets_;

    bool simpleProblem_;

    int firstLocalOffset_, lastLocalOffset_;

    std::vector<int> eqnNumbers_;

    bool newInitData_;
    bool initCompleteAlreadyCalled_;

    std::string name_;
    std::string dbgprefix_;
    bool checkSharedIDs_;
  }; // class fei::VectorSpace

  inline std::vector<int>& VectorSpace::getEqnNumbers()
    {
      return( eqnNumbers_ );
    }

  inline const std::vector<int>& VectorSpace::getEqnNumbers() const
    {
      return( eqnNumbers_ );
    }

  inline snl_fei::PointBlockMap* VectorSpace::getPointBlockMap()
    {
      return( ptBlkMap_ );
    }

  inline const snl_fei::PointBlockMap* VectorSpace::getPointBlockMap() const
    {
      return( ptBlkMap_ );
    }

} // namespace fei

#endif // _fei_VectorSpace_hpp_
