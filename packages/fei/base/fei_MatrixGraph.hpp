/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_MatrixGraph_hpp_
#define _fei_MatrixGraph_hpp_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_Reducer.hpp>
#include <snl_fei_Constraint.hpp>
#include <fei_Record.hpp>
#include <fei_SparseRowGraph.hpp>

#include <vector>

namespace fei {
  class ConnectivityBlock;
  class Pattern;
  class SparseRowGraph;
  /** alias for constraint type */
  typedef snl_fei::Constraint<fei::Record<int>*> ConstraintType;

/** A container for the data that defines connectivity, and which will
    ultimately be used to generate a matrix graph.
*/

class MatrixGraph {
 public:
  /** MatrixGraph Factory interface */
  class Factory {
  public:
    /** Usual virtual destructor */
    virtual ~Factory(){}

    /** Produce an instance of a MatrixGraph. Either or both of columnSpace
        and name may be NULL. If columnSpace is NULL, it will be assumed that
        the structure to be created/defined is symmetric. i.e., columnSpace
        will be assumed to be identically equal to rowSpace. */
    virtual fei::SharedPtr<fei::MatrixGraph>
                 createMatrixGraph(fei::SharedPtr<fei::VectorSpace> rowSpace,
                                   fei::SharedPtr<fei::VectorSpace> columnSpace,
                                   const char* name) = 0;
  };

  enum { REDUCED_INDICES   = 0,
         UNREDUCED_INDICES = 1,
         BLOCK_ENTRY_GRAPH = 2,
         POINT_ENTRY_GRAPH = 3};

  /** Destructor. */
  virtual ~MatrixGraph(){}

  /** Set parameters from a ParameterSet object.
      Currently two parameters are recognized:
      "debugOutput 'path'" where 'path' is the path to the location where
      debug-log files will be produced.<br>
      "name 'string'" where 'string' is an identifier that will be used in
      debug-log file-names.
   */
  virtual void setParameters(const fei::ParameterSet& params) = 0;

  /** Provide a VectorSpace to be used for looking up indices, field-masks,
      etc., for the row-space. If no column-VectorSpace is provided, it will
      be assumed that the column-space equals the row-space.

      @return error-code 0 if successful
  */
  virtual void setRowSpace(fei::SharedPtr<fei::VectorSpace> rowSpace) = 0;

  /** Obtain the VectorSpace that corresponds to the row-space for this
      MatrixGraph object.
  */
  virtual fei::SharedPtr<fei::VectorSpace> getRowSpace() = 0;

  /** Provide a VectorSpace to be used for looking up indices, field-masks,
      etc., for the column-space. If no column-VectorSpace is provided, it
      will be assumed that the column-space equals the row-space.

      @return error-code 0 if successful
  */
  virtual void setColumnSpace(fei::SharedPtr<fei::VectorSpace> columnSpace) = 0;

  /** Obtain the VectorSpace that corresponds to the column-space for this
      MatrixGraph object.
  */
  virtual fei::SharedPtr<fei::VectorSpace> getColSpace() = 0;

  /** Define a pattern to use for subsequent blocked-contributions. Examples
      include element-contributions. Return an int patternID that can be used
      to reference this pattern in calls to initConnectivityBlock, etc<br>

      This is the simplest of the pattern-definition methods. IMPORTANT NOTE:
      this method does not associate a field with the identifiers. Only use
      this method for problems where you explicitly don't want or need to
      associate fields with identifiers. Examples would include problems
      where only a single scalar field exists across the entire mesh and thus
      doesn't need to be explicitly referenced. Other cases where this might
      be used is for non finite-element problems that don't have
      identifier/field pairs.

      @param numIDs Input. number of identifiers per pattern 'instance'.
      @param idType Input. Specifies which type of identifiers are associated
      with instances of this pattern. Must be one of the idTypes defined for a
      VectorSpace that is associated with this MatrixGraph. idTypes are
      defined via the method VectorSpace::defineIDTypes().

      @return patternID
    */
  virtual int definePattern(int numIDs,
                     int idType) = 0;

    /** Define a pattern to use for subsequent blocked-contributions. Examples
      include element-contributions. Return an int patternID that can be used
      to reference this pattern in calls to initConnectivityBlock, etc<br>

        This is the simplest of the 3 pattern-definition methods that
        associate fields with identifiers (there is one pattern-definition
        method above that allows for specifying a pattern of identifiers that
        don't have associated fields). This method defines
        patterns for contributions where a single field is associated with each
        identifier in a list of identifiers, and all the identifiers in the list
        are of the same type.<br>

        @param numIDs Input. number of identifiers per pattern 'instance'.
        @param idType Input. Specifies which type of identifiers are associated
        with instances of this pattern. Must be one of the idTypes defined for a
        VectorSpace that is associated with this MatrixGraph. idTypes are
        defined via the method VectorSpace::defineIDTypes().
        @param fieldID Input. field-identifier for the single field that is to
        reside at each identifier.

        @return patternID
    */
  virtual int definePattern(int numIDs,
                    int idType,
                    int fieldID) = 0;

    /** Define a pattern to use for subsequent blocked-contributions. Examples
        include element-contributions. Return an int patternID that can be used
        to reference this pattern in calls to initConnectivityBlock, etc.<br>

        This is the 'middle' of the pattern-definition methods, in terms of
        the complexity of pattern that can be defined. This method
        defines patterns for contributions where the identifiers are all of the
        same type, but an arbitrary list of fields can be associated with each
        identifier. <br>

        @param numIDs Input. number of identifiers per pattern 'instance'.
        @param idType Input. Specifies which type of identifiers are associated
        with instances of this pattern. Must be one of the idTypes defined for
        a VectorSpace that is associated with this MatrixGraph. idTypes are
        defined via the method VectorSpace::defineIDTypes().
        @param numFieldsPerID Input. List of length numIDs. i-th entry ives the
        number of fields to be associated with the i-th identifier in a
        contribution.
        @param fieldIDs Input. Packed list of length sum(numFieldsPerID[i]).
        Contains the fieldIDs to be associated with the identifiers for a
        contribution.
        @return patternID Input. Identifier to be used later when referring to
        this pattern.
    */
   virtual int definePattern(int numIDs,
                     int idType,
                     const int* numFieldsPerID,
                     const int* fieldIDs) = 0;

    /** Define a pattern to use for subsequent blocked-contributions. Examples
        include element-contributions. Return an int patternID that can be used
        to reference this pattern in calls to initConnectivityBlock, etc.<br>

        This is the most general of the pattern-definition methods. This
        method defines a pattern consisting of a mixture of identifier-types,
        with each identifier having an arbitrary list of associated fields.<br>

        @param numIDs Input. number of identifiers per pattern 'instance'.
        @param idTypes Input. List of length numIDs. Specifies the type of each
        identifier to be contributed for instances of this pattern. Each of the
        idTypes must be one of the idTypes defined for a VectorSpace that is
        associated with this MatrixGraph. idTypes are defined via the method
        VectorSpace::defineIDTypes().
        @param numFieldsPerID Input. List of length numIDs. i-th entry gives the
        number of fields to be associated with the i-th identifier in a
        contribution.
        @param fieldIDs Input. Packed list of length sum(numFieldsPerID[i]).
        Contains the fieldIDs to be associated with the identifiers for a
        contribution.
        @return patternID Input. Identifier to be used later when referring to
        this pattern.
    */
   virtual int definePattern(int numIDs,
                     const int* idTypes,
                     const int* numFieldsPerID,
                     const int* fieldIDs) = 0;

    /** Initialize a block of connectivity contributions. An example
        is a block of elements which share a common layout of nodes/fields per
        element.<br>
        This method accepts only one pattern-id, implying that connectivities in
        this block describe a symmetric structure. See the other overloading of
        this method for the non-symmetric case.

        @param blockID Input. User-specified identifier for this block. Will
        generally be required to be non-negative.
        @param numConnectivityLists Input. Number of connectivity-lists that
        will be supplied for this block.
        @param patternID Input. Descriptor for the connectivities to be
        provided. Must be a pattern that was previously defined via
        definePattern().

        @param diagonal Optional argument, defaults to false. If specified as true,
        each connectivity list will only contribute diagonal entries to the graph.
        This is used if the connectivity-block represents a collection of lumped-
        mass submatrix contributions, or something similar.

        @return error-code 0 if successful
    */
   virtual int initConnectivityBlock(int blockID,
                             int numConnectivityLists,
                             int patternID,
                             bool diagonal=false) = 0;

    /** Initialize a block of connectivity contributions. An example
        is a block of elements which share a common layout of nodes/fields per
        element.<br>
        This method accepts only one pattern-id, implying that connectivities in
        this block describe a symmetric structure. See the other overloading of
        this method for the non-symmetric case.

        @param blockID Input. User-specified identifier for this block. Will
        generally be required to be non-negative.
        @param numConnectivityLists Input. Number of connectivity-lists that
        will be supplied for this block.
        @param patternID Input. Descriptor for the connectivities to be
        provided. Must be a pattern that was previously defined via
        definePattern().

        @param diagonal Optional argument, defaults to false. If specified as true,
        each connectivity list will only contribute diagonal entries to the graph.
        This is used if the connectivity-block represents a collection of lumped-
        mass submatrix contributions, or something similar.

        @return identifier for the new connectivity-block.
    */
   virtual int initConnectivityBlock(int numConnectivityLists,
                                     int patternID,
                                     bool diagonal=false) = 0;

    /** Initialize a block of connectivity contributions. An example
        is a block of elements which share a common layout of nodes/fields per
        element.<br>
        This method accepts two pattern-ids, implying that connectivities in
        this block describe a non-symmetric structure. See the other overloading
        of this method for the symmetric case.

        @param blockID Input. User-specified identifier for this block. Will
        generally be required to be non-negative.
        @param numConnectivityLists Input. Number of connectivity-lists that
        will be supplied for this block.
        @param rowPatternID Input. Descriptor for the row-connectivities to be
        provided.
        Must be a pattern that was previously defined via definePattern().
        @param colPatternID Input. Descriptor for the column-connectivities to
        be provided.
        Must be a pattern that was previously defined via definePattern().

        @return error-code 0 if successful
    */
   virtual int initConnectivityBlock(int blockID,
                             int numConnectivityLists,
                             int rowPatternID,
                             int colPatternID) = 0;

    /** Make a contribution to the MatrixGraph's connectivity. Examples would
        include element-node connectivity lists, etc.

        @param blockID Input. Must correspond to a blockID that was previously
        used in a call to initConnectivityBlock().
        @param connectivityID Input. Identifier for this connectivity list.
        May be an element-identifier, etc.
        @param connectedIdentifiers Input. List of the identifiers that form
        this connectivity list.

        @return error-code 0 if successful
     */
   virtual int initConnectivity(int blockID,
                        int connectivityID,
                        const int* connectedIdentifiers) = 0;

    /** Make a contribution to the MatrixGraph's connectivity. This overloading
        of initConnectivity() provides for structurally non-symmetric entries.

        @param blockID Input. Must correspond to a blockID that was previously
        used in a call to initConnectivityBlock().
        @param connectivityID Input. Identifier for this connectivity list.
        May be an element-identifier, etc.
        @param rowConnectedIdentifiers Input. List of the identifiers that form
        the connectivity list for the row-space.
        @param colConnectedIdentifiers Input. List of the identifiers that form
        the connectivity list for the column-space.

        @return error-code 0 if successful
     */
   virtual int initConnectivity(int blockID,
                        int connectivityID,
                        const int* rowConnectedIdentifiers,
                        const int* colConnectedIdentifiers) = 0;

    /** Make a contribution to the MatrixGraph's connectivity. This overloading
        of initConnectivity() assumes structurally symmetric entries.

        @param patternID Input. Must correspond to a Pattern ID that was
        previously used in a call to definePattern().

        @param connectedIdentifiers Input. List of the identifiers that form
        the connectivity list for the row-space (and the column-space, since this
        is a structurally symmetric contribution).

        @return error-code 0 if successful
     */
   virtual int initConnectivity(int patternID,
                        const int* connectedIdentifiers) = 0;

    /** Make a contribution to the MatrixGraph's connectivity. This overloading
        of initConnectivity() provides for structurally non-symmetric entries.

        @param rowPatternID Input. Must correspond to a Pattern ID that was
        previously used in a call to definePattern().
        @param rowConnectedIdentifiers Input. List of the identifiers that form
        the connectivity list for the row-space.
        @param colPatternID Input. Must correspond to a Pattern ID that was
        previously used in a call to definePattern().
        @param colConnectedIdentifiers Input. List of the identifiers that form
        the connectivity list for the column-space.

        @return error-code 0 if successful
     */
   virtual int initConnectivity(int rowPatternID,
                        const int* rowConnectedIdentifiers,
                        int colPatternID,
                        const int* colConnectedIdentifiers) = 0;

    /** Initialize a set of arbitrary positions in the graph by providing
        data in a "raw" or "purely algebraic" format similar to what might be
        used with a standard sparse CSR (compressed sparse row) matrix.

        @param idType identifier-type
        @param numRows Number of rows, length of the following 'rowIDs' list.
        @param rowIDs List of length 'numRows', specifying identifiers in
        the row-space.
        @param rowOffsets List of length numRows+1, giving offsets into the
        'packedColumnIDs' list at which each row begins. i.e., the column IDs
        for rowIDs[i] are packedColumnIDs[rowOffsets[i]...rowOffsets[i+1]-1].

        @param packedColumnIDs Packed list of length rowOffsets[numRows],
        containing the column IDs.
    */
   virtual int initConnectivity(int idType,
                        int numRows,
                        const int* rowIDs,
                        const int* rowOffsets,
                        const int* packedColumnIDs) = 0;

    /** Initialize a set of arbitrary positions in the graph by providing
        data in a "raw" or "purely algebraic" format similar to what might be
        used with a standard sparse CSR (compressed sparse row) matrix. Also
        specify a fieldID to be associated with these graph positions.

        @param idType identifier-type
        @param fieldID field-identifier
        @param numRows Number of rows, length of the following 'rowIDs' list.
        @param rowIDs List of length 'numRows', specifying identifiers in
        the row-space.
        @param rowOffsets List of length numRows+1, giving offsets into the
        'packedColumnIDs' list at which each row begins. i.e., the column IDs
        for rowIDs[i] are packedColumnIDs[rowOffsets[i]...rowOffsets[i+1]-1].

        @param packedColumnIDs Packed list of length rowOffsets[numRows],
        containing the column IDs.
    */
   virtual int initConnectivity(int idType,
                        int fieldID,
                        int numRows,
                        const int* rowIDs,
                        const int* rowOffsets,
                        const int* packedColumnIDs) = 0;

    /** Initialize a set of arbitrary positions in the graph by providing
        data in a "raw" or "purely algebraic" format similar to what might be
        used with a standard sparse CSR (compressed sparse row) matrix.

        @param idType identifier-type
        @param numRows Number of rows, length of the following 'rowIDs' list.
        @param rowIDs List of length 'numRows', specifying identifiers in
        the row-space.
        @param rowLengths List of length numRows, giving the number of column
        IDs for each row ID.

        @param columnIDs C-style table (list of lists) containing the column
        IDs. Number of rows is numRows, length of i-th row is rowLengths[i].
    */
   virtual int initConnectivity(int idType,
                        int numRows,
                        const int* rowIDs,
                        const int* rowLengths,
                        const int*const* columnIDs) = 0;

    /** Initialize a lagrange-multiplier constraint.
    */
   virtual int initLagrangeConstraint(int constraintID,
                              int constraintIDType,
                              int numIDs,
                              const int* idTypes,
                              const int* IDs,
                              const int* fieldIDs) = 0;

    /** Initialize a penalty constraint.
    */
   virtual int initPenaltyConstraint(int constraintID,
                             int constraintIDType,
                             int numIDs,
                             const int* idTypes,
                             const int* IDs,
                             const int* fieldIDs) = 0;

    /** Initialize a slave constraint. (Note to self: document the parameters.)
     */
   virtual int initSlaveConstraint(int numIDs,
                           const int* idTypes,
                           const int* IDs,
                           const int* fieldIDs,
                           int offsetOfSlave,
                           int offsetIntoSlaveField,
                           const double* weights,
                           double rhsValue) = 0;

   virtual bool newSlaveData() = 0;

   /** Query whether a given mesh object has one or more slave DOFs.
   */
   virtual bool hasSlaveDof(int ID, int idType) = 0;

    /** Signal the MatrixGraph object that initialization is complete.
        At this point the MatrixGraph implementation performs internal
        synchronizations etc. This is a collective method.
    */
   virtual int initComplete() = 0;

    /** Generate a sparse row-based graph from structural data that has been
        accumulated. Don't use this until after initComplete() has been called.

        @param locallyOwnedRows Those rows of a matrix that would be owned by the
        local processor.

        @param blockEntryGraph Specifies whether the graph should be constructed
        on a block-entry or point-entry basis. If there is only 1 scalar DOF at
        each mesh-object, then a block-entry graph is the same as a point-entry
        graph.
    */
   virtual fei::SharedPtr<fei::SparseRowGraph>
     createGraph(bool blockEntryGraph,
                 bool localRowGraph_includeSharedRows=false) = 0;

    /** Query whether the specified MatrixGraph is structurally equivalent to
        this MatrixGraph.
    */
   virtual int compareStructure(const fei::MatrixGraph& matrixGraph,
                        bool& equivalent) const = 0;

    /** Query how many connectivity blocks have been initialized. */
   virtual int getNumConnectivityBlocks() const = 0;

    /** Query for the container of connectivity-blocks. */
   virtual std::map<int,fei::ConnectivityBlock*>& getConnectivityBlocks() = 0;
    /** Query for the list of connectivity-block-IDs. */
   virtual int getConnectivityBlockIDs(std::vector<int>& blockIDs) const = 0;

   /** Query how many IDs are in each connectivity list in the specified
       connectivity block. */
   virtual int getNumIDsPerConnectivityList(int blockID) const = 0;

    /** Query how many scatter-indices are associated with each connectivity
        list for a given connectivity-block.
     */
   virtual int getConnectivityNumIndices(int blockID) const = 0;

    /** Query how many scatter-indices are associated with each connectivity
        list for a given connectivity-block,
        in both the row-dimension and the column-dimension.
     */
   virtual int getConnectivityNumIndices(int blockID,
                                 int& numRowIndices,
                                 int& numColIndices) = 0;

    /** Obtain the scatter-indices associated with a connectivity list.
     */
   virtual int getConnectivityIndices(int blockID,
                              int connectivityID,
                              int indicesAllocLen,
                              int* indices,
                              int& numIndices) = 0;

    /** Obtain the scatter-indices for both the row- and column-dimension,
        associated with a connectivity list.
    */
   virtual int getConnectivityIndices(int blockID,
                              int connectivityID,
                              int rowIndicesAllocLen,
                              int* rowIndices,
                              int& numRowIndices,
                              int colIndicesAllocLen,
                              int* colIndices,
                              int& numColIndices) = 0;

   /** Query associated with Pattern rather than connectivity-block.
    */
   virtual int getPatternNumIndices(int patternID,
                            int& numIndices) = 0;

   /** Query associated with Pattern rather than connectivity-block.
    */
   virtual int getPatternIndices(int patternID,
                         const int* IDs,
                         std::vector<int>& indices) = 0;

   /** Query number of local lagrange constraints */
   virtual int getLocalNumLagrangeConstraints() const = 0;

   /** Query number of slave-constraints
    */
   virtual int getGlobalNumSlaveConstraints() const = 0;

   /** Won't typically be of
       interest to application users of fei:: methods.
    */
   virtual ConstraintType* getLagrangeConstraint(int constraintID) = 0;

   /** Won't typically be of
       interest to application users of fei:: methods.
    */
   virtual std::map<int, ConstraintType* >& getLagrangeConstraints() = 0;

   /** Won't typically be of
       interest to application users of fei:: methods.
    */
   virtual ConstraintType* getPenaltyConstraint(int constraintID) = 0;

   /** Won't typically be of
       interest to application users of fei:: methods.
    */
   virtual ConstraintType* getSlaveConstraint(int constraintID) = 0;

   /** Won't typically be of
       interest to application users of fei:: methods.
    */
   virtual int getConstraintConnectivityIndices(ConstraintType* cr,
                                        std::vector<int>& globalIndices) = 0;

   /** Won't typically be of
       interest to application users of fei:: methods.
    */
   virtual const fei::ConnectivityBlock* getConnectivityBlock(int blockID) const = 0;

   /** Won't typically be of
       interest to application users of fei:: methods.
    */
   virtual fei::ConnectivityBlock* getConnectivityBlock(int blockID) = 0;

   /** Utility method. */
   virtual void setIndicesMode(int mode) = 0;

   /** Utility method. */
   virtual fei::SharedPtr<fei::FillableMat> getSlaveDependencyMatrix() = 0;

   /** Retrieve pointer to specified Pattern object.
       If specified pattern is not found, return NULL.
    */
   virtual fei::Pattern* getPattern(int patternID) = 0;

   /** power-users only */
   virtual int createSlaveMatrices() = 0;

   /** Query for the object that manages the equation-space reductions
     associated with removing slave constraints from the linear system.
   */
   virtual fei::SharedPtr<fei::Reducer> getReducer() = 0;

   /** query for shared-but-not-owned graph rows */
   virtual fei::SharedPtr<fei::SparseRowGraph> getRemotelyOwnedGraphRows() = 0;

   /** fill a vector with eqn-numbers of constrained ids */
   virtual void getConstrainedIndices(std::vector<int>& crindices) const = 0;
};//class MatrixGraph
}//namespace fei

#endif

