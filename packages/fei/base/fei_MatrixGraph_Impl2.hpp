/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_MatrixGraph_Impl2_hpp_
#define _fei_MatrixGraph_Impl2_hpp_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_Reducer.hpp>
#include <fei_Graph.hpp>
#include <snl_fei_Constraint.hpp>
#include <fei_Record.hpp>
#include <fei_Logger.hpp>
#include <fei_SparseRowGraph.hpp>
#include <fei_MatrixGraph.hpp>

#include <vector>
#include <map>

namespace fei {
  class ConnectivityBlock;
  class Pattern;

/** A container for the data that defines connectivity, and which will
    ultimately be used to generate a matrix graph.
*/
class MatrixGraph_Impl2 : public fei::MatrixGraph, private fei::Logger {
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
                                   const char* name);
  };

  /** Constructor.
      @param rowSpace
      @param colSpace
      @param name
  */
  MatrixGraph_Impl2(fei::SharedPtr<fei::VectorSpace> rowSpace,
              fei::SharedPtr<fei::VectorSpace> colSpace,
              const char* name = NULL);

  /** Destructor. */
  virtual ~MatrixGraph_Impl2();

  /** Set parameters from a ParameterSet object.
      Currently two parameters are recognized:
      "debugOutput 'path'" where 'path' is the path to the location where
      debug-log files will be produced.<br>
      "name 'string'" where 'string' is an identifier that will be used in
      debug-log file-names.
   */
  void setParameters(const fei::ParameterSet& params);

  /** Provide a VectorSpace to be used for looking up indices, field-masks,
      etc., for the row-space. If no column-VectorSpace is provided, it will
      be assumed that the column-space equals the row-space.

      @return error-code 0 if successful
  */
  void setRowSpace(fei::SharedPtr<fei::VectorSpace> rowSpace);

  /** Obtain the VectorSpace that corresponds to the row-space for this
      MatrixGraph object.
  */
  fei::SharedPtr<fei::VectorSpace> getRowSpace();

  /** Provide a VectorSpace to be used for looking up indices, field-masks,
      etc., for the column-space. If no column-VectorSpace is provided, it
      will be assumed that the column-space equals the row-space.

      @return error-code 0 if successful
  */
  void setColumnSpace(fei::SharedPtr<fei::VectorSpace> columnSpace);

  /** Obtain the VectorSpace that corresponds to the column-space for this
      MatrixGraph object.
  */
  fei::SharedPtr<fei::VectorSpace> getColSpace();

  /** Define a pattern to use for subsequent blocked-contributions. Examples
      include element-contributions. returns patternID.<br>

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
  int definePattern(int numIDs,
                    int idType);

    /** Define a pattern to use for subsequent blocked-contributions. Examples
        include element-contributions. returns patternID.<br>

        This is the simplest of the 3 pattern-definition methods that
        associate fields with identifiers (there is one pattern-definition
        method above that allows for specifying a pattern of identifiers that
        don't have associated fields). This method defines
        patterns for contributions where a single field is associated with each
        identifier in a list of identifiers, and all the identifiers in the list
        are of the same type.<br>

        @param numIDs Input. number of identifiers per pattern 'instance'.
        @param idType Input. Specifies which type of identifiers are associated with instances of this pattern. Must be one of the idTypes defined for a
        VectorSpace that is associated with this MatrixGraph. idTypes are
        defined via the method VectorSpace::defineIDTypes().
        @param fieldID Input. field-identifier for the single field that is to
        reside at each identifier.
        @return patternID Identifier to be used later when referring to
        this pattern.
    */
  int definePattern(int numIDs,
                    int idType,
                    int fieldID);

    /** Define a pattern to use for subsequent blocked-contributions. Examples
        include element-contributions. returns patternID<br>

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
        @return patternID Identifier to be used later when referring to
        this pattern.
    */
   int definePattern(int numIDs,
                     int idType,
                     const int* numFieldsPerID,
                     const int* fieldIDs);

    /** Define a pattern to use for subsequent blocked-contributions. Examples
        include element-contributions.<br>

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
        @return patternID Identifier to be used later when referring to
        this pattern.
    */
   int definePattern(int numIDs,
                     const int* idTypes,
                     const int* numFieldsPerID,
                     const int* fieldIDs);

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
   int initConnectivityBlock(int blockID,
                             int numConnectivityLists,
                             int patternID,
                             bool diagonal=false);

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
   int initConnectivityBlock(int numConnectivityLists,
                             int patternID,
                             bool diagonal=false);

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
   int initConnectivityBlock(int blockID,
                             int numConnectivityLists,
                             int rowPatternID,
                             int colPatternID);

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
   int initConnectivity(int blockID,
                        int connectivityID,
                        const int* connectedIdentifiers);

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
   int initConnectivity(int blockID,
                        int connectivityID,
                        const int* rowConnectedIdentifiers,
                        const int* colConnectedIdentifiers);

    /** Make a contribution to the MatrixGraph's connectivity. This overloading
        of initConnectivity() assumes structurally symmetric entries.

        @param patternID Input. Must correspond to a Pattern ID that was
        previously used in a call to definePattern().

        @param connectedIdentifiers Input. List of the identifiers that form
        the connectivity list for the row-space.

        @return error-code 0 if successful
     */
   int initConnectivity(int patternID,
                        const int* connectedIdentifiers);

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
   int initConnectivity(int rowPatternID,
                        const int* rowConnectedIdentifiers,
                        int colPatternID,
                        const int* colConnectedIdentifiers);

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
   int initConnectivity(int idType,
                        int numRows,
                        const int* rowIDs,
                        const int* rowOffsets,
                        const int* packedColumnIDs);

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
   int initConnectivity(int idType,
                        int fieldID,
                        int numRows,
                        const int* rowIDs,
                        const int* rowOffsets,
                        const int* packedColumnIDs);

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
   int initConnectivity(int idType,
                        int numRows,
                        const int* rowIDs,
                        const int* rowLengths,
                        const int*const* columnIDs);

    /** Initialize a lagrange-multiplier constraint.
    */
   int initLagrangeConstraint(int constraintID,
                              int constraintIDType,
                              int numIDs,
                              const int* idTypes,
                              const int* IDs,
                              const int* fieldIDs);

    /** Initialize a penalty constraint.
    */
   int initPenaltyConstraint(int constraintID,
                             int constraintIDType,
                             int numIDs,
                             const int* idTypes,
                             const int* IDs,
                             const int* fieldIDs);

    /** Initialize a slave constraint. (Note to self: document the parameters.)
     */
   int initSlaveConstraint(int numIDs,
                           const int* idTypes,
                           const int* IDs,
                           const int* fieldIDs,
                           int offsetOfSlave,
                           int offsetIntoSlaveField,
                           const double* weights,
                           double rhsValue);

   bool newSlaveData();

   /** Query whether a given mesh object has one or more slave DOFs.
   */
   bool hasSlaveDof(int ID, int idType);

    /** Signal the MatrixGraph object that initialization is complete.
        At this point the MatrixGraph implementation performs internal
        synchronizations etc. This is a collective method.
    */
   int initComplete();

    /** Generate a sparse row-based graph from structural data that has been
        accumulated. Don't use this until after initComplete() has been called.

        @param locallyOwnedRows Those rows that are owned by the
        local processor.

        @param blockEntryGraph Specifies whether the graph should be constructed
        on a block-entry or point-entry basis. If there is only 1 scalar DOF at
        each mesh-object, then a block-entry graph is the same as a point-entry
        graph.
    */
   fei::SharedPtr<fei::SparseRowGraph>
     createGraph(bool blockEntryGraph,
                 bool localRowGraph_includeSharedRows=false);

    /** Query whether the specified MatrixGraph is structurally equivalent to
        this MatrixGraph.
    */
   int compareStructure(const fei::MatrixGraph& matrixGraph,
                        bool& equivalent) const;

    /** Query how many connectivity blocks have been initialized. */
   int getNumConnectivityBlocks() const;

    /** Query for the container of connectivity-blocks. */
   std::map<int,fei::ConnectivityBlock*>& getConnectivityBlocks();

    /** Query for the list of connectivity-block-IDs. */
   int getConnectivityBlockIDs(std::vector<int>& blockIDs) const;

   /** Query how many IDs are in each connectivity list in the specified
       connectivity block. */
   int getNumIDsPerConnectivityList(int blockID) const;

    /** Query how many scatter-indices are associated with each connectivity
        list for a given connectivity-block.
     */
   int getConnectivityNumIndices(int blockID) const;

    /** Query how many scatter-indices are associated with each connectivity
        list for a given connectivity-block,
        in both the row-dimension and the column-dimension.
     */
   int getConnectivityNumIndices(int blockID,
                                 int& numRowIndices,
                                 int& numColIndices);

    /** Obtain the scatter-indices associated with a connectivity list.
     */
   int getConnectivityIndices(int blockID,
                              int connectivityID,
                              int indicesAllocLen,
                              int* indices,
                              int& numIndices);

    /** Obtain the scatter-indices for both the row- and column-dimension,
        associated with a connectivity list.
    */
   int getConnectivityIndices(int blockID,
                              int connectivityID,
                              int rowIndicesAllocLen,
                              int* rowIndices,
                              int& numRowIndices,
                              int colIndicesAllocLen,
                              int* colIndices,
                              int& numColIndices);

   /** Query associated with Pattern rather than connectivity-block.
    */
   int getPatternNumIndices(int patternID,
                            int& numIndices);

   /** Given a Pattern and list of IDs, fill output vector with associated indices.
    */
   int getPatternIndices(int patternID,
                         const int* IDs,
                         std::vector<int>& indices);

   /** Query number of local lagrange constraints */
   int getLocalNumLagrangeConstraints() const;

   /** Query number of slave-constraints
    */
   int getGlobalNumSlaveConstraints() const;

   /** Won't typically be of
       interest to application users of fei:: methods.
    */
   ConstraintType* getLagrangeConstraint(int constraintID);

   /** Won't typically be of
       interest to application users of fei:: methods.
    */
   std::map<int, ConstraintType* >& getLagrangeConstraints();

   /** Won't typically be of
       interest to application users of fei:: methods.
    */
   ConstraintType* getPenaltyConstraint(int constraintID);

   /** Won't typically be of
       interest to application users of fei:: methods.
    */
   ConstraintType* getSlaveConstraint(int constraintID);

   /** Won't typically be of
       interest to application users of fei:: methods.
    */
   int getConstraintConnectivityIndices(ConstraintType* cr,
                                        std::vector<int>& globalIndices);

   /** Won't typically be of
       interest to application users of fei:: methods.
    */
   const fei::ConnectivityBlock* getConnectivityBlock(int blockID) const;

   /** Won't typically be of
       interest to application users of fei:: methods.
    */
   fei::ConnectivityBlock* getConnectivityBlock(int blockID);

   /** Utility method. */
   void setIndicesMode(int mode);

   /** Utility method. */
   fei::SharedPtr<FillableMat> getSlaveDependencyMatrix();

   /** Retrieve pointer to specified Pattern object.
       If specified pattern is not found, return NULL.
    */
   fei::Pattern* getPattern(int patternID);

   /** power-users only */
   int createSlaveMatrices();

   /** query for the equation-reduction manager. */
   fei::SharedPtr<fei::Reducer> getReducer();

   /** query for shared-but-not-owned graph rows */
   fei::SharedPtr<fei::SparseRowGraph> getRemotelyOwnedGraphRows();

   /** fill a vector with eqn-numbers of constrained ids */
   void getConstrainedIndices(std::vector<int>& crindices) const;

 private:
   int createAlgebraicGraph(bool blockEntryGraph,
                            fei::Graph* graph,
                            bool gatherFromOverlap);

   int addBlockToGraph_multiField_symmetric(fei::Graph* graph,
                                            fei::ConnectivityBlock* cblock);
   int addBlockToGraph_multiField_nonsymmetric(fei::Graph* graph,
                                            fei::ConnectivityBlock* cblock);
   int addBlockToGraph_singleField_symmetric(fei::Graph* graph,
                                             fei::ConnectivityBlock* cblock);
   int addBlockToGraph_singleField_nonsymmetric(fei::Graph* graph,
                                             fei::ConnectivityBlock* cblock);
   int addBlockToGraph_noField_symmetric(fei::Graph* graph,
                                         fei::ConnectivityBlock* cblock);
   int addBlockToGraph_sparse(fei::Graph* graph,
                              fei::ConnectivityBlock* cblock);

   int addPattern(fei::Pattern* pattern);

   int getConnectivityIndices_multiField(const snl_fei::RecordCollection*const* recordCollections,
                                         int* records, int numRecords,
                                         const int* numFieldsPerID,
                                         const int* fieldIDs,
                                         const int* fieldSizes,
                                         int indicesAllocLen,
                                         int* indices,
                                         int& numIndices);

   int getConnectivityIndices_singleField(const snl_fei::RecordCollection*const* recordCollections,
                                          int* records, int numRecords,
                                          int fieldID, int fieldSize,
                                          int indicesAllocLen,
                                          int* indices,
                                          int& numIndices);

   int getConnectivityIndices_noField(const snl_fei::RecordCollection*const* recordCollections,
                                      int* records,
                                      int numRecords,
                                      int indicesAllocLen,
                                      int* indices,
                                      int& numIndices);

   int getConnectivityRecords(fei::VectorSpace* vecSpace,
                              int idType,
                              int numIDs,
                              const int* IDs,
                              int* records);

   int getConnectivityRecords(fei::VectorSpace* vecSpace,
                              int idType,
                              int fieldID,
                              int numIDs,
                              const int* IDs,
                              int* records);

   int getConnectivityRecords(fei::Pattern* pattern,
                              fei::VectorSpace* solnSpace,
                              const int* connectedIdentifiers,
                              int* recordList);

   int exchangeBlkEqnSizes(fei::Graph* graph);

   int addLagrangeConstraintsToGraph(fei::Graph* graph);

   int addPenaltyConstraintsToGraph(fei::Graph* graph);

   void setName(const char* name);

 private:
   int localProc_, numProcs_;

   MPI_Comm comm_;

   fei::SharedPtr<fei::VectorSpace> rowSpace_;
   fei::SharedPtr<fei::VectorSpace> colSpace_;
   bool haveRowSpace_;
   bool haveColSpace_;
   bool symmetric_;

   fei::SharedPtr<fei::SparseRowGraph> remotelyOwnedGraphRows_;

   bool simpleProblem_;
   bool blockEntryGraph_;

   std::map<int,fei::Pattern*> patterns_;

   std::map<int,fei::ConnectivityBlock*> connectivityBlocks_;
   int arbitraryBlockCounter_;

   std::vector<fei::ConnectivityBlock*> sparseBlocks_;

   std::map<int, ConstraintType* >
     lagrangeConstraints_, penaltyConstraints_, slaveConstraints_;

   bool ptEqualBlk_;

   bool newSlaveData_;

   int localNumSlaves_;
   int globalNumSlaves_;

   fei::SharedPtr<FillableMat> D_;
   fei::SharedPtr<CSVec> g_;
   bool g_nonzero_;

   fei::SharedPtr<fei::Reducer> reducer_;

   std::string name_;
   std::string dbgprefix_;

   std::vector<int> tmpIntArray1_, tmpIntArray2_;

   int* vspcEqnPtr_;

   std::set<int> constrained_indices_;

   bool includeAllSlaveConstraints_;
};//class MatrixGraph_Impl2

 inline fei::SharedPtr<fei::VectorSpace> MatrixGraph_Impl2::getRowSpace()
   {
     return(rowSpace_);
   }

 inline fei::SharedPtr<fei::VectorSpace> MatrixGraph_Impl2::getColSpace()
   {
     return(colSpace_);
   }

 inline std::map<int,fei::ConnectivityBlock*>& MatrixGraph_Impl2::getConnectivityBlocks()
   {
     return(connectivityBlocks_);
   }

 inline int MatrixGraph_Impl2::getGlobalNumSlaveConstraints() const
   {
     return( globalNumSlaves_ );
   }

 inline std::map<int, ConstraintType* >& MatrixGraph_Impl2::getLagrangeConstraints()
   {
     return( lagrangeConstraints_ );
   }
}//namespace fei

#endif

