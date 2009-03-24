#ifndef _fei_FiniteElementData_hpp_
#define _fei_FiniteElementData_hpp_

#include <fei_Lookup.hpp>

/** This interface is used to pass finite-element data through to a
  solver library, from the FEI implementation. This is primarily intended
  to be an internal interface in the FEI implementation, acting as an
  abstraction between the FEI and underlying solvers. The FEI layer is the
  'data-source', and this FiniteElementData interface is the 'pipe' through
  which data passes to the solver that's assembling an underlying set of
  data-structures.
  This interface is an alternative to the LinearSystemCore and ESI interfaces,
  which are the other major 'pipes' for passing data through to solver
  libraries. The difference is, LinearSystemCore and ESI interfaces accept
  almost entirely algebraic data, e.g., equation-numbered matrix/vector
  contributions, whereas this FiniteElementData interface accepts data in terms
  of nodeNumbers and degree-of-freedom offsets. Equation numbers don't appear
  in this interface.

  Currently (as of August 2001) the only implementation of this interface is
  for the FETI solver. FETI, as well as the FETI implementation of
  FiniteElementData, are written by Kendall Pierson.

  Note: This interface may eventually be either replaced or combined with
  Charles Tong's FiniteElementGrid interface.

  Note2: Three types of "ID"s appear in this interface: elemBlockID, elemID
   and nodeNumbers. elemBlockID's and elemID's are not-necessarily globally
   unique. In certain problems where it matters (e.g., multi-physics, in
   parallel) element-blocks may span processor boundaries, in which case the
   implementer of FiniteElementData may assume that elemBlockID's are globally
   unique. In most cases it won't matter. Elements are never shared between
   processors, and so elemID's don't need to be globally unique. They may be
   assumed to be locally zero-based numbers. nodeNumbers are globally unique,
   and are also globally contiguous and zero-based.
*/

class FiniteElementData {
 public:

  virtual ~FiniteElementData() {};

  /** For setting argc/argv style parameters.
     @param numParams Number of strings in the params argument
     @param params A list of strings which will usually contain space-separated
        key-value pairs. Example: "debugOutput /usr/users/me/work_dir"
  */
  virtual int parameters(int numParams, char** params) = 0;


  /** Supply the FiniteElementData implementation with an object (created and
      owned by the FEI layer) that can be used to obtain various information
      about problem layout, shared finite-element nodes, etc.
      For details, see the documentation for the Lookup interface.
      @param lookup Input. Reference to an implementation of the
      Lookup interface
  */
  virtual int setLookup(Lookup& lookup) = 0;


  /** For describing the general structure of the finite-element problem that
      is to be assembled.
      @param numElemBlocks Number of element-blocks. An element-block is a
      collection of elements having the same number of nodes and the same
      pattern of num-degrees-of-freedom-per-node.
      @param numElemsPerBlock List of length numElemBlocks.
      @param numNodesPerElem List of length numElemBlocks.
      @param elemMatrixSizePerBlock List of length numElemBlocks, i-th entry
        gives the number of scalar degrees-of-freedom in element-matrices for
        the i-th element-block.
      @param totalNumNodes Total number of local nodes.
      @param numSharedNodes Number of nodes that are shared with other
      processors.
      @return error-code 0 if successful
  */
  virtual int describeStructure(int numElemBlocks,
                                const int* numElemsPerBlock,
                                const int* numNodesPerElem,
                                const int* elemMatrixSizePerBlock,
                                int totalNumNodes,
                                int numSharedNodes,
                                int numMultCRs) = 0;

   /** For passing element-connectivity arrays.
       @param elemBlockID Identifier for the element-block that this element
        belongs to.
        @param elemID Locally zero-based identifier for this element.
        @param numNodes Number of nodes for this element.
        @param nodeNumbers List of length numNodes.
        @param numDofPerNode List of length numNodes.
   */
   virtual int setConnectivity(int elemBlockID,
                               int elemID,
                               int numNodes,
                               const int* nodeNumbers,
                               const int* numDofPerNode) = 0;

   /** For passing element-stiffness arrays.
    @param elemBlockID Identifier for the element-block that these elements
        belong to.
    @param elemID Locally zero-based identifier for this element.
    @param numNodes Number of nodes on this element.
    @param nodeNumbers List of length numNodes
    @param dofPerNode List of length numNodes.
    @param coefs C-style table (list of pointers). Each row of the table is of
    length sum(dofPerNode[i]), and that is also the number of rows.
   */
   virtual int setElemMatrix(int elemBlockID,
                             int elemID,
                             int numNodes,
                             const int* nodeNumbers,
                             const int* dofPerNode,
                             const double *const * coefs) = 0;

   /** For passing element-load vectors.
    @param elemBlockID Identifier for the element-block that this element
        belongs to.
    @param elemID Locally zero-based identifier for this element.
    @param numNodes Number of nodes on this element.
    @param nodeNumbers
    @param dofPerNode
    @param coefs Packed list, length sum(dofPerNode[i]).
   */
   virtual int setElemVector(int elemBlockID,
                             int elemID,
                             int numNodes,
                             const int* nodeNumbers,
                             const int* dofPerNode,
                             const double* coefs) = 0;

   virtual int setDirichletBCs(int numBCs,
                               const int* nodeNumbers,
                               const int* dofOffsets,
                               const double* values) = 0;

   virtual int sumIntoMatrix(int numRowNodes,
                             const int* rowNodeNumbers,
                             const int* rowDofOffsets,
                             const int* numColNodesPerRow,
                             const int* colNodeNumbers,
                             const int* colDofOffsets,
                             const double* coefs) = 0;

   virtual int sumIntoRHSVector(int numNodes,
                             const int* nodeNumbers,
                             const int* dofOffsets,
                             const double* coefs) = 0;

   virtual int putIntoRHSVector(int numNodes,
                             const int* nodeNumbers,
                             const int* dofOffsets,
                             const double* coefs) = 0;

   /** Function called to signal to the FiniteElementData implementation that
       data-loading is complete and any synchronization or other final
       operations may be performed now. This is a collective function, must be
       called by all processors.
       @return error-code 0 if successful
   */
   virtual int loadComplete() = 0;


   /** Function called to request the launching of the linear solver.
    @param solveStatus Output, should indicate the status of the solve. A
    successful solve is usually indicated by a value of 0.
    @param iterations Output, how many iterations were performed.
    @return error-code, 0 if convergence tolerance was achieved within the
     specified maximum number of iterations. If error return is non-zero, the
    calling application will be expected to check solveStatus, and consult the
    solver-library's documentation to figure out exactly what happened.
   */
   virtual int launchSolver(int& solveStatus, int& iterations) = 0;

   /** Function to signal that all coefficient values should be zero'd in
       preparation for a new assemble/solve phase. This is to handle cases
       where no structural information is changing from one time-step to the
       next, but new coefficient values need to be loaded for the next solve.
   */
   virtual int reset() = 0;

   /** Function to signal that lagrange multiplier constraints should be
       deleted.
   */
   virtual int deleteConstraints() = 0;

   virtual int getSolnEntry(int nodeNumber,
                            int dofOffset,
                            double& value) = 0;

   /** Function for requesting the solution value of a Lagrange Multiplier.
    */
   virtual int getMultiplierSoln(int CRID, double& lagrangeMultiplier) = 0;

   /** Pass nodal data that probably doesn't mean anything to the FEI
     implementation, but may mean something to the linear solver. Examples:
     geometric coordinates, nullspace data, etc.
    @param fieldID Identifier for the field that describes this data. Lists of
       field identifiers and field sizes defined for the finite-element problem
       may be obtained from the Lookup interface that is supplied by the
       FEI implementation.
    @param nodeNumbers List of nodes for which data is being supplied.
    @param numNodes
    @param data List of length numNodes * (size of field 'fieldID')
   */
   virtual int putNodalFieldData(int fieldID,
                                 int fieldSize,
                                 int numNodes,
                                 const int* nodeNumbers,
                                 const double* coefs) = 0;

   virtual int setMultiplierCR(int CRID,
                               int numNodes,
                               const int* nodeNumbers,
                               const int* dofOffsets,
                               const double* coefWeights,
                               double rhsValue) = 0;

   virtual int setPenaltyCR(int CRID,
                            int numNodes,
                            const int* nodeNumbers,
                            const int* dofOffsets,
                            const double* coefWeights,
                            double penaltyValue,
                            double rhsValue) = 0;
};

#endif

