#ifndef _ESI_Broker_h_
#define _ESI_Broker_h_

/**
  This interface "brokers" instances of ESI interfaces.

It is intended to play the role of a "factory", producing interface instances
on request. The interfaces produced are intended to be abstract ESI interfaces,
but of course can be any named interface supported by the particular ESI_Broker
implementation.

The ESI effort hasn't produced an official factory specification, so this
particular interface was concocted by me (Alan Williams) for use in connecting
ESI implementations to the Sandia FEI implementation. The FEI is a linear-system
assembly library, and as such may have relatively specialized needs which may
be reflected in this interface, making it less general than others may need.

So, in the context of the FEI's needs, this interface works as follows.
It primarily manages the linear system which is being assembled by the FEI 
implementation: i.e., 1 to several matrices, 1 to several vector(s), and
the tools to solve a linear system: solver, preconditioner, etc.
The implementation of these interfaces is library-specific -- i.e., concrete
implementations (such as ISIS_ESI or Trilinos classes) lie beneath the ESI
interfaces being managed or brokered by ESI_Broker instances.

FEI-specific backwards compatibility note: This interface is replacing the old
LinearSystemCore interface used by the FEI implementation (that interface was
used to interact with non-ESI solver libraries). But in recognition of the fact
that many solver libraries (and thus ESI_Broker implementations) will not yet
have ESI interface implementations, this interface can also be implemented as
simply a shell for LinearSystemCore. i.e., the function 'hasESIinterfaces()'
may return false, and instead this object may only successfully answer requests
for either a LinearSystemCore interface, or a FiniteElementData interface.
See the already-provided adapter-interface LSC_Broker, which may be used with
any correct implementation of LinearSystemCore or FiniteElementData. i.e., a
user of FEI_Implementation may first instantiate their LinearSystemCore
implementation, use that as a constructor argument in instantiating LSC_Broker,
then pass THAT as the constructor argument to FEI_Implementation.<br>
Example:<p>
LinearSystemCore* isis_lsc = new ISIS_LinSysCore(...);<br>
ESI_Broker* broker = new LSC_Broker(isis_lsc);<br>
FEI* fei = new FEI_Implementation(broker, ...);<br>
<p>

The general ESI approach is to always inherit from ESI_Object, which has
a function called 'getInterface' through which the implementing object
can be asked to supply particular interfaces for itself. e.g., a matrix
object can be asked for its row-oriented write-access interface, or its
Operator interface, etc. ESI_Broker is somewhat different than a normal
esi::Object, in the sense that an ESI_Broker implementation manages/constructs
several objects of different fundamental types. Thus, having this (ESI_Broker)
object provide an implementation of esi::Object's 'getInterface' is ambiguous
because an ESI_Broker instance isn't simply providing interfaces or "views"
of itself, but rather is handing out instances of interfaces to separate
objects. Those interface instances are described by two strings which specify
both the type and an instance-name of the particular interface. The ESI_Broker
implementation keeps a pointer to the requested interfaces, and retains the
responsibility for deleting them when the ESI_Broker is destroyed. The reason
for this is that an esi object produced by this broker may be used by more than
one scope of the calling code. As an example, an application code may
instantiate the ESI_Broker implementation, hand it to the FEI implementation,
and assemble finite-element data into the FEI. That finite-element data ends
up in an ESI matrix. So during the course of the assembly, the FEI
implementation has requested and used matrix interfaces for depositing data
into. Once the assembly is finished, the application code may then request
that same matrix and pass it to an eigensolver, or separate linear
solver, etc. Thus, since neither the FEI implementation nor the application
code can always be sure when the matrix should be destroyed, that
responsibility lies with the ESI_Broker object.

Now a few comments on the usage of this object, for constructing vectors and
matrices. There is an assumed calling sequence which must be followed in order
to create a vector or matrix instance. Vector and matrix instances are requested
via the getInterfaceInstance function, but those requests will not be successful
unless one or two other functions have already been called.<br>
Specifically, requests for vectors can not be filled until after the method
ESI_Broker::setGlobalOffsets has been called to establish an index space that
the vector belongs to. In other words, to establish the size and distribution
of the vector.<br>
Requests for matrices can not be filled until after both
ESI_Broker::setGlobalOffsets and ESI_Broker::setMatrixStructure have been
called, the latter method providing the graph of the matrix, which may or may
not be required when the underlying concrete matrix object is constructed.<p>

Important note: these comments refer to ESI interfaces using the prefix 'ESI_',
which will soon be incorrect. The ESI interface effort is adopting C++ 
namespaces, and that prefix will be replaced with 'esi::'.<br>
e.g., ESI_Vector will be replaced with esi::Vector.

ESI_Broker's primary query function is:

   getInterfaceInstance -- for getting an interface to some particular object
   which resides in this ESI_Broker implementation.
  
  Example calls that may be made to the getInterfaceInstance function:
 
  getInterfaceInstance("x_Vector", "ESI_Vector",...);
  
  getInterfaceInstance("b_Vector_0", "ESI_Vector",...);
  
  getInterfaceInstance("A_Matrix", "ESI_MatrixRowWriteAccess",...);
  
  getInterfaceInstance("A_Matrix", "ESI_MatrixRowReadAccess",...);

  getInterfaceInstance("solver", "ESI_Solver",...);

  getInterfaceInstance("precond", "ESI_Preconditioner",...);

  getInterfaceInstance("FE_Data", "FiniteElementData",...);
 
  whenever the return-value 'err' is non-zero, it indicates that this
  implementation of ESI_Broker is unable to provide the interface that
  was requested.


  ESI_Broker implementers may also expect the FEI implementation to make the
  following call to the setInterfaceInstance function:

  setInterfaceInstance("lookup", "Lookup",...);

*/

class ESI_Broker { // : public virtual ESI_Object {
 public:

   /** Destructor. Destroys all internally generated objects and allocated
       memory. Does NOT destroy any objects set via 'setInterfaceInstance'.
   */
   virtual ~ESI_Broker() {};

   virtual int setMPIComm(MPI_Comm comm) = 0;

   /** Backwards compatibility function. Allows calling code to query
     whether this object is a 'real' ESI interface factory/broker/manager.
     If this function returns false, an FEI implementation will then request,
     via getInterfaceInstance, a LinearSystemCore and/or a FiniteElementData
     interface. One of these should be successful. If neither of those is
     successful, then execution will probably be terminated.

     @return true if ESI interfaces may be requested from this object.
   */
   virtual bool hasESIinterfaces() = 0;


   /** Obtain a name for the underlying library.
      @param name E.g., ISIS-ESI, Trilinos, etc. This string is owned
     (allocated) by this object, the caller should not delete it.
   */
   virtual int getLibraryName(char*& libName) = 0;

   /** Pass parameter strings to this object, which will in turn pass them on to
      any internal objects that take parameter strings. 

      @param numParams Number of parameters being passed.

      @param paramStrings Array of parameters. Each parameter should be a null-
        terminated string. Usually these strings will be key-value pairs,
        separated by a space. Example: "debugOutput /usr/users/me/work_dir"
   */
   virtual int parameters(int numParams, char** paramStrings) = 0;

   /** Create a pointer to a clone (new instantiation) of this object. The
       caller will be responsible for destroying the new instantiation.

       @param esi_broker Output. Pointer to the new object.
    */
   virtual int clone(ESI_Broker*& esi_broker) = 0;

  /** Supply ESI_Broker with global offset information for the problem
      being assembled. This is basically parallel decomposition information for
      three quantities: finite-element-nodes, point-entry equations and
      block-entry equations. Particular implementations may ignore one or more
      of these lists. e.g., a purely algebraic solver library doesn't care
      about finite-element-nodes, and a library that only has point-entry 
      matrices doesn't care about block-entry equation decompositions.
      @param len Input. Length of the following list arguments. This will be
      numProcs+1
      @param nodeOffsets The FEI implementation assigns a global 0-based
          numbering to the finite-element nodes in the problem. Each processor
          is given ownership of a contiguous subset of these node-numbers.
          nodeOffsets[i] gives the first local node-number for processor i.
          nodeOffsets[len-1] gives the total number of nodes.
      @param eqnOffsets eqnOffsets[i] gives the first local equation number for
         processor i, eqnOffsets[len-1] gives the global number of equations.
      @param blkEqnOffsets Contains the same kind of information as eqnOffsets,
          but for 'block-equations'. A block-equation contains all of the
          point-equations present at a finite-element node. Special case: if
          this problem contains Lagrange Multiplier constraints, they will
          be equations that don't correspond to any node, and there will only be
          one of these equations mapped to a block-equation.
  */
   virtual int setGlobalOffsets(int len, int* nodeOffsets,
                        int* eqnOffsets, int* blkEqnOffsets) = 0;

   /** Specify which global equation numbers correspond to Lagrange Multiplier
      equations. This function won't be called if there are no Lagrange
      Multiplier constraints in the problem. If this function is called, it is
      guaranteed to be called after 'setGlobalOffsets' and before
      'setMatrixStructure'. The primary purpose of this function is to give
      ESI_Broker implementers the opportunity to deal with constraints in a
      special way, rather than assembling everything into one matrix. If the
      problem being assembled does have Lagrange constraints, then the FEI
      implementation will request an ESI_MatrixRowWriteAccess interface for
      "C_Matrix" from ESI_Broker. If that is not available, then the FEI
      implementation will request "A_Matrix" and assemble everything into that.
     @param numCRs number of constraint relations
     @param numNodesPerCR number of constrained node in each constraint relation
     @param nodeNumbers Table of constrained nodes. 'numCRs' rows, with the i-th
        row being of length numNodesPerCR[i].
     @param eqnNumbers Table, same dimensions as 'nodeNumbers'. These are the
        global 0-based matrix column indices of the constraint coefficients.
     @param multiplierEqnNumbers Equation numbers that the Lagrange Multipliers
        correspond to.
   */
   virtual int setMultCREqns(int numCRs, 
                             const int* numNodesPerCR,
                             const int* const* nodeNumbers,
                             const int* const* eqnNumbers,
                             const int* fieldIDs,
                             const int* multiplierEqnNumbers) = 0;

   /** Supply ESI_Broker with information defining the structure of the
       sparse matrix to be assembled. Implementers of ESI_Broker may safely
       assume that this function will not be called until after the
       function 'setGlobalOffsets' has been called. Using the information 
       provided via setGlobalOffsets, the number-of-local-equations can be
       trivially calculated. After setMatrixStructure has been called, there
       should be enough information to instantiate internal linear-algebra
       entities, such as vectors, matrix, etc.
     @param ptColIndices Table, with num-local-eqns rows, and the i-th row is of
             length ptRowLengths[i].
     @param ptRowLengths
     @param blkColIndices Table, with num-local-blkEqns rows, and the i-th row
             is of length blkRowLengths[i].
     @param blkRowLengths
     @param ptRowsPerBlkRow The i-th local block-equation corresponds to
          ptRowsPerBlkRow[i] point-equations.
   */
   virtual int setMatrixStructure(int** ptColIndices,
                                   int* ptRrowLengths,
                                   int** blkColIndices,
                                   int* blkRowLengths,
                                   int* ptRowsPerBlkRow) = 0;


   /** Request an interface from ESI_Broker. This will be called by the
     FEI implementation to request, for example, a matrix interface into which
     matrix-data (coefficients, indices) will be passed. All purely algebraic
     data will be assembled directly into matrix/vector interfaces by the FEI
     implementation. The behavior of this function should be to check whether
     an object already exists with the requested instanceName. If so, then just
     return the appropriate interface for that object. If not, then create the
     appropriate concrete class instance, associate it with this instanceName,
     and finally return the requested interface for it.
     @param instanceName Instance name. E.g., "b_Vector" or "x_Vector" if the
        interface being requested is "ESI_Vector".
     @param interface Interface name. E.g., "ESI_Solver".
   */
   virtual int getInterfaceInstance(const char* instanceName,
				    const char* interface,
				    void*& objPtr) = 0;

   /** Set an interface into ESI_Broker. The FEI implementation may use this
     function to give ESI_Broker an object instantiated somewhere else. An
     example would be passing in an already-instantiated preconditioner if, for
     some reason, it was desired to mix components from another library or
     simply from another instance of this library.
    @param instanceName Same as argument to 'getInterfaceInstance'.
    @param interface Same as argument to 'getInterfaceInstance'.
   */
   virtual int setInterfaceInstance(const char* instanceName,
				    const char* interface, void* objPtr) = 0;
};

#endif

