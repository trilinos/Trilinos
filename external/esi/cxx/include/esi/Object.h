#ifndef __ESI_Object_h
#define __ESI_Object_h

/** ESI (Equation Solver Interfaces) is an open standard
    API for both plug-n-play and mix-&-match linear solvers.

    A brief description of the ESI architecture follows.
    Detailed discussions will be found in the report (REFERENCE HERE).
    The API architecture is, to the maximum extent possible,
    a flat set of abstract interfaces, which are compatible
    to the maximum extent possible with 'Ports' in the
    HPC component standard evolving in the DOE Common
    Component Architecture forum (http://www.cca-forum.org).
    To allow delegation to be used in implementations, a
    base Object interface supplying elementary introspections
    is included, though in future work the ESI forum
    may replace this Object interface with common interfaces
    yet to be specified by the CCA or by the Babel/SIDL project
    in the CASC division of LLNL.

    The Equation Solver Interfaces (ESI) forum is a
    group of university, industry, and DOE Laboratory research 
    and applications programmers interested in creating
    standards for scalable parallel linear solvers,
    principally iterative solvers. The standards
    are intended to capture both plug-n-play library
    interchange and, at a deeper level, to allow the mixing
    of libraries to allow experiments in applications,
    particularly experiments involving preconditioners
    and hybrid high performance matrix implementations.

    The ESI forum maintains a mailing list for interface
    discussions: if-forum@z.ca.sandia.gov, aweb voting site 
    (courtesy LLNL/CASC) at https://www-casc.llnl.gov/quorum/
    and a download and development site at 
    http://www.eterascale.com/esi/.
*/
namespace esi {

//! The core ESI interface so we can support both dynamic_cast and delegation.
/**
   This is the ESI "abstract object", from which everything else in ESI is
   derived. It contains the getInterface function, which is used
   to obtain different "windows"/"views" on implementations.
   It features reference counting as outlined below.
*/
class Object 
{
 public:

  /** Default destructor.  */
  virtual ~Object( void ) {};

  /**
   * Obtaining another 'view' on an object.
   * This may be done with dynamic_cast<>() or with delegation internally. 
   * Implementations may make all their interfaces and implementation
   * classes available via this method or only their ESI interfaces.
   * 
   * @param name A string containing the c++ class name of the
   *        view desired.
   * @param iface A pointer to be filled up (or set to 0) with a
   *        pointer to the interface requested. iface should be
   *        cast properly on either input or output, e.g. input:
   *        \verbatim
   *        esi::Object *iface;
   *        iface = 0;
   *        o->getInterface("esi::Object",static_cast<void *>(iface));
   *        \endverbatim
    *
   * @return A number < 0 if the name requested is not supported;
   * in this case iface will be 0/NULL.
   */
  virtual ErrorCode getInterface( const char * name, 
                                  void * & iface ) = 0;

  /**
   * Returns a container of the names of interfaces supported by this.
   * The caller should, when finished, destroy it.
   * In typical code:
     \verbatim
        esi::Argv *list;
        SomeConcreteArgvImpl argvimpl;
        // perhaps the argv class from the reference implementation.
        esi::Argv *list = dynamic_cast<esi::Argv *>(&argvimpl);
        getRunTimeModelsSupported(list);
        // do something
     \endverbatim
   * @param list Output container. Will have zero or more elements inside.
   */
  virtual ErrorCode getInterfacesSupported( Argv * list ) = 0;

  /**
   *  Function to add a "communicator" (RunTimeModel) to the object
   *  possibly from application level code or possibly required
   *  during construction. The underlying implementation must allow
   *  arbitrarily many distinct runtime models to be set.
   *
   *  @param name A C string, for which the only "well known value"
   *              is "MPI", to register the given void pointer by.
   *  @param comm An unknown object to be associated with this. 
   *         If comm is 0, any already stored pointer with the
   *         same name is forgotten and 0 is returned.
   *         If comm is non-0 and the name is already in use
   *         storing a _different_ pointer, an error is returned and
   *         comm is ignored. Reset the pointer for a given name to 0
   *         before assigning a new non-0 value.
   *         The typical case is expected to be 
            \verbatim
             MPI_Comm myComm; // ...
             err = setRunTimeModel("MPI",static_cast<void *>(&myComm));
            \endverbatim
   *         We cannot assume MPI_Comm is a pointer, so &myComm may be
   *         in fact a pointer to a pointer in some MPI implementations.
   *         NOTE: comm may in fact be any kind of pointer, except in
   *         the case of name=="MPI".
   *
   *  @return 0 if ok, -1 if the name given is already in use.
   */
  virtual ErrorCode setRunTimeModel( const char * name, 
                                     void * comm ) = 0;

  /**
   * Retrieves the named comm from the underlying object
   * if the name is known.
   * On input, comm must have the value 0.
   * @param name The C string name of the run time model desired. "MPI"
   *             corresponds to the mpi communicator on which the
   *             object lives. Other values may be implemented as well.
   * @param comm A pointer that will be the run time model requested 
   *             on successful returns.
   *        The typical case is expected to be
           \verbatim
            esi::Object *o;
            MPI_Comm *commPointer;
            o->getRunTimeModel("MPI",static_cast<void *>(commPointer));
           \endverbatim
   *
   * @return 0 if successful, -1 if the name given is not known.
   */
  virtual ErrorCode getRunTimeModel( const char * name, 
                                     void * & comm )  = 0;

  /**
   * Returns a container of the names of RTModels known to the object
   * by setRTModel() or construction.
   * The caller should, when finished, destroy it.
   * In typical code:
     \verbatim
        esi::Argv *list;
        getRTModelsSupported(list);
        // do something
        delete list; 
     \endverbatim
   * @param list Output container. Will have zero or more elements inside.
   */
  virtual ErrorCode getRunTimeModelsSupported( Argv * list ) = 0;

  /** Increment reference count.
   *  Design notes:
     \verbatim
      One cannot memory-safely use ESI objects
      of one library (e.g., ISIS_vector) as parts when
      constructing objects from another library unless
      we standardize a reference counting or cloning mechanism 
      for all objects. We choose to do reference counting. Those
      wishing to wrap the reference counting mechanism with
      smart pointers are at liberty to do so.
     \endverbatim

     What follows are examples of the problem to be solved
     and guidelines on how one must use the functions 
     addReference/deleteReference to solve it.

     \verbatim
      As a specific example of the problem:
        We want either cloning or reference counting
        on the ESI IndexSpace if SNL_Vector has the
        constructor SNL_Vector(esi::IndexSpace *is);
    
        If the SNL_Vector stores the value 'is' (most memory 
        efficient), then 'is' needs to be reference counted or 
        deleting 'is' leaves a dangling pointer in the vector.
        
        If the SNL_Vector queries 'is' to get its data
        and then builds an SNL_IndexSpace from that data,
        it is robust but inefficient.
        
        If the SNL_Vector can clone 'is', it is more robust
        (reduced chance of copy error) and still inefficient.
        
        If SNL_Vector were instead SNL_Solver, it could
        never even consider cloning a MatrixData given during 
        construction. This would rule out cloning as an option, 
        or else require that Solvers be essentially state-data 
        free(highly unlikely).
     \endverbatim

     \verbatim
     Guidelines for memory management of ESI-derived objects.
     1)  delete is never explicitly called on an ESI interface,
         instead deleteReference is called any place that delete
         would have been appropriate.
     2)  ESI objects returned from functions (including getInterface)
         already have their reference counts adjusted accordingly, so
         only deleteReference is used.
         That is, if you fetch a pointer, ptr, to an ESI object by any
         method, you call ptr->deleteReference() when done with
         the pointer.
         In the particular case of getInterface, each time an interface
         is handed out, the reference count goes up.
     3)  When the reference count on an object drops to 0, it
         deletes itself from within the deleteReference function.
     4)  Constructors of objects set the reference count to 1.
     5)  When an ESI interface is passed to a function, that function
         will call addReference if it caches a pointer to the interface
         that persists past the return of the function.
     \endverbatim

     \verbatim
     Examples:
     See the reference implementations (REFERENCE HERE).
     \endverbatim
  */
  virtual ErrorCode addReference( void ) = 0;

  /** Decrement reference count, and call delete if it reaches 0. */
  virtual ErrorCode deleteReference( void ) = 0;

}; // esi::Object
}; // esi namespace
#endif // __ESI_Object_h
