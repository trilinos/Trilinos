#ifndef __ESI_Properties_h_seen
#define __ESI_Properties_h_seen

namespace esi {

//! Abstract algorithm control parameters and string named data io interface.
/** This is a generic properties input/output interface which may 
    be supported at the discretion of the concrete class implementer.
    Since it is discretionary, these functions are not required in or
    directly accessible from the base classes such as Object, Vector,
    or Solver. The example below shows how properties are accessed.

    Solver and preconditioner authors are free to use this interface
    both for input of control parameters and for output of profiling
    data and even result data where those, including arrays.

    This is a Closed hash -- one value (and value type) per key.
    Setting a new value and type for a key which has a value of 
    another type already is an error.

    Some properties/control parameters (such as iterationLimit
    which must be present in all ESI_Solvers) may be represented
    redundantly in this <em>optional<\em>Properties itnerface.


    As an example:
    If you have a solver that cares about the properties of
    a Matrix or Operator, simply :
      esi::Properties *p;
      err = mtx->getInterface("esi::Properties", (void *)p);
      bool is_psym = false;
      if (err >= 0) {  // <0 means no properties interface p.
        err2 = p->getBool("esi.P-Symmetric",psym); // ignore err2 is ok
        // if esi.P-Symmetric unknown, nothing psym is returned unchanged.
      }

    Implementation notes:

    Since this depends on Ordinal and Scalar, it
    cannot be in the Object interface.  Since some
    operators/matrices/solvers/vectors/pcs may have no special
    properties, the interface is not required (by inheritance)
    into those classes.

    We will create a reference implementation of this interface which
    can be inherited or used as an internal data member (delegated to)
    to provide this interface. 

    Probably most implementers will support it by delegation
    or inheritance of the reference-- no bonus points in the math
    world for implementing hashes. Those implementers that want to
    back this interface by an existing resource mechanism
    are of course free to do so.
*/
template<class Scalar, class Ordinal>
class Properties {
public:
  /** Defines a new string property or changes the value of an existing one.
      Because the string arguments are copied, it is okay to call this
      with string literal arguments.
      @param key The property name, which will be copied.
      @param value The property value, which will be copied.
                   A 0 value will cause the key to remain
       with a NULL value. Use 'unset' to delete keys.
      @return 0 if ok, -2 if not ok due to a type conflict,
      and -1 for all other reason.
  */
  virtual ErrorCode setString( const char * key, 
                               const char * value ) = 0;

  /** Defines a new Scalar property or changes its value.
      @param key The property name, which will be copied.
      @param value The property value to be set.
      @return 0 if ok, -2 if not ok due to a type conflict,
      and -1 for all other reason.
  */
  virtual ErrorCode setScalar( const char * key, 
                               Scalar floatingPointProp ) = 0;

  /** Defines a new integral property or changes its value.
      @param key The property name, which will be copied.
      @param value The property value to be set.
      @return 0 if ok, -2 if not ok due to a type conflict,
      and -1 for all other reason.
  */
  virtual ErrorCode setOrdinal( const char * key, 
                                Ordinal integerProp ) = 0;

  /** Defines a new boolean property or changes its value.
      @param key The property name, which will be copied.
      @param value The property value to be set.
      @return 0 if ok, -2 if not ok due to a type conflict,
      and -1 for all other reason.
  */
  virtual ErrorCode setBool( const char * key, 
                             bool boolProp ) = 0;

  /** Defines a new arbitrary object property, or substitutes
      a new object for the currently stored one.
      @param key The property name, which will be copied.
      @param value If the void pointer given should be reference 
      counted, that is the caller's responsibility.  The pointer 
      will be stored until the next set or unset on the same key.
      The void pointer may be to an object, an array, or anything else.
      @return 0 if ok, -2 if not ok due to a type conflict,
      and -1 for all other reasons.
  */
  virtual ErrorCode setPointer( const char * key, 
                                void * objProp ) = 0;

  /**  Get a C string property value by C string name.
       @param key The name of the property desired.
       @param value output of the value requested.
       @return 0 if value exists. -2 if it exists and is of 
       the wrong type. -1 does not exist or any other error.
  */
  virtual ErrorCode getString( const char * propName, 
                               const char * & value ) = 0;

  /**  Get a Scalar property value by C string name.
       @param key The string name of the property desired.
       @param value output of the value requested.
       @return 0 if value exists. -2 if it exists and is of 
       the wrong type. -1 does not exist or any other error.
  */
  virtual ErrorCode getScalar( const char * propName, 
                               Scalar & value ) = 0;

  /**  Get an Ordinal property value by C string name.
       @param key The string name of the property desired.
       @param value output of the value requested.
       @return 0 if value exists. -2 if it exists and is of 
       the wrong type. -1 does not exist or any other error.
  */
  virtual ErrorCode getOrdinal( const char * propName, 
                                Ordinal & value ) = 0;

  /**  Get a bool property value by C string name.
       @param key The string name of the property desired.
       @param value output of the value requested.
       @return 0 if value exists. -2 if it exists and is of 
       the wrong type. -1 does not exist or any other error.
   */
  virtual ErrorCode getBool( const char * propName, 
                             bool & value ) = 0;

  /**  Get a pointer property value by C string name.
       @param key The string name of the property desired.
       @param value output of the value requested.
       @return 0 if value exists. -2 if it exists and is of 
       the wrong type. -1 does not exist or any other error.
   */
  virtual ErrorCode getPointer( const char * propName, 
                                void * & value ) = 0;

  /** Remove a key from the properties. All errors are silently 
      ignored.
      @param pointerOut If the key corresponds to a void pointer 
      property, the pointer is returned in pointerOut,
      otherwise pointerOut is 0.
  */
  virtual void unset( const char * propName, 
                      void * & pointerOut ) = 0;

  /** Produce a list of all the keys in an abstract container.
      @param keylist input/output Argv that should be supplied 
      empty (containing no keys). The names of keys in the 
      Properties are appended.
  */
  virtual void getKeys( Argv * keylist ) = 0;

  /** Produce a list of the keys in the main underlying object
      that the object "supports", i.e. sets or reads during computations.
      Other keys may be present for use by third parties that have
      nothing to do with "the algorithm", but these keys will not appear
      in the getKeysSupported list.
      @param keylist input/output Argv that should be supplied empty
      (containing no keys). The names of supported keys in the Properties 
      are appended.
  */
  virtual void getKeysSupported( Argv * keylist ) = 0;

};   // esi::Properties
};   // esi namespace
#endif // __ESI_Properties_h_seen
