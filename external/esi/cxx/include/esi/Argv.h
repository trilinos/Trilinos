#ifndef __ESI_Argv_h
#define __ESI_Argv_h

namespace esi {

//! Simplest string container interface, ala C's argv.
/** A little iterator interface interface around arrays of 
    strings, whereever tbey may be fetched from an ESI 
    interface. 
    
    In particular:
    Exists because we want a simple, abstract alternative
    to the STL. Implementors are free to use the
    STL underneath this interface if they so choose.
   \verbatim
    An example of its use:

    SomeConcreteArgvImpl argvimpl;
    // perhaps the argv class from the reference implementation.
    esi::Argv *argvi = dynamic_cast<esi::Argv *>(&argvimpl);
    esi::ErrorCode err;
    err = object->getRunTimeModelsSupported(argvi);
    if (!err) {
      const char *arg;
      for (i = 0; (arg = argvi->get(i)) != 0; i++) {
        do.something(arg);
      }
      delete argvi;
    }
    argvi = 0;

    And inside the getRunTimeModelsSupported call, a set of calls on
    argvi->appendArg(rtmodelName[i]); occurs, filling in the answers.
   \endverbatim
   This pattern keeps creation/management of all the memory associated with
   the Argv in the scope of the caller or internal to the Argv itself.
*/
class Argv {

 public:

  /** Default destructor. */
  virtual ~Argv( void ) {};

  /** @return the ith element of argv (not to be modified) if exists.
   *       returns 0 otherwise.
   */
  virtual const char * get( int index ) = 0;

  /** @return the current number of strings. 
   * This number obviously changes when appendArg
   */ 
  virtual int getArgCount( void ) = 0; 

  /** Copy the arg given into the next slot in the argv. 
   */
  virtual int appendArg( const char *arg ) = 0;

};     // esi::Argv class
};     // esi namespace
#endif // __ESI_Argv_h
