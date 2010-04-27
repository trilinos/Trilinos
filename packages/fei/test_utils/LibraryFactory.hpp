/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _LibraryFactory_hpp_
#define _LibraryFactory_hpp_

#include <fei_macros.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_LibraryWrapper.hpp>

#include <fei_Factory.hpp>

namespace fei {

  /** Create an instance of LibraryWrapper. Throws std::runtime_error
      if the input libraryName is not recognized. (names are case-sensitive)

      @param libraryName Input name of the solver library that is to be used.
      Valid values of this parameter are:<br>
      <ul>
      <li>Aztec
      <li>FETI
      </ul>

      @return shared-pointer holding newly-created LibraryWrapper instance.
   */
  fei::SharedPtr<LibraryWrapper> create_LibraryWrapper(MPI_Comm comm,
						       const char* libraryName);

  /** Create an instance of the fei::Factory interface. Throws std::runtime_error
      if the input libraryName is not recognized. (names are case-sensitive)

      @param libraryName Input name of solver library, same valid values as for
      'create_LibraryWrapper' above, as well as allowing "Trilinos".

      @return shared-pointer holding newly-created fei::Factory instance.
  */
  fei::SharedPtr<fei::Factory> create_fei_Factory(MPI_Comm comm,
						  const char* libraryName);
}//namespace fei

#endif // _LibraryFactory_hpp_

