#ifndef PANZER_PARAMETER_LIBRARY_ACCEPTOR_HPP
#define PANZER_PARAMETER_LIBRARY_ACCEPTOR_HPP

#include "Panzer_ParameterLibrary.hpp"
#include "Teuchos_RCP.hpp"

namespace panzer {

  /** \brief Pure Virtual base class for accepting the parameter library

      This class is used to retrieve the parameter library from an
      object.
  */
  class ParameterLibraryAcceptor {

  public:

    virtual ~ParameterLibraryAcceptor() {}

    virtual Teuchos::RCP<panzer::ParamLib> getParameterLibrary() const = 0;

  };

}

#endif
