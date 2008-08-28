#ifndef I_O_CONT_FILE_UTILS_H
#define I_O_CONT_FILE_UTILS_H

#include "Teuchos_ParameterList.hpp"

// Write the header of the continuation file
bool WriteHeaderToContFile( const std::string & fileName,
    const Teuchos::ParameterList & fileParams );

// Add a line to the continuation file
bool UpdateContFile( const std::string & fileName,
    const int & idStep,
    const Teuchos::ParameterList & fileParams );

// Restart a continuation file from idStep
bool RestartContFile( const std::string & fileName, const int & idStep );

// Removes the rather annoying "[unused]" flag in the parameters
bool TouchContFileParameters( Teuchos::ParameterList & fileParams );

#endif
