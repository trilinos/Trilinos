//Muemex callbacks
//Brian Kelley

#ifndef MUEMEX_CALLBACKS_H
#define MUEMEX_CALLBACKS_H

#include "Teuchos_ParameterList.hpp"
#include <string>
#include <complex>
#include <stdexcept>
#include <vector>
#include "muemexTypes.h"
#include "mex.h"

#if !defined(HAVE_MUELU_MATLAB) || !defined(HAVE_MUELU_EPETRA) || !defined(HAVE_MUELU_TPETRA)
#error "Muemex callbacks require MATLAB, Epetra and Tpetra."
#else

namespace MuemexCallback
{
	//The two callback functions that MueLu can call to run anything in MATLAB
	void callMatlabNoArgs(std::string function);
	std::vector<Teuchos::RCP<MuemexArg>> callMatlab(std::string function, int numOutputs, std::vector<Teuchos::RCP<MuemexArg>> args);
}

#endif //HAVE_MUELU_MATLAB
#endif //MUEMEX_CALLBACKS_H guard
