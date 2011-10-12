#ifndef MUELU_HPP
#define MUELU_HPP

// This header is to simplify examples or user drivers.

#include "MueLu_ConfigDefs.hpp"

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Xpetra
#include <Xpetra_Parameters.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_VectorFactory.hpp>

// Gallery
#define XPETRA_ENABLED // == Gallery have to be build with the support of Xpetra matrices.
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>

// MueLu
#include "MueLu_FactoryManager.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_Utilities.hpp"

#endif
