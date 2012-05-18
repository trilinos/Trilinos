// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov),
//                    Denis Ridzal  (dridzal@sandia.gov),
//                    Kara Peterson (kjpeter@sandia.gov).
//
// ************************************************************************
// @HEADER

#ifndef __TrilinosCouplings_IntrepidPoissonExampleHelpers_hpp
#define __TrilinosCouplings_IntrepidPoissonExampleHelpers_hpp

/// \file TrilinosCouplings_IntrepidPoissonExampleHelpers.hpp
/// \brief Helper functions for Poisson test problem with Intrepid + Pamgen.
///
/// This directory contains two versions of a Poisson test problem
/// that uses Intrepid to generate a finite element discretization on
/// a Pamgen-generated mesh.  One version uses Epetra, and the other
/// uses Tpetra.  This header file contains utility functions that are
/// useful for both versions.

#include <Teuchos_CommandLineProcessor.hpp>
#include <string>

namespace TrilinosCouplings {
namespace IntrepidPoissonExample {

/// \brief Make a Pamgen mesh specification for the Poisson test problem.
///
/// Pamgen accepts mesh descriptions as human-readable strings in a
/// Pamgen-specific mesh description language.  This function creates
/// a mesh description for the Poisson test problem in this directory.
/// The mesh description returned by this function describes a
/// rectangular prism with cubic elements, with nx elements along the
/// x dimension, ny elements along the y dimension, and nz elements
/// along the z dimension.
///
/// \pre nx > 0 && ny > 0 && nz > 0
std::string
makeMeshInput (const int nx, const int ny, const int nz);

/// \brief Set default command-line argument values for Poisson test problem.
///
/// nx, ny, and nz mean the same thing as the input arguments of
/// makeMeshInput().
///
/// \param nx [out] Number of elements along the x dimension.
/// \param ny [out] Number of elements along the y dimension.
/// \param nz [out] Number of elements along the z dimension.
/// \param xmlInputParamsFile [out] Name of XML file encoding
///   an input ParameterList for the Poisson test problem.
/// \param verbose [out] Whether to print verbose status output.
/// \param debug [out] Whether to print debugging output.
void
setCommandLineArgumentDefaults (int& nx,
                                int& ny,
                                int& nz,
                                std::string& xmlInputParamsFile,
                                bool& verbose,
                                bool& debug);

void
setUpCommandLineArguments (Teuchos::CommandLineProcessor& cmdp,
                           int& nx,
                           int& ny,
                           int& nz,
                           std::string& xmlInputParamsFile,
                           bool& verbose,
                           bool& debug);

/// \brief Parse and partially validate the command-line arguments.
///
/// \param printedHelp [out] Whether the --help option was specified
///   on the command line.
/// \param argc [in] Same as the argument to main().
/// \param argv [in/out] Same as the argument to main().
///
/// All other arguments are the same as
/// setCommandLineArgumentDefaults().
void
parseCommandLineArguments (Teuchos::CommandLineProcessor& cmdp,
                           bool& printedHelp,
                           int argc,
                           char* argv[],
                           int& nx,
                           int& ny,
                           int& nz,
                           std::string& xmlInputParamsFile,
                           bool& verbose,
                           bool& debug);

} // namespace IntrepidPoissonExample
} // namespace TrilinosCouplings

#endif // __TrilinosCouplings_IntrepidPoissonExampleHelpers_hpp
