// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
/// \namespace IntrepidPoissonExample
/// \brief Intrepid Poisson test problem example: common functionality.
///
/// The Intrepid Poisson test problem uses Pamgen to construct a 3-D
/// mesh (a simple rectangular prism with hex elements) in parallel,
/// Sacado automatic differentiation to construct a right-hand side of
/// the PDE corresponding to a given exact solution, and Intrepid to
/// build a discretization.
///
/// We provide two variants of the Intrepid Poisson test: one that
/// fills Epetra objects, and one that fills Tpetra objects.  The two
/// variants do exactly the same things otherwise, so you can use them
/// to compare the performance of Epetra and Tpetra fill.
///
/// This namespace contains functions which both variants can share,
/// because they do not depend on Epetra or Tpetra.  You can include
/// this file in your main() test driver in order to set up and read
/// in command-line arguments, and prepare the Pamgen mesh
/// specification.
namespace IntrepidPoissonExample {

/// \brief Get off-diagonal value for material tensor.
///
/// You can use this value to control the iteration count.  See the
/// documentation of setMaterialTensorOffDiagonalValue() below for
/// examples.
///
/// This value has to be a double, because we read it from the command
/// line; Teuchos' facility for this stores floating-point values as
/// double.
double getMaterialTensorOffDiagonalValue ();

/// \brief Set off-diagonal value for material tensor.
///
/// You can use this value to control the iteration count.  The
/// iteration counts below are for Belos' GMRES with no
/// preconditioning, using the default problem size.
///
/// newVal = -5/4: 209 iterations (CG breaks!)
/// newVal = -1/2: 47 iterations
/// newVal = 0: 40 iterations (CG works)
/// newVal = 1/2: 46 iterations
/// newVal = 3/4: 47 iterations
/// newVal = 1: 59 iterations
/// newVal = 5/4: 183 iterations
/// newVal = 3/2: 491 iterations
/// newVal = 2: 939 iterations (CG breaks!)
void setMaterialTensorOffDiagonalValue (const double newVal);


/// \brief Set parameters for anisotropic diffusion
/// 
/// Q = Rz(theta_z) Ry(theta_y) Rx(theta_x)
///   Rotate in x, then y , then z
///
/// D = diag(str_x,str_y, str_z)
///   Diagonal diffusion strength pre-rotation
///
/// Diffusion Tensor = Q D Q^T
void setDiffusionRotationAndStrength(const std::vector<double>& diff_rotation_angle,  const std::vector<double>& diff_strength);


/// \brief Gets the diffusion Tensor
const std::vector<double> & getDiffusionMatrix();


/// \brief Gets the 2D diffusion Tensor (for 2D problems)
const std::vector<double> & getDiffusionMatrix2D();


/// \brief Use the diffusion tensor rather than the off-diagonal values
bool useDiffusionMatrix();


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
/// \param solverName [out] Name of Belos solver; anything that
///   Belos::SolverFactory understands is valid here.
/// \param verbose [out] Whether to print verbose status output.
/// \param debug [out] Whether to print debugging output.
void
setCommandLineArgumentDefaults (int& nx,
                                int& ny,
                                int& nz,
                                std::string& xmlInputParamsFile,
                                std::string& solverName,
                                bool& verbose,
                                bool& debug);

/// \brief Prepare for parsing command-line arguments.
///
/// This sets up command-line options for the given arguments, which
/// correspond to the arguments of setCommandLineArgumentDefaults.
/// This function reads in the default values of the arguments on
/// input.  When the command-line arguments are parsed (by
/// parseCommandLineArguments()), their values will be overwritten
/// with the values specified on the command line.
void
setUpCommandLineArguments (Teuchos::CommandLineProcessor& cmdp,
                           int& nx,
                           int& ny,
                           int& nz,
                           std::string& xmlInputParamsFile,
                           std::string& solverName,
                           double& tol,
                           int& maxNumIters,
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
                           std::string& solverName,
                           bool& verbose,
                           bool& debug);

} // namespace IntrepidPoissonExample
} // namespace TrilinosCouplings

#endif // __TrilinosCouplings_IntrepidPoissonExampleHelpers_hpp
