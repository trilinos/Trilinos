/*@HEADER
// ***********************************************************************
//
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef IFPACK_DYNAMIC_FACTORY_H
#define IFPACK_DYNAMIC_FACTORY_H

#if defined(Ifpack_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Ifpack package is deprecated"
#endif
#endif

#include <ostream>
#include <string>
#include <map>
#include <algorithm>

#include "Ifpack_ConfigDefs.h"
#include "Ifpack_Preconditioner.h"
#include "Teuchos_iostream_helpers.hpp"
#include "Ifpack_AdditiveSchwarz.h"


#ifdef HAVE_HYPRE
#include "Ifpack_Hypre.h"
#endif


//! Ifpack_DynamicFactory
/*!
  TODO: write class description

  \author Radu Popescu <radu.popescu@epfl.ch>
  \date 2013-05-02
*/

class Ifpack_DynamicFactory {
public:
  // The prototype of the preconditioner builder function
  typedef Ifpack_Preconditioner* (*builderFunction)(Epetra_RowMatrix*, int, bool, bool);

  /** \brief Creates an instance of Ifpack_Preconditioner given the std::string
   * name of the preconditioner type (can fail with bad input).
   *
   * \param PrecType (In) - String name of preconditioner type to be created.
   *
   * \param Matrix (In) - Matrix used to define the preconditioner
   *
   * \param overlap (In) - specified overlap, defaulted to 0.
   *
   * Returns <tt>0</tt> if the preconditioner with that input name does not
   * exist.  Otherwise, return a newly created preconditioner object.  Note
   * that the client is responsible for calling <tt>delete</tt> on the
   * returned object once it is finished using it!
   */
  Ifpack_Preconditioner* Create(const std::string PrecType,
                                Epetra_RowMatrix* Matrix,
                                const int overlap = 0,
                                bool overrideSerialDefault = false);

  // Static methods
  /** \brief Initializes the static data of the Ifpac_DynamicFactory class
   *
   * Returns true if initialization succeeded, otherwise false or
   * FILE_NOT_FOUND
   */
  static bool Initialize();

  /** \brief Register a new preconditioner with the factory
   *
   * \param PrecName - String name of the new preconditioner
   *
   * \param PrecBuilder - function pointer to the builder function
   *
   * Returns <tt>0</tt> if ok, otherwise <tt>1</tt>
   */
  static int RegisterPreconditioner(const std::string PrecName,
                                                                    builderFunction PrecBuilder);

  // Static methods
  /** \brief Prints the current list of registered preconditioners
   *
   */
  static void Print(std::ostream& os = std::cout);

  // Templated build function
  template <typename PrecType, bool StandAlone>
  static Ifpack_Preconditioner* buildPreconditioner(Epetra_RowMatrix* Matrix,
                                                                                                    int Overlap,
                                                                                                    bool Serial,
                                                                                                    bool OverrideSerialDefault);

private:
  static std::map<std::string, builderFunction> PreconditionerMap_;
  static int NumPreconditioners_;
  static bool Initialized_;
};

// Templated build function
template <typename PrecType, bool StandAlone>
Ifpack_Preconditioner*
Ifpack_DynamicFactory::buildPreconditioner(Epetra_RowMatrix* Matrix,
                                                                                   int Overlap,
                                                                                   bool Serial,
                                                                                   bool OverrideSerialDefault)
{
        if (StandAlone || (Serial && !OverrideSerialDefault)) {
                return new PrecType(Matrix);
        } else {
                return new Ifpack_AdditiveSchwarz<PrecType>(Matrix, Overlap);
        }
}

#endif // IFPACK_DYNAMIC_FACTORY_H
