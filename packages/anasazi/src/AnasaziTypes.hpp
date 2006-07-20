// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef ANASAZI_TYPES_HPP
#define ANASAZI_TYPES_HPP

#include "AnasaziConfigDefs.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ScalarTraits.hpp"

/*! \file AnasaziTypes.hpp
  \brief Types and exceptions used within Anasazi solvers and interfaces.
*/

namespace Anasazi {
  //@{ \name Anasazi Exceptions

  /** \brief An exception class parent to all Anasazi exceptions.
   */
  class AnasaziError : public std::logic_error { 
    public: AnasaziError(const std::string& what_arg) : std::logic_error(what_arg) {} 
  };

  //@}

}

namespace Anasazi {

  //@{ \name Anasazi Structures

  /*! \struct Eigensolution
   *  \brief Struct for storing an eigenproblem solution.
   */
  template <class ScalarType, class MV>
  struct Eigensolution {
    Teuchos::RefCountPtr<MV> Evecs, Espace;
    std::vector<int>         index;
    std::vector<typename Teuchos::ScalarTraits<ScalarType>::magnitudeType>  Evals;
    int numVecs;

    Eigensolution() : Evecs(),Espace(),index(0),Evals(0),numVecs(0) {}
  };

  //@}

  //@{ \name Anasazi Enumerations

  /*!  \enum ReturnType    
       \brief Enumerated type used to pass back information from a solver manager.
  */
  enum ReturnType 
  {
    Converged,       /*!< The solver manager computed the requested eigenvalues. */
    Unconverged      /*!< This solver manager did not compute all of the requested eigenvalues. */
  };


  /*!  \enum ConjType
   *    
   *    \brief Enumerated types used to specify conjugation arguments.
   */
  enum ConjType 
  {
    NO_CONJ,      /*!< Not conjugated */
    CONJ          /*!< Conjugated */
  };


  /*!  \enum TestStatus
       \brief Enumerated type used to pass back information from a StatusTest
  */
  enum TestStatus
  {
    Passed,       /*!< The solver passed the test */
    Failed,       /*!< The solver failed the test */
    Undefined     /*!< The test has not been evaluated on the solver */ 
  };


  /*! \enum MsgType
      \brief Enumerated list of available message types recognized by the eigensolvers.
  */
  enum MsgType 
  {
    Errors = 0,                 /*!< Errors [ always printed ] */
    Warnings = 0x1,             /*!< Internal warnings */
    IterationDetails = 0x2,     /*!< Approximate eigenvalues, errors */
    OrthoDetails = 0x4,         /*!< Orthogonalization/orthonormalization details */
    FinalSummary = 0x8,         /*!< Final computational summary */
    TimingDetails = 0x10,       /*!< Timing details */
    Debug = 0x20                /*!< Debugging information */
  };

  //@}

} // end of namespace Anasazi
#endif
// end of file AnasaziTypes.hpp
