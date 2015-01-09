// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_SAVEEIGENDATA_ABSTRACTSTRATEGY_H
#define LOCA_SAVEEIGENDATA_ABSTRACTSTRATEGY_H

#include "Teuchos_RCP.hpp"

#include "NOX_Abstract_Group.H" // for ReturnType

namespace LOCA {

  //! Collection strategies to save eigenvector/value data
  namespace SaveEigenData {

    //! Abstract interface class strategies to save eigenvector/value data
    /*!
     * AbstractStrategy defines an abstract interface for saving eigenvectors
     * and eigenvalues that are computed at each continuation step.  This
     * is important because this data is often useful for restarting
     * continuations near bifurcation points and gives the user flexibility
     * in how this data is stored.
     *
     * The interface currently defines one pure virtual method,
     * save(), to save any eigenvectors or values as specified by the user.
     * Derived classes should implement this method for a particular strategy
     * to save this data, which is usually highly application code dependent.
     * Constructors for derived classes should be of the form:
     *
     * \code
     * class Derived : public AbstractStrategy {
     * public:
     *   Derived(
     *   const Teuchos::RCP<LOCA::GlobalData>& global_data,
     *   const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& eigenParams);
     *   ...
     * };
     * \endcode
     *
     * where \c global_data is the LOCA global data object, \c topParams is
     * the parsed top-level parameter list, and \c eigenParams
     * is a parameter list of eigensolver parameters.  This list should also
     * specify which and how many eigenvectors/values to save as defined
     * by the strategy.
     *
     * This class and its children follow the Strategy pattern as defined
     * in Erich Gamma, et al. "Design Patterns:  Elements of Reusable
     * Object-Oriented Software." Addison Wesley, Boston, MA, 1995.
     */
    class AbstractStrategy {

    public:

      //! Constructor
      AbstractStrategy() {}

      //! Destructor
      virtual ~AbstractStrategy() {}

      //! Save eigenvalues/eigenvectors
      /*!
       * \param evals_r [out] Real eigenvalues
       * \param evals_i [out] Imaginary eigenvalues
       * \param evecs_r [out] Real eigenvectors
       * \param evecs_i [out] Imaginary eigenvectors
       * \returns ReturnType code indicating success or failure
       */
      virtual NOX::Abstract::Group::ReturnType
      save(Teuchos::RCP< std::vector<double> >& evals_r,
       Teuchos::RCP< std::vector<double> >& evals_i,
       Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_r,
       Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_i) = 0;

    private:

      //! Private to prohibit copying
      AbstractStrategy(const AbstractStrategy&);

      //! Private to prohibit copying
      AbstractStrategy& operator = (const AbstractStrategy&);

    }; // Class AbstractStrategy

  } // Namespace SaveEigenData

} // Namespace LOCA

#endif // LOCA_SAVEEIGENDATA_ABSTRACTSTRATEGY_H
