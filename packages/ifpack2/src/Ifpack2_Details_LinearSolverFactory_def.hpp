/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

/// \file   Ifpack2_Details_LinearSolverFactory_def.hpp
/// \author Mark Hoemmen
/// \brief  Definition of Ifpack2::Details::LinearSolverFactory.

#ifndef IFPACK2_DETAILS_LINEARSOLVERFACTORY_DEF_HPP
#define IFPACK2_DETAILS_LINEARSOLVERFACTORY_DEF_HPP

#include <Ifpack2_Details_LinearSolver.hpp>
#include <Trilinos_Details_LinearSolverFactory.hpp>
#include <type_traits> // std::is_same

namespace Ifpack2 {
namespace Details {

template<class SC, class LO, class GO, class NT>
Teuchos::RCP<typename LinearSolverFactory<SC, LO, GO, NT>::solver_type>
LinearSolverFactory<SC, LO, GO, NT>::
getLinearSolver (const std::string& solverName)
{
  return Teuchos::rcp (new Ifpack2::Details::LinearSolver<SC, LO, GO, NT> (solverName));
}

template<class SC, class LO, class GO, class NT>
void
LinearSolverFactory<SC, LO, GO, NT>::
registerLinearSolverFactory ()
{
  typedef Tpetra::MultiVector<SC, LO, GO, NT> MV;
  typedef Tpetra::Operator<SC, LO, GO, NT> OP;
  typedef typename MV::mag_type mag_type;
  typedef Trilinos::Details::LinearSolverFactory<MV, OP, mag_type> factory_base_type;
  typedef Ifpack2::Details::LinearSolverFactory<SC, LO, GO, NT> factory_impl_type;

#ifdef HAVE_TEUCHOSCORE_CXX11
  typedef std::shared_ptr<factory_base_type> base_ptr_type;
  typedef std::shared_ptr<factory_impl_type> impl_ptr_type;
#else
  typedef Teuchos::RCP<factory_base_type> base_ptr_type;
  typedef Teuchos::RCP<factory_impl_type> impl_ptr_type;
#endif // HAVE_TEUCHOSCORE_CXX11

  impl_ptr_type factory (new factory_impl_type ());
  base_ptr_type factoryBase = factory; // implicit cast to base class

  TEUCHOS_TEST_FOR_EXCEPTION
    (factoryBase.get () == NULL, std::logic_error, "Factory is null!  This "
     "should never happen!  Please report this bug to the Ifpack2 developers.");

#ifdef HAVE_TEUCHOS_DEBUG
  {
    using std::cerr;
    using std::endl;
    using Teuchos::TypeNameTraits;
    cerr << "Registering Ifpack2 LinearSolverFactory for"
         << " SC = " << TypeNameTraits<SC>::name ()
         << ", LO = " << TypeNameTraits<LO>::name ()
         << ", GO = " << TypeNameTraits<GO>::name ()
         << ", NT = " << TypeNameTraits<NT>::name ()
         << ", and mag_type = " << TypeNameTraits<mag_type>::name ()
         << endl;
  }
#endif // HAVE_TEUCHOS_DEBUG
  Trilinos::Details::registerLinearSolverFactory<MV, OP, mag_type> ("Ifpack2", factoryBase);
}

} // namespace Details
} // namespace Ifpack2

// Do explicit instantiation of Ifpack2::Details::LinearSolverFactory,
// for Tpetra objects, with the given Tpetra template parameters.
#define IFPACK2_DETAILS_LINEARSOLVERFACTORY_INSTANT( SC, LO, GO, NT ) \
  template class Ifpack2::Details::LinearSolverFactory<SC, LO, GO, NT>;

#endif // IFPACK2_DETAILS_LINEARSOLVERFACTORY_DEF_HPP
