// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOSSOLVERFACTORY_TPETRA_HPP
#define BELOSSOLVERFACTORY_TPETRA_HPP

#include "Belos_Details_Tpetra_registerSolverFactory.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosMultiVecTraits_Tpetra.hpp"
#include "BelosOperatorTraits_Tpetra.hpp"

namespace Belos {

template<class SC, class MV, class OP, class DM = DefaultDenseMatrix<int, SC>>
class TpetraSolverFactory : public Impl::SolverFactoryParent<SC, MV, OP, DM>
{
  public:
    TpetraSolverFactory() {
      Details::Tpetra::registerSolverFactory();
    };
};

namespace Impl {

template<class SC, class LO, class GO, class NT, class DM>
class SolverFactorySelector<SC,Tpetra::MultiVector<SC, LO, GO, NT>,Tpetra::Operator<SC, LO, GO, NT>,DM> {
  public:
    typedef TpetraSolverFactory<SC,Tpetra::MultiVector<SC, LO, GO, NT>,Tpetra::Operator<SC, LO, GO, NT>, DM> type;
};

} // namespace Impl
} // namespace Belos

// Every consumer of this header (this DLL's own solver-registrar .cpp files,
// and any downstream executable or shared library) must import the single,
// DLL-exported instantiation of SolverFactoryParent defined in
// BelosSolverFactory_Tpetra_ETI.cpp, rather than silently creating its own
// private one. SolverFactoryParent's registry (get_solverManagers()) is a
// function-local static inside a template static method - vague linkage -
// and Windows PE does not merge vague-linkage template statics across DLL/EXE
// boundaries the way ELF does, so without this extern template declaration
// every consumer would get its own, empty, unsynchronized copy of the solver
// registry.
#include "BelosTpetra_DLLExportMacro.h"
#include "TpetraCore_ETIHelperMacros.h"

TPETRA_ETI_MANGLING_TYPEDEFS()

#define BELOS_TPETRA_SOLVERFACTORYPARENT_EXTERN_INSTANT( SC, LO, GO, NT ) \
  extern template class BELOSTPETRA_LIB_DLL_EXPORT \
    Belos::Impl::SolverFactoryParent<SC, ::Tpetra::MultiVector<SC,LO,GO,NT>, \
      ::Tpetra::Operator<SC,LO,GO,NT>, ::Belos::DefaultDenseMatrix<int,SC> >;

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( BELOS_TPETRA_SOLVERFACTORYPARENT_EXTERN_INSTANT )

#undef BELOS_TPETRA_SOLVERFACTORYPARENT_EXTERN_INSTANT

#endif // BELOSSOLVERFACTORY_TPETRA_HPP
