// We need both the decl and def files here, whether or not ETI is on.
#include "Ifpack2_Details_LinearSolverFactory_decl.hpp"
#include "Ifpack2_Details_LinearSolverFactory_def.hpp"
// We need this whether or not ETI is on, in order to define typedefs
// for making Tpetra's macros work.
#include "TpetraCore_ETIHelperMacros.h"

// Define typedefs that make the Tpetra macros work.
TPETRA_ETI_MANGLING_TYPEDEFS()

#ifdef HAVE_IFPACK2_EXPLICIT_INSTANTIATION

// Macro for doing explicit instantiation of
// Ifpack2::Details::LinearSolverFactory, for Tpetra objects, with
// given Tpetra template parameters (Scalar, LocalOrdinal,
// GlobalOrdinal, Node).
#define LCLINST(SC, LO, GO, NT) \
  template class Ifpack2::Details::LinearSolverFactory<SC, LO, GO, NT>;

// Do explicit instantiation of Ifpack2::Details::LinearSolverFactory, for
// Tpetra objects, for all combinations of Tpetra template parameters
// for which Tpetra does explicit template instantiation (ETI).
TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( LCLINST )

#endif // HAVE_IFPACK2_EXPLICIT_INSTANTIATION

//
// Register Ifpack2's solver factory/ies with the central registry.
// Do this whether or not ETI is on.
//

namespace { // (anonymous)

// Inside a macro, ## means "join up the characters".
// For example, if LO=int, "_##LO##_" becomes "_int_".

#define IFPACK2_DETAILS_REGISTER(SC, LO, GO, NT) \
  Ifpack2::Details::RegisterLinearSolverFactory<SC, LO, GO, NT> \
    registerer_Tpetra_##SC##_##LO##_##GO##_##NT ;

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( IFPACK2_DETAILS_REGISTER )

} // namespace (anonymous)
