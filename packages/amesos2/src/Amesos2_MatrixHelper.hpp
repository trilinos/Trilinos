#ifndef AMESOS2_MATRIXHELPER_HPP
#define AMESOS2_MATRIXHELPER_HPP

#include "Amesos2_MatrixAdapter.hpp"
#include "Amesos2_MultiVecAdapter.hpp"

namespace Amesos {

/** \brief convert Matrices and MultiVectors into the appropriate format for a
 * third-party solver.
 *
 * The particular functions that must be implemented, and the signatures in
 * each will vary depending on the third-party solver's needs.  This templated
 * \c struct will just provide a central location for functions which deal
 * with Matrix and MultiVector conversions.
 *
 * \tparam ConcreteSolver The Amesos::Solver instance for which the functions hold.
 */
template <template <typename,typename> class ConcreteSolver>
struct MatrixHelper 
{};

} // end namespace Amesos

#endif  // AMESOS2_MATRIXHELPER_HPP
