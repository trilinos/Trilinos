#ifndef XPETRA_REGION_UTILS_IMPL_HPP_
#define XPETRA_REGION_UTILS_IMPL_HPP_

// Teuchos
#include <Teuchos_Comm.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_RCP.hpp>

// Xpetra
#include "Xpetra_RegionUtils_decl.hpp"

template<class GO>
bool Xpetra::compareRegions(const std::tuple<GO,GO>& lhs,
    const std::tuple<GO,GO>& rhs)
{
  // First, we prioritize the sorting according to the region label.
  // If the region is the same, then the sorting looks at the composite node index.
  if (std::get<1>(lhs) < std::get<1>(rhs))
    return true;
  else if (std::get<1>(lhs) == std::get<1>(rhs))
    return std::get<0>(lhs) < std::get<0>(rhs);
  else
    return false;
}

template<class GO>
bool Xpetra::compareNodes(const std::tuple<GO,GO>& lhs,
    const std::tuple<GO,GO>& rhs)
{
  return std::get<0>(lhs) < std::get<0>(rhs);
}

#endif /* #ifndef XPETRA_REGION_UTILS_IMPL_HPP_ */
