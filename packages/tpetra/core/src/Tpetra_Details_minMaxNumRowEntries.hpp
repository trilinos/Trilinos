#ifndef TPETRA_DETAILS_MINMAXNUMROWENTRIES_HPP
#define TPETRA_DETAILS_MINMAXNUMROWENTRIES_HPP

#include "TpetraCore_config.h"
#include "Kokkos_ArithTraits.hpp"
#include "Kokkos_Core.hpp"

namespace Tpetra {
namespace Details {

/// \brief Given the row offsets array of a local sparse graph or
///   matrix, return the min and max of the number of entries in each
///   local row.
template<class ViewType, class LO>
std::pair<LO, LO>
minMaxNumRowEntries(const char kernelLabel[],
                    const ViewType& ptr,
                    const LO lclNumRows);

namespace Impl {

template<class ViewType, class LO>
class MinMaxNumRowEntries {
private:
  using result_type = Kokkos::pair<LO, LO>;
  typename ViewType::const_type ptr_;

public:
  MinMaxNumRowEntries(const ViewType& ptr) : ptr_(ptr) {}

  KOKKOS_INLINE_FUNCTION void
  init(result_type& dst) const
  {
    dst.first = Kokkos::ArithTraits<LO>::max();
    dst.second = Kokkos::ArithTraits<LO>::min();
  }

  KOKKOS_INLINE_FUNCTION void
  join(volatile result_type& dst,
       const volatile result_type& src) const
  {
    dst.first = (src.first < dst.first) ? src.first : dst.first;
    dst.second = (src.second > dst.second) ? src.second : dst.second;
  }

  KOKKOS_INLINE_FUNCTION void
  operator() (const LO lclRow, result_type& minMax) const
  {
    // LO (the local ordinal type) must be able to express the number
    // of entries in any one row of the sparse graph or matrix.  The
    // total number of entries over all the rows on a process need not
    // fit in LO, though.
    const LO numEnt = static_cast<LO>(ptr_(lclRow+1) - ptr_(lclRow));
    minMax.first = numEnt < minMax.first ? numEnt : minMax.first;
    minMax.second = numEnt > minMax.second ? numEnt : minMax.second;
  }
};

} // namespace Impl

template<class ViewType, class LO>
std::pair<LO, LO>
minMaxNumRowEntries(const char kernelLabel[],
                    const ViewType& ptr,
                    const LO lclNumRows)
{
  using execution_space = typename ViewType::execution_space;
  using range_type = Kokkos::RangePolicy<execution_space, LO>;
  using result_type = Kokkos::pair<LO, LO>;

  if (lclNumRows == 0) {
    return {LO(0), LO(0)};
  }
  else {
    result_type minMaxNumEnt{Kokkos::ArithTraits<LO>::max(),
                             Kokkos::ArithTraits<LO>::min()};
    Kokkos::parallel_reduce(
      kernelLabel, range_type(0, lclNumRows),
      Impl::MinMaxNumRowEntries<ViewType, LO>(ptr),
      minMaxNumEnt);
    return {minMaxNumEnt.first, minMaxNumEnt.second};
  }
}

} // namespace Details
} // namespace Tpetra

#endif // TPETRA_DETAILS_MINMAXNUMROWENTRIES_HPP
