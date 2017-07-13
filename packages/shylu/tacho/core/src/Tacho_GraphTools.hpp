#ifndef __TACHO_GRAPH_TOOLS_HPP__
#define __TACHO_GRAPH_TOOLS_HPP__

/// \file Tacho_GraphTools.hpp
/// \brief Interface to scotch reordering
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

  template<typename OrdinalType, 
           typename SizeType = OrdinalType,
           typename SpaceType = void>
  class GraphTools {
  public:
    typedef OrdinalType ordinal_type;
    typedef SizeType    size_type;
    typedef SpaceType   space_type;

    typedef Kokkos::View<ordinal_type*,space_type> ordinal_type_array;
    typedef Kokkos::View<size_type*,   space_type> size_type_array;

    template<typename CrsMatBaseType>
    KOKKOS_INLINE_FUNCTION
    static void getGraph(size_type_array &rptr,
                         ordinal_type_array &cidx,
                         const CrsMatBaseType A) {
      // parallel version is not effective here
      size_type jj = 0;
      for (ordinal_type i=0;i<A.NumRows();++i) {
        const size_type jbegin = A.RowPtrBegin(i), jend = A.RowPtrEnd(i);
        rptr[i] = jj;
        for (size_type j=jbegin;j<jend;++j) {
          const auto aj = A.Col(j);
          if (i != aj)
            cidx[jj++] = aj;
        }
      }
      rptr[A.NumRows()] = jj;
    }
  };

}

#endif
