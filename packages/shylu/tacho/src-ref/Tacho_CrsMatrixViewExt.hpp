#ifndef __TACHO_CRS_MATRIX_VIEW_EXT_HPP__
#define __TACHO_CRS_MATRIX_VIEW_EXT_HPP__


/// \file Tacho_CrsMatrixViewExt.hpp
/// \brief Extended matrix view that has nested dense block.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho {

  template<typename VT, typename OT, typename ST, typename SpT>
  class DenseMatrixBase;

  template<typename CrsMatBaseType>
  class CrsMatrixViewExt : public CrsMatrixView<CrsMatBaseType> {
  public:
    typedef typename CrsMatBaseType::space_type    space_type;

    typedef typename CrsMatBaseType::value_type    value_type;
    typedef typename CrsMatBaseType::ordinal_type  ordinal_type;
    typedef typename CrsMatBaseType::size_type     size_type;

    typedef DenseMatrixBase<value_type,ordinal_type,size_type,space_type> flat_mat_base_type;
    // typedef DenseMatrixBase<TaskView<DenseMatrixView<flat_mat_base_type> >,
    //                         ordinal_type,size_type,space_type> hier_mat_base_type;

  private:

    flat_mat_base_type _A;
    // hier_mat_base_type _H;

  public:

    KOKKOS_INLINE_FUNCTION
    flat_mat_base_type& Flat() { return _A; }

    KOKKOS_INLINE_FUNCTION
    flat_mat_base_type  Flat() const { return _A; }

    KOKKOS_INLINE_FUNCTION
    void copyToFlat() { 
      const ordinal_type m = this->NumRows();
      for (ordinal_type i=0;i<m;++i) {
        const auto &row = this->RowView(i);
        const size_type nnz = row.NumNonZeros();
        for (ordinal_type j=0;j<nnz;++j)
          _A.Value(i, row.Col(j)) = row.Value(j);
      }
    }

    KOKKOS_INLINE_FUNCTION
    void copyToView() { 
      const ordinal_type m = this->NumRows();
      for (ordinal_type i=0;i<m;++i) {
        auto &row = this->RowView(i);
        const size_type nnz = row.NumNonZeros();
        for (ordinal_type j=0;j<nnz;++j)
          row.Value(j) = _A.Value(i, row.Col(j));
      }
    }

    KOKKOS_INLINE_FUNCTION
    CrsMatrixViewExt()
      : CrsMatrixView<CrsMatBaseType>(), _A()//, _H()
    { }

    template<typename MT>
    KOKKOS_INLINE_FUNCTION
    CrsMatrixViewExt(const CrsMatrixViewExt<MT> &b)
      : CrsMatrixView<MT>(b), _A(b._A)//, _H(b._H)
    { }

    KOKKOS_INLINE_FUNCTION
    CrsMatrixViewExt(const CrsMatBaseType &b)
      : CrsMatrixView<CrsMatBaseType>(b), _A()//, _H()
    { }

    KOKKOS_INLINE_FUNCTION
    CrsMatrixViewExt(const CrsMatBaseType &b,
                     const ordinal_type offm, const ordinal_type m,
                     const ordinal_type offn, const ordinal_type n)
      : CrsMatrixView<CrsMatBaseType>(b, offm, m, offn, n), _A()//, _H()
    { }

    KOKKOS_INLINE_FUNCTION
    ~CrsMatrixViewExt() = default;
  };
}

#endif
