#pragma once
#ifndef __CRS_MATRIX_VIEW_EXT_HPP__
#define __CRS_MATRIX_VIEW_EXT_HPP__

/// \file crs_matrix_view_ext.hpp
/// \brief Extended matrix view that has nested dense block.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

namespace Tacho {

  template<typename ValueType,
           typename OrdinalType,
           typename SizeType,
           typename SpaceType,
           typename MemoryTraits>
  class DenseMatrixBase;

  template<typename DenseMatBaseType>
  class DenseMatrixView;

  template<typename CrsMatViewType,
           typename DenseMatViewType>
  class CrsMatrixViewExt : public CrsMatViewType {
  public:
    typedef typename CrsMatViewType::space_type    space_type;
    typedef typename CrsMatViewType::memory_traits memory_traits;

    typedef typename CrsMatViewType::value_type    value_type;
    typedef typename CrsMatViewType::ordinal_type  ordinal_type;
    typedef typename CrsMatViewType::size_type     size_type;

    typedef typename DenseMatViewType::mat_base_type dense_mat_base_type;
    typedef DenseMatViewType dense_mat_view_type;

  private:
    dense_mat_base_type _A;

  public:
    bool hasDenseMatBase() const {
      return (_A.NumRows() && _A.NumCols());
    }

    bool isDenseMatBaseValid() const {
      return (_A.NumRows() >= this->NumRows() && _A.NumCols() >= this->NumCols());
    }

    void createDenseMatBase() {
      if (hasDenseMatBase() && isDenseMatBaseValid()) ;
      else
        _A = dense_mat_base_type("NestedDenseBase", this->NumRows(), this->NumCols());
    }

    dense_mat_base_type* DenseBaseObject() { return &_A; }

    int copyToDenseMatBase() {
      int r_val = 0;
      if (hasDenseMatBase() && isDenseMatBaseValid())  {
        const ordinal_type nrows = this->NumRows();
        for (ordinal_type i=0;i<nrows;++i) {
          auto row = this->RowView(i);
          const ordinal_type nnz = row.NumNonZeros();
          for (ordinal_type j=0;j<nnz;++j) 
            _A.Value(i, row.Col(j)) = row.Value(j);
        }
      } else {
        r_val = -1;
      }
      return r_val;
    }

    int copyToCrsMatrixView() {
      int r_val = 0;
      if (hasDenseMatBase() && isDenseMatBaseValid())  {
        const ordinal_type nrows = this->NumRows();
        for (ordinal_type i=0;i<nrows;++i) {
          auto row = this->RowView(i);
          const ordinal_type nnz = row.NumNonZeros();
          for (ordinal_type j=0;j<nnz;++j)
            row.Value(j) = _A.Value(i, row.Col(j));
        }
      } else {
        r_val = -1;
      }
      return r_val;
    }

    CrsMatrixViewExt()
      : CrsMatViewType(), _A()
    { }

    CrsMatrixViewExt(const CrsMatrixViewExt &b)
      : CrsMatViewType(b), _A(b._A)
    { }

    CrsMatrixViewExt(typename CrsMatViewType::mat_base_type *b)
      : CrsMatViewType(b), _A()
    { }

    CrsMatrixViewExt(typename CrsMatViewType::mat_base_type *b,
                 const ordinal_type offm, const ordinal_type m,
                 const ordinal_type offn, const ordinal_type n)
      : CrsMatViewType(b, offm, m, offn, n), _A()
    { }

  };
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
  namespace Impl {

    //  The Kokkos::View allocation will by default assign each allocated datum to zero.
    //  This is not the required initialization behavior when
    //  non-trivial objects are used within a Kokkos::View.
    //  Create a partial specialization of the Kokkos::Impl::AViewDefaultConstruct
    //  to replace the assignment initialization with placement new initialization.
    //
    //  This work-around is necessary until a TBD design refactorization of Kokkos::View.

    template< class ExecSpace , typename T1, typename T2 >
    struct ViewDefaultConstruct< ExecSpace , Tacho::CrsMatrixViewExt<T1,T2> , true >
    {
      typedef Tacho::CrsMatrixViewExt<T1,T2> type ;
      type * const m_ptr ;

      KOKKOS_FORCEINLINE_FUNCTION
      void operator()( const typename ExecSpace::size_type& i ) const
      { new(m_ptr+i) type(); }

      ViewDefaultConstruct( type * pointer , size_t capacity )
        : m_ptr( pointer )
      {
        Kokkos::RangePolicy< ExecSpace > range( 0 , capacity );
        parallel_for( range , *this );
        ExecSpace::fence();
      }
    };

  } // namespace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif
