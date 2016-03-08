#ifndef __TACHO_DENSE_MATRIX_BASE_HPP__
#define __TACHO_DENSE_MATRIX_BASE_HPP__

/// \file Tacho_DenseMatrixBase.hpp
/// \brief dense matrix base object interfaces 
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho { 

  /// \class DenseMatrixBase
  /// \breif Dense matrix base object using Kokkos view and subview
  template<typename ValueType,
           typename OrdinalType, 
           typename SizeType = OrdinalType,
           typename SpaceType = void>
  class DenseMatrixBase {
  public:
    typedef ValueType    value_type;
    typedef OrdinalType  ordinal_type;
    typedef SizeType     size_type;
    typedef SpaceType    space_type;

    // 1D view has some advantange in terms of workspace allocation
    typedef Kokkos::View<value_type*,space_type>   value_type_array;
    typedef Kokkos::View<ordinal_type*,space_type> ordinal_type_array;

    template<typename, typename, typename, typename>
    friend class DenseMatrixBase ;
    
  private:                          
    char               _label[Util::LabelSize]; //!< object label
    
    ordinal_type       _m;         //!< # of rows
    ordinal_type       _n;         //!< # of cols

    ordinal_type       _rs;        //!< row stride    
    ordinal_type       _cs;        //!< column stride
    
    value_type_array   _a;         //!< values
    
  protected:

    /// Callable: Device (o), KokkosFunctors (x)

    KOKKOS_INLINE_FUNCTION    
    void createInternalArrays(const ordinal_type m, 
                              const ordinal_type n,
                              const ordinal_type rs,
                              const ordinal_type cs,
                              const bool align = true) {
      // compute necessary storage and adjust strides
      const size_type value_type_size = sizeof(value_type);
      size_type size = 0;
      if (m == 0 || n == 0) { // empty
        _rs = 1;
        _cs = 1; 
        size = 0;
      } else if (rs == 1 || rs < cs) { // column major storage
        _rs = rs;
        _cs = align ? Util::alignDimension<space_type>(cs, value_type_size) : cs;
        size = _cs*n;
      } else if (cs == 1 || rs > cs) { // row major storage
        _rs = align ? Util::alignDimension<space_type>(rs, value_type_size) : rs;
        _cs = cs;
        size = m*_rs;
      } else {              // general storage
        _rs = rs;
        _cs = align ? Util::alignDimension<space_type>(cs, value_type_size) : cs;
        size = m*n*cs;
      }

      _m = m;
      _n = n;

      // grow buffer dimension
      if (static_cast<size_type>(_a.dimension_0()) < size) {
        char label[Util::LabelSize*2];
        strcat(label, _label); 
        strcat(label, "::ValueArray");
        _a = value_type_array(label, size);
      } else {
        // otherwise initialize it
        Kokkos::Impl::ViewFill<value_type_array>(_a, value_type());
      }
    }

    KOKKOS_INLINE_FUNCTION    
    void createInternalArrays(const ordinal_type m, 
                              const ordinal_type n) {
      createInternalArrays(m, n, 1, m);
    }

  public:

    /// Callable: Device (o), KokkosFunctors (o)
    
    KOKKOS_INLINE_FUNCTION    
    void setLabel(const char *label) { 
      strncpy(_label, label, Util::min(strlen(label)+1, Util::LabelSize));
    }
    
    KOKKOS_INLINE_FUNCTION
    char* Label() const { return _label; }
    
    KOKKOS_INLINE_FUNCTION
    ordinal_type NumRows() const { return _m; }
    
    KOKKOS_INLINE_FUNCTION
    ordinal_type NumCols() const { return _n; }
    
    KOKKOS_INLINE_FUNCTION
    ordinal_type ColStride() const { return _cs; }

    KOKKOS_INLINE_FUNCTION
    ordinal_type RowStride() const { return _rs; }

    KOKKOS_INLINE_FUNCTION
    value_type& Value(const ordinal_type i, 
                      const ordinal_type j) { 
      return _a[i*_rs + j*_cs]; 
    }

    KOKKOS_INLINE_FUNCTION
    value_type Value(const ordinal_type i, 
                     const ordinal_type j) const { 
      return _a[i*_rs + j*_cs]; 
    }

    // KOKKOS_INLINE_FUNCTION
    // value_type Value(const Kokkos::pair<ordinal_type,ordinal_type> idx) { 
    //   return _a[idx.first*_rs + idx.second*_cs];
    // }

    // KOKKOS_INLINE_FUNCTION
    // value_type Value(const Kokkos::pair<ordinal_type,ordinal_type> idx) const { 
    //   return _a[idx.first*_rs + idx.second*_cs];
    // }

    KOKKOS_INLINE_FUNCTION
    value_type* ValuePtr() const { return &_a[0]; }

    KOKKOS_INLINE_FUNCTION
    value_type_array ValueArray() const { return _a; }
    
    /// Callable: Device (o), KokkosFunctors (x)

    /// \brief Default constructor.
    KOKKOS_INLINE_FUNCTION
    DenseMatrixBase() 
      : _m(0),
        _n(0),
        _rs(0),
        _cs(0),
        _a()
    { 
      setLabel("DenseMatrixBase"); 
    }

    /// \brief Constructor with label
    //KOKKOS_INLINE_FUNCTION
    DenseMatrixBase(const char *label) 
      : _m(0),
        _n(0),
        _rs(0),
        _cs(0),
        _a()
    { 
      setLabel(label); 
    }

    /// \brief Copy constructor (shallow copy)
    /// Question : what would happen when SpT is different ? Does Kokkos can filter it out ? UVM allows this ?
    template<typename VT,
             typename OT,
             typename ST,
             typename SpT>
    KOKKOS_INLINE_FUNCTION
    DenseMatrixBase(const DenseMatrixBase<VT,OT,ST,SpT> &b) 
      : _m(b._m),
        _n(b._n),
        _rs(b._rs),
        _cs(b._cs),
        _a(b._a) 
    { 
      setLabel(b._label); 
    }
    
    /// \brief Constructor to allocate internal data structures.
    ///        By default, it uses the column oriented format.
    //KOKKOS_INLINE_FUNCTION
    DenseMatrixBase(const char *label,
                    const ordinal_type m, 
                    const ordinal_type n)
      : _m(m),
        _n(n)
    { 
      setLabel(label); 
      createInternalArrays(m, n);
    }
      
    /// \brief Constructor to attach external arrays to the matrix.
    ///        This is advanced constructor and trust user inputs
    ///        Todo :: later get an input of 2D container
    KOKKOS_INLINE_FUNCTION
    DenseMatrixBase(const char *label,
                    const ordinal_type m, 
                    const ordinal_type n,
                    const ordinal_type cs,
                    const ordinal_type rs,
                    const value_type_array &a) 
      : _m(m),
        _n(n),
        _rs(rs == -1 ? 1 : rs),
        _cs(cs == -1 ? m : cs),
        _a(a) 
    {
      setLabel(label); 
    }

    template<typename SpT>
    KOKKOS_INLINE_FUNCTION
    void 
    createConfTo(const DenseMatrixBase<value_type,ordinal_type,size_type,SpT> &b) {
      createInternalArrays(b._m, b._n, b._rs, b._cs, false);
    }
    
    /// \brief deep copy of matrix b
    /// Callable: Device (o), KokkosFunctors (x), Blocking (o)
    template<typename SpT>
    KOKKOS_INLINE_FUNCTION
    void 
    mirror(const DenseMatrixBase<value_type,ordinal_type,size_type,SpT> &b) {
      if (Kokkos::Impl::is_same<SpT,space_type>::value) {
        // when the space is same, everything is shallow copy 
        // setLabel(b._label);
        _m  = b._m;
        _n  = b._n;
        _rs = b._rs;
        _cs = b._cs;
        _a  = b._a;
      } else {
        // when the space is different, perform deep copy
        createInternalArrays(b.NumRows(), b.NumCols(), b.RowStride(), b.ColStride(), false);
        
        const auto range 
          = Kokkos::pair<ordinal_type,ordinal_type>(0, Util::min(_a.dimension_0(), b._a.dimension_0())); 
        
        space_type::execution_space::fence();      
        Kokkos::deep_copy(Kokkos::subview(_a, range), Kokkos::subview(b._a, range));
        space_type::execution_space::fence();
      }
    }

    /// \brief elementwise copy of matrix b
    /// Callable: Device (o), KokkosFunctors (o), Blocking (o)
    KOKKOS_INLINE_FUNCTION
    void
    copy(const DenseMatrixBase &b,
         const ordinal_type_array &ip = ordinal_type_array(),
         const ordinal_type_array &jp = ordinal_type_array()) { 
      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Static> > range_policy;
      
      const auto ip_dim = ip.dimension_0();
      const auto jp_dim = jp.dimension_0();

      space_type::execution_space::fence();      

      // loop fusion assumes column major format. 
      if (ip_dim && jp_dim) { 
        // row/col permutation
        Kokkos::parallel_for( range_policy(0, b._n), 
                              [&](const ordinal_type j) 
                              {
#pragma unroll 
                                for (auto i=0;i<b._m;++i)
                                  this->Value(i, j) = b.Value(ip(i), jp(j));
                              } );
      } else if (ip_dim) {
        // row permutation
        Kokkos::parallel_for( range_policy(0, b._n), 
                              [&](const ordinal_type j) 
                              {
#pragma unroll 
                                for (auto i=0;i<b._m;++i)
                                  this->Value(i, j) = b.Value(ip(i), j);
                              } );
      } else if (jp_dim) {
        // col permutation
        Kokkos::parallel_for( range_policy(0, b._n), 
                              [&](const ordinal_type j) 
                              {
                                const ordinal_type jj = jp(j);
#pragma unroll 
                                for (auto i=0;i<b._m;++i)
                                  this->Value(i, j) = b.Value(i, jj);
                              } );
      } else {
        // no permutation
        Kokkos::parallel_for( range_policy(0, b._n), 
                              [&](const ordinal_type j) 
                              {
#pragma unroll
                                for (auto i=0;i<b._m;++i)
                                  this->Value(i,j) = b.Value(i,j);
                              } );

        // // fused loop
        // Kokkos::View<Kokkos::pair<ordinal_type,ordinal_type>,space_type> idx;
        // Kokkos::parallel_for( range_policy(0, b._m*b._n), 
        //                       [=](const ordinal_type k) 
        //                       {
        //                         Util::unrollIndex(idx(), k, b._m);
        //                         this->Value(idx().first, idx().second) = b.Value(idx().first, idx().second);
        //                       } );
      }

      space_type::execution_space::fence();
    }

    /// \brief elementwise copy of lower/upper triangular of matrix b
    /// Callable: Device (o), KokkosFunctors (o)
    template<typename VT,
             typename OT,
             typename ST>
    KOKKOS_INLINE_FUNCTION
    void
    copy(const int uplo, 
         const DenseMatrixBase<VT,OT,ST,space_type> &b) { 

      typedef Kokkos::RangePolicy<space_type,Kokkos::Schedule<Kokkos::Dynamic> > range_policy;

      space_type::execution_space::fence();

      switch (uplo) {
      case Uplo::Lower: {
        Kokkos::parallel_for( range_policy(0, b._n), 
                              [&](const ordinal_type j) 
                              { 
#pragma unroll 
                                for (ordinal_type i=j;i<b._m;++i) 
                                  this->Value(i, j) = b.Value(i, j);
                              } );
        break;
      }
      case Uplo::Upper: {
        Kokkos::parallel_for( range_policy(0, b._n), 
                              [&](const ordinal_type j) 
                              { 
#pragma unroll 
                                for (ordinal_type i=0;i<(j+1);++i) 
                                  this->Value(i, j) = b.Value(i, j);
                              } );
        break;
      }
      }

      space_type::execution_space::fence();
    }

    /// \brief debugging output
    /// Callable: Device (x), KokkosFunctors (x)
    
    std::ostream& showMe(std::ostream &os) const {
      std::streamsize prec = os.precision();
      os.precision(8);
      os << std::scientific;

      os << " -- " << _label << " -- " << std::endl
         << "    # of Rows              = " << _m << std::endl
         << "    # of Cols              = " << _n << std::endl
         << "    Col Stride             = " << _cs << std::endl
         << "    Row Stride             = " << _rs << std::endl
         << std::endl
         << "    ValueArray dimensions  = " << _a.dimension_0() << std::endl
         << std::endl;
      
      const int w = 4;
      if (_a.size()) {
        for (ordinal_type i=0;i<_m;++i) {
          for (ordinal_type j=0;j<_n;++j) {
            const value_type val = this->Value(i,j);
            os << std::setw(w) << val << "  ";
          }
          os << std::endl;
        }
      }
      
      os.unsetf(std::ios::scientific);
      os.precision(prec);
      
      return os;
    }

    template<typename VT, typename OT, typename ST, typename SpT>
    friend std::ostream& std::operator<<(std::ostream &os, const DenseMatrixBase<VT,OT,ST,SpT> &self) {
      return self.showMe(os);
    }
  };
  
}

#endif
