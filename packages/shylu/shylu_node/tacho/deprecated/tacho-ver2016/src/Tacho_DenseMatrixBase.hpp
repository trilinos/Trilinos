#ifndef __TACHO_DENSE_MATRIX_BASE_HPP__
#define __TACHO_DENSE_MATRIX_BASE_HPP__

/// \file Tacho_DenseMatrixBase.hpp
/// \brief dense matrix base object.
/// \author Kyungjoo Kim (kyukim@sandia.gov)

#include "Tacho_Util.hpp"

namespace Tacho { 

  /// \class DenseMatrixBase
  /// \breif Dense matrix base object using Kokkos view. 
  ///        The base object hold actual matrix storage (1D) and 
  ///        is responsible for mirroring to a device.
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

    // range type
    template<typename T> using range_type = Kokkos::pair<T,T>;

    template<typename, typename, typename, typename>
    friend class DenseMatrixBase;

    //friend class DenseMatrixTools;

  private:                          
    char               _label[Util::LabelSize]; //!< object label
    
    ordinal_type       _m;         //!< # of rows
    ordinal_type       _n;         //!< # of cols

    ordinal_type       _rs;        //!< row stride    
    ordinal_type       _cs;        //!< column stride
    
    value_type_array   _a;         //!< values
    
  protected:

    /// Properties: 
    /// - Compile with Device (o), 
    /// - Callable in KokkosFunctors (x)
    inline
    void createInternalArrays(const ordinal_type m, 
                              const ordinal_type n,
                              const ordinal_type rs,
                              const ordinal_type cs) {
      // compute necessary storage and adjust strides
      size_type size = 0;
      if (m == 0 || n == 0) { // empty
        _rs = 1;
        _cs = 1; 
        size = 0;
      } else if (rs == 1 || rs < cs) { // column major storage
        _rs = rs;
        _cs = cs;
        size = _cs*n;
      } else if (cs == 1 || rs > cs) { // row major storage
        _rs = rs;
        _cs = cs;
        size = m*_rs;
      } else {              // general storage
        _rs = rs;
        _cs = cs;
        size = m*n*cs;
      }

      _m = m;
      _n = n;

      // grow buffer dimension
      if (static_cast<size_type>(_a.extent(0)) < size) {
        _a = value_type_array("DenseMatrixBase::ValueArray", size);
      } else {
        // otherwise initialize it
        Kokkos::Impl::ViewFill<value_type_array>(_a, value_type());
      }
    }

    /// Properties: 
    /// - Compile with Device (o), 
    /// - Callable in KokkosFunctors (x)
    inline
    void createInternalArrays(const ordinal_type m, 
                              const ordinal_type n) {
      createInternalArrays(m, n, 1, m);
    }

  public:

    /// Interface functions
    /// ------------------------------------------------------------------
    /// Properties: 
    /// - Compile with Device (o), 
    /// - Callable in KokkosFunctors (o)

    KOKKOS_INLINE_FUNCTION
    void setExternalMatrix(const ordinal_type m, 
                           const ordinal_type n,
                           const ordinal_type cs,
                           const ordinal_type rs,
                           const value_type_array &a) {
      _m = m;
      _n = n; 
      _rs = (rs == -1 ? 1 : rs);
      _cs = (cs == -1 ? m : cs);
      _a = a;
    }

    KOKKOS_INLINE_FUNCTION    
    bool isValueArrayNull() const {
      return !_a.extent(0);
    }

    KOKKOS_INLINE_FUNCTION    
    void setLabel(const char *label) { 
      int i = 0 ;
      for ( ; i < int(Util::LabelSize - 1) && label[i] ; ++i )
        _label[i] = label[i] ;
      _label[i] = 0 ;
    }
    
    KOKKOS_INLINE_FUNCTION
    const char* Label() const { return _label; }
    
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
                      const ordinal_type j) const { 
      return _a[i*_rs + j*_cs]; 
    }

    // KOKKOS_INLINE_FUNCTION
    // value_type Value(const ordinal_type i, 
    //                 const ordinal_type j) const { 
    //  return _a[i*_rs + j*_cs]; 
    // }

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
    value_type_array Values() const { return _a; }

    /// ------------------------------------------------------------------
    
    /// Constructors
    /// ------------------------------------------------------------------
    /// Properties: 
    /// - Compile with Device (o), 
    /// - Callable in KokkosFunctors 
    ///   - Default and copy constructors are allowed in KokkosFunctors.
    ///   - Creating internal workspace is not allowed in KokkosFunctors.

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
    KOKKOS_INLINE_FUNCTION
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
    KOKKOS_INLINE_FUNCTION
    DenseMatrixBase(const DenseMatrixBase &b) 
      : _m(b._m),
        _n(b._n),
        _rs(b._rs),
        _cs(b._cs),
        _a(b._a) 
    { 
      // need static assert to evaluate space type
      setLabel(b._label); 
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

    /// \brief Constructor to allocate internal data structures.
    ///        By default, it uses the column oriented format.
    DenseMatrixBase(const char *label,
                    const ordinal_type m, 
                    const ordinal_type n)
      : _m(m),
        _n(n)
    { 
      setLabel(label); 
      createInternalArrays(m, n);
    }

    /// ------------------------------------------------------------------

    /// Destructor
    /// ------------------------------------------------------------------
    /// Properties:
    /// - Compile with Device (o),
    /// - Callable in KokkosFunctors (x)
    KOKKOS_INLINE_FUNCTION
    ~DenseMatrixBase() = default;

    /// Create and mirror
    /// ------------------------------------------------------------------
    /// Properties: 
    /// - Compile with Device (o), 
    /// - Callable in KokkosFunctors (x) 

    inline
    void 
    create(const ordinal_type m, 
           const ordinal_type n) {
      createInternalArrays(m, n);
    }
    
    template<typename SpT>
    inline
    void 
    createConfTo(const DenseMatrixBase<value_type,ordinal_type,size_type,SpT> &b) {
      createInternalArrays(b._m, b._n, b._rs, b._cs);
    }
    
    /// \brief deep copy of matrix b
    /// Callable: Device (o), KokkosFunctors (x), Blocking (o)
    template<typename SpT>
    inline
    typename std::enable_if< Kokkos::Impl::is_same<SpT,space_type>::value >::type
    mirror(const DenseMatrixBase<value_type,ordinal_type,size_type,SpT> &b) {
        // when the space is same, everything is shallow copy 
        // setLabel(b._label);
        _m  = b._m;
        _n  = b._n;
        _rs = b._rs;
        _cs = b._cs;
        _a  = b._a;
    }

    template<typename SpT>
    inline
    typename std::enable_if< ! Kokkos::Impl::is_same<SpT,space_type>::value >::type
    mirror(const DenseMatrixBase<value_type,ordinal_type,size_type,SpT> &b) {
        // when the space is different, perform deep copy
        createInternalArrays(b._m, b._n, b._rs, b._cs);
        
        const auto range = range_type<ordinal_type>(0, Util::min(_a.extent(0), b._a.extent(0))); 
        
        space_type::execution_space::fence();      
        Kokkos::deep_copy(Kokkos::subview(_a, range), Kokkos::subview(b._a, range));
        space_type::execution_space::fence();
    }

    /// ------------------------------------------------------------------


    /// Print out 
    /// ------------------------------------------------------------------
    /// Properties: 
    /// - Compile with Device (x), 
    /// - Callable in KokkosFunctors (x) 

    /// \brief print out to stream
    std::ostream& showMe(std::ostream &os) const {
      std::streamsize prec = os.precision();
      os.precision(16);
      os << std::scientific;

      os << " -- " << _label << " -- " << std::endl
         << "    # of Rows              = " << _m << std::endl
         << "    # of Cols              = " << _n << std::endl
         << "    Col Stride             = " << _cs << std::endl
         << "    Row Stride             = " << _rs << std::endl
         << std::endl
         << "    ValueArray dimensions  = " << _a.extent(0) << std::endl
         << std::endl;

      const int w = 4;
      if (_a.size()) {
        for (ordinal_type i=0;i<_m;++i) {
          for (ordinal_type j=0;j<_n;++j) {
            const value_type val = this->Value(i,j);
            os << std::setw(w) << std::showpos << val << std::noshowpos << "  ";
          }
          os << std::endl;
        }
      }

      os.unsetf(std::ios::scientific);
      os.precision(prec);
      
      return os;
    }

    /// ------------------------------------------------------------------
  };


}

#endif
