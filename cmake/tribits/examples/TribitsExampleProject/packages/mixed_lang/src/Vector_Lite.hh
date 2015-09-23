//---------------------------------------------------------------------------//
// Vector_Lite.hh
//---------------------------------------------------------------------------//

#ifndef mixed_language_src_Vector_Lite_hh
#define mixed_language_src_Vector_Lite_hh

#include <iostream>
#include <numeric>
#include <cstddef>

namespace tribits_mixed
{

//---------------------------------------------------------------------------//

template <class T, size_t N>
class Vector_Lite
{
  public:
    //@{
    //! Typedefs.
    typedef T                 value_type;
    typedef value_type*       pointer;
    typedef const T*          const_pointer;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef ptrdiff_t         difference_type;
    typedef size_t            size_type;
    typedef pointer           iterator;
    typedef const_pointer     const_iterator;
    //@}
    
  private:
    // >>> DATA

    T d_U[N];
    
  public:
    // Constructor based on a scalar value.
    inline explicit Vector_Lite(const T &u = T());

    // Constructor for N = 2.
    inline Vector_Lite(const T &u0, const T &u1);

    // Constructor for N = 3.
    inline Vector_Lite(const T &u0, const T &u1, const T &u2);

    // Constructor for N = 4.
    inline Vector_Lite(const T &u0, const T &u1, const T &u2, const T &u3);

    // Constructor for N = 5.
    inline Vector_Lite(const T &u0, const T &u1, const T &u2,
                       const T &u3, const T &u4);

    // Fill in from C array.
    inline void fill(const T u[N]);

    //! Destructor.
    ~Vector_Lite(void) { }

    // >>> MANIPULATORS

    // Assignment to a another Vector_Lite.
    inline Vector_Lite &operator=(const Vector_Lite &rhs);

    // Assignment to a scalar.
    inline Vector_Lite &operator=(const T &rhs);

    // Comparisons to another Vector_Lite.
    inline bool operator==(const Vector_Lite &a) const;
    inline bool operator<(const Vector_Lite &a) const;

    // Basic arithmetic operations, vector right-hand side.
    inline Vector_Lite &operator+=(const Vector_Lite &a);
    inline Vector_Lite &operator-=(const Vector_Lite &a);
    inline Vector_Lite &operator*=(const Vector_Lite &a);
    inline Vector_Lite &operator/=(const Vector_Lite &a);

    // Basic arithmetic operations, scalar right-hand side.
    inline Vector_Lite &operator+=(const T &a);
    inline Vector_Lite &operator-=(const T &a);
    inline Vector_Lite &operator*=(const T &a);
    inline Vector_Lite &operator/=(const T &a);
    
    //! Returns true if \a i is a valid array index.
    bool valid_index(const size_type i) const
    {
        return i < N;
    }

    // >>> ACCESSORS

    //! Indexing using ().
    reference operator()(const size_type i)
    {
        return d_U[i];
    }
    
    //! Const indexing using ().
    const_reference operator()(const size_type i) const
    {
        return d_U[i];
    }
    
    //! Indexing using [].
    reference operator[](const size_type i)
    {
        return d_U[i];
    }
    
    //! const indexing using [].
    const_reference operator[](const size_type i) const
    {
        return d_U[i];
    }

    // >>> ITERATOR SUPPORT
    
    //! Iterator begin.
    iterator begin() { return d_U; }

    //! Const iterator begin.
    const_iterator begin() const { return d_U; }

    //! Iterator end.
    iterator end() { return d_U + N; }

    //! Const iterator end.
    const_iterator end() const { return d_U + N; }

    //! Number of elements (\a N); for STL support.
    size_type size() const { return N; }

    //! Max number of elements (\a N); for STL support.
    size_type max_size() const { return N; }

    //! True if \a N = 0; for STL support.
    bool empty() const { return N == 0; }
};

//---------------------------------------------------------------------------//
// MEMBER FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor based on a scalar value.
 *
 * Initializes all values to \a u.  Note that this ctor also acts as the
 * default ctor.
 * 
 * \param u  Scalar value.
 */
template <class T, size_t N>
Vector_Lite<T, N>::Vector_Lite(const T &u)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] = u;
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for N = 2.
 * 
 * \param u0 1st element.    
 * \param u1 2nd element.
 */
template <class T, size_t N>
Vector_Lite<T, N>::Vector_Lite(const T &u0,
                               const T &u1)
{
    d_U[0] = u0;
    d_U[1] = u1;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for N = 3.
 * 
 * \param u0 1st element.    
 * \param u1 2nd element.
 * \param u2 3rd element.
 */
template <class T, size_t N>
Vector_Lite<T, N>::Vector_Lite(const T &u0,
                               const T &u1,
                               const T &u2)
{
    d_U[0] = u0;
    d_U[1] = u1;
    d_U[2] = u2;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for N = 4.
 * 
 * \param u0 1st element.    
 * \param u1 2nd element.
 * \param u2 3rd element.
 * \param u3 4th element.
 */
template <class T, size_t N>
Vector_Lite<T, N>::Vector_Lite(const T &u0,
                               const T &u1,
                               const T &u2,
                               const T &u3)
{
    d_U[0] = u0;
    d_U[1] = u1;
    d_U[2] = u2;
    d_U[3] = u3;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Constructor for N = 5.
 * 
 * \param u0 1st element.    
 * \param u1 2nd element.
 * \param u2 3rd element.
 * \param u3 4th element.
 * \param u4 5th element.
 */
template <class T, size_t N>
Vector_Lite<T, N>::Vector_Lite(const T &u0,
                               const T &u1,
                               const T &u2,
                               const T &u3,
                               const T &u4)
{
    d_U[0] = u0;
    d_U[1] = u1;
    d_U[2] = u2;
    d_U[3] = u3;
    d_U[4] = u4;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Fill in from C array.
 * 
 * \param u c-style pointer to array of size N.
 *
 */
template <class T, size_t N>
void Vector_Lite<T, N>::fill(const T u[N])
{
    for (size_t i=0; i<N; ++i) 
    {
        d_U[i] = u[i];
    }
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment to another Vector_Lite.
 */
template <class T, size_t N>
Vector_Lite<T, N>& Vector_Lite<T, N>::operator=(const Vector_Lite &rhs)
{
    // return if objects are the same
    if (&rhs == this)
        return *this;

    // otherwise do member assignment
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] = rhs.d_U[i];
    }
	
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Assignment to a scalar.
 */
template <class T, size_t N>
Vector_Lite<T, N>& Vector_Lite<T, N>::operator=(const T &rhs)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] = rhs;
    }
	
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Comparison to another Vector_Lite.
 */
template <class T, size_t N>
bool Vector_Lite<T, N>::operator==(const Vector_Lite<T, N> &a) const
{
    for ( size_type i = 0; i < N; ++i )
    {
        if ( d_U[i] != a.d_U[i] ) return false;
    }
    
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Comparison to another Vector_Lite.
 */
template <class T, size_t N>
bool Vector_Lite<T, N>::operator<(const Vector_Lite<T, N> &a) const
{
    for ( size_type i = 0; i < N; ++i )
    {
        if ( ! ( d_U[i] <  a.d_U[i] ) ) return false;
    }
    
    return true;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Element-wise not-equals.
 */
template<class T, size_t N>
inline bool operator!=(Vector_Lite<T,N> const & lhs, 
                       Vector_Lite<T,N> const & rhs)
{
    return !(lhs == rhs);
}

//---------------------------------------------------------------------------//
/*!
 * \brief Element-wise greater than. 
 */
template<class T, size_t N>
inline bool operator>(Vector_Lite<T,N> const & lhs, 
                      Vector_Lite<T,N> const & rhs)
{
    return rhs < lhs;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Element-wise less-than-or-equal.
 */
template<class T, size_t N>
inline bool operator<=(Vector_Lite<T,N> const & lhs, 
                       Vector_Lite<T,N> const & rhs)
{
    return !(lhs > rhs);
}

//---------------------------------------------------------------------------//
/*!
 * \brief  Element-wise greater-than-or-equal.
 */
template<class T, size_t N> 
inline bool operator>=(Vector_Lite<T,N> const & lhs, 
                       Vector_Lite<T,N> const & rhs)
{
    return !(rhs > lhs);
}

//---------------------------------------------------------------------------//
// BASIC ARITHMETIC MEMBER FUNCTIONS, VECTOR RIGHT-HAND SIDE
//---------------------------------------------------------------------------//
/*!
 * \brief Support for +=, Vector_Lite right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator+=(const Vector_Lite<T, N> &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] += a.d_U[i];
    }
    
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Support for -=, Vector_Lite right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator-=(const Vector_Lite<T, N> &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] -= a.d_U[i];
    }
    
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Support for *=, Vector_Lite right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator*=(const Vector_Lite<T, N> &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] *= a.d_U[i];
    }
    
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Support for /=, Vector_Lite right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator/=(const Vector_Lite<T, N> &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] /= a.d_U[i];
    }
    
    return *this;
}

//---------------------------------------------------------------------------//
// BASIC ARITHMETIC MEMBER FUNCTIONS, SCALAR RIGHT-HAND SIDE
//---------------------------------------------------------------------------//
/*!
 * \brief Support for +=, scalar right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator+=(const T &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] += a;
    }
    
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Support for -=, scalar right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator-=(const T &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] -= a;
    }
    
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Support for *=, scalar right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator*=(const T &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] *= a;
    }
    
    return *this;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Support for /=, scalar right-hand side.
 */
template <class T, size_t N>
Vector_Lite<T,N>& Vector_Lite<T, N>::operator/=(const T &a)
{
    for ( size_type i = 0; i < N; ++i )
    {
        d_U[i] /= a;
    }
    
    return *this;
}

//---------------------------------------------------------------------------//
// FREE FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Computes the inner product between two vectors.
 * 
 * \param a 1st vector.
 * \param b 2nd vector.
 */
template <class T, size_t N>
T inner_product(const Vector_Lite<T, N> &a,
                const Vector_Lite<T, N> &b)
{
    return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

//---------------------------------------------------------------------------//
// GLOBAL OPERATORS
//---------------------------------------------------------------------------//
/*!
 * \brief \a a + \a b, element by element.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator+(const Vector_Lite<T, N> &a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(a) += b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a - \a b, element by element.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator-(const Vector_Lite<T, N> &a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(a) -= b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a * \a b, element by element.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator*(const Vector_Lite<T, N> &a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(a) *= b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a / \a b, element by element.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator/(const Vector_Lite<T, N> &a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(a) /= b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a b added to all elements of \a a.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator+(const Vector_Lite<T, N> &a,
                                         const T                  b)
{
    return Vector_Lite<T, N>(a) += b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a added to all elements of \a b.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator+(const T                  a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(b) += a;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a b subracted from all elements of \a a.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator-(const Vector_Lite<T, N> &a,
                                         const T                  b)
{
    return Vector_Lite<T, N>(a) -= b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a subtracted from all elements of \a b.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator-(const T                  a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(b) -= a;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a b multiplied with all elements of \a a.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator*(const Vector_Lite<T, N> &a,
                                         const T                  b)
{
    return Vector_Lite<T, N>(a) *= b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a a multiplied with all elements of \a b.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator*(const T                  a,
                                         const Vector_Lite<T, N> &b)
{
    return Vector_Lite<T, N>(b) *= a;
}

//---------------------------------------------------------------------------//
/*!
 * \brief \a b divided into all elements of \a a.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator/(const Vector_Lite<T, N> &a,
                                         const T                  b)
{
    return Vector_Lite<T, N>(a) /= b;
}

//---------------------------------------------------------------------------//
/*!
 * \brief negates all of \a a.
 */
template <class T, size_t N>
inline const Vector_Lite<T, N> operator-(const Vector_Lite<T, N> &a)
{
    Vector_Lite<T, N> neg(a);
    
    for (size_t i = 0; i < N; ++i )
        neg(i) = -a(i);
    
    return neg;
}

//---------------------------------------------------------------------------//
// STREAM OPEATORS
//---------------------------------------------------------------------------//
/*!
 * \brief Write the elements of \a a to stream \a os.
 */
template <class T, size_t N>
std::ostream &operator<<(std::ostream            &os,
                         const Vector_Lite<T, N> &a)
{
    os << a(0);
    
    for (size_t i = 1; i < N; ++i) 
    {
        os << " " << a(i);
    }
    
    return os;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Read elements into \a a from \a is.
 */
template <class T, size_t N>
std::istream &operator>>(std::istream      &is,
                         Vector_Lite<T, N> &a)
{
    for (size_t i = 0; i < N; ++i) 
    {
        is >> a(i);
    }
    
    return is;
}

} // namespace tribits_mixed

#endif

//---------------------------------------------------------------------------//
// end of Vector_Lite.hh
//---------------------------------------------------------------------------//

