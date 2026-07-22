// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_MP_VECTOR_MASKTRAITS_HPP
#define STOKHOS_MP_VECTOR_MASKTRAITS_HPP

#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include <iostream>
#include <cmath>

template <typename T>
struct EnsembleTraits_m {
  static const std::size_t size = 1;
  typedef T value_type;
  KOKKOS_INLINE_FUNCTION
  static const value_type& coeff(const T& x, int i) { return x; }
  KOKKOS_INLINE_FUNCTION
  static value_type& coeff(T& x, int i) { return x; }
};

template <typename S>
struct EnsembleTraits_m< Sacado::MP::Vector<S> > {
  static const std::size_t size = S::static_size ? S::static_size : 1;
  typedef typename S::value_type value_type;
  KOKKOS_INLINE_FUNCTION
  static const value_type& coeff(const Sacado::MP::Vector<S>& x, int i) {
    return x.fastAccessCoeff(i);
  }
  KOKKOS_INLINE_FUNCTION
  static value_type& coeff(Sacado::MP::Vector<S>& x, int i) {
    return x.fastAccessCoeff(i);
  }
};

template<typename scalar> class MaskedAssign
{
private:
    scalar &data;
    bool m;

public:
    KOKKOS_INLINE_FUNCTION MaskedAssign(scalar &data_, bool m_) : data(data_), m(m_) {};

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator = (const scalar & KOKKOS_RESTRICT s)
    {
        if(m)
            data = s;

        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator = (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();

        if(m)
            data = st_array[0];
        else
            data = st_array[1];
        
        return *this;
    }


    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator += (const scalar & KOKKOS_RESTRICT s)
    {
        if(m)
            data += s;

        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator += (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();

        if(m)
            data = st_array[0]+st_array[1];
        else
            data = st_array[2];
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator -= (const scalar & KOKKOS_RESTRICT s)
    {
        if(m)
            data -= s;

        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator -= (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();

        if(m)
            data = st_array[0]-st_array[1];
        else
            data = st_array[2];
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator *= (const scalar & KOKKOS_RESTRICT s)
    {
        if(m)
            data *= s;

        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator *= (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();

        if(m)
            data = st_array[0]*st_array[1];
        else
            data = st_array[2];
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator /= (const scalar & KOKKOS_RESTRICT s)
    {
        if(m)
            data /= s;

        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator /= (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();

        if(m)
            data = st_array[0]/st_array[1];
        else
            data = st_array[2];
        
        return *this;
    }
};

template<typename scalar> class Mask;

template<typename S> class MaskedAssign< Sacado::MP::Vector<S> >
{
    typedef Sacado::MP::Vector<S> scalar;
private:
    static const std::size_t size = EnsembleTraits_m<scalar>::size;
    scalar &data;
    Mask<scalar> m;

public:
    KOKKOS_INLINE_FUNCTION MaskedAssign(scalar &data_, Mask<scalar> m_) : data(data_), m(m_) {};

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator = (const scalar & KOKKOS_RESTRICT s)
    {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<size; ++i)
            if(m.get(i))
                data[i] = s[i];
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator = (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<size; ++i)
            if(m.get(i))
                data[i] = st_array[0][i];
            else
                data[i] = st_array[1][i];
        
        return *this;
    }


    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator += (const scalar & KOKKOS_RESTRICT s)
    {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<size; ++i)
            if(m.get(i))
                data[i] += s[i];
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator += (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<size; ++i)
            if(m.get(i))
                data[i] = st_array[0][i]+st_array[1][i];
            else
                data[i] = st_array[2][i];
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator -= (const scalar & KOKKOS_RESTRICT s)
    {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<size; ++i)
            if(m.get(i))
                data[i] -= s[i];
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator -= (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<size; ++i)
            if(m.get(i))
                data[i] = st_array[0][i]-st_array[1][i];
            else
                data[i] = st_array[2][i];
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator *= (const scalar & KOKKOS_RESTRICT s)
    {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<size; ++i)
            if(m.get(i))
                data[i] *= s[i];
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator *= (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<size; ++i)
            if(m.get(i))
                data[i] = st_array[0][i]*st_array[1][i];
            else
                data[i] = st_array[2][i];
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator /= (const scalar & KOKKOS_RESTRICT s)
    {
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<size; ++i)
            if(m.get(i))
                data[i] /= s[i];
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator /= (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for(std::size_t i=0; i<size; ++i)
            if(m.get(i))
                data[i] = st_array[0][i]/st_array[1][i];
            else
                data[i] = st_array[2][i];
        
        return *this;
    }
};

template <typename ordinal_t, typename value_t, typename device_t> class MaskedAssign< Sacado::MP::Vector<Stokhos::DynamicStorage<ordinal_t,value_t,device_t>> >
{
    typedef Sacado::MP::Vector<Stokhos::DynamicStorage<ordinal_t,value_t,device_t>> scalar;
private:
    static const std::size_t size = EnsembleTraits_m<scalar>::size;
    scalar &data;
    Mask<scalar> m;

public:
    KOKKOS_INLINE_FUNCTION MaskedAssign(scalar &data_, Mask<scalar> m_) : data(data_), m(m_) {};

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator = (const scalar & KOKKOS_RESTRICT s)
    {
        if(m.get(0))
            data = s;
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator = (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();
        
        if(m.get(0))
            data = st_array[0];
        else
            data = st_array[1];
        
        return *this;
    }


    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator += (const scalar & KOKKOS_RESTRICT s)
    {
        if(m.get(0))
            data += s;
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator += (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();

        if(m.get(0))
            data = st_array[0]+st_array[1];
        else
            data = st_array[2];
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator -= (const scalar & KOKKOS_RESTRICT s)
    {
        if(m.get(0))
            data -= s;
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator -= (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();

        if(m.get(0))
            data = st_array[0]-st_array[1];
        else
            data = st_array[2];
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator *= (const scalar & KOKKOS_RESTRICT s)
    {
        if(m.get(0))
            data *= s;
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator *= (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();

        if(m.get(0))
            data = st_array[0]*st_array[1];
        else
            data = st_array[2];
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator /= (const scalar & KOKKOS_RESTRICT s)
    {
        if(m.get(0))
            data /= s;
        
        return *this;
    }

    KOKKOS_INLINE_FUNCTION MaskedAssign<scalar>& operator /= (const std::initializer_list<scalar> & KOKKOS_RESTRICT st)
    {
        auto st_array = st.begin();

        if(m.get(0))
            data = st_array[0]/st_array[1];
        else
            data = st_array[2];
        
        return *this;
    }
};

template<typename scalar> class Mask
{
public:
    static const std::size_t size = EnsembleTraits_m<scalar>::size;

private:
    bool data[size];

public:
    KOKKOS_INLINE_FUNCTION Mask(){
        for(std::size_t i=0; i<size; ++i)
            this->set(i,false);
    }

    KOKKOS_INLINE_FUNCTION Mask(bool a){
        for(std::size_t i=0; i<size; ++i)
            this->set(i,a);
    }

    KOKKOS_INLINE_FUNCTION Mask(const Mask &a){
        for(std::size_t i=0; i<size; ++i)
            data[i] = a.data[i];
    }

    KOKKOS_INLINE_FUNCTION Mask& operator=(const Mask &a){
        for(std::size_t i=0; i<size; ++i)
            this->data[i] = a.data[i];
        return *this;
    }

    KOKKOS_INLINE_FUNCTION std::size_t getSize() const {return size;}

    KOKKOS_INLINE_FUNCTION bool operator> (double v)
    {
        double sum = 0;
        for(std::size_t i=0; i<size; ++i)
            sum = sum + this->get(i);

        return sum > v*size;
    }

    KOKKOS_INLINE_FUNCTION bool operator< (double v)
    {
        double sum = 0;
        for(std::size_t i=0; i<size; ++i)
            sum = sum + this->get(i);

        return sum < v*size;
    }

    KOKKOS_INLINE_FUNCTION bool operator>= (double v)
    {
        double sum = 0;
        for(std::size_t i=0; i<size; ++i)
            sum = sum + this->get(i);

        return sum >= v*size;
    }

    KOKKOS_INLINE_FUNCTION bool operator<= (double v)
    {
        double sum = 0;
        for(std::size_t i=0; i<size; ++i)
            sum = sum + this->get(i);

        return sum <= v*size;
    }

    KOKKOS_INLINE_FUNCTION bool operator== (double v)
    {
        double sum = 0;
        for(std::size_t i=0; i<size; ++i)
            sum = sum + this->get(i);

        return sum == v*size;
    }

    KOKKOS_INLINE_FUNCTION bool operator!= (double v)
    {
        double sum = 0;
        for(std::size_t i=0; i<size; ++i)
            sum = sum + this->get(i);

        return sum != v*size;
    }

    KOKKOS_INLINE_FUNCTION bool operator== (const Mask<scalar> &m2)
    {
        bool all = true;
        for (std::size_t i = 0; i < size; ++i) {
            all = all && (this->get(i) == m2.get(i));
        }
        return all;
    }

    KOKKOS_INLINE_FUNCTION bool operator!= (const Mask<scalar> &m2)
    {
        return !(this==m2);
    }
  
    KOKKOS_INLINE_FUNCTION Mask<scalar> operator&& (const Mask<scalar> &m2)
    {
        Mask<scalar> m3;
        for(std::size_t i=0; i<size; ++i)
            m3.set(i,(this->get(i) && m2.get(i)));
        return m3;
    }
      
    KOKKOS_INLINE_FUNCTION Mask<scalar> operator|| (const Mask<scalar> &m2)
    {
        Mask<scalar> m3;
        for(std::size_t i=0; i<size; ++i)
            m3.set(i,(this->get(i) || m2.get(i)));

        return m3;
    }

    KOKKOS_INLINE_FUNCTION Mask<scalar> operator&& (bool m2)
    {
        Mask<scalar> m3;
        for(std::size_t i=0; i<size; ++i)
            m3.set(i,(this->get(i) && m2));
        return m3;
    }

    KOKKOS_INLINE_FUNCTION Mask<scalar> operator|| (bool m2)
    {
        Mask<scalar> m3;
        for(std::size_t i=0; i<size; ++i)
            m3.set(i,(this->get(i) || m2));

        return m3;
    }

    KOKKOS_INLINE_FUNCTION Mask<scalar> operator+ (const Mask<scalar> &m2)
    {
        Mask<scalar> m3;
        for(std::size_t i=0; i<size; ++i)
            m3.set(i,(this->get(i) + m2.get(i)));

        return m3;
    }

    KOKKOS_INLINE_FUNCTION Mask<scalar> operator- (const Mask<scalar> &m2)
    {
        Mask<scalar> m3;
        for(std::size_t i=0; i<size; ++i)
            m3.set(i,(this->get(i) - m2.get(i)));

        return m3;
    }

    KOKKOS_INLINE_FUNCTION scalar operator* (const scalar &v)
    {
        typedef EnsembleTraits_m<scalar> ET;
        scalar v2;
        for(std::size_t i=0; i<size; ++i)
            ET::coeff(v2,i) = ET::coeff(v,i)*this->get(i);

        return v2;
    }
    
    KOKKOS_INLINE_FUNCTION bool get (int i) const
    {
        return this->data[i];
    }

    KOKKOS_INLINE_FUNCTION void set (int i, bool b)
    {
        this->data[i] = b;
    }
    
    KOKKOS_INLINE_FUNCTION Mask<scalar> operator! ()
    {
        Mask<scalar> m2;
        for(std::size_t i=0; i<size; ++i)
            m2.set(i,!(this->get(i)));
        return m2;
    }
    
    KOKKOS_INLINE_FUNCTION operator bool() const
    {
        return this->get(0);
    }

    KOKKOS_INLINE_FUNCTION operator double() const
    {
        double sum = 0;
        for(std::size_t i=0; i<size; ++i)
            sum = sum + this->get(i);

        return sum/size;
    }
};

template<typename scalar> std::ostream &operator<<(std::ostream &os, const Mask<scalar>& m) {
    os << "[ ";
    for(std::size_t i=0; i<m.getSize(); ++i)
        os << m.get(i) << " ";
    return os << "]";
}

template<typename S> KOKKOS_INLINE_FUNCTION Sacado::MP::Vector<S> operator* (const Sacado::MP::Vector<S> &a1, const Mask<Sacado::MP::Vector<S>> &m)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S>> ET;
    Sacado::MP::Vector<S> mul;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for(std::size_t i=0; i<ET::size; ++i){
        ET::coeff(mul,i) = ET::coeff(a1,i)*m.get(i);
    }
    return mul;
}

template<typename S> KOKKOS_INLINE_FUNCTION Sacado::MP::Vector<S> operator* (const typename S::value_type &a1, const Mask<Sacado::MP::Vector<S>> &m)
{
    Sacado::MP::Vector<S> mul;
    typedef EnsembleTraits_m<Sacado::MP::Vector<S>> ET;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for(std::size_t i=0; i<ET::size; ++i){
        ET::coeff(mul,i) = m.get(i)*a1;
    }
    return mul;
}

template<typename S> KOKKOS_INLINE_FUNCTION Sacado::MP::Vector<S> operator* (const Mask<Sacado::MP::Vector<S>> &m, const typename S::value_type &a1)
{
    Sacado::MP::Vector<S> mul;
    typedef EnsembleTraits_m<Sacado::MP::Vector<S>> ET;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for(std::size_t i=0; i<ET::size; ++i){
        ET::coeff(mul,i) = m.get(i)*a1;
    }
    return mul;
}

template<typename scalar> KOKKOS_INLINE_FUNCTION MaskedAssign<scalar> mask_assign(bool b, scalar *s)
{
    Mask<scalar> m = Mask<scalar>(b);
    MaskedAssign<scalar> maskedAssign = MaskedAssign<scalar>(*s,m);
    return maskedAssign;
}

template<typename scalar> KOKKOS_INLINE_FUNCTION MaskedAssign<scalar> mask_assign(Mask<scalar> m, scalar *s)
{
    MaskedAssign<scalar> maskedAssign = MaskedAssign<scalar>(*s,m);
    return maskedAssign;
}

template<typename scalar> KOKKOS_INLINE_FUNCTION MaskedAssign<scalar> mask_assign(bool b, scalar &s)
{
    Mask<scalar> m = Mask<scalar>(b);
    MaskedAssign<scalar> maskedAssign = MaskedAssign<scalar>(s,m);
    return maskedAssign;
}

template<typename scalar> KOKKOS_INLINE_FUNCTION MaskedAssign<scalar> mask_assign(Mask<scalar> m, scalar &s)
{
    MaskedAssign<scalar> maskedAssign = MaskedAssign<scalar>(s,m);
    return maskedAssign;
}

namespace Sacado {
    namespace MP {
        template <typename S> KOKKOS_INLINE_FUNCTION Vector<S> copysign(const Vector<S> &a1, const Vector<S> &a2)
        {
            typedef EnsembleTraits_m< Vector<S> > ET;
            
            Vector<S> a_out;
            
            using std::copysign;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
            for(std::size_t i=0; i<ET::size; ++i){
                ET::coeff(a_out,i) = copysign(ET::coeff(a1,i),ET::coeff(a2,i));
            }
            
            return a_out;
        }
    }
}


template<typename S> KOKKOS_INLINE_FUNCTION Mask<Sacado::MP::Vector<S> > signbit_v(const Sacado::MP::Vector<S> &a1)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S> > ET;
    using std::signbit;
    
    Mask<Sacado::MP::Vector<S> > mask;
#ifdef STOKHOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef STOKHOS_HAVE_PRAGMA_VECTOR_ALIGNED
#pragma vector aligned
#endif
#ifdef STOKHOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for(std::size_t i=0; i<ET::size; ++i)
        mask.set(i, signbit(ET::coeff(a1,i)));
    return mask;
}

// Relation operations for vector:

#define OPNAME ==
#include "Stokhos_MP_Vector_MaskTraits_vector_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME !=
#include "Stokhos_MP_Vector_MaskTraits_vector_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME >
#include "Stokhos_MP_Vector_MaskTraits_vector_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME >=
#include "Stokhos_MP_Vector_MaskTraits_vector_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME <
#include "Stokhos_MP_Vector_MaskTraits_vector_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME <=
#include "Stokhos_MP_Vector_MaskTraits_vector_relops_tmpl.hpp"
#undef OPNAME

// Relation operations for expressions:

#define OPNAME ==
#include "Stokhos_MP_Vector_MaskTraits_expr_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME !=
#include "Stokhos_MP_Vector_MaskTraits_expr_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME <
#include "Stokhos_MP_Vector_MaskTraits_expr_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME >
#include "Stokhos_MP_Vector_MaskTraits_expr_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME <=
#include "Stokhos_MP_Vector_MaskTraits_expr_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME >=
#include "Stokhos_MP_Vector_MaskTraits_expr_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME <<=
#include "Stokhos_MP_Vector_MaskTraits_expr_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME >>=
#include "Stokhos_MP_Vector_MaskTraits_expr_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME &
#include "Stokhos_MP_Vector_MaskTraits_expr_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME |
#include "Stokhos_MP_Vector_MaskTraits_expr_relops_tmpl.hpp"
#undef OPNAME

#if STOKHOS_USE_MP_VECTOR_SFS_SPEC

// Relation operations for static fixed storage:

#define OPNAME ==
#include "Stokhos_MP_Vector_MaskTraits_sfs_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME !=
#include "Stokhos_MP_Vector_MaskTraits_sfs_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME >
#include "Stokhos_MP_Vector_MaskTraits_sfs_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME >=
#include "Stokhos_MP_Vector_MaskTraits_sfs_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME <
#include "Stokhos_MP_Vector_MaskTraits_sfs_relops_tmpl.hpp"
#undef OPNAME

#define OPNAME <=
#include "Stokhos_MP_Vector_MaskTraits_sfs_relops_tmpl.hpp"
#undef OPNAME

#endif

namespace MaskLogic{

    template<typename T> KOKKOS_INLINE_FUNCTION bool OR(Mask<T> m){
        return (((double) m)!=0.);
    }

    KOKKOS_INLINE_FUNCTION bool OR(bool m){
        return m;
    }

    template<typename T> KOKKOS_INLINE_FUNCTION bool XOR(Mask<T> m){
        return (((double) m)==1./m.getSize());
    }

    KOKKOS_INLINE_FUNCTION bool XOR(bool m){
        return m;
    }

    template<typename T> KOKKOS_INLINE_FUNCTION bool AND(Mask<T> m){
        return (((double) m)==1.);
    }

    KOKKOS_INLINE_FUNCTION bool AND(bool m){
        return m;
    }

}

#endif // STOKHOS_MP_VECTOR_MASKTRAITS_HPP
