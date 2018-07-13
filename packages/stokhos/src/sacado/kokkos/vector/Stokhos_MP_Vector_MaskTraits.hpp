// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_MP_VECTOR_MASKTRAITS_HPP
#define STOKHOS_MP_VECTOR_MASKTRAITS_HPP

#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include <iostream>
#include <cmath>

template <typename T>
struct EnsembleTraits_m {
  static const int size = 1;
  typedef T value_type;
  static const value_type& coeff(const T& x, int i) { return x; }
  static value_type& coeff(T& x, int i) { return x; }
};

template <typename S>
struct EnsembleTraits_m< Sacado::MP::Vector<S> > {
  static const int size = S::static_size;
  typedef typename S::value_type value_type;
  static const value_type& coeff(const Sacado::MP::Vector<S>& x, int i) {
    return x.fastAccessCoeff(i);
  }
  static value_type& coeff(Sacado::MP::Vector<S>& x, int i) {
    return x.fastAccessCoeff(i);
  }
};

template<typename scalar> class Mask 
{
    private: 
        static const int size = EnsembleTraits_m<scalar>::size;
        bool data[size];

    public:
        Mask(){
            for(int i=0; i<size; ++i)
                data[i]=false;
        }

        Mask(bool a){
            for(int i=0; i<size; ++i)
                data[i]=a;
        }

        int getSize() const {return size;}

        bool operator> (double v)
        {
            double sum = 0;
            for(int i=0; i<size; ++i)
                sum = sum + this->data[i];

            return sum > v*size;
        }

        bool operator< (double v)
        {
            double sum = 0;
            for(int i=0; i<size; ++i)
                sum = sum + this->data[i];

            return sum < v*size;
        }

        bool operator>= (double v)
        {
            double sum = 0;
            for(int i=0; i<size; ++i)
                sum = sum + this->data[i];
            
            return sum >= v*size;
        }
    
        bool operator<= (double v)
        {
            double sum = 0;
            for(int i=0; i<size; ++i)
                sum = sum + this->data[i];
            
            return sum <= v*size;
        }

        bool operator== (double v)
        {
            double sum = 0;
            for(int i=0; i<size; ++i)
                sum = sum + this->data[i];

            return sum == v*size;
        } 

        bool operator!= (double v)
        {
            double sum = 0;
            for(int i=0; i<size; ++i)
                sum = sum + this->data[i];

            return sum != v*size;
        }
    
        bool operator== (const Mask<scalar> &m2)
        {
            bool all = true;
            for (int i = 0; i < size; ++i) {
                all && (this->data[i] == m2.data[i]);
            }
            return all;
        }
    
        bool operator!= (const Mask<scalar> &m2)
        {
            return !(this==m2);
        }

        Mask<scalar> operator&& (const Mask<scalar> &m2)
        {
            Mask<scalar> m3;
            for(int i=0; i<size; ++i)
                m3.data[i] = (this->data[i] && m2.data[i]);

            return m3;
        } 

        Mask<scalar> operator|| (const Mask<scalar> &m2)
        {
            Mask<scalar> m3;
            for(int i=0; i<size; ++i)
                m3.data[i] = (this->data[i] || m2.data[i]);

            return m3;
        }
    
        Mask<scalar> operator&& (bool m2)
        {
            Mask<scalar> m3;
            for(int i=0; i<size; ++i)
                m3.data[i] = (this->data[i] && m2);
            
            return m3;
        }
    
        Mask<scalar> operator|| (bool m2)
        {
            Mask<scalar> m3;
            for(int i=0; i<size; ++i)
                m3.data[i] = (this->data[i] || m2);
            
            return m3;
        }
    
        Mask<scalar> operator+ (const Mask<scalar> &m2)
        {
            Mask<scalar> m3;
            for(int i=0; i<size; ++i)
                m3.data[i] = (this->data[i] + m2.data[i]);

            return m3;
        } 

        Mask<scalar> operator- (const Mask<scalar> &m2)
        {
            Mask<scalar> m3;
            for(int i=0; i<size; ++i)
                m3.data[i] = (this->data[i] - m2.data[i]);

            return m3;
        }

        scalar operator* (const scalar &v)
        {
            typedef EnsembleTraits_m<scalar> ET; 
            scalar v2;
            for(int i=0; i<size; ++i)
                ET::coeff(v2,i) = ET::coeff(v,i)*this->data[i];

            return v2;            
        }

        bool operator[] (int i) const
        {
            return this->data[i];
        }

        bool & operator[] (int i)
        {
            return this->data[i];
        }

        Mask<scalar> operator! ()
        {
            Mask<scalar> m2;
            for(int i=0; i<size; ++i)
                m2.data[i] = !(this->data[i]);

            return m2;            
        }

        operator bool() const
        {
            return this->data[0];                            
        }   

        operator double() const
        {
            double sum = 0;
            for(int i=0; i<size; ++i)
                sum = sum + this->data[i];

            return sum/size;    
        }                               
};

template<typename scalar> std::ostream &operator<<(std::ostream &os, const Mask<scalar>& m) {
    os << "[ ";
    for(int i=0; i<m.getSize(); ++i)
        os << m[i] << " ";
    return os << "]";
}

template<typename S>  Sacado::MP::Vector<S> operator* (const Sacado::MP::Vector<S> &a1, const Mask<Sacado::MP::Vector<S>> &m)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S>> ET;
    Sacado::MP::Vector<S> mul;
    for(int i=0; i<ET::size; ++i){
        ET::coeff(mul,i) = ET::coeff(a1,i)*m[i];
    }
    return mul;
}

template<typename S>  Sacado::MP::Vector<S> operator* (const typename S::value_type &a1, const Mask<Sacado::MP::Vector<S>> &m)
{
    Sacado::MP::Vector<S> mul;
    typedef EnsembleTraits_m<Sacado::MP::Vector<S>> ET;
    for(int i=0; i<ET::size; ++i){
        ET::coeff(mul,i) = m[i]*a1;
    }
    return mul;
}

template<typename S>  Sacado::MP::Vector<S> operator* (const Mask<Sacado::MP::Vector<S>> &m, const typename S::value_type &a1)
{
    Sacado::MP::Vector<S> mul;
    typedef EnsembleTraits_m<Sacado::MP::Vector<S>> ET;
    for(int i=0; i<ET::size; ++i){
        ET::coeff(mul,i) = m[i]*a1;
    }
    return mul;
}

namespace Sacado {
    namespace MP {
        template <typename S> Vector<S> copysign(const Vector<S> &a1, const Vector<S> &a2)
        {
            typedef EnsembleTraits_m< Vector<S> > ET;
            
            Vector<S> a_out;
            
            using std::copysign;
            for(int i=0; i<ET::size; ++i){
                ET::coeff(a_out,i) = copysign(ET::coeff(a1,i),ET::coeff(a2,i));
            }
            
            return a_out;
        }
    }
}


template<typename S> Mask<Sacado::MP::Vector<S> > signbit_v(const Sacado::MP::Vector<S> &a1)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S> > ET;
    using std::signbit;
    
    Mask<Sacado::MP::Vector<S> > mask;
    for(int i=0; i<ET::size; ++i)
        mask[i] = signbit(ET::coeff(a1,i));
    return mask;
}

// Vector - vector comparisons

template<typename S> Mask<Sacado::MP::Vector<S> > operator> (const Sacado::MP::Vector<S> &a1, const Sacado::MP::Vector<S> &a2)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S> > ET;
    
    Mask<Sacado::MP::Vector<S> > mask;
    for(int i=0; i<ET::size; ++i)
        mask[i] = ET::coeff(a1,i) > ET::coeff(a2,i);
    return mask;
}

template<typename S> Mask<Sacado::MP::Vector<S> > operator>= (const Sacado::MP::Vector<S> &a1, const Sacado::MP::Vector<S> &a2)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S> > ET;
    
    Mask<Sacado::MP::Vector<S> > mask;
    for(int i=0; i<ET::size; ++i)
        mask[i] = ET::coeff(a1,i) >= ET::coeff(a2,i);
    return mask;
}

template<typename S> Mask<Sacado::MP::Vector<S> > operator< (const Sacado::MP::Vector<S> &a1, const Sacado::MP::Vector<S> &a2)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S> > ET;
    
    Mask<Sacado::MP::Vector<S> > mask;
    for(int i=0; i<ET::size; ++i)
        mask[i] = ET::coeff(a1,i) < ET::coeff(a2,i);
    return mask;
}

template<typename S> Mask<Sacado::MP::Vector<S> > operator<= (const Sacado::MP::Vector<S> &a1, const Sacado::MP::Vector<S> &a2)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S> > ET;
    
    Mask<Sacado::MP::Vector<S> > mask;
    for(int i=0; i<ET::size; ++i)
        mask[i] = ET::coeff(a1,i) <= ET::coeff(a2,i);
    return mask;
}

// Vector - scalar comparisons

template<typename S> Mask<Sacado::MP::Vector<S> > operator> (const Sacado::MP::Vector<S> &a1, const typename S::value_type &a2)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S> > ET;
    
    Mask<Sacado::MP::Vector<S> > mask;
    for(int i=0; i<ET::size; ++i)
        mask[i] = ET::coeff(a1,i) > a2;
    return mask;
}

template<typename S> Mask<Sacado::MP::Vector<S> > operator>= (const Sacado::MP::Vector<S> &a1, const typename S::value_type &a2)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S> > ET;
    
    Mask<Sacado::MP::Vector<S> > mask;
    for(int i=0; i<ET::size; ++i)
        mask[i] = ET::coeff(a1,i) >= a2;
    return mask;
}

template<typename S> Mask<Sacado::MP::Vector<S> > operator< (const Sacado::MP::Vector<S> &a1, const typename S::value_type &a2)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S> > ET;
    
    Mask<Sacado::MP::Vector<S> > mask;
    for(int i=0; i<ET::size; ++i)
        mask[i] = ET::coeff(a1,i) < a2;
    return mask;
}

template<typename S> Mask<Sacado::MP::Vector<S> > operator<= (const Sacado::MP::Vector<S> &a1, const typename S::value_type &a2)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S> > ET;
    
    Mask<Sacado::MP::Vector<S> > mask;
    for(int i=0; i<ET::size; ++i)
        mask[i] = ET::coeff(a1,i) <= a2;
    return mask;
}

// Scalar -vector comparisons

template<typename S> Mask<Sacado::MP::Vector<S> > operator> (const typename S::value_type &a2, Sacado::MP::Vector<S> &a1)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S> > ET;
    
    Mask<Sacado::MP::Vector<S> > mask;
    for(int i=0; i<ET::size; ++i)
        mask[i] = a2 > ET::coeff(a1,i);
    return mask;
}

template<typename S> Mask<Sacado::MP::Vector<S> > operator>= (const typename S::value_type &a2, const Sacado::MP::Vector<S> &a1)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S> > ET;
    
    Mask<Sacado::MP::Vector<S> > mask;
    for(int i=0; i<ET::size; ++i)
        mask[i] = a2 >= ET::coeff(a1,i);
    return mask;
}

template<typename S> Mask<Sacado::MP::Vector<S> > operator< (const typename S::value_type &a2, const Sacado::MP::Vector<S> &a1)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S> > ET;
    
    Mask<Sacado::MP::Vector<S> > mask;
    for(int i=0; i<ET::size; ++i)
        mask[i] = a2 < ET::coeff(a1,i);
    return mask;
}

template<typename S> Mask<Sacado::MP::Vector<S> > operator<= (const typename S::value_type &a2, const Sacado::MP::Vector<S> &a1)
{
    typedef EnsembleTraits_m<Sacado::MP::Vector<S> > ET;
    
    Mask<Sacado::MP::Vector<S> > mask;
    for(int i=0; i<ET::size; ++i)
        mask[i] = a2 <= ET::coeff(a1,i);
    return mask;
}

namespace MaskLogic{

class OR
{
private:
    bool value;
    
public:
    OR(){
        value = false;
    }

    OR(bool value_){
        value = value_;
    }

    template<typename T> OR(T m){
        value = (((double) m)!=0.);
    }
    
    operator bool() const
    {
        return value;
    }
};

class XOR
{
private:
    bool value;
    
public:
    XOR(){
        value = false;
    }

    XOR(bool value_){
        value = value_;
    }

    template<typename T> XOR(T m){
        value = (((double) m)==1./m.getSize());
    }
    
    operator bool() const
    {
        return value;
    }
};

class AND
{
private:
    bool value;
    
public:
    AND(){
        value = false;
    }

    AND(bool value_){
        value = value_;
    }
   
    template<typename T> AND(T m){
        value = (((double) m)==1.);
    }
    
    operator bool() const
    {
        return value;
    }
};
}
#endif // STOKHOS_MP_VECTOR_MASKTRAITS_HPP
