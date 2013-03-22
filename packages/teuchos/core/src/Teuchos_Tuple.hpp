// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_TUPLE_HPP
#define TEUCHOS_TUPLE_HPP


#include "Teuchos_ArrayView.hpp"


namespace Teuchos {


/** \brief Statically sized simple array (tuple) class.
 *
 * ToDo: Finish documentation!
 *
 * \section Teuchos_Tuple_DesignDiscussion_sec Design Discussion
 *
 * This class derives from ArrayView and therefore inherits all of the
 * features line element access and iterators along with full runtime
 * checking.  It does so at the expense of a little extra data (an extra
 * pointer and an extra integer).  However, this overhead will be small for
 * any reasonable size of N.
 *
 * \ingroup teuchos_mem_mng_grp
 */
template<typename T, int N>
class Tuple : public ArrayView<T> {
public:
  
	/** \brief Default construct raw storage.
	 */
	inline Tuple();
  
	/** \brief Copy constructor
	 */
	Tuple( const Tuple<T,N> &t );
  
	/** \brief Copy constructor
	 */
	Tuple<T,N>& operator=( const Tuple<T,N> &t );

private:

	T array_[N];

};


/** \brief Create a Tuple<T,1>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,1> tuple(const T& a);

                      
/** \brief Create a Tuple<T,2>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,2> tuple(const T& a, const T& b);


/** \brief Create a Tuple<T,3>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,3> tuple(const T& a, const T& b, const T& c);


/** \brief Create a Tuple<T,4>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,4> tuple(const T& a, const T& b, const T& c, const T& d);


/** \brief Create a Tuple<T,5>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,5> tuple(const T& a, const T& b, const T& c, const T& d, const T& e);


/** \brief Create a Tuple<T,6>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,6> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f);


/** \brief Create a Tuple<T,7>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,7> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g);


/** \brief Create a Tuple<T,8>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,8> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h);


/** \brief Create a Tuple<T,9>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,9> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i);


/** \brief Create a Tuple<T,10>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,10> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i, const T& j);


/** \brief Create a Tuple<T,11>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,11> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i, const T& j, const T& k);


/** \brief Create a Tuple<T,12>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,12> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i, const T& j, const T& k,
  const T& l);


/** \brief Create a Tuple<T,13>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,13> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i, const T& j, const T& k,
  const T& l, const T& m);


/** \brief Create a Tuple<T,14>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,14> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i, const T& j, const T& k,
  const T& l, const T& m, const T& n);


/** \brief Create a Tuple<T,15>.
 *
 * \relates Tuple
 */
template<typename T> inline
Tuple<T,15> tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i, const T& j, const T& k,
  const T& l, const T& m, const T& n, const T& o);


//
// Implementations
//


template<typename T, int N> inline
Tuple<T,N>::Tuple()
:ArrayView<T>() // To get rid of warnings!
{
  ArrayView<T>::operator=(ArrayView<T>(&array_[0],N));
}


template<typename T, int N>
Tuple<T,N>::Tuple( const Tuple<T,N> &t )
:ArrayView<T>() // To get rid of warnings!
{
  for( int i = 0; i < N; ++i )
    array_[i] = t[i];
  // Above, this loop with static N should allow the compiler to unroll this
  // entire loop!
  ArrayView<T>::operator=(ArrayView<T>(&array_[0],N));
}


template<typename T, int N>
Tuple<T,N>& Tuple<T,N>::operator=( const Tuple<T,N> &t )
{
  for( int i = 0; i < N; ++i )
    array_[i] = t[i];
  // Above, this loop with static N should allow the compiler to unroll this
  // entire loop!
  return *this;
}


} // end namespace Teuchos


//
// Nonmember function implementations
//


template<typename T> inline
Teuchos::Tuple<T,1>
Teuchos::tuple(const T& a)
{
  Tuple<T,1> rtn;
  rtn[0] = a;
  return rtn;
}


template<typename T> inline
Teuchos::Tuple<T,2>
Teuchos::tuple(const T& a, const T& b)
{
  Tuple<T,2> rtn;
  rtn[0] = a;
  rtn[1] = b;
  return rtn;
}


template<typename T> inline
Teuchos::Tuple<T,3>
Teuchos::tuple(const T& a, const T& b, const T& c)
{
  Tuple<T,3> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  return rtn;
}


template<typename T> inline
Teuchos::Tuple<T,4>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d)
{
  Tuple<T,4> rtn;;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  return rtn;
}


template<typename T> inline
Teuchos::Tuple<T,5>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e)
{
  Tuple<T,5> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  return rtn;
}


template<typename T> inline
Teuchos::Tuple<T,6>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f)
{
  Tuple<T,6> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  return rtn;
}


template<typename T> inline
Teuchos::Tuple<T,7>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g)
{
  Tuple<T,7> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  return rtn;
}


template<typename T> inline
Teuchos::Tuple<T,8>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h)
{
  Tuple<T,8> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  rtn[7] = h;
  return rtn;
}


template<typename T> inline
Teuchos::Tuple<T,9>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i)
{
  Tuple<T,9> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  rtn[7] = h;
  rtn[8] = i;
  return rtn;
}


template<typename T> inline
Teuchos::Tuple<T,10>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i, const T& j)
{
  Tuple<T,10> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  rtn[7] = h;
  rtn[8] = i;
  rtn[9] = j;
  return rtn;
}


template<typename T> inline
Teuchos::Tuple<T,11>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i, const T& j, const T& k)
{
  Tuple<T,11> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  rtn[7] = h;
  rtn[8] = i;
  rtn[9] = j;
  rtn[10] = k;
  return rtn;
}

template<typename T> inline
Teuchos::Tuple<T,12>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i, const T& j, const T& k,
  const T& l)
{
  Tuple<T,12> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  rtn[7] = h;
  rtn[8] = i;
  rtn[9] = j;
  rtn[10] = k;
  rtn[11] = l;
  return rtn;
}

template<typename T> inline
Teuchos::Tuple<T,13>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i, const T& j, const T& k,
  const T& l, const T& m)
{
  Tuple<T,13> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  rtn[7] = h;
  rtn[8] = i;
  rtn[9] = j;
  rtn[10] = k;
  rtn[11] = l;
  rtn[12] = m;
  return rtn;
}


template<typename T> inline
Teuchos::Tuple<T,14>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i, const T& j, const T& k,
  const T& l, const T& m, const T& n)
{
  Tuple<T,14> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  rtn[7] = h;
  rtn[8] = i;
  rtn[9] = j;
  rtn[10] = k;
  rtn[11] = l;
  rtn[12] = m;
  rtn[13] = n;
  return rtn;
}


template<typename T> inline
Teuchos::Tuple<T,15>
Teuchos::tuple(const T& a, const T& b, const T& c, const T& d, const T& e,
  const T& f, const T& g, const T& h, const T& i, const T& j, const T& k,
  const T& l, const T& m, const T& n, const T& o)
{
  Tuple<T,15> rtn;
  rtn[0] = a;
  rtn[1] = b;
  rtn[2] = c;
  rtn[3] = d;
  rtn[4] = e;
  rtn[5] = f;
  rtn[6] = g;
  rtn[7] = h;
  rtn[8] = i;
  rtn[9] = j;
  rtn[10] = k;
  rtn[11] = l;
  rtn[12] = m;
  rtn[13] = n;
  rtn[14] = o;
  return rtn;
}


#endif	// TEUCHOS_TUPLE_HPP
