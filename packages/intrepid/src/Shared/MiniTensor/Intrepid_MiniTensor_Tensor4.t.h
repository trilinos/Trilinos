// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
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
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER

#if !defined(Intrepid_MiniTensor_Tensor4_t_h)
#define Intrepid_MiniTensor_Tensor4_t_h

namespace Intrepid {

namespace {

template< typename T, Index N>
inline
void ones_in_ikjl(Tensor4<T, N> & A)
{

  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          if (i == k && j == l) {
            A(i,j,k,l) = 1;
          }
        }
      }
    }
  }

  return;
}

template< typename T, Index N>
inline
void ones_in_iljk(Tensor4<T, N> & A)
{

  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          if (i == l && j == k) {
            A(i,j,k,l) = 1;
          }
        }
      }
    }
  }

  return;
}

template< typename T, Index N>
inline
void ones_in_ijkl(Tensor4<T, N> & A)
{

  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          if (i == j && k == l) {
            A(i,j,k,l) = 1;
          }
        }
      }
    }
  }

  return;
}

} // anonymous namespace

//
// 4th-order identity I1
// \return \f$ \delta_{ik} \delta_{jl} \f$ such that \f$ A = I_1 A \f$
//
template<typename T, Index N>
const Tensor4<T, N>
identity_1()
{
  Tensor4<T, N> I(N, ZEROS);
  ones_in_ikjl(I);
  return I;
}

template<typename T>
const Tensor4<T, DYNAMIC>
identity_1(Index const dimension)
{
  Tensor4<T, DYNAMIC> I(dimension, ZEROS);
  ones_in_ikjl(I);
  return I;
}

template<typename T, Index N>
const Tensor4<T, N>
identity_1(Index const dimension)
{
  if (N != DYNAMIC) assert(dimension == N);

  Tensor4<T, N> I(dimension, ZEROS);
  ones_in_ikjl(I);
  return I;
}

//
// 4th-order identity I2
// \return \f$ \delta_{il} \delta_{jk} \f$ such that \f$ A^T = I_2 A \f$
//
template<typename T, Index N>
const Tensor4<T, N>
identity_2()
{
  Tensor4<T, N> I(N, ZEROS);
  ones_in_iljk(I);
  return I;
}

template<typename T>
const Tensor4<T, DYNAMIC>
identity_2(Index const dimension)
{
  Tensor4<T, DYNAMIC> I(dimension, ZEROS);
  ones_in_iljk(I);
  return I;
}

template<typename T, Index N>
const Tensor4<T, N>
identity_2(Index const dimension)
{
  if (N != DYNAMIC) assert(dimension == N);

  Tensor4<T, N> I(dimension, ZEROS);
  ones_in_iljk(I);
  return I;
}

//
// 4th-order identity I3
// \return \f$ \delta_{ij} \delta_{kl} \f$ such that \f$ I_A I = I_3 A \f$
//
template<typename T, Index N>
const Tensor4<T, N>
identity_3()
{
  Tensor4<T, N> I(N, ZEROS);
  ones_in_ijkl(I);
  return I;
}

template<typename T>
const Tensor4<T, DYNAMIC>
identity_3(Index const dimension)
{
  Tensor4<T, DYNAMIC> I(dimension, ZEROS);
  ones_in_ijkl(I);
  return I;
}

template<typename T, Index N>
const Tensor4<T, N>
identity_3(Index const dimension)
{
  if (N != DYNAMIC) assert(dimension == N);

  Tensor4<T, N> I(dimension, ZEROS);
  ones_in_ijkl(I);
  return I;
}

//
// 4th-order tensor transpose
// per Holzapfel 1.157
//
template<typename T, Index N>
Tensor4<T, N>
transpose(Tensor4<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  Tensor4<T, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          B(i, j, k, l) = A(k, l, i, j);
        }
      }
    }
  }

  return B;
}

//
// 4th-order tensor vector dot product
// \param A 4th-order tensor
// \param u vector
// \return 3rd-order tensor \f$ B = A \cdot u := B_{ijk}=A_{ijkp} u_{p} \f$
//
template<typename S, typename T, Index N>
Tensor3<typename Promote<S, T>::type, N>
dot(Tensor4<T, N> const & A, Vector<S, N> const & u)
{
  Index const
  dimension= A.get_dimension();

  assert(u.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(i,j,k,p) * u(p);
        }
        B(i,j,k) = s;
      }
    }
  }

  return B;
}

//
// vector 4th-order tensor dot product
// \param A 4th-order tensor
// \param u vector
// \return 3rd-order tensor \f$ u dot A \f$ as \f$ B_{ijk}=u_{p} A_{pijk} \f$
//
template<typename S, typename T, Index N>
Tensor3<typename Promote<S, T>::type, N>
dot(Vector<S, N> const & u, Tensor4<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += u(p) * A(p,i,j,k);
        }
        B(i,j,k) = s;
      }
    }
  }

  return B;
}

//
// 4th-order tensor vector dot2 product
// \param A 4th-order tensor
// \param u vector
// \return 3rd-order tensor \f$ B = A \cdot u := B_{ijk} = A_{ijpk} u_{p} \f$
//
template<typename S, typename T, Index N>
Tensor3<typename Promote<S, T>::type, N>
dot2(Tensor4<T, N> const & A, Vector<S, N> const & u)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += A(i,j,p,k) * u(p);
        }
        B(i,j,k) = s;
      }
    }
  }

  return B;
}

//
// vector 4th-order tensor dot2 product
// \param A 4th-order tensor
// \param u vector
// \return 3rd-order tensor \f$ u dot2 A \f$ as \f$ B_{ijk}=u_{p} A_{ipjk} \f$
//
template<typename S, typename T, Index N>
Tensor3<typename Promote<S, T>::type, N>
dot2(Vector<S, N> const & u, Tensor4<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  assert(u.get_dimension() == dimension);

  Tensor3<typename Promote<S, T>::type, N>
  B(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index p = 0; p < dimension; ++p) {
          s += u(p) * A(i,p,j,k);
        }
        B(i,j,k) = s;
      }
    }
  }

  return B;
}

//
// \return 2nd-order tensor \f$ C = A : B := C_{ij} = A_{ijpq} B_{pq} \f$
//
template<typename S, typename T, Index N>
Tensor<typename Promote<S, T>::type, N>
dotdot(Tensor4<T, N> const & A, Tensor<S, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        for (Index q = 0; q < dimension; ++q) {
          s += A(i,j,p,q) * B(p,q);
        }
      }
      C(i,j) = s;
    }
  }

  return C;
}

//
// \return 2nd-order tensor \f$ C = B : A := C_{ij} = B_{pq} A_{pqij} \f$
//
template<typename S, typename T, Index N>
Tensor<typename Promote<S, T>::type, N>
dotdot(Tensor<S, N> const & B, Tensor4<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {

      typename Promote<S, T>::type
      s = 0.0;

      for (Index p = 0; p < dimension; ++p) {
        for (Index q = 0; q < dimension; ++q) {
          s += B(p,q) * A(p,q,i,j);
        }
      }
      C(i,j) = s;
    }
  }

  return C;
}

//
// \return \f$ C = A : B := C_{ijkl} = A_{ijpq} B{pqkl} \f$
//
template<typename S, typename T, Index N>
Tensor4<typename Promote<S, T>::type, N>
dotdot(Tensor4<S, N> const & A, Tensor4<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            for (Index q = 0; q < dimension; ++q) {
              s += A(i,j,p,q) * B(p,q,k,l);
            }
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \otimes B := C_{ijkl} = A_{ij} B_{kl} \f$
//
template<typename S, typename T, Index N>
Tensor4<typename Promote<S, T>::type, N>
tensor(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          C(i,j,k,l) = A(i,j) * B(k,l);
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C_{ijkl} = A_{ik} B_{jl} \f$
//
template<typename S, typename T, Index N>
Tensor4<typename Promote<S, T>::type, N>
tensor2(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          C(i,j,k,l) = A(i,k) * B(j,l);
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C_{ijkl} = A_{il} B_{kj} \f$
//
template<typename S, typename T, Index N>
Tensor4<typename Promote<S, T>::type, N>
tensor3(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          C(i,j,k,l) = A(i,l) * B(k,j);
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B := C_{ijkl} = A_{ijkp} B_{pl} \f$
//
template<typename S, typename T, Index N>
Tensor4<typename Promote<S, T>::type, N>
dot(Tensor4<T, N> const & A, Tensor<S, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(i,j,k,p) * B(p,l);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B^T := C_{ijkl} = A_{ijkp} B_{lp} \f$
//
template<typename S, typename T, Index N>
Tensor4<typename Promote<S, T>::type, N>
dot_t(Tensor4<T, N> const & A, Tensor<S, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(i,j,k,p) * B(l,p);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B := C_{ijkl} = A_{ip} B_{pjkl} \f$
//
template<typename S, typename T, Index N>
Tensor4<typename Promote<S, T>::type, N>
dot(Tensor<S, N> const & A, Tensor4<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(i,p) * B(p,j,k,l);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A^T \cdot B := C_{ijkl} = A_{pi} B_{pjkl} \f$
//
template<typename S, typename T, Index N>
Tensor4<typename Promote<S, T>::type, N>
t_dot(Tensor<S, N> const & A, Tensor4<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(p,i) * B(p,j,k,l);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B := C_{ijkl} = A_{ijpl} B_{pk} \f$
//
template<typename S, typename T, Index N>
Tensor4<typename Promote<S, T>::type, N>
dot2(Tensor4<T, N> const & A, Tensor<S, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(i,j,p,l) * B(p,k);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B^T := C_{ijkl} = A_{ijpl} B_{kp} \f$
//
template<typename S, typename T, Index N>
Tensor4<typename Promote<S, T>::type, N>
dot2_t(Tensor4<T, N> const & A, Tensor<S, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(i,j,p,l) * B(k,p);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A \cdot B := C_{ijkl} = A_{jp} B_{ipkl} \f$
//
template<typename S, typename T, Index N>
Tensor4<typename Promote<S, T>::type, N>
dot2(Tensor<S, N> const & A, Tensor4<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(j,p) * B(i,p,k,l);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C = A^T \cdot B := C_{ijkl} = A_{pj} B_{ipkl} \f$
//
template<typename S, typename T, Index N>
Tensor4<typename Promote<S, T>::type, N>
t_dot2(Tensor<S, N> const & A, Tensor4<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index p = 0; p < dimension; ++p) {
            s += A(p,j) * B(i,p,k,l);
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// odot operator useful for \f$ \frac{\partial A^{-1}}{\partial A} \f$
// see Holzapfel eqn 6.165
// \param A 2nd-order tensor
// \param B 2nd-order tensor
// \return \f$ A \odot B \f$ which is
// \f$ C_{ijkl} = \frac{1}{2}(A_{ik} B_{jl} + A_{il} B_{jk}) \f$
//
template<typename S, typename T, Index N>
Tensor4<typename Promote<S, T>::type, N>
odot(Tensor<S, N> const & A, Tensor<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          C(i,j,k,l) = 0.5 * (A(i,k) * B(j,l) + A(i,l) * B(j,k));
        }
      }
    }
  }

  return C;
}

//
// \return \f$ C'_{i'j'k'l'} = Q_{i'i} Q_{j'j} Q_{k'k} Q_{l'l} C_{ijkl}  \f$
//
template<typename S, typename T, Index N>
Tensor4<typename Promote<S, T>::type, N>
kronecker(Tensor<S, N> const & A, Tensor4<T, N> const & B)
{
  Index const
  dimension = A.get_dimension();

  assert(B.get_dimension() == dimension);

  Tensor4<typename Promote<S, T>::type, N>
  C(dimension);

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          // we assume A is the direction cosine matrix
          // have s_{i}\be_{i} = s_{i'}\be_{i'}, want s_{i}
          // s_{i} = A_{i',i}s_{i'}
          for (Index p = 0; p < dimension; ++p) {
            for (Index q = 0; q < dimension; ++q) {
              for (Index m = 0; m < dimension; ++m) {
                for (Index n = 0; n < dimension; ++n) {
                  s += A(p,i) * A(q,j) * A(m,k) * A(n,l) * B(p,q,m,n);
                }
              }
            }
          }
          C(i,j,k,l) = s;
        }
      }
    }
  }

  return C;
}

//
// 4th-order input
// \param A 4th-order tensor
// \param is input stream
// \return is input stream
//
template<typename T, Index N>
std::istream &
operator>>(std::istream & is, Tensor4<T, N> & A)
{
  Index const
  dimension = A.get_dimension();

  for (Index i = 0; i < dimension; ++i) {
    for (Index j = 0; j < dimension; ++j) {
      for (Index k = 0; k < dimension; ++k) {
        for (Index l = 0; l < dimension; ++l) {
          is >> A(i,j,k,l);
        }
      }
    }
  }

  return is;
}

//
// 4th-order output
// \param A 4th-order tensor
// \param os output stream
// \return os output stream
//
template<typename T, Index N>
std::ostream &
operator<<(std::ostream & os, Tensor4<T, N> const & A)
{
  Index const
  dimension = A.get_dimension();

  if (dimension == 0) {
    return os;
  }

  for (Index i = 0; i < dimension; ++i) {

    for (Index j = 0; j < dimension; ++j) {

      for (Index k = 0; k < dimension; ++k) {

        os << std::scientific << std::setprecision(16) << A(i,j,k,0);

        for (Index l = 1; l < dimension; ++l) {

          os << "," << std::scientific  << std::setprecision(16) << A(i,j,k,l);
        }

        os << std::endl;

      }

      os << std::endl;
      os << std::endl;

    }

    os << std::endl;

  }

  return os;
}

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Tensor4_t_h
