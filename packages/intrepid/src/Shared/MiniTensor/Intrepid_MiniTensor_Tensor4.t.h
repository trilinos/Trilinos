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

  //
  // set dimension
  //
  //
  template<typename T>
  void
  Tensor4<T>::set_dimension(Index const N)
  {
    if (N == dimension) return;

    if (e != NULL) {
      delete [] e;
    }

    Index const
    number_components = N * N * N * N;

    e = new T[number_components];

    dimension = N;

    return;
  }

  //
  // R^N 4th-order tensor default constructor
  //
  template<typename T>
  Tensor4<T>::Tensor4() :
    dimension(0),
    e(NULL)
  {
    return;
  }

  //
  // R^N 4th-order tensor constructor with NaNs
  //
  template<typename T>
  Tensor4<T>::Tensor4(Index const N) :
    dimension(0),
    e(NULL)
  {
    set_dimension(N);

    Index const
    number_components = N * N * N * N;

    for (Index i = 0; i < number_components; ++i) {
      e[i] = not_a_number<T>();
    }

    return;
  }

  //
  // R^N 4th-order tensor constructor with a scalar
  // \param s all components set to this scalar
  //
  template<typename T>
  Tensor4<T>::Tensor4(Index const N, T const & s) :
    dimension(0),
    e(NULL)
  {
    set_dimension(N);

    Index const
    number_components = N * N * N * N;

    for (Index i = 0; i < number_components; ++i) {
      e[i] = s;
    }

    return;
  }

  //
  // R^N copy constructor
  // 4th-order tensor constructor with 4th-order tensor
  // \param A from which components are copied
  //
  template<typename T>
  Tensor4<T>::Tensor4(Tensor4<T> const & A) :
    dimension(0),
    e(NULL)
  {
    Index const
    N = A.get_dimension();

    set_dimension(N);

    Index const
    number_components = N * N * N * N;

    for (Index i = 0; i < number_components; ++i) {
      e[i] = A.e[i];
    }

    return;
  }

  //
  // R^N 4th-order tensor simple destructor
  //
  template<typename T>
  Tensor4<T>::~Tensor4()
  {
    if (e != NULL) {
      delete [] e;
    }
    return;
  }

  //
  // R^N 4th-order tensor copy assignment
  //
  template<typename T>
  Tensor4<T> &
  Tensor4<T>::operator=(Tensor4<T> const & A)
  {
    if (this != &A) {
      Index const
      N = A.get_dimension();

      set_dimension(N);

      Index const
      number_components = N * N * N * N;

      for (Index i = 0; i < number_components; ++i) {
        e[i] = A.e[i];
      }

    }

    return *this;
  }

  //
  // 4th-order tensor increment
  // \param A added to this tensor
  //
  template<typename T>
  Tensor4<T> &
  Tensor4<T>::operator+=(Tensor4<T> const & A)
  {
    Index const
    N = get_dimension();

    assert(A.get_dimension() == N);

    Index const
    number_components = N * N * N * N;

    for (Index i = 0; i < number_components; ++i) {
      e[i] += A.e[i];
    }

    return *this;
  }

  //
  // 4th-order tensor decrement
  // \param A substracted from this tensor
  //
  template<typename T>
  Tensor4<T> &
  Tensor4<T>::operator-=(Tensor4<T> const & A)
  {
    Index const
    N = get_dimension();

    assert(A.get_dimension() == N);

    Index const
    number_components = N * N * N * N;

    for (Index i = 0; i < number_components; ++i) {
      e[i] -= A.e[i];
    }

    return *this;
  }

  //
  // R^N fill 4th-order tensor with zeros
  //
  template<typename T>
  void
  Tensor4<T>::clear()
  {
    Index const
    N = get_dimension();

    Index const
    number_components = N * N * N * N;

    for (Index i = 0; i < number_components; ++i) {
      e[i] = 0.0;;
    }

    return;
  }

  //
  // 4th-order tensor addition
  // \param A 4th-order tensor
  // \param B 4th-order tensor
  // \return \f$ A + B \f$
  //
  template<typename S, typename T>
  Tensor4<typename Promote<S, T>::type>
  operator+(Tensor4<S> const & A, Tensor4<T> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor4<typename Promote<S, T>::type>
    C(N);


    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {
            C(i,j,k,l) = A(i,j,k,l) + B(i,j,k,l);
          }
        }
      }
    }

    return C;
  }

  //
  // 4th-order tensor substraction
  // \param A 4th-order tensor
  // \param B 4th-order tensor
  // \return \f$ A - B \f$
  //
  template<typename S, typename T>
  Tensor4<typename Promote<S, T>::type>
  operator-(Tensor4<S> const & A, Tensor4<T> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor4<typename Promote<S, T>::type>
    C(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {
            C(i,j,k,l) = A(i,j,k,l) - B(i,j,k,l);
          }
        }
      }
    }

    return C;
  }

  //
  // 4th-order tensor minus
  // \return \f$ -A \f$
  //
  template<typename T>
  Tensor4<T>
  operator-(Tensor4<T> const & A)
  {
    Index const
    N = A.get_dimension();

    Tensor4<T> S(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {
            S(i,j,k,l) = - A(i,j,k,l);
          }
        }
      }
    }

    return S;
  }

  //
  // 4th-order equality
  // Tested by components
  //
  template<typename T>
  inline bool
  operator==(Tensor4<T> const & A, Tensor4<T> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {
            if (A(i,j,k,l) != B(i,j,k,l)) {
              return false;
            }
          }
        }
      }
    }

    return true;
  }

  //
  // 4th-order inequality
  // Tested by components
  //
  template<typename T>
  inline bool
  operator!=(Tensor4<T> const & A, Tensor4<T> const & B)
  {
    return !(A==B);
  }

  //
  // Scalar 4th-order tensor product
  // \param s scalar
  // \param A 4th-order tensor
  // \return \f$ s A \f$
  //
  template<typename S, typename T>
  typename lazy_disable_if< order_1234<S>, apply_tensor4< Promote<S,T> > >::type
  operator*(S const & s, Tensor4<T> const & A)
  {
    Index const
    N = A.get_dimension();

    Tensor4<typename Promote<S, T>::type>
    B(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {
            B(i,j,k,l) = s * A(i,j,k,l);
          }
        }
      }
    }

    return B;
  }

  //
  // 4th-order tensor scalar product
  // \param A 4th-order tensor
  // \param s scalar
  // \return \f$ s A \f$
  //
  template<typename S, typename T>
  typename lazy_disable_if< order_1234<S>, apply_tensor4< Promote<S,T> > >::type
  operator*(Tensor4<T> const & A, S const & s)
  {
    return s * A;
  }

  //
  // 4th-order tensor scalar division
  // \param A 4th-order tensor
  // \param s scalar
  // \return \f$ s A \f$
  //
  template<typename S, typename T>
  Tensor4<typename Promote<S, T>::type>
  operator/(Tensor4<T> const & A, S const & s)
  {
    Index const
    N = A.get_dimension();

    Tensor4<typename Promote<S, T>::type>
    B(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {
            B(i,j,k,l) = A(i,j,k,l) / s;
          }
        }
      }
    }

    return B;
  }

  //
  // 4th-order identity I1
  // \return \f$ \delta_{ik} \delta_{jl} \f$ such that \f$ A = I_1 A \f$
  //
  template<typename T>
  const Tensor4<T>
  identity_1(Index const N)
  {
    Tensor4<T> I(N, T(0.0));

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {
            if (i == k && j == l) {
              I(i,j,k,l) = 1.0;
            }
          }
        }
      }
    }

    return I;
  }

  //
  // 4th-order identity I2
  // \return \f$ \delta_{il} \delta_{jk} \f$ such that \f$ A^T = I_2 A \f$
  //
  template<typename T>
  const Tensor4<T>
  identity_2(Index const N)
  {
    Tensor4<T> I(N, T(0.0));

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {
            if (i == l && j == k) {
              I(i,j,k,l) = 1.0;
            }
          }
        }
      }
    }

    return I;
  }

  //
  // 4th-order identity I3
  // \return \f$ \delta_{ij} \delta_{kl} \f$ such that \f$ I_A I = I_3 A \f$
  //
  template<typename T>
  const Tensor4<T>
  identity_3(Index const N)
  {
    Tensor4<T> I(N, T(0.0));

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {
            if (i == j && k == l) {
              I(i,j,k,l) = 1.0;
            }
          }
        }
      }
    }

    return I;
  }

  //
  // 4th-order tensor vector dot product
  // \param A 4th-order tensor
  // \param u vector
  // \return 3rd-order tensor \f$ A dot u \f$ as \f$ B_{ijk}=A_{ijkl}u_{l} \f$
  //
  template<typename S, typename T>
  Tensor3<typename Promote<S, T>::type>
  dot(Tensor4<T> const & A, Vector<S> const & u)
  {
    Index const
    N = A.get_dimension();

    assert(u.get_dimension() == N);

    Tensor3<typename Promote<S, T>::type>
    B(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index l = 0; l < N; ++l) {
            s += A(i,j,k,l) * u(l);
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
  // \return 3rd-order tensor \f$ u dot A \f$ as \f$ B_{jkl}=u_{i}A_{ijkl} \f$
  //
  template<typename S, typename T>
  Tensor3<typename Promote<S, T>::type>
  dot(Vector<S> const & u, Tensor4<T> const & A)
  {
    Index const
    N = A.get_dimension();

    assert(u.get_dimension() == N);

    Tensor3<typename Promote<S, T>::type>
    B(N);

    for (Index j = 0; j < N; ++j) {
      for (Index k = 0; k < N; ++k) {
        for (Index l = 0; l < N; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index i = 0; i < N; ++i) {
            s += u(i) * A(i,j,k,l);
          }
          B(j,k,l) = s;
        }
      }
    }

    return B;
  }

  //
  // 4th-order tensor vector dot2 product
  // \param A 4th-order tensor
  // \param u vector
  // \return 3rd-order tensor \f$ A dot2 u \f$ as \f$ B_{ijl}=A_{ijkl}u_{k} \f$
  //
  template<typename S, typename T>
  Tensor3<typename Promote<S, T>::type>
  dot2(Tensor4<T> const & A, Vector<S> const & u)
  {
    Index const
    N = A.get_dimension();

    assert(u.get_dimension() == N);

    Tensor3<typename Promote<S, T>::type>
    B(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index l = 0; l < N; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index k = 0; k < N; ++k) {
            s += A(i,j,k,l) * u(k);
          }
          B(i,j,l) = s;
        }
      }
    }

    return B;
  }

  //
  // vector 4th-order tensor dot2 product
  // \param A 4th-order tensor
  // \param u vector
  // \return 3rd-order tensor \f$ u dot2 A \f$ as \f$ B_{ikl}=u_{j}A_{ijkl} \f$
  //
  template<typename S, typename T>
  Tensor3<typename Promote<S, T>::type>
  dot2(Vector<S> const & u, Tensor4<T> const & A)
  {
    Index const
    N = A.get_dimension();

    assert(u.get_dimension() == N);

    Tensor3<typename Promote<S, T>::type>
    B(N);

    for (Index i = 0; i < N; ++i) {
      for (Index k = 0; k < N; ++k) {
        for (Index l = 0; l < N; ++l) {

          typename Promote<S, T>::type
          s = 0.0;

          for (Index j = 0; j < N; ++j) {
            s += u(j) * A(i,j,k,l);
          }
          B(i,k,l) = s;
        }
      }
    }

    return B;
  }

  //
  // \return 2nd-order tensor \f$ C = A : B := C_{ij}=A_{ijkl}B_{kl} \f$
  //
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  dotdot(Tensor4<T> const & A, Tensor<S> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor<typename Promote<S, T>::type>
    C(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {
            s += A(i,j,k,l) * B(k,l);
          }
        }
        C(i,j) = s;
      }
    }

    return C;
  }

  //
  // \return 2nd-order tensor \f$ C = B : A := C_{kl} = A_{ijkl} B_{ij} \f$
  //
  template<typename S, typename T>
  Tensor<typename Promote<S, T>::type>
  dotdot(Tensor<S> const & B, Tensor4<T> const & A)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor<typename Promote<S, T>::type>
    C(N);

    for (Index k = 0; k < N; ++k) {
      for (Index l = 0; l < N; ++l) {

        typename Promote<S, T>::type
        s = 0.0;

        for (Index i = 0; i < N; ++i) {
          for (Index j = 0; j < N; ++j) {
            s += A(i,j,k,l) * B(i,j);
          }
        }
        C(k,l) = s;
      }
    }

    return C;
  }

  //
  // \return \f$ C = A : B := C_{ijkl} = A_{ijmn} B{mnkl} \f$
  //
  template<typename S, typename T>
  Tensor4<typename Promote<S, T>::type>
  dotdot(Tensor4<S> const & A, Tensor4<T> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor4<typename Promote<S, T>::type>
    C(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {

            typename Promote<S, T>::type
            s = 0.0;

            for (Index m = 0; m < N; ++m) {
              for (Index n = 0; n < N; ++n) {
                s += A(i,j,m,n) * B(m,n,k,l);
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
  template<typename S, typename T>
  Tensor4<typename Promote<S, T>::type>
  tensor(Tensor<S> const & A, Tensor<T> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor4<typename Promote<S, T>::type>
    C(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {
            C(i,j,k,l) = A(i,j) * B(k,l);
          }
        }
      }
    }

    return C;
  }

  //
  // \return \f$ C = A \cdot B := C_{ijkl} = A_{ijkp} B_{pl} \f$
  //
  template<typename S, typename T>
  Tensor4<typename Promote<S, T>::type>
  dot(Tensor4<T> const & A, Tensor<S> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor4<typename Promote<S, T>::type>
    C(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {

            typename Promote<S, T>::type
            s = 0.0;

            for (Index p = 0; p < N; ++p) {
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
  template<typename S, typename T>
  Tensor4<typename Promote<S, T>::type>
  dot_t(Tensor4<T> const & A, Tensor<S> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor4<typename Promote<S, T>::type>
    C(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {

            typename Promote<S, T>::type
            s = 0.0;

            for (Index p = 0; p < N; ++p) {
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
  template<typename S, typename T>
  Tensor4<typename Promote<S, T>::type>
  dot(Tensor<S> const & A, Tensor4<T> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor4<typename Promote<S, T>::type>
    C(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {

            typename Promote<S, T>::type
            s = 0.0;

            for (Index p = 0; p < N; ++p) {
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
  template<typename S, typename T>
  Tensor4<typename Promote<S, T>::type>
  t_dot(Tensor<S> const & A, Tensor4<T> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor4<typename Promote<S, T>::type>
    C(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {

            typename Promote<S, T>::type
            s = 0.0;

            for (Index p = 0; p < N; ++p) {
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
  // \return \f$ C = A \cdot B := C_{ijkl} = A_{ijpk} B_{pl} \f$
  //
  template<typename S, typename T>
  Tensor4<typename Promote<S, T>::type>
  dot2(Tensor4<T> const & A, Tensor<S> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor4<typename Promote<S, T>::type>
    C(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {

            typename Promote<S, T>::type
            s = 0.0;

            for (Index p = 0; p < N; ++p) {
              s += A(i,j,p,k) * B(p,l);
            }
            C(i,j,k,l) = s;
          }
        }
      }
    }

    return C;
  }

  //
  // \return \f$ C = A \cdot B^T := C_{ijkl} = A_{ijpk} B_{lp} \f$
  //
  template<typename S, typename T>
  Tensor4<typename Promote<S, T>::type>
  dot2_t(Tensor4<T> const & A, Tensor<S> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor4<typename Promote<S, T>::type>
    C(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {

            typename Promote<S, T>::type
            s = 0.0;

            for (Index p = 0; p < N; ++p) {
              s += A(i,j,p,k) * B(l,p);
            }
            C(i,j,k,l) = s;
          }
        }
      }
    }

    return C;
  }

  //
  // \return \f$ C = A \cdot B := C_{ijkl} = A_{ip} B_{jpkl} \f$
  //
  template<typename S, typename T>
  Tensor4<typename Promote<S, T>::type>
  dot2(Tensor<S> const & A, Tensor4<T> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor4<typename Promote<S, T>::type>
    C(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {

            typename Promote<S, T>::type
            s = 0.0;

            for (Index p = 0; p < N; ++p) {
              s += A(i,p) * B(j,p,k,l);
            }
            C(i,j,k,l) = s;
          }
        }
      }
    }

    return C;
  }

  //
  // \return \f$ C = A^T \cdot B := C_{ijkl} = A_{pi} B_{jpkl} \f$
  //
  template<typename S, typename T>
  Tensor4<typename Promote<S, T>::type>
  t_dot2(Tensor<S> const & A, Tensor4<T> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor4<typename Promote<S, T>::type>
    C(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {

            typename Promote<S, T>::type
            s = 0.0;

            for (Index p = 0; p < N; ++p) {
              s += A(p,i) * B(j,p,k,l);
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
  template<typename S, typename T>
  Tensor4<typename Promote<S, T>::type>
  odot(Tensor<S> const & A, Tensor<T> const & B)
  {
    Index const
    N = A.get_dimension();

    assert(B.get_dimension() == N);

    Tensor4<typename Promote<S, T>::type>
    C(N);

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {
            C(i,j,k,l) = 0.5 * (A(i,k) * B(j,l) + A(i,l) * B(j,k));
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
  template<typename T>
  std::istream &
  operator>>(std::istream & is, Tensor4<T> & A)
  {
    Index const
    N = A.get_dimension();

    for (Index i = 0; i < N; ++i) {
      for (Index j = 0; j < N; ++j) {
        for (Index k = 0; k < N; ++k) {
          for (Index l = 0; l < N; ++l) {
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
  template<typename T>
  std::ostream &
  operator<<(std::ostream & os, Tensor4<T> const & A)
  {
    Index const
    N = A.get_dimension();

    if (N == 0) {
      return os;
    }

    for (Index i = 0; i < N; ++i) {

      for (Index j = 0; j < N; ++j) {

        for (Index k = 0; k < N; ++k) {

          os << std::scientific << A(i,j,k,0);

          for (Index l = 1; l < N; ++l) {

            os << std::scientific << "," << A(i,j,k,l);
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
