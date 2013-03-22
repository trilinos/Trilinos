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

#if !defined(Intrepid_MiniTensor_Mechanics_t_h)
#define Intrepid_MiniTensor_Mechanics_t_h

namespace Intrepid {

  //
  // Push forward covariant vector
  // \param \f$ F, u \f$
  // \return \f$ F^{-T} u \f$
  //
  template<typename T>
  Vector<T>
  push_forward_covariant(Tensor<T> const & F, Vector<T> const & u)
  {
    Index const
    N = F.get_dimension();

    Vector<T>
    v(N);

    T const
    J = det(F);

    assert(J > 0.0);

    switch (N) {

      default:
        std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
        std::cerr << std::endl;
        std::cerr << "Supports only 2D and 3D. Found dimension: " << N;
        std::cerr << std::endl;
        exit(1);
        break;

      case 3:
        v(0) = ((-F(1,2)*F(2,1) + F(1,1)*F(2,2)) * u(0) +
                ( F(1,2)*F(2,0) - F(1,0)*F(2,2)) * u(1) +
                (-F(1,1)*F(2,0) + F(1,0)*F(2,1)) * u(2)) / J;

        v(1) = (( F(0,2)*F(2,1) - F(0,1)*F(2,2)) * u(0) +
                (-F(0,2)*F(2,0) + F(0,0)*F(2,2)) * u(1) +
                ( F(0,1)*F(2,0) - F(0,0)*F(2,1)) * u(2)) / J;

        v(2) = ((-F(0,2)*F(1,1) + F(0,1)*F(1,2)) * u(0) +
                ( F(0,2)*F(1,0) - F(0,0)*F(1,2)) * u(1) +
                (-F(0,1)*F(1,0) + F(0,0)*F(1,1)) * u(2)) / J;

        break;

      case 2:
        v(0) = ( F(1,1) * u(0) - F(1,0) * u(1)) / J;
        v(1) = (-F(0,1) * u(0) + F(0,0) * u(1)) / J;
        break;

    }

    return v;
  }

  //
  // Pull back covariant vector
  // \param \f$ F, v \f$
  // \return \f$ F^T v \f$
  //
  template<typename T>
  Vector<T>
  pull_back_covariant(Tensor<T> const & F, Vector<T> const & u)
  {
    Index const
    N = F.get_dimension();

    Vector<T>
    v(N);

    switch (N) {

      default:
        std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
        std::cerr << std::endl;
        std::cerr << "Supports only 2D and 3D. Found dimension: " << N;
        std::cerr << std::endl;
        exit(1);
        break;

      case 3:
        v(0) = F(0,0) * u(0) + F(1,0) * u(1) + F(2,0) * u(2);
        v(1) = F(0,1) * u(0) + F(1,1) * u(1) + F(2,1) * u(2);
        v(2) = F(0,2) * u(0) + F(1,2) * u(1) + F(2,2) * u(2);

        break;

      case 2:
        v(0) = F(0,0) * u(0) + F(1,0) * u(1);
        v(1) = F(0,1) * u(0) + F(1,1) * u(1);

        break;

    }

    return v;
  }

  //
  // Push forward contravariant vector
  // \param \f$ F, u \f$
  // \return \f$ F u \f$
  //
  template<typename T>
  Vector<T>
  push_forward_contravariant(Tensor<T> const & F, Vector<T> const & u)
  {
    Index const
    N = F.get_dimension();

    Vector<T>
    v(N);

    switch (N) {

      default:
        std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
        std::cerr << std::endl;
        std::cerr << "Supports only 2D and 3D. Found dimension: " << N;
        std::cerr << std::endl;
        exit(1);
        break;

      case 3:
        v(0) = F(0,0) * u(0) + F(0,1) * u(1) + F(0,2) * u(2);
        v(1) = F(1,0) * u(0) + F(1,1) * u(1) + F(1,2) * u(2);
        v(2) = F(2,0) * u(0) + F(2,1) * u(1) + F(2,2) * u(2);

        break;

      case 2:
        v(0) = F(0,0) * u(0) + F(0,1) * u(1);
        v(1) = F(1,0) * u(0) + F(1,1) * u(1);

        break;

    }

    return v;
  }

  //
  // Pull back contravariant vector
  // \param \f$ F, u \f$
  // \return \f$ F^{-1} u \f$
  //
  template<typename T>
  Vector<T>
  pull_back_contravariant(Tensor<T> const & F, Vector<T> const & u)
  {
    Index const
    N = F.get_dimension();

    Vector<T>
    v(N);

    T const
    J = det(F);

    assert(J > 0.0);

    switch (N) {

      default:
        std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
        std::cerr << std::endl;
        std::cerr << "Supports only 2D and 3D. Found dimension: " << N;
        std::cerr << std::endl;
        exit(1);
        break;

      case 3:
        v(0) = ((-F(1,2)*F(2,1) + F(1,1)*F(2,2)) * u(0) +
                ( F(0,2)*F(2,1) - F(0,1)*F(2,2)) * u(1) +
                (-F(0,2)*F(1,1) + F(0,1)*F(1,2)) * u(2)) / J;

        v(1) = (( F(1,2)*F(2,0) - F(1,0)*F(2,2)) * u(0) +
                (-F(0,2)*F(2,0) + F(0,0)*F(2,2)) * u(1) +
                ( F(0,2)*F(1,0) - F(0,0)*F(1,2)) * u(2)) / J;

        v(2) = ((-F(1,1)*F(2,0) + F(1,0)*F(2,1)) * u(0) +
                ( F(0,1)*F(2,0) - F(0,0)*F(2,1)) * u(1) +
                (-F(0,1)*F(1,0) + F(0,0)*F(1,1)) * u(2)) / J;

        break;

      case 2:
        v(0) = ( F(1,1) * u(0) - F(0,1) * u(1)) / J;
        v(1) = (-F(1,0) * u(0) + F(0,0) * u(1)) / J;
        break;

    }

    return v;
  }

  //
  // Push forward covariant tensor
  // \param \f$ F, A \f$
  // \return \f$ F^{-T} A F^{-1} \f$
  //
  template<typename T>
  Tensor<T>
  push_forward_covariant(Tensor<T> const & F, Tensor<T> const & A)
  {
    Index const
    N = F.get_dimension();

    Tensor<T>
    G(N);

    T const
    J = det(F);

    assert(J > 0.0);

    switch (N) {

      default:
        std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
        std::cerr << std::endl;
        std::cerr << "Supports only 2D and 3D. Found dimension: " << N;
        std::cerr << std::endl;
        exit(1);
        break;

      case 3:
        G(0,0) = (-F(1,2)*F(2,1) + F(1,1)*F(2,2)) / J;
        G(0,1) = ( F(0,2)*F(2,1) - F(0,1)*F(2,2)) / J;
        G(0,2) = (-F(0,2)*F(1,1) + F(0,1)*F(1,2)) / J;

        G(1,0) = ( F(1,2)*F(2,0) - F(1,0)*F(2,2)) / J;
        G(1,1) = (-F(0,2)*F(2,0) + F(0,0)*F(2,2)) / J;
        G(1,2) = ( F(0,2)*F(1,0) - F(0,0)*F(1,2)) / J;

        G(2,0) = (-F(1,1)*F(2,0) + F(1,0)*F(2,1)) / J;
        G(2,1) = ( F(0,1)*F(2,0) - F(0,0)*F(2,1)) / J;
        G(2,2) = (-F(0,1)*F(1,0) + F(0,0)*F(1,1)) / J;
        break;

      case 2:
        G(0,0) =  F(1,1) / J;
        G(0,1) = -F(0,1) / J;

        G(1,0) = -F(1,0) / J;
        G(1,1) =  F(0,0) / J;
        break;

    }

    return t_dot(G, dot(A, G));
  }

  //
  // Pull back covariant tensor
  // \param \f$ F, A \f$
  // \return \f$ F^T A F\f$
  //
  template<typename T>
  Tensor<T>
  pull_back_covariant(Tensor<T> const & F, Tensor<T> const & A)
  {
    return t_dot(F, dot(A, F));
  }

  //
  // Push forward contravariant tensor
  // \param \f$ F, A \f$
  // \return \f$ F A F^T \f$
  //
  template<typename T>
  Tensor<T>
  push_forward_contravariant(Tensor<T> const & F, Tensor<T> const & A)
  {
    return dot_t(dot(F, A), F);
  }

  //
  // Pull back contravariant tensor
  // \param \f$ F, A \f$
  // \return \f$ F^{-1} A F^{-T} \f$
  //
  template<typename T>
  Tensor<T>
  pull_back_contravariant(Tensor<T> const & F, Tensor<T> const & A)
  {
    Index const
    N = F.get_dimension();

    Tensor<T>
    G(N);

    T const
    J = det(F);

    assert(J > 0.0);

    switch (N) {

      default:
        std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
        std::cerr << std::endl;
        std::cerr << "Supports only 2D and 3D. Found dimension: " << N;
        std::cerr << std::endl;
        exit(1);
        break;

      case 3:
        G(0,0) = (-F(1,2)*F(2,1) + F(1,1)*F(2,2)) / J;
        G(0,1) = ( F(0,2)*F(2,1) - F(0,1)*F(2,2)) / J;
        G(0,2) = (-F(0,2)*F(1,1) + F(0,1)*F(1,2)) / J;

        G(1,0) = ( F(1,2)*F(2,0) - F(1,0)*F(2,2)) / J;
        G(1,1) = (-F(0,2)*F(2,0) + F(0,0)*F(2,2)) / J;
        G(1,2) = ( F(0,2)*F(1,0) - F(0,0)*F(1,2)) / J;

        G(2,0) = (-F(1,1)*F(2,0) + F(1,0)*F(2,1)) / J;
        G(2,1) = ( F(0,1)*F(2,0) - F(0,0)*F(2,1)) / J;
        G(2,2) = (-F(0,1)*F(1,0) + F(0,0)*F(1,1)) / J;
        break;

      case 2:
        G(0,0) =  F(1,1) / J;
        G(0,1) = -F(0,1) / J;

        G(1,0) = -F(1,0) / J;
        G(1,1) =  F(0,0) / J;
        break;

    }

    return dot_t(dot(G, A), G);
  }

  //
  // Piola transformation for vector
  // \param \f$ F, u \f$
  // \return \f$ \det F F^{-1} u \f$
  //
  template<typename T>
  Vector<T>
  piola(Tensor<T> const & F, Vector<T> const & u)
  {
    Index const
    N = F.get_dimension();

    Vector<T>
    v(N);

    switch (N) {

      default:
        std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
        std::cerr << std::endl;
        std::cerr << "Supports only 2D and 3D. Found dimension: " << N;
        std::cerr << std::endl;
        exit(1);
        break;

      case 3:
        v(0) = ((-F(1,2)*F(2,1) + F(1,1)*F(2,2)) * u(0) +
                ( F(0,2)*F(2,1) - F(0,1)*F(2,2)) * u(1) +
                (-F(0,2)*F(1,1) + F(0,1)*F(1,2)) * u(2));

        v(1) = (( F(1,2)*F(2,0) - F(1,0)*F(2,2)) * u(0) +
                (-F(0,2)*F(2,0) + F(0,0)*F(2,2)) * u(1) +
                ( F(0,2)*F(1,0) - F(0,0)*F(1,2)) * u(2));

        v(2) = ((-F(1,1)*F(2,0) + F(1,0)*F(2,1)) * u(0) +
                ( F(0,1)*F(2,0) - F(0,0)*F(2,1)) * u(1) +
                (-F(0,1)*F(1,0) + F(0,0)*F(1,1)) * u(2));

        break;

      case 2:
        v(0) = ( F(1,1) * u(0) - F(0,1) * u(1));
        v(1) = (-F(1,0) * u(0) + F(0,0) * u(1));
        break;

    }

    return v;
  }

  //
  // Inverse Piola transformation for vector
  // \param \f$ F, u \f$
  // \return \f$ (\det F)^{-1} F u \f$
  //
  template<typename T>
  Vector<T>
  piola_inverse(Tensor<T> const & F, Vector<T> const & u)
  {
    Index const
    N = F.get_dimension();

    Vector<T>
    v(N);

    T const
    J = det(F);

    assert(J > 0.0);

    switch (N) {

      default:
        std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
        std::cerr << std::endl;
        std::cerr << "Supports only 2D and 3D. Found dimension: " << N;
        std::cerr << std::endl;
        exit(1);
        break;

      case 3:
        v(0) = (F(0,0) * u(0) + F(0,1) * u(1) + F(0,2) * u(2)) / J;
        v(1) = (F(1,0) * u(0) + F(1,1) * u(1) + F(1,2) * u(2)) / J;
        v(2) = (F(2,0) * u(0) + F(2,1) * u(1) + F(2,2) * u(2)) / J;

        break;

      case 2:
        v(0) = (F(0,0) * u(0) + F(0,1) * u(1)) / J;
        v(1) = (F(1,0) * u(0) + F(1,1) * u(1)) / J;

        break;

    }

    return v;
  }

  //
  // Piola transformation for tensor, applied on second
  // index. Useful for transforming Cauchy stress to 1PK stress.
  // \param \f$ F, \sigma \f$
  // \return \f$ \det F \sigma F^{-T} \f$
  //
  template<typename T>
  Tensor<T>
  piola(Tensor<T> const & F, Tensor<T> const & sigma)
  {
    Index const
    N = F.get_dimension();

    Tensor<T>
    G(N);

    switch (N) {

      default:
        std::cerr << "ERROR: " << __PRETTY_FUNCTION__;
        std::cerr << std::endl;
        std::cerr << "Supports only 2D and 3D. Found dimension: " << N;
        std::cerr << std::endl;
        exit(1);
        break;

      case 3:
        G(0,0) = (-F(1,2)*F(2,1) + F(1,1)*F(2,2));
        G(0,1) = ( F(0,2)*F(2,1) - F(0,1)*F(2,2));
        G(0,2) = (-F(0,2)*F(1,1) + F(0,1)*F(1,2));

        G(1,0) = ( F(1,2)*F(2,0) - F(1,0)*F(2,2));
        G(1,1) = (-F(0,2)*F(2,0) + F(0,0)*F(2,2));
        G(1,2) = ( F(0,2)*F(1,0) - F(0,0)*F(1,2));

        G(2,0) = (-F(1,1)*F(2,0) + F(1,0)*F(2,1));
        G(2,1) = ( F(0,1)*F(2,0) - F(0,0)*F(2,1));
        G(2,2) = (-F(0,1)*F(1,0) + F(0,0)*F(1,1));
        break;

      case 2:
        G(0,0) =  F(1,1);
        G(0,1) = -F(0,1);

        G(1,0) = -F(1,0);
        G(1,1) =  F(0,0);
        break;

    }

    return dot_t(sigma, G);
  }

  //
  // Inverse Piola transformation for tensor, applied on second
  // index. Useful for transforming 1PK stress to Cauchy stress.
  // \param \f$ F, P \f$
  // \return \f$ (\det F)^{-1} P F^T \f$
  //
  template<typename T>
  Tensor<T>
  piola_inverse(Tensor<T> const & F, Tensor<T> const & P)
  {
    T const
    J = det(F);

    assert(J > 0.0);

    return dot_t(P, F) / J;
  }

} // namespace Intrepid

#endif // Intrepid_MiniTensor_Mechanics_t_h




