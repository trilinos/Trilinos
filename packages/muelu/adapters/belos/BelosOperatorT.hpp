// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef BELOS_OPERATORT_HPP
#define BELOS_OPERATORT_HPP

#include "BelosOperatorTraits.hpp"

namespace Belos {

  // TODO: this file should certainly be moved to Belos.

  // Base class for the Belos OP template type. Similar to Belos::Operator<> but this one deals with any kind of vector (not only Belos::MultiVec as the Belos::Operator interface)
  template <class MV>
  class OperatorT {

  public:

    //! @name Constructor/Destructor
    //@{

    //! Default constructor
    OperatorT() {};

    //! Destructor.
    virtual ~OperatorT() {};
    //@}

    //! @name Operator application method
    //@{

    /*! \brief This routine takes the Belos::MultiVec \c x and applies the operator
      to it resulting in the Belos::MultiVec \c y, which is returned.
      \note It is expected that any problem with applying this operator to \c x will be
      indicated by an std::exception being thrown.
    */
    virtual void Apply ( const MV & x, MV & y, ETrans trans=NOTRANS ) const = 0;
  };

  /// \brief Specialization of OperatorTraits for OperatorT.
  ///
  /// This is a partial template specialization of
  /// Belos::OperatorTraits class using the Belos::OperatorT
  /// abstract interface. Any class that inherits
  /// from Belos::OperatorT will be accepted by the Belos templated
  /// solvers, due to this specialization of Belos::OperatorTraits.
  template <class ScalarType, class MV>
  class OperatorTraits<ScalarType, MV, OperatorT<MV> >
  {

  public:
    //! Specialization of Apply() for OperatorT.
    static void Apply (const OperatorT<MV>& Op,
                       const MV& x,
                       MV& y, ETrans trans=NOTRANS) {
      Op.Apply (x, y, trans);
    }
  };

} // namespace Belos

#endif // BELOS_OPERATORT_HPP
