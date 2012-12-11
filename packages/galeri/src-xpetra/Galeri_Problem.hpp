// @HEADER
//
// ***********************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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

#ifndef GALERI_PROBLEM_HPP
#define GALERI_PROBLEM_HPP

#include <Teuchos_RCP.hpp>

#include "Galeri_ConfigDefs.h"

namespace Galeri {

  namespace Xpetra {

    template<typename Map, typename Matrix, typename MultiVector>
    class Problem : public Teuchos::Describable {
    public:
      Problem(Teuchos::ParameterList& list) : list_(list) { };
      Problem(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : list_(list) { Map_ = map; };
      virtual ~Problem() { }

      virtual Teuchos::RCP<Matrix>      BuildMatrix() = 0;
      virtual Teuchos::RCP<MultiVector> BuildNullspace() {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Nullspace is not implemented for this problem");
      }

      // Get methods
      Teuchos::RCP<const Map>         getMap()       const { return Map_; }
      Teuchos::RCP<const Matrix>      getMatrix()    const { return A_; }
      Teuchos::RCP<const MultiVector> getNullspace() const { return Nullspace_; }

      // Set methods
      // Teuchos::RCP<const Matrix> setCoords(const Teuchos::RCP<const Coords>& coords) { Coords_ = coords; }
      Teuchos::RCP<const Map> setMap   (const Teuchos::RCP<const Map>& map)       { Map_ = map; }

    protected:
      Teuchos::ParameterList&                list_;
      Teuchos::RCP<const Map>                Map_;
      Teuchos::RCP<Matrix>                   A_;
      Teuchos::RCP<MultiVector>              Nullspace_;
    };

  } // namespace Xpetra

} // namespace Galeri

#endif // GALERI_PROBLEM_HPP
