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

#ifndef GALERI_PROBLEM_HELMHOLTZ_HPP
#define GALERI_PROBLEM_HELMHOLTZ_HPP

#include <Teuchos_RCP.hpp>

#include "Galeri_ConfigDefs.h"
#include "Galeri_Problem.hpp"

namespace Galeri {

  namespace Xpetra {

    template<typename Map, typename Matrix, typename MultiVector>
    class Problem_Helmholtz : public Problem<Map,Matrix,MultiVector> {
    public:
      Problem_Helmholtz(Teuchos::ParameterList& list)                                     : Problem<Map,Matrix,MultiVector>(list) { }
      Problem_Helmholtz(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) : Problem<Map,Matrix,MultiVector>(list, map) { }
      virtual ~Problem_Helmholtz() { }

      virtual std::pair< Teuchos::RCP<Matrix>, Teuchos::RCP<Matrix> > BuildMatrices() = 0;

      // Get methods
      Teuchos::RCP<const Matrix>      getStiff()     const { return K_; }
      Teuchos::RCP<const Matrix>      getMass()      const { return M_; }

    protected:
      Teuchos::RCP<Matrix> K_, M_;
    };

  } // namespace Xpetra

} // namespace Galeri

#endif // GALERI_PROBLEM_HELMHOLTZ_HPP
