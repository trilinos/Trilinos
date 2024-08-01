// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
