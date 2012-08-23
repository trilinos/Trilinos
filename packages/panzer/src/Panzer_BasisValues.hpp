// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_BASIS_VALUES_DECL_HPP
#define PANZER_BASIS_VALUES_DECL_HPP

#include "Teuchos_RCP.hpp"

#include "Intrepid_Basis.hpp"
#include "Panzer_BasisIRLayout.hpp"

namespace panzer {

  template<typename Scalar,typename Array>
  struct BasisValues { 
    static const Array dummyArray;    
 
    //! Sizes/allocates memory for arrays
    template <typename ArrayFactory>
    void setupArrays(const Teuchos::RCP<const panzer::BasisIRLayout>& basis,
                     const ArrayFactory & af);

    void evaluateValues(const Array& cub_points,
			const Array& jac,
			const Array& jac_det,
			const Array& jac_inv);

    void evaluateValues(const Array& cub_points,
			const Array& jac,
			const Array& jac_det,
			const Array& jac_inv,
			const Array& node_coordinates);

    void evaluateValues(const Array& cub_points,
			const Array& jac,
			const Array& jac_det,
			const Array& jac_inv,
			const Array& weighted_measure,
			const Array& node_coordinates);

    PureBasis::EElementSpace getElementSpace() const; 

    Array basis_ref;           // <BASIS,IP>
    Array basis;               // <Cell,BASIS,IP>
    Array grad_basis_ref;      // <BASIS,IP,Dim>
    Array grad_basis;          // <Cell,BASIS,IP,Dim>
    Array curl_basis_ref;      // <BASIS,IP,Dim> (dimension dependent)
    Array curl_basis;          // <Cell,BASIS,IP,Dim> (dimension dependent)
    Array weighted_basis;      // <Cell,BASIS,IP>
    Array weighted_grad_basis; // <Cell,BASIS,IP,Dim>
    Array weighted_curl_basis; // <Cell,BASIS,IP,Dim> (dimension dependent)

    /** Carterisan coordinates for basis coefficients

        NOTE: This quantity is not always available.  Certain bases
        may not have a corresponding coordiante value
    */
    Array basis_coordinates_ref;     // <Cell,BASIS>

    /** Carterisan coordinates for basis coefficients

        NOTE: This quantity is not always available.  Certain bases
        may not have a corresponding coordiante value
    */
    Array basis_coordinates;         // <Cell,BASIS,Dim>

    Teuchos::RCP<const panzer::BasisIRLayout> basis_layout;
    
    Teuchos::RCP<Intrepid::Basis<Scalar,Array> > intrepid_basis;
  };

} // namespace panzer

#include "Panzer_BasisValues_impl.hpp"

#endif
