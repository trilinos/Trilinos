// $Id$ 
// $Source$ 
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

#ifndef STOKHOS_BASIS_FACTORY_HPP
#define STOKHOS_BASIS_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"

namespace Stokhos {

  //! Factory for building multivariate orthogonal polynomial bases.
  template <typename ordinal_type, typename value_type>
  class BasisFactory {
  public:

    //! Constructor
    BasisFactory() {};

    //! Destructor
    virtual ~BasisFactory() {};

    //! Generate multivariate basis
    static Teuchos::RCP<const Stokhos::OrthogPolyBasis<ordinal_type, value_type> >
    create(Teuchos::ParameterList& sgParams);

  protected:

    //! Generate 1-D basis
    static Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<ordinal_type, value_type> > create1DBasis(Teuchos::ParameterList& params);

  private:

    // Prohibit copying
    BasisFactory(const BasisFactory&);

    // Prohibit Assignment
    BasisFactory& operator=(const BasisFactory& b);

  }; // class BasisFactory

} // Namespace Stokhos

// Include template definitions
#include "Stokhos_BasisFactoryImp.hpp"

#endif // STOKHOS_BASIS_FACTORY_HPP
