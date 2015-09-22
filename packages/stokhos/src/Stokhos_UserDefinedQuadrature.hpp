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

#ifndef STOKHOS_USERDEFINEDQUADRATURE
#define STOKHOS_USERDEFINEDQUADRATURE

#include "Stokhos_Quadrature.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Teuchos_RCP.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class UserDefinedQuadrature : public Quadrature<ordinal_type,value_type> {
  public:

    //! Constructor
    UserDefinedQuadrature(
     const Teuchos::RCP<const OrthogPolyBasis<ordinal_type,value_type> >& basis,
     const Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > >& points,
     const Teuchos::RCP<const Teuchos::Array<value_type> >& weights);

    //! Constructor
    UserDefinedQuadrature(
     const Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > >& points,
     const Teuchos::RCP<const Teuchos::Array<value_type> >& weights,
     const Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > >& values);

    //! Destructor
    virtual ~UserDefinedQuadrature() {}

    //! Get number of quadrature points
    virtual ordinal_type size() const { return quad_weights->size(); }

    //! Get quadrature points
    virtual const Teuchos::Array< Teuchos::Array<value_type> >& 
    getQuadPoints() const;

    //! Get quadrature weights
    virtual const Teuchos::Array<value_type>& 
    getQuadWeights() const;

    //! Get values of basis at quadrature points
    virtual const Teuchos::Array< Teuchos::Array<value_type> > & 
    getBasisAtQuadPoints() const;

    //! Print quadrature data
    virtual std::ostream& print(std::ostream& os) const;

  private:

    // Prohibit copying
    UserDefinedQuadrature(const UserDefinedQuadrature&);

    // Prohibit Assignment
    UserDefinedQuadrature& operator=(const UserDefinedQuadrature& b);

  protected:

    //! Quadrature points
    Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > > quad_points;

    //! Quadrature weights
    Teuchos::RCP<const Teuchos::Array<value_type> > quad_weights;

    //! Quadrature values
    Teuchos::RCP<const Teuchos::Array< Teuchos::Array<value_type> > > quad_values;

  }; // class UserDefinedQuadrature

} // namespace Stokhos

// Include template definitions
#include "Stokhos_UserDefinedQuadratureImp.hpp"

#endif // STOKHOS_USERDEFINEDQUADRATURE
