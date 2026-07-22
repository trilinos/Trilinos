// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_REDUCED_PCE_BASIS_HPP
#define STOKHOS_REDUCED_PCE_BASIS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Array.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_Quadrature.hpp"

namespace Stokhos {

  /*! 
   * \brief Abstract base class for reduced basis strategies built from 
   * polynomial chaos expansions in some other basis.
   */
  template <typename ordinal_type, typename value_type>
  class ReducedPCEBasis : 
    public virtual OrthogPolyBasis<ordinal_type,value_type> {
  public:

    //! Default constructor
    ReducedPCEBasis() {}

    //! Destructor
    virtual ~ReducedPCEBasis() {}

    //! \name ReducedBasis virtual methods
    //@{

    //! Transform coefficients to original basis from this basis
    virtual void 
    transformToOriginalBasis(const value_type *in, 
			     value_type *out,
			     ordinal_type ncol = 1, 
			     bool transpose = false) const = 0;

    //! Transform coefficients from original basis to this basis
    virtual void 
    transformFromOriginalBasis(const value_type *in, 
			       value_type *out,
			       ordinal_type ncol = 1, 
			       bool transpose = false) const = 0;

    //! Get reduced quadrature object
    virtual Teuchos::RCP<const Stokhos::Quadrature<ordinal_type, value_type> >
    getReducedQuadrature() const = 0;

    //@}

  private:

    // Prohibit copying
    ReducedPCEBasis(const ReducedPCEBasis&);

    // Prohibit Assignment
    ReducedPCEBasis& operator=(const ReducedPCEBasis&);

  }; // class ReducedPCEBasis

} // Namespace Stokhos

#endif
