// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_ORTHOGPOLYEXPANSIONBASE_HPP
#define STOKHOS_ORTHOGPOLYEXPANSIONBASE_HPP

#include "Stokhos_OrthogPolyExpansion.hpp"
#include "Stokhos_DivisionExpansionStrategy.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Stokhos {

  //! Base class for consolidating common expansion implementations
  /*!
   * Implements a few of the OrthogPolyExpansion virtual methods that are 
   * common to several implementations so those implementations are duplicated.
   */
  template <typename ordinal_type, typename value_type, typename node_type> 
  class OrthogPolyExpansionBase : 
    public OrthogPolyExpansion<ordinal_type, value_type, node_type> {
  public:

    //! Constructor
    OrthogPolyExpansionBase(
      const Teuchos::RCP<const OrthogPolyBasis<ordinal_type, value_type> >& basis,
      const Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> >& Cijk,
      const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);

    //! Destructor
    virtual ~OrthogPolyExpansionBase() {}

    //! Get expansion size
    ordinal_type size() const { return sz; }

    //! Get basis
    Teuchos::RCP< const OrthogPolyBasis<ordinal_type, value_type> > 
    getBasis() const {return basis; }

    //! Get triple product
    virtual Teuchos::RCP<const Sparse3Tensor<ordinal_type, value_type> >
    getTripleProduct() const { return Cijk; }
 
    // Operations
    void unaryMinus(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);

    void plusEqual(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& x);
    void minusEqual(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& x);
    void timesEqual(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& x);
    void divideEqual(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const value_type& x);

    void plusEqual(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const OrthogPolyApprox<ordinal_type, value_type, node_type>& x);
    void minusEqual(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const OrthogPolyApprox<ordinal_type, value_type, node_type>& x);
    void timesEqual(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const OrthogPolyApprox<ordinal_type, value_type, node_type>& x);
    void divideEqual(
      OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
      const OrthogPolyApprox<ordinal_type, value_type, node_type>& x);

    void plus(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
              const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void plus(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
              const value_type& a, 
              const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void plus(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
              const value_type& b);
    void minus(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void minus(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const value_type& a, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void minus(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
               const value_type& b);
    void times(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void times(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const value_type& a, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void times(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
               const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
               const value_type& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
                const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
                const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
                const value_type& a, 
                const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void divide(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
                const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
                const value_type& b);

    
    void abs(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void fabs(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
              const OrthogPolyApprox<ordinal_type, value_type, node_type>& a);
    void max(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void max(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const value_type& a, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void max(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
             const value_type& b);
    void min(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& a,
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void min(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const value_type& a, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& b);
    void min(OrthogPolyApprox<ordinal_type, value_type, node_type>& c, 
             const OrthogPolyApprox<ordinal_type, value_type, node_type>& a, 
             const value_type& b);

  private:

    // Prohibit copying
    OrthogPolyExpansionBase(const OrthogPolyExpansionBase&);

    // Prohibit Assignment
    OrthogPolyExpansionBase& operator=(const OrthogPolyExpansionBase& b);

  protected:

     //! Basis
    Teuchos::RCP<const OrthogPolyBasis<ordinal_type, value_type> > basis;

    //! Short-hand for Cijk
    typedef Stokhos::Sparse3Tensor<ordinal_type, value_type> Cijk_type;

    //! Triple-product tensor
    Teuchos::RCP<const Stokhos::Sparse3Tensor<ordinal_type, value_type> > Cijk;

    //! Parameter list
    Teuchos::RCP<Teuchos::ParameterList> params;

    //! Division expansion strategy
    Teuchos::RCP<Stokhos::DivisionExpansionStrategy<ordinal_type,value_type,node_type> > division_strategy;

    //! Expansions size
    ordinal_type sz;
    
  }; // class OrthogPolyExpansionBase

} // namespace Stokhos

#include "Stokhos_OrthogPolyExpansionBaseImp.hpp"

#endif // STOKHOS_ORTHOGPOLYEXPANSIONBASE_HPP
