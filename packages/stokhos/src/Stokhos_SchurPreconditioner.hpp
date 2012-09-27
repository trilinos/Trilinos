

#ifndef STOKHOS_SCHURPRECONDITIONER_HPP
#define STOKHOS_SCHURPRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_Operator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class SchurPreconditioner : 
    public Stokhos::Operator<ordinal_type, value_type> {
  public:

    //! Constructor 
    SchurPreconditioner(
      const Teuchos::SerialDenseMatrix<ordinal_type,value_type> & K, 
      const ordinal_type p, const ordinal_type m, const ordinal_type diag);
    
    //! Destructor
    virtual ~SchurPreconditioner(); 
    
    virtual ordinal_type ApplyInverse(
      const Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Input, 
      Teuchos::SerialDenseMatrix<ordinal_type, value_type>& Result, 
      ordinal_type prec_iters) const;
      
  protected:
    ordinal_type fact(ordinal_type n) const;
    ordinal_type size(ordinal_type n, ordinal_type m) const;

    const Teuchos::SerialDenseMatrix<ordinal_type,value_type> & K;
    const ordinal_type p;
    const ordinal_type m;
    const ordinal_type diag;

  }; // class SchurPreconditioner

} // namespace Stokhos

#include "Stokhos_SchurPreconditionerImp.hpp"

#endif // STOKHOS_SCHURPRECONDITIONER_HPP

