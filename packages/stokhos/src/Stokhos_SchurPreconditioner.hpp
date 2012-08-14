

#ifndef STOKHOS_SCHURPRECONDITIONER_HPP
#define STOKHOS_SCHURPRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_Operator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {


  class SchurPreconditioner : public Stokhos::Operator {

  public:

    //! Constructor 
    SchurPreconditioner(const Teuchos::SerialDenseMatrix<int,double> & K, const int p, const int m, const int diag);
    
  
    //! Destructor
    virtual ~SchurPreconditioner(); 
    
  
    virtual int ApplyInverse(const Teuchos::SerialDenseMatrix<int,double>& Input,
                             Teuchos::SerialDenseMatrix<int,double>& Output, int prec_iters) const;
      
  protected:
    const Teuchos::SerialDenseMatrix<int,double> & K;
    const int p;
    const int m;
    const int diag;

  }; // class SchurPreconditioner

} // namespace Stokhos

#endif // STOKHOS_SCHURPRECONDITIONER_HPP

