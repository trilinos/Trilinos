

#ifndef STOKHOS_DIAGPRECONDITIONER_HPP
#define STOKHOS_DIAGPRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_Operator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {


  class DiagPreconditioner : public Stokhos::Operator {

  public:

    //! Constructor 
    DiagPreconditioner(const Teuchos::SerialDenseMatrix<int,double> & A);
    
  
    //! Destructor
    virtual ~DiagPreconditioner(); 
    
  
    virtual int ApplyInverse(const Teuchos::SerialDenseMatrix<int,double>& Input,
                             Teuchos::SerialDenseMatrix<int,double>& Output, int m) const;
   
  protected:
    const Teuchos::SerialDenseMatrix<int,double> & A;
  }; // class DiagPreconditioner

} // namespace Stokhos

#endif // STOKHOS_DIAGPRECONDITIONER_HPP

