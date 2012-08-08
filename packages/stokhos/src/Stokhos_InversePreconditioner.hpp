

#ifndef STOKHOS_INVERSEPRECONDITIONER_HPP
#define STOKHOS_INVERSEPRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_Operator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {


  class InversePreconditioner : public Stokhos::Operator {

  public:

    //! Constructor 
    InversePreconditioner(const Teuchos::SerialDenseMatrix<int,double> & A);
    
  
    //! Destructor
    virtual ~InversePreconditioner(); 
    
  
    virtual int ApplyInverse(const Teuchos::SerialDenseMatrix<int,double>& Input,
                             Teuchos::SerialDenseMatrix<int,double>& Output, int m) const;
   
  protected:
    const Teuchos::SerialDenseMatrix<int,double> & A;
  }; // class InversePreconditioner

} // namespace Stokhos

#endif // STOKHOS_INVERSEPRECONDITIONER_HPP

