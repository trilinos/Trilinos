

#ifndef STOKHOS_JACOBIPRECONDITIONER_HPP
#define STOKHOS_JACOBIPRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_Operator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {


  class JacobiPreconditioner : public Stokhos::Operator {

  public:

    //! Constructor 
    JacobiPreconditioner(const Teuchos::SerialDenseMatrix<int,double> & A);
    
  
    //! Destructor
    virtual ~JacobiPreconditioner(); 
    
  
    virtual int ApplyInverse(const Teuchos::SerialDenseMatrix<int,double>& Input,
                             Teuchos::SerialDenseMatrix<int,double>& Output, int m) const;
   
  protected:
    const Teuchos::SerialDenseMatrix<int,double> & A;
  }; // class JacobiPreconditioner

} // namespace Stokhos

#endif // STOKHOS_JACOBIPRECONDITIONER_HPP

