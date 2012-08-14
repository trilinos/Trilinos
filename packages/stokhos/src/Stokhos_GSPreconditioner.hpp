

#ifndef STOKHOS_GSPRECONDITIONER_HPP
#define STOKHOS_GSPRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_Operator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {


  class GSPreconditioner : public Stokhos::Operator {

  public:

    //! Constructor 
    GSPreconditioner(const Teuchos::SerialDenseMatrix<int,double> & A, const int sym);
    
  
    //! Destructor
    virtual ~GSPreconditioner(); 
    
  
    virtual int ApplyInverse(const Teuchos::SerialDenseMatrix<int,double>& Input,
                             Teuchos::SerialDenseMatrix<int,double>& Output, int m) const;
   
  protected:
    const Teuchos::SerialDenseMatrix<int,double> & A;

    const int sym;
  }; // class GSPreconditioner

} // namespace Stokhos

#endif // STOKHOS_GSPRECONDITIONER_HPP

