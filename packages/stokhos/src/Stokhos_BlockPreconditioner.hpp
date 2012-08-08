

#ifndef STOKHOS_BLOCKPRECONDITIONER_HPP
#define STOKHOS_BLOCKPRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Stokhos_Operator.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

namespace Stokhos {


  class BlockPreconditioner : public Stokhos::Operator {

  public:

    //! Constructor 
    BlockPreconditioner(const Teuchos::SerialDenseMatrix<int,double> & K, const int p, const int m);
    
  
    //! Destructor
    virtual ~BlockPreconditioner(); 
    
  
    virtual int ApplyInverse(const Teuchos::SerialDenseMatrix<int,double>& Input,
                             Teuchos::SerialDenseMatrix<int,double>& Output, int m) const;
   
  protected:
    const Teuchos::SerialDenseMatrix<int,double> & K;
    const int p;
    const int m;

  }; // class BlockPreconditioner

} // namespace Stokhos

#endif // STOKHOS_BLOCKPRECONDITIONER_HPP

