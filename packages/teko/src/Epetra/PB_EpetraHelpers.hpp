#ifndef __PB_EpetraHelpers_hpp__
#define __PB_EpetraHelpers_hpp__

// stl includes
#include <string>

// Epetra includes
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_VectorBase.hpp"
#include "Thyra_DefaultSpmdMultiVector.hpp"

namespace PB {

namespace Epetra {

/** \brief Builds an epetra operator from a block 2x2 matrix.
  *
  * Builds an epetra operator from a block 2x2 matrix.
  * <em>*** The calling user is responsible for deleting the resulting Epetra_Operator! ***</em>
  */
Epetra_Operator * block2x2(const Epetra_Operator * sub00,const Epetra_Operator * sub01,
                           const Epetra_Operator * sub10,const Epetra_Operator * sub11,
                           const std::string & str="ANYM");

/** \brief Swaps the Apply/ApplyInverse to ApplyInverse/Apply.
  *
  * Swaps the Apply/ApplyInverse operations to ApplyInverse/Apply.
  * <em>*** The calling user is responsible for deleting the resulting Epetra_Operator! ***</em>
  */
Epetra_Operator * mechanicalInverse(const Epetra_Operator * inverse);

/** \brief Fill a Thyra vector with the contents of an epetra vector. This prevents the
  *
  * Fill a Thyra vector with the contents of an epetra vector. This prevents the need
  * to reallocate memory using a create_MultiVector routine. It also allows an aritrary
  * Thyra vector to be filled.
  *
  * \param[in,out] spmdMV Multi-vector to be filled.
  * \param[in]     mv     Epetra multi-vector to be used in filling the Thyra vector.
  */    
void fillDefaultSpmdMultiVector(Teuchos::RCP<Thyra::DefaultSpmdMultiVector<double> > & spmdMV,
                                Teuchos::RCP<Epetra_MultiVector> & epetraMV);

const Teuchos::RCP<const Thyra::LinearOpBase<double> > thyraDiagOp(const Teuchos::RCP<const Epetra_Vector> & ev,const Epetra_Map & map,const std::string & lbl="ANYM");

} // end namespace Epetra
} // end namespace PB

#endif
