#include "NS/PB_StaticLSCStrategy.hpp"

#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_EpetraLinearOp.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"

#include "Teuchos_Time.hpp"

// Teko includes
#include "PB_Utilities.hpp"

using Teuchos::RCP;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcp_const_cast;

namespace Teko {
namespace NS {

   // Staiblized constructor
StaticLSCStrategy::StaticLSCStrategy(const LinearOp & invF,
                                     const LinearOp & invBQBtmC,
                                     const LinearOp & invD,
                                     const LinearOp & invMass)
   : invF_(invF), invBQBtmC_(invBQBtmC), invD_(invD), invMass_(invMass)
{ }
 
   // Stable constructor
StaticLSCStrategy::StaticLSCStrategy(const LinearOp & invF,
                                     const LinearOp & invBQBtmC,
                                     const LinearOp & invMass)
   : invF_(invF), invBQBtmC_(invBQBtmC), invD_(Teuchos::null), invMass_(invMass)
{ }

} // end namespace NS
} // end namespace Teko
