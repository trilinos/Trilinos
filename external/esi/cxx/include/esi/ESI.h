#ifndef __ESI_h_seen
#define __ESI_h_seen

/* ESI.h for C++
 * A library header for the equation solver interfaces.
 *
 */

#define ESI_MAJOR_VERSION    1
#define ESI_MINOR_VERSION    0
#define ESI_PATCH_LEVEL      0

#include "../esi/basicTypes.h"
#include "../esi/ordinalTraits.h"
#include "../esi/scalarTraits.h"

// core ESI interfaces
#include "../esi/Argv.h"

#include "../esi/Object.h"

#include "../esi/IndexSpace.h"

#include "../esi/Vector.h"
#include "../esi/VectorReplaceAccess.h"

#include "../esi/Operator.h"
#include "../esi/OperatorTranspose.h"

#include "../esi/MatrixData.h"
#include "../esi/MatrixRowReadAccess.h"
#include "../esi/MatrixRowWriteAccess.h"
#include "../esi/MatrixRowPointerAccess.h"

#include "../esi/Preconditioner.h"
#include "../esi/PreconditionerTranspose.h"

#include "../esi/Solver.h"
#include "../esi/SolverIterative.h"

#endif /* __ESI_h_seen */
