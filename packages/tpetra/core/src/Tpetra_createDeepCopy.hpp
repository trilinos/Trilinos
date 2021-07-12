#ifndef TPETRA_CREATEDEEPCOPY_HPP
#define TPETRA_CREATEDEEPCOPY_HPP

#ifdef TPETRA_ENABLE_DEPRECATED_CODE

#include "TpetraCore_config.h"

namespace Tpetra {

template<class OutputType, class InputType>
OutputType 
TPETRA_DEPRECATED
createDeepCopy (const InputType& in);

} // namespace Tpetra

#endif // TPETRA_ENABLE_DEPRECATED_CODE

#endif // TPETRA_CREATEDEEPCOPY_HPP
