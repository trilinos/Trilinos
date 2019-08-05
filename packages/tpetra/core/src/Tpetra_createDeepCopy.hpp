#ifndef TPETRA_CREATEDEEPCOPY_HPP
#define TPETRA_CREATEDEEPCOPY_HPP

#include "TpetraCore_config.h"

namespace Tpetra {

template<class OutputType, class InputType>
OutputType createDeepCopy (const InputType& in);

} // namespace Tpetra

#endif // TPETRA_CREATEDEEPCOPY_HPP
