// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_Details_MpiCommRequest.hpp>

namespace Teuchos {
namespace Details {

MpiCommRequest::
MpiCommRequest (MPI_Request rawMpiRequest,
		const ArrayRCP<const char>& buffer) :
  MpiCommRequestBase<int> (rawMpiRequest),
  buffer_ (buffer)
{}

MpiCommRequest::~MpiCommRequest () {}

RCP<MpiCommRequest>
mpiCommRequest (MPI_Request rawMpiRequest,
		const ArrayRCP<const char>& buffer)
{
  return rcp (new MpiCommRequest (rawMpiRequest, buffer));
}

} // namespace Details
} // namespace Teuchos

