// ////////////////////////////////////////////////////
// ref_count_ptr.cpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ThrowException.hpp"

void Teuchos::PrivateUtilityPack::assert_not_null(const void *ptr)
{
	THROW_EXCEPTION(
		!ptr, std::logic_error
		,"ref_count_ptr<...>::assert_not_null() : You can not "
		" call operator->() or operator*() if get() == 0" );
}
