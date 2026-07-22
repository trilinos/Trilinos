// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_StringInputSource.hpp"
#include "Teuchos_StringInputStream.hpp"

using namespace Teuchos;


StringInputSource::StringInputSource(const std::string& text)
	: XMLInputSource(), text_(text)
{;}

RCP<XMLInputStream> StringInputSource::stream() const
{
	return rcp(new StringInputStream(text_), true);
}

