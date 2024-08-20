// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_FileInputStream.hpp"

using namespace Teuchos;


FileInputSource::FileInputSource(const std::string& filename)
	: XMLInputSource(), filename_(filename)
{;}

RCP<XMLInputStream> FileInputSource::stream() const
{
	return rcp(new FileInputStream(filename_), true);
}

