// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef PYTRILINOS_DOMI_UTIL_H
#define PYTRILINOS_DOMI_UTIL_H

// Include the PyTrilinos Distributed Array Protocol header
#include "PyTrilinos_DAP.hpp"

// Include Domi headers
#include "Domi_MDMap.hpp"

namespace PyTrilinos
{

////////////////////////////////////////////////////////////////////////

// Given a 'distarray' object returned by the DistArray protocol,
// convert to an RCP of a Domi MDComm.
Teuchos::RCP< const Domi::MDComm >
convertToMDComm(const Domi::TeuchosCommRCP teuchosComm,
                const DistArrayProtocol & distarray);

////////////////////////////////////////////////////////////////////////

// Given a 'distarray' object returned by the DistArray protocol,
// convert to an RCP of a Domi MDMap.
Teuchos::RCP< const Domi::MDMap<> >
convertToMDMap(const Domi::TeuchosCommRCP teuchosComm,
               const DistArrayProtocol & distarray);

////////////////////////////////////////////////////////////////////////

}
