//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   <pkg>/Release.cc
 * \author <user>
 * \date   <date>
 * \brief  Release function implementation for <pkg> library
 * \note   Copyright (C) 2008 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
// $Id: Release.cc,v 1.3 2008/01/02 17:18:47 9te Exp $
//---------------------------------------------------------------------------//

#include "Release.hh"

namespace <namespace>
{

using std::string;

/*!  
 * \return string of the release number
 *
 * Function definition for Release, define the local version number for
 * this library in the form <spkg>-\#_\#_\# in pkg_release variable 
 */
const string release()
{
    string pkg_release = "<spkg>(nemesis-<start>#_#_#)";
    return pkg_release;
}

}  // end of <namespace>

//---------------------------------------------------------------------------//
//                             end of Release.cc
//---------------------------------------------------------------------------//
