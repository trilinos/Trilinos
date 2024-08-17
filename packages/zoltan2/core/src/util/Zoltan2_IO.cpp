// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_IO.cpp
 *  \brief Definition of methods to assist in file input/output.
 *
 *   \todo write a solution and its model to an exodus file
 */

#include <Zoltan2_IO.hpp>

namespace Zoltan2{

void addNumberToFileName(int number, std::string fname, std::string &newf)
{
  std::ostringstream id;
  id << "_";
  id.width(6);
  id.fill('0');
  id << number;

  std::ostringstream localFileName;
  std::string::size_type loc = fname.find('.');

  if (loc == std::string::npos)
    localFileName << fname << id.str();
  else
    localFileName << fname.substr(0, loc) << id.str() << fname.substr(loc);

  newf = localFileName.str();
}

} // namespace Zoltan2

