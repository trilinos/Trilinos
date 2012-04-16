// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_IO.cpp
 *  \brief Definition of methods to assist in file input/output.
 *
 *   \todo write a solution and its model to an exodus file
 */

#include <Zoltan2_IO.hpp>

namespace Zoltan2{

/*! \brief Helper method to add number to a file name.
 *    \param number the number (such as process rank) to add
 *    \param fname  the file name to modify
 *    \param newf   on return newf is fname with the rank added to the name.
 *
 *   If fname has no dot in it, then the rank is added to the end of the name.
 *   Otherwise the rank is added before the first dot in fname.
 */

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

