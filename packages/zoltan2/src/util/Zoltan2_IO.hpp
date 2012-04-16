// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_IO.hpp
 *  \brief Declaration of methods to assist in file input/output.
 */

#ifndef ZOLTAN2_IO_HPP
#define ZOLTAN2_IO_HPP

#include <string>
#include <ostream>
#include <sstream>

namespace Zoltan2{

/*! \brief Helper method to add number to a file name.
 *    \param number the number (such as process rank) to add
 *    \param fname  the file name to modify
 *    \param newf   on return newf is fname with the rank added to the name.
 *
 *   If fname has no dot in it, then the rank is added to the end of the name.
 *   Otherwise the rank is added before the first dot in fname.
 */

void addNumberToFileName(int number, std::string fname, std::string &newf);

} // namespace Zoltan2

#endif
