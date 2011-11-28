// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#ifndef _ZOLTAN2_MEMORY_HPP_
#define _ZOLTAN2_MEMORY_HPP_

/*! \file Zoltan2_Memory.hpp

  \brief Memory related declarations.
*/

#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Ptr.hpp>
#include <Zoltan2_config.h>

namespace Zoltan2{

#ifdef HAVE_MALLINFO
// Note: Calling Zoltan's meminfo helpers works better
// on Linux nodes.
int getAllocatedMemory();
int mallocCount(std::string label);
int getMallocCount(std::string label);
void printMallocCount();
void eraseMallocCount();
#endif


} //namespace Zoltan2

#endif


