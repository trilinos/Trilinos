/*--------------------------------------------------------------------*/
/*    Copyright 2008 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// -*- Mode: c++ -*-
#ifndef SIERRA_Ioss_SurfaceSplit_h
#define SIERRA_Ioss_SurfaceSplit_h
namespace Ioss {
  enum SurfaceSplitType {
    SPLIT_INVALID = -1,
    SPLIT_BY_TOPOLOGIES = 1,
    SPLIT_BY_ELEMENT_BLOCK = 2,
    SPLIT_BY_DONT_SPLIT = 3,
    SPLIT_LAST_ENTRY = 4
  };

  inline SurfaceSplitType int_to_surface_split(int split_int) {
    SurfaceSplitType split_type = Ioss::SPLIT_INVALID;
    if (split_int == 1) split_type = Ioss::SPLIT_BY_TOPOLOGIES;
    if (split_int == 2) split_type = Ioss::SPLIT_BY_ELEMENT_BLOCK;
    if (split_int == 3) split_type = Ioss::SPLIT_BY_DONT_SPLIT;
    return split_type;
  }
}
#endif
