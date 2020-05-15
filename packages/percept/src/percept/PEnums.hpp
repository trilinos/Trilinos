// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef PEnums_hpp
#define PEnums_hpp

namespace percept {

  enum NodeRegistryFieldData {
    NR_FIELD_OWNING_ELEMENT_ID     = 0,
    NR_FIELD_OWNING_ELEMENT_RANK   = 1,
    NR_FIELD_MARK                  = 2,
    NR_FIELD_OWNING_SUBDIM_RANK    = 3,
    NR_FIELD_OWNING_SUBDIM_ORDINAL = 4,
    NR_FIELD_OWNING_SUBDIM_SIZE    = 5,
    NUM_NR_FIELD_SLOTS
  };


}

#endif
