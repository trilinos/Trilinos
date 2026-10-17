// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_PHASE_AKRI_PHASETAGREGISTRY_HPP_
#define KRINO_KRINO_PHASE_AKRI_PHASETAGREGISTRY_HPP_

#include <Akri_PhaseTag.hpp>

#include <set>

namespace krino {

// Holds one copy of each distinct PhaseTag, so that the many objects that share a phase (for example
// the elements of a CDFEM mesh) can each store a pointer to it instead of their own copy.  A
// registered tag never moves or changes, so its address stays valid for the life of the registry.
class PhaseTagRegistry
{
public:
  const PhaseTag & get(const PhaseTag & phase) { return *myPhases.insert(phase).first; }
  size_t size() const { return myPhases.size(); }

  static const PhaseTag & empty_phase();

private:
  std::set<PhaseTag> myPhases;
};

} // namespace krino

#endif /* KRINO_KRINO_PHASE_AKRI_PHASETAGREGISTRY_HPP_ */
