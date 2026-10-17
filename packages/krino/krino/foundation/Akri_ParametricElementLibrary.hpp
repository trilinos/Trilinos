// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef Akri_ParametricElementLibrary_h
#define Akri_ParametricElementLibrary_h

#include <Akri_ParametricElement.hpp>

#include <memory>
#include <vector>

namespace krino { class FieldRef; }
namespace stk { namespace mesh { class Bucket; } }
namespace stk { struct topology; }

namespace krino {

// Owns the cache of ParametricElements, one per stk::topology, lazily constructed on
// first request.  The storage is a member rather than static, so its lifetime (and
// the lifetime of the Kokkos Views held by the Basis objects inside each
// ParametricElement) is controlled by the owner.
//
// NOT THREAD SAFE: the query methods mutate their caches on a miss without any
// synchronization.  All current callers use them from serial host setup code.  Do
// not call concurrently on the same library instance.
class ParametricElementLibrary
{
public:
  ParametricElementLibrary();
  ~ParametricElementLibrary();

  ParametricElementLibrary(const ParametricElementLibrary &) = delete;
  ParametricElementLibrary & operator=(const ParametricElementLibrary &) = delete;

  // Returns the cached ParametricElement for this topology, creating it on first use.
  // The returned reference is stable for the lifetime of this library: the cache
  // owns each element through a unique_ptr and never relocates it.
  const ParametricElement & get_parametric_element(stk::topology topology) const;

  // Returns the ParametricElement for the topology of the given nodal field on the bucket.
  const ParametricElement & get_parametric_element(const stk::mesh::Bucket & bucket, const FieldRef field) const;

  // Resolve the topology of a nodal field on a bucket (memoized internally).
  stk::topology get_field_topology(const stk::mesh::Bucket & b, const FieldRef field) const;

  // Destroy all cached ParametricElements.  Needed (at least with ParametricElementIntrepid)
  // to release the Kokkos Views held by Basis objects before Kokkos::finalize().
  void clear();

private:
  // mutable: the query methods are logically const (memoization) but fill on miss.
  mutable std::vector<std::unique_ptr<ParametricElement>> myParametricElements;
};

}  // namespace krino

#endif // Akri_ParametricElementLibrary_h
