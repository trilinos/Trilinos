// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_INCLUDE_AKRI_SNAP_H_
#define KRINO_INCLUDE_AKRI_SNAP_H_
#include <Akri_FieldRef.hpp>
#include <Akri_NodeToCapturedDomains.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Part.hpp>

namespace krino
{
class InterfaceGeometry;
class QualityMetric;

NodeToCapturedDomainsMap snap_as_much_as_possible_while_maintaining_quality(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const FieldSet & interpolationFields,
    const InterfaceGeometry & geometry,
    const bool globalIDsAreParallelConsistent);

double determine_quality(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const QualityMetric &qualityMetric);

}



#endif /* KRINO_INCLUDE_AKRI_SNAP_H_ */
