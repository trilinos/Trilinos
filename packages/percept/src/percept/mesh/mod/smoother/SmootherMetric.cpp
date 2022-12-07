// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/Percept.hpp>
#include <percept/MeshType.hpp>

#include <percept/mesh/mod/smoother/SmootherMetric.hpp>

namespace percept {
//madbrew: putting everything into header file for kokkkos inlining
template<>
SmootherMetricUntangleImpl<STKMesh>::
SmootherMetricUntangleImpl(PerceptMesh *eMesh) : SmootherMetricImpl<STKMesh>(eMesh) {
	m_beta_mult = 0.05;
}

template<>
double SmootherMetricUntangleImpl<STKMesh>::
metric(typename STKMesh::MTElement element, bool& valid)
{
	valid = true;

	JacobianUtilImpl<STKMesh> jacA, jacW;

	double A_ = 0.0, W_ = 0.0; // current and reference detJ
	jacA(A_, *Base::m_eMesh, element, Base::m_coord_field_current, Base::m_topology_data);
	jacW(W_, *Base::m_eMesh, element, Base::m_coord_field_original, Base::m_topology_data);
	double val_untangle=0.0;

	for (int i=0; i < jacA.m_num_nodes; i++)
	{
		double detAi = jacA.m_detJ[i];
		double detWi = jacW.m_detJ[i];

		if (detAi <= 0.) valid = false;

		double vv = m_beta_mult*detWi - detAi;
		vv = std::max(vv, 0.0);
		val_untangle += vv*vv;
	}
	return val_untangle;
}

} // namespace percept
