// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef EvaluateGregoryPatch_hpp
#define EvaluateGregoryPatch_hpp

#include <map>
#include <vector>

#include <percept/PerceptMesh.hpp>

/** Evaluate Gregory patches
 */
namespace percept {

class EvaluateGregoryPatch
{
public:
  EvaluateGregoryPatch(PerceptMesh& eMesh, bool debug = false)
    : m_eMesh(eMesh), m_debug(debug)
  {
  }

  // return true if error - this is a pure-STK interface
  static bool
  evaluateGregoryPatch(const stk::mesh::BulkData& bulk, const double *uv, stk::mesh::Entity face, double *xyz,
           const GregoryControlPointsType *sideset_field = 0, const GregoryControlPointsType *shell_field=0);

  static bool
  normalGregoryPatch(const stk::mesh::BulkData& bulk, const double *uv, stk::mesh::Entity face, double *normal,
            const GregoryControlPointsType *sideset_field = 0, const GregoryControlPointsType *shell_field = 0);

  // return true if error
  bool
  evaluateGregoryPatch(const double *uv, stk::mesh::Entity face, double *xyz);

  bool
  normalGregoryPatch(const double *uv, stk::mesh::Entity face, double *normal);

  // return true if error
  bool
  findClosestPoint(const double *input_xyz, stk::mesh::Entity face, double *closest_xyz, double *found_uv=0,
                   bool linearOnly=false, stk::mesh::Entity node = stk::mesh::Entity());

  static void evaluateLinearQuad(PerceptMesh& eMesh, const double *uv, stk::mesh::Entity face, double *xyz, double G[3][2], double H[3][2][2]);

  const std::string& error_message() { return m_errorMsg; }

private:

  void project_to_face_boundary(stk::mesh::Entity face, double *xyz_in, double *uv) const;

  // found_uv is also the initial guess
  template<class F, class G, class H>
  bool
  findClosestPointInternal(const double *input_xyz, stk::mesh::Entity face, double *closest_xyz, double *found_uv);

  PerceptMesh& m_eMesh;
  bool m_debug;
  std::string m_errorMsg;
};

}

#endif
