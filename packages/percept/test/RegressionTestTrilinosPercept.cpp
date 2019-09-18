// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <gtest/gtest.h>

#include <percept/PerceptMesh.hpp>
#include <adapt/UniformRefinerPattern.hpp>
#include <adapt/UniformRefiner.hpp>

TEST(UniformRefineNoGeometry, BDot)
{
  using namespace percept;

  PerceptMesh eMesh(3u);
  const std::string fileroot="BDot";
  eMesh.open(fileroot+".g");
  
  Tet4_Tet4_8 break_pattern(eMesh);
  
  eMesh.commit();

  UniformRefiner breaker(eMesh, break_pattern);
  breaker.doBreak();

#ifdef HAVE_ACIS

  //std::string geomfile = fileroot+".sat";
  //std::string assocfile = fileroot+".m2g";

#endif // HAVE_ACIS

  eMesh.save_as(fileroot+"_R1.g");
}
