// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Ioss_Bar2.h>
#include <Ioss_Bar3.h>
#include <Ioss_Edge2.h>
#include <Ioss_Edge2D2.h>
#include <Ioss_Edge2D3.h>
#include <Ioss_Edge3.h>
#include <Ioss_Hex20.h>
#include <Ioss_Hex27.h>
#include <Ioss_Hex8.h>
#include <Ioss_Initializer.h>
#include <Ioss_Node.h>
#include <Ioss_Pyramid5.h>
#include <Ioss_Pyramid13.h>
#include <Ioss_Pyramid14.h>
#include <Ioss_Quad4.h>
#include <Ioss_Quad8.h>
#include <Ioss_Quad9.h>
#include <Ioss_Shell4.h>
#include <Ioss_Shell8.h>
#include <Ioss_Shell9.h>
#include <Ioss_ShellLine2D2.h>
#include <Ioss_ShellLine2D3.h>
#include <Ioss_Sphere.h>
#include <Ioss_Super.h>
#include <Ioss_Tet4.h>
  // Disable for now.  Might be needed for adagio/salinas coupling...
//#include <Ioss_Tet7.h>
#include <Ioss_Tet8.h>
#include <Ioss_Tet10.h>
#include <Ioss_Tet11.h>
#include <Ioss_Tri3.h>
#include <Ioss_Tri4.h>
// Disable for now.  Might be needed for adagio/salinas coupling...
//#include <Ioss_Tri4a.h>
#include <Ioss_Tri6.h>
#include <Ioss_TriShell3.h>
#include <Ioss_TriShell4.h>
#include <Ioss_TriShell6.h>
#include <Ioss_Unknown.h>
#include <Ioss_Wedge6.h>
#include <Ioss_Wedge15.h>
#include <Ioss_Wedge18.h>

Ioss::Initializer::Initializer()
{
  // List all storage types here with a call to their factory method.
  // This is Used to get the linker to pull in all needed libraries.
  Ioss::Sphere::factory();
  
  Ioss::Edge2D2::factory();
  Ioss::Edge2D3::factory();
  Ioss::Edge2::factory();
  Ioss::Edge3::factory();

  Ioss::Bar2::factory();
  Ioss::Bar3::factory();
  Ioss::ShellLine2D2::factory();
  Ioss::ShellLine2D3::factory();

  Ioss::Hex8::factory();
  Ioss::Hex20::factory();
  Ioss::Hex27::factory();

  Ioss::Node::factory();

  Ioss::Pyramid5::factory();
  Ioss::Pyramid13::factory();
  Ioss::Pyramid14::factory();
  
  Ioss::Quad4::factory();
  Ioss::Quad8::factory();
  Ioss::Quad9::factory();

  Ioss::Shell4::factory();
  Ioss::Shell8::factory();
  Ioss::Shell9::factory();

  Ioss::Tet4::factory();

  // Disable for now.  Might be needed for adagio/salinas coupling...
  //Ioss::Tet7::factory();
  Ioss::Tet8::factory();
  Ioss::Tet10::factory();
  Ioss::Tet11::factory();

  Ioss::Tri3::factory();
  Ioss::Tri4::factory();
  // Disable for now.  Might be needed for adagio/salinas coupling...
  //Ioss::Tri4a::factory();
  Ioss::Tri6::factory();

  Ioss::TriShell3::factory();
  Ioss::TriShell4::factory();
  Ioss::TriShell6::factory();

  Ioss::Unknown::factory();

  Ioss::Wedge6::factory();
  Ioss::Wedge15::factory();
  Ioss::Wedge18::factory();

  Ioss::Super::factory();
}
