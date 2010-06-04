/*--------------------------------------------------------------------*/
/*    Copyright 2000,2002,2009 Sandia Corporation.                         */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_Initializer.h>

#include <Ioss_Bar2.h>
#include <Ioss_Bar3.h>
#include <Ioss_ShellLine2D2.h>
#include <Ioss_ShellLine2D3.h>
#include <Ioss_Edge2.h>
#include <Ioss_Edge3.h>
#include <Ioss_Edge2D2.h>
#include <Ioss_Edge2D3.h>
#include <Ioss_ElementVariableType.h>
#include <Ioss_Hex8.h>
#include <Ioss_Hex20.h>
#include <Ioss_Hex27.h>
#include <Ioss_Node.h>
#include <Ioss_Pyramid5.h>
#include <Ioss_Pyramid13.h>
#include <Ioss_Quad4.h>
#include <Ioss_Quad8.h>
#include <Ioss_Quad9.h>
#include <Ioss_Shell4.h>
#include <Ioss_Shell8.h>
#include <Ioss_Shell9.h>
#include <Ioss_Sphere.h>
#include <Ioss_TriShell3.h>
#include <Ioss_TriShell6.h>
#include <Ioss_Tri3.h>
#include <Ioss_Tri4.h>
#include <Ioss_Tri6.h>
#include <Ioss_Unknown.h>
#include <Ioss_Wedge6.h>
#include <Ioss_Wedge15.h>
#include <Ioss_Tet4.h>
#include <Ioss_Tet8.h>
#include <Ioss_Tet10.h>
#include <Ioss_Super.h>

Ioss::Initializer::Initializer()
{
  // List all storage types here with a call to their factory method.
  // This is Used to get the linker to pull in all needed libraries.
  Ioss::Sphere::factory();
  
  Ioss::Edge2::factory();
  Ioss::Edge3::factory();
  Ioss::Edge2D2::factory();
  Ioss::Edge2D3::factory();

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
  
  Ioss::Quad4::factory();
  Ioss::Quad8::factory();
  Ioss::Quad9::factory();

  Ioss::Shell4::factory();
  Ioss::Shell8::factory();
  Ioss::Shell9::factory();

  Ioss::Tet4::factory();
  Ioss::Tet8::factory();
  Ioss::Tet10::factory();

  Ioss::Tri3::factory();
  Ioss::Tri4::factory();
  Ioss::Tri6::factory();

  Ioss::TriShell3::factory();
  Ioss::TriShell6::factory();

  Ioss::Unknown::factory();

  Ioss::Wedge6::factory();
  Ioss::Wedge15::factory();

  Ioss::Super::factory();
}
