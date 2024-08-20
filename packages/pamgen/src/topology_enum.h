// @HEADER
// *****************************************************************************
//                     Pamgen Package
//
// Copyright 2004 NTESS and the Pamgen contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef topology_enumH
#define topology_enumH
  enum Topo_Loc      {MINUS_I=0,PLUS_I,MINUS_J,PLUS_J,MINUS_K,PLUS_K,
		      EDGE0,EDGE1,EDGE2,EDGE3,EDGE4,EDGE5,EDGE6,EDGE7,EDGE8,EDGE9,
		      EDGE10,EDGE11,VERTEX0,VERTEX1,VERTEX2,VERTEX3,VERTEX4,VERTEX5,
		      VERTEX6,VERTEX7,NUM_TOPO_CONNECTIONS,Z_AXIS,PROCESSOR_NODE=100,ALL_NODES=200};
#endif
