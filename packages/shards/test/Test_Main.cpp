// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

extern void test_shards_array();
extern void test_shards_cell_topology();

int main( int , char ** )
{
  test_shards_array();
  test_shards_cell_topology();
}

