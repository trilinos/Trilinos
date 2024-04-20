/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <Ionit_Initializer.h>     // for Initializer
#include <stddef.h>                     // for size_t, nullptr
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD
#include "stk_io/DatabasePurpose.hpp"
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include <stk_unit_test_utils/FaceCreationTestUtils.hpp>

namespace
{

//The Magical Alphabet of Hexes, Shells & Sidesets
//
// A = hex (quad in 2D) in block A
// B = hex (quad in 2D) in block B
// e = shell in block E
// f = shell in block F
// r = bar in block E
// m = beam in block E
// t = truss in block E
// L = sideset associated with the side on the left
// R = "        "           "   "  "     "   " right
// D = sideset containing 2 sides, one associated to left and one to right
// X = sideset associated with all sides on this surface
// J = two hexes in block A connected to the same 8 nodes
// Z = degenerate hex in block Z
// Y = degenerate hex in block Y
// T = tet in block T
// P = hex that is partially coincident with block A on a face
// g = degenerate shell (usually attached to face of tet T)
// .e = the language of our Patron Saint Exodus
//
// RR = pronounced like a pirate
// RRR = roll the R

const SideTestUtil::TestCaseData exposedBoundaryTestCases =
{
  /* filename, max#procs, #side,  sideset */
  {"A.e",         1,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}}},
  {"AL.e",        1,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}}},

  {"Re.e",        1,       2,     {{1, 0}, {1, 1}}},
  {"ReL.e",       1,       2,     {{1, 0}, {1, 1}}},
  {"e.e",         1,       2,     {{1, 0}, {1, 1}}},
  {"eL.e",        1,       2,     {{1, 0}, {1, 1}}},

  {"Ae.e",        2,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}}},
  {"ALe.e",       2,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}}},
  {"ARe.e",       2,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}}},
  {"ADe.e",       2,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}}},


  {"AA.e",        2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AB.e",        2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ADA.e",       2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ADB.e",       2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ADDA.e",      2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ADDB.e",      2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ADReA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ADeDA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ADeLA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ADeRA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALA.e",       2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALB.e",       2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALRA.e",      2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALRB.e",      2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALReA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALeDA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALeLA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALeRA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ARA.e",       2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ARB.e",       2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ARReA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AReDA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AReLA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AReRA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AeA.e",       3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},

  {"ADReB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
  {"ADeDB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
  {"ADeLB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
  {"ADeRB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
  {"ALReB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
  {"ALeDB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
  {"ALeLB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
  {"ALeRB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
  {"ARReB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
  {"AReDB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
  {"AReLB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
  {"AReRB.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},
  {"AeB.e",       3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},

  {"AB_doubleKissing.e",  2,  8,  {{1, 0}, {1, 3}, {1, 4}, {1, 5}, {2, 1}, {2, 2}, {2, 4}, {2, 5}}},
  {"Tg.e",        2,      6,      {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}}},
  {"ZY.e",        2,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},

  {"2D_A.e",      1,       4,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}}},
  {"2D_AL.e",     1,       4,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}}},
  {"2D_AA.e",     2,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},
  {"2D_AB.e",     2,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},
  {"2D_ALA.e",    2,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},
  {"2D_ALB.e",    2,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},
  {"2D_ALRA.e",   2,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},
  {"2D_ALRB.e",   2,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},
  {"2D_ADA.e",    2,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},
  {"2D_ADB.e",    2,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},

  {"2D_AeA.e",    3,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},
  {"2D_AeB.e",    3,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},

  {"2D_AtA.e",    3,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},
  {"2D_AtB.e",    3,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},
  {"2D_ArA.e",    3,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},
  {"2D_ArB.e",    3,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},
  {"2D_AmA.e",    3,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}},
  {"2D_AmB.e",    3,       6,     {{1, 1}, {1, 2}, {1, 3}, {2, 0}, {2, 1}, {2, 3}}}
};

const SideTestUtil::TestCaseData exposedBoundaryCoincidentElementTestCases =
{
  {"Aef.e",       3,       6,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {3, 0}}},
  {"ALeDfRA.e",   4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALeDfRB.e",   4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALeXfRA.e",   4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALefRA.e",    4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ARefLA.e",    4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AeDfA.e",     4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AefA.e",      4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AefB.e",      4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {4, 0}, {4, 1}, {4, 2}, {4, 3}, {4, 5}}},

  {"ef.e",        2,      2,      {{1, 0}, {1, 1}, {2, 0}, {2, 1}}},
  {"eff.e",       3,      2,      {{1, 0}, {1, 1}, {2, 0}, {2, 1}, {3,0}, {3,1}}},
  {"ALJ.e",       3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}, {3, 0}, {3, 1}, {3, 2}, {3, 3}, {3, 5}}},

  {"AP.e",        2,      11,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 4}, {2, 5}}},
};

const SideTestUtil::TestCaseData createExposedBoundaryForOneBlockTestCases =
{
  /* filename, max#procs, #side,  sideset */
  {"AB.e",        2,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ADB.e",       2,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ADDB.e",      2,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ADReB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ADeDB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ADeLB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ADeRB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ALB.e",       2,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ALRB.e",      2,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ALReB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ALeDB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ALeLB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ALeRB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ARB.e",       2,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ARReB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"AReDB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"AReLB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"AReRB.e",     3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"AeB.e",       3,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"Ae.e",        2,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},

  {"ALReA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALeDA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALeLA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALeRA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ARReA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AReDA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AReLA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AReRA.e",     3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AeA.e",       3,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},

  {"AB_doubleKissing.e",  2,  4,  {{1, 0}, {1, 3}, {1, 4}, {1, 5}}},
  {"Tg.e",        2,      4,      {{1, 0}, {1, 1}, {1, 2}, {1, 3}}},
  {"ZY.e",        2,      5,      {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},

  {"2D_AB.e",     2,       3,     {{1, 1}, {1, 2}, {1, 3}}},
  {"2D_ALB.e",    2,       3,     {{1, 1}, {1, 2}, {1, 3}}},
  {"2D_ALRB.e",   2,       3,     {{1, 1}, {1, 2}, {1, 3}}},
  {"2D_ADB.e",    2,       3,     {{1, 1}, {1, 2}, {1, 3}}},
  {"2D_AeB.e",    3,       3,     {{1, 1}, {1, 2}, {1, 3}}},
  {"2D_AtB.e",    3,       3,     {{1, 1}, {1, 2}, {1, 3}}},
  {"2D_ArB.e",    3,       3,     {{1, 1}, {1, 2}, {1, 3}}},
  {"2D_AmB.e",    3,       3,     {{1, 1}, {1, 2}, {1, 3}}}
};

const SideTestUtil::TestCaseData createExposedBoundaryForOneBlockCoincidentElementTestCases =
{
  {"AefB.e",      4,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ALeDfRB.e",   4,       5,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}}},
  {"ALeDfRA.e",   4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALeXfRA.e",   4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ALefRA.e",    4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"ARefLA.e",    4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AeDfA.e",     4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},
  {"AefA.e",      4,      10,     {{1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4}, {2, 0}, {2, 1}, {2, 2}, {2, 3}, {2, 5}}},

  {"ef.e",        2,      2,      {{1, 0}, {1, 1}, {2, 0}, {2, 1}}},
  {"eff.e",       3,      2,      {{1, 0}, {1, 1}, {2, 0}, {2, 1}, {3,0}, {3,1}}},
};

const SideTestUtil::TestCaseData createExposedBoundaryForDegenerateElementTestCases =
{
  {"quadDegenTriSandwich.g",  3,  7, {}},
  {"hexDegenWedgeSandwich.g", 3,  13, {}}
};

class SkinnedMesh : public SideTestUtil::SideCreationTester
{
public:
  SkinnedMesh() : SideTestUtil::SideCreationTester(MPI_COMM_WORLD) {}
protected:
  virtual void test_side_creation(stk::mesh::BulkData& bulkData,
                                  const SideTestUtil::TestCase& testCase)
  {
    stk::mesh::Part& skinnedPart = SideTestUtil::run_skin_mesh(bulkData, get_things_to_skin(bulkData));
    expect_exposed_sides_connected_as_specified_in_test_case(bulkData, testCase, get_things_to_skin(bulkData), skinnedPart);
  }

  virtual stk::mesh::Selector get_things_to_skin(const stk::mesh::BulkData& bulkData)
  {
    return bulkData.mesh_meta_data().universal_part();
  }
};

class SkinSingleBlock : public SkinnedMesh
{
protected:
  virtual stk::mesh::Selector get_things_to_skin(const stk::mesh::BulkData& bulkData)
  {
    return *bulkData.mesh_meta_data().get_part("block_1");
  }
};

TEST(ExposedBlockBoundaryTest, run_all_test_cases_aura)
{
  SkinnedMesh().run_all_test_cases(exposedBoundaryTestCases, stk::mesh::BulkData::AUTO_AURA);
}

TEST(ExposedBlockBoundaryTest, run_all_test_cases_no_aura)
{
  SkinnedMesh().run_all_test_cases(exposedBoundaryTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST(ExposedBlockBoundaryTest, run_some_degenerate_cases_aura)
{
  SkinnedMesh().run_all_test_cases(createExposedBoundaryForDegenerateElementTestCases, stk::mesh::BulkData::AUTO_AURA);
}

TEST(ExposedBlockBoundaryTest, run_some_degenerate_cases_no_aura)
{
  SkinnedMesh().run_all_test_cases(createExposedBoundaryForDegenerateElementTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST(CreateExposedBoundaryForSingleBlockTest, run_all_test_cases_aura)
{
  SkinSingleBlock().run_all_test_cases(createExposedBoundaryForOneBlockTestCases, stk::mesh::BulkData::AUTO_AURA);
}

TEST(CreateExposedBoundaryForSingleBlockTest, run_all_test_cases_no_aura)
{
  SkinSingleBlock().run_all_test_cases(createExposedBoundaryForOneBlockTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}

// now using graph to create faces during mesh read, split coincident test cases fail
TEST(ExposedBlockBoundaryTest, DISABLED_run_coincident_element_test_cases_no_aura)
{
  SkinnedMesh().run_all_test_cases(exposedBoundaryCoincidentElementTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST(CreateExposedBoundaryForSingleBlockTest, DISABLED_run_coincident_element_test_cases_no_aura)
{
  SkinSingleBlock().run_all_test_cases(createExposedBoundaryForOneBlockCoincidentElementTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST(ExposedBlockBoundaryTest, np1_run_coincident_element_test_cases_no_aura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
    SkinnedMesh().run_all_test_cases(exposedBoundaryCoincidentElementTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}
TEST(CreateExposedBoundaryForSingleBlockTest, np1_run_coincident_element_test_cases_no_aura)
{
  if(stk::parallel_machine_size(MPI_COMM_WORLD) == 1)
    SkinSingleBlock().run_all_test_cases(createExposedBoundaryForOneBlockCoincidentElementTestCases, stk::mesh::BulkData::NO_AUTO_AURA);
}
} //namespace
