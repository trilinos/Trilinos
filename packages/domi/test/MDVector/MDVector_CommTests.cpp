/*
// @HEADER
// ***********************************************************************
//
//     Domi: Multi-dimensional Distributed Linear Algebra Services
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER
*/

// Teuchos includes
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_Utils.hpp"
#include "Domi_MDVector.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace
{

using std::string;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::Tuple;
using Teuchos::tuple;
typedef Domi::dim_type dim_type;
typedef Domi::size_type size_type;
using Domi::splitStringOfIntsWithCommas;

string dims      = "20,16";
string commDims  = "-1";
int    commPad   = 1;
string commPads  = "";
int    bndryPad  = 0;
string bndryPads = "";
string periodic  = "";
bool   verbose   = false;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption("dims"     , &dims,
                "Comma-separated global dimensions of Field");
  clp.setOption("commDims" , &commDims,
                "Comma-separated number of processors along each axis");
  clp.setOption("commPad"  , &commPad,
                "CommPad size along every axis");
  clp.setOption("commPads" , &commPads,
                "Comma-separated list of commPad sizes along each axis");
  clp.setOption("bndryPad" , &bndryPad,
                "CBndryPad size along every axis");
  clp.setOption("bndryPads", &bndryPads,
                "Comma-separated list of bndryPad sizes on each axis");
  clp.setOption("periodic" , &periodic,
                "Comma-separated list of axis periodicity flags (use 0,1)");
  clp.setOption("verbose"  , "quiet"       , &verbose,
                "Verbose or quiet output");
}

//
// Templated Unit Test
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, mdVectorComm, Sca )
{
  // Construct the communicator
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();
  int pid = comm->getRank();

  // Convert the command-line arguments into usable arrays
  Array< dim_type > dimVals;
  Array< int > commDimVals;
  Array< int > commPadVals;
  Array< int > bndryPadVals;
  Array< int > periodicFlags;
  dimVals       = splitStringOfIntsWithCommas(dims     );
  commDimVals   = splitStringOfIntsWithCommas(commDims );
  commPadVals   = splitStringOfIntsWithCommas(commPads );
  bndryPadVals  = splitStringOfIntsWithCommas(bndryPads);
  periodicFlags = splitStringOfIntsWithCommas(periodic );

  // Check for 2D or 3D
  int numDims = dimVals.size();
  if (numDims == 3 && dimVals[2] == 1)
  {
    numDims = 2;
    dimVals.pop_back();
  }
  TEST_ASSERT(numDims >= 1 && numDims <= 3);

  // Print the arrays that will be passed to the MDMap constructor
  if (verbose && pid == 0)
    cout << endl
         << "dims:      " << dimVals       << endl
         << "commDims:  " << commDimVals   << endl
         << "commPad:   " << commPad       << endl
         << "commPads:  " << commPadVals   << endl
         << "bndryPad:  " << bndryPad      << endl
         << "bndryPads: " << bndryPadVals  << endl
         << "periodic:  " << periodicFlags << endl;

  // Construct the MDVector ParameterList
  Teuchos::ParameterList plist;
  if (! commDimVals.empty())
    plist.set("comm dimensions", commDimVals);
  if (! periodicFlags.empty())
    plist.set("periodic", periodicFlags);
  plist.set("dimensions", dimVals);
  if (bndryPad)
    plist.set("boundary pad size", bndryPad);
  if (! bndryPadVals.empty())
    plist.set("boundary pad sizes", bndryPadVals);
  if (commPad)
    plist.set("communication pad size", commPad);
  if (! commPadVals.empty())
    plist.set("communication pad sizes", commPadVals);
  if (verbose && pid == 0)
    cout << endl << "MDVector constructor ParameterList =" << endl << plist
         << endl;

  // Construct the MDVector and extract the MDArrayView and MDMap
  Domi::MDVector< Sca >    mdVector(comm, plist);
  Domi::MDArrayView< Sca > mdArray = mdVector.getDataNonConst();
  Teuchos::RCP< const Domi::MDMap<> > mdMap = mdVector.getMDMap();

  // Reconstruct the periodicity flags so that we can legally check
  // periodicity along each axis
  periodicFlags.resize(numDims);
  for (int axis = 0; axis < numDims; ++axis)
    periodicFlags[axis] = mdMap->isPeriodic(axis) ? 1 : 0;

  // Initilize with -1 everywhere
  mdVector.putScalar(-1.0);

  // Assign each owned element the value of its global ID
  if (numDims == 1)
  {
    Domi::Slice iBounds = mdVector.getLocalBounds(0);
    for (int i = iBounds.start(); i < iBounds.stop(); ++i)
    {
      int lid = mdMap->getLocalID(tuple(i));
      mdArray(i) = (Sca) mdMap->getGlobalID(lid);
    }
  }
  else if (numDims == 2)
  {
    Domi::Slice iBounds = mdVector.getLocalBounds(0);
    Domi::Slice jBounds = mdVector.getLocalBounds(1);
    for (int j = jBounds.start(); j < jBounds.stop(); ++j)
    {
      for (int i = iBounds.start(); i < iBounds.stop(); ++i)
      {
        int lid = mdMap->getLocalID(tuple(i,j));
        mdArray(i,j) = (Sca) mdMap->getGlobalID(lid);
      }
    }
  }
  else  // 3D
  {
    Domi::Slice iBounds = mdVector.getLocalBounds(0);
    Domi::Slice jBounds = mdVector.getLocalBounds(1);
    Domi::Slice kBounds = mdVector.getLocalBounds(2);
    for (int k = kBounds.start(); k < kBounds.stop(); ++k)
    {
      for (int j = jBounds.start(); j < jBounds.stop(); ++j)
      {
        for (int i = iBounds.start(); i < iBounds.stop(); ++i)
        {
          int lid = mdMap->getLocalID(tuple(i,j,k));
          mdArray(i,j,k) = (Sca) mdMap->getGlobalID(lid);
        }
      }
    }
  }

  // Update the communication padding values.  After returning
  // from this method, communication padding points should now
  // have values that correspond to their global IDs.  Boundary
  // padding points should still be equal to -1.
  mdVector.updateCommPad();

  // Check the owned, non-padding values.  Interior data should be
  // unchanged
  if (verbose)
    cout << pid << ": checking interior data" << endl;
  if (numDims == 1)
  {
    Domi::Slice iBounds = mdVector.getLocalBounds(0);
    for (int i = iBounds.start(); i < iBounds.stop(); ++i)
    {
      int lid = mdMap->getLocalID(tuple(i));
      Sca gid = (Sca) mdMap->getGlobalID(lid);
      if (verbose)
        cout << pid << ": mdVector(" << i << ") = "
             << mdArray(i) << " (should be " << gid << ")"
             << endl;
      TEST_EQUALITY(mdArray(i), gid);
    }
  }
  else if (numDims == 2)
  {
    Domi::Slice iBounds = mdVector.getLocalBounds(0);
    Domi::Slice jBounds = mdVector.getLocalBounds(1);
    for (int j = jBounds.start(); j < jBounds.stop(); ++j)
    {
      for (int i = iBounds.start(); i < iBounds.stop(); ++i)
      {
        int lid = mdMap->getLocalID(tuple(i,j));
        Sca gid = (Sca) mdMap->getGlobalID(lid);
        if (verbose)
          cout << pid << ": mdVector(" << i << "," << j << ") = "
               << mdArray(i,j) << " (should be " << gid << ")"
               << endl;
        TEST_EQUALITY(mdArray(i,j), gid);
      }
    }
  }
  else  // 3D
  {
    Domi::Slice iBounds = mdVector.getLocalBounds(0);
    Domi::Slice jBounds = mdVector.getLocalBounds(1);
    Domi::Slice kBounds = mdVector.getLocalBounds(2);
    for (int k = kBounds.start(); k < kBounds.stop(); ++k)
    {
      for (int j = jBounds.start(); j < jBounds.stop(); ++j)
      {
        for (int i = iBounds.start(); i < iBounds.stop(); ++i)
        {
          int lid = mdMap->getLocalID(tuple(i,j,k));
          Sca gid = (Sca) mdMap->getGlobalID(lid);
          if (verbose)
            cout << pid << ": mdVector(" << i << "," << j << "," << k
                 << ") = " << mdArray(i,j,k) << " (should be " << gid << ")"
                 << endl;
          TEST_EQUALITY(mdArray(i,j,k), gid);
        }
      }
    }
  }

  // Check the lower padding along axis 0
  if (verbose)
    cout << pid << ": checking lower padding along axis 0" << endl;
  Domi::MDArrayView< const Sca > view = mdVector.getLowerPadData(0);
  if (numDims == 1)
  {
    for (int i = 0; i < view.dimension(0); ++i)
    {
      int lid = mdMap->getLocalID(tuple(i));
      Sca gid = (Sca) mdMap->getGlobalID(lid);
      if (mdMap->isPeriodic(0))
      {
        Teuchos::Array< dim_type > globalIndex =
          mdMap->getGlobalIndex((size_type)gid);
        if (globalIndex[0] < mdMap->getCommPadSize(0))
          gid += mdMap->getGlobalDim(0);
      }
      if (mdMap->isBndryPad(tuple(i)))
        gid = -1;
      if (verbose)
        cout << pid << ": mdVector(" << i << ") = "
             << mdArray(i) << " (should be " << gid << ")"
             << endl;
      TEST_EQUALITY(view(i), gid);
    }
  }
  else if (numDims == 2)
  {
    for (int j = 0; j < view.dimension(1); ++j)
    {
      for (int i = 0; i < view.dimension(0); ++i)
      {
        int lid = mdMap->getLocalID(tuple(i,j));
        Sca gid = (Sca) mdMap->getGlobalID(lid);
        if (mdMap->isPeriodic(0))
        {
          Teuchos::Array< dim_type > globalIndex =
            mdMap->getGlobalIndex((size_type)gid);
          if (globalIndex[0] < mdMap->getCommPadSize(0))
            gid += mdMap->getGlobalDim(0);
        }
        if (mdMap->isPeriodic(1))
        {
          Teuchos::Array< dim_type > globalIndex =
            mdMap->getGlobalIndex((size_type)gid);
          if (globalIndex[1] < mdMap->getGlobalBounds(1).start())
            gid += mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1);
          if (globalIndex[1] >= mdMap->getGlobalBounds(1).stop())
            gid -= mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1);
        }
        if (mdMap->isBndryPad(tuple(i,j)))
          gid = -1;
        if (verbose)
          cout << pid << ": mdVector(" << i << "," << j << ") = "
               << mdArray(i,j) << " (should be " << gid << ")"
               << endl;
        TEST_EQUALITY(view(i,j), gid);
      }
    }
  }
  else  // 3D
  {
    for (int k = 0; k < view.dimension(2); ++k)
    {
      for (int j = 0; j < view.dimension(1); ++j)
      {
        for (int i = 0; i < view.dimension(0); ++i)
        {
          int lid = mdMap->getLocalID(tuple(i,j,k));
          Sca gid = (Sca) mdMap->getGlobalID(lid);
          if (mdMap->isPeriodic(0))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[0] < mdMap->getCommPadSize(0))
              gid += mdMap->getGlobalDim(0);
          }
          if (mdMap->isPeriodic(1))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[1] < mdMap->getCommPadSize(1))
              gid += mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1);
            if (globalIndex[1] >= mdMap->getGlobalBounds(1).stop())
              gid -= mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1);
          }
          if (mdMap->isPeriodic(2))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[2] < mdMap->getCommPadSize(2))
              gid += mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1) *
                mdMap->getGlobalDim(2);
            if (globalIndex[2] >= mdMap->getGlobalBounds(2).stop())
              gid -= mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1) *
                mdMap->getGlobalDim(2);
          }
          if (mdMap->isBndryPad(tuple(i,j,k)))
            gid = -1;
          if (verbose)
            cout << pid << ": mdVector(" << i << "," << j << "," << k
                 << ") = " << view(i,j,k) << " (should be " << gid << ")"
                 << endl;
          TEST_EQUALITY(view(i,j,k), gid);
        }
      }
    }
  }

  // Check the upper padding along axis 0
  if (verbose)
    cout << pid << ": checking upper padding along axis 0" << endl;
  view = mdVector.getUpperPadData(0);
  if (numDims == 1)
  {
    int ioff = mdMap->getLocalBounds(0).stop();
    for (int i = 0; i < view.dimension(0); ++i)
    {
      int lid = mdMap->getLocalID(tuple(i+ioff));
      Sca gid = (Sca) mdMap->getGlobalID(lid);
      if (mdMap->isPeriodic(0))
      {
        Teuchos::Array< dim_type > globalIndex =
          mdMap->getGlobalIndex((size_type)gid);
        if (globalIndex[0] >= mdMap->getGlobalBounds(0).stop())
          gid -= mdMap->getGlobalDim(0);
      }
      if (mdMap->isBndryPad(tuple(i+ioff)))
        gid = -1;
      if (verbose)
        cout << pid << ": mdVector(" << i+ioff << ") = "
             << view(i) << " (should be " << gid << ")" << endl;
      TEST_EQUALITY(view(i), gid);
    }
  }
  else if (numDims == 2)
  {
    int ioff = mdMap->getLocalBounds(0).stop();
    for (int j = 0; j < view.dimension(1); ++j)
    {
      for (int i = 0; i < view.dimension(0); ++i)
      {
        int lid = mdMap->getLocalID(tuple(i+ioff,j));
        Sca gid = (Sca) mdMap->getGlobalID(lid);
        if (mdMap->isPeriodic(0))
        {
          Teuchos::Array< dim_type > globalIndex =
            mdMap->getGlobalIndex((size_type)gid);
          if (globalIndex[0] >= mdMap->getGlobalBounds(0).stop())
            gid -= mdMap->getGlobalDim(0);
        }
        if (mdMap->isPeriodic(1))
        {
          Teuchos::Array< dim_type > globalIndex =
            mdMap->getGlobalIndex((size_type)gid);
          if (globalIndex[1] < mdMap->getGlobalBounds(1).start())
            gid += mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1);
          if (globalIndex[1] >= mdMap->getGlobalBounds(1).stop())
            gid -= mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1);
        }
        if (mdMap->isBndryPad(tuple(i+ioff,j)))
          gid = -1;
        if (verbose)
          cout << pid << ": mdVector(" << i+ioff << "," << j << ") = "
               << view(i,j) << " (should be " << gid << ")" << endl;
        TEST_EQUALITY(view(i,j), gid);
      }
    }
  }
  else  // 3D
  {
    int ioff = mdMap->getLocalBounds(0).stop();
    for (int k = 0; k < view.dimension(2); ++k)
    {
      for (int j = 0; j < view.dimension(1); ++j)
      {
        for (int i = 0; i < view.dimension(0); ++i)
        {
          int lid = mdMap->getLocalID(tuple(i+ioff,j,k));
          Sca gid = (Sca) mdMap->getGlobalID(lid);
          if (mdMap->isPeriodic(0))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[0] >= mdMap->getGlobalBounds(0).stop())
              gid -= mdMap->getGlobalDim(0);
          }
          if (mdMap->isPeriodic(1))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[1] < mdMap->getGlobalBounds(1).start())
              gid += mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1);
            if (globalIndex[1] >= mdMap->getGlobalBounds(1).stop())
              gid -= mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1);
          }
          if (mdMap->isPeriodic(2))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[2] < mdMap->getGlobalBounds(2).start())
              gid += mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1) *
                mdMap->getGlobalDim(2);
            if (globalIndex[2] >= mdMap->getGlobalBounds(2).stop())
              gid -= mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1) *
                mdMap->getGlobalDim(2);
          }
          if (mdMap->isBndryPad(tuple(i+ioff,j,k)))
            gid = -1;
          if (verbose)
            cout << pid << ": mdVector(" << i+ioff << "," << j << "," << k
                 << ") = " << view(i,j,k) << " (should be " << gid << ")"
                 << endl;
          TEST_EQUALITY(view(i,j,k), gid);
        }
      }
    }
  }

  // Check the lower padding along axis 1
  if (numDims > 1)
  {
    if (verbose)
      cout << pid << ": checking lower padding along axis 1" << endl;
    view = mdVector.getLowerPadData(1);
  }
  if (numDims == 2)
  {
    for (int j = 0; j < view.dimension(1); ++j)
    {
      for (int i = 0; i < view.dimension(0); ++i)
      {
        int lid = mdMap->getLocalID(tuple(i,j));
        Sca gid = (Sca) mdMap->getGlobalID(lid);
        if (mdMap->isPeriodic(0))
        {
          Teuchos::Array< dim_type > globalIndex =
            mdMap->getGlobalIndex((size_type)gid);
          if (globalIndex[0] < mdMap->getCommPadSize(0))
            gid += mdMap->getGlobalDim(0);
          if (globalIndex[0] >= mdMap->getGlobalBounds(0).stop())
            gid -= mdMap->getGlobalDim(0);
        }
        if (mdMap->isPeriodic(1))
        {
          Teuchos::Array< dim_type > globalIndex =
            mdMap->getGlobalIndex((size_type)gid);
          if (globalIndex[1] < mdMap->getCommPadSize(1))
            gid += mdMap->getGlobalDim(0)*mdMap->getGlobalDim(1);
        }
        if (mdMap->isBndryPad(tuple(i,j)))
          gid = -1;
        if (verbose)
          cout << pid << ": mdVector(" << i << "," << j << ") = "
               << mdArray(i,j) << " (should be " << gid << ")"
               << endl;
        TEST_EQUALITY(view(i,j), gid);
      }
    }
  }
  else if (numDims == 3)
  {
    for (int k = 0; k < view.dimension(2); ++k)
    {
      for (int j = 0; j < view.dimension(1); ++j)
      {
        for (int i = 0; i < view.dimension(0); ++i)
        {
          int lid = mdMap->getLocalID(tuple(i,j,k));
          Sca gid = (Sca) mdMap->getGlobalID(lid);
          if (mdMap->isPeriodic(0))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[0] < mdMap->getCommPadSize(0))
              gid += mdMap->getGlobalDim(0);
            if (globalIndex[0] >= mdMap->getGlobalBounds(0).stop())
              gid -= mdMap->getGlobalDim(0);
          }
          if (mdMap->isPeriodic(1))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[1] < mdMap->getCommPadSize(1))
              gid += mdMap->getGlobalDim(0)*mdMap->getGlobalDim(1);
          }
          if (mdMap->isPeriodic(2))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[2] < mdMap->getCommPadSize(2))
              gid += mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1) *
                mdMap->getGlobalDim(2);
            if (globalIndex[2] >= mdMap->getGlobalBounds(2).stop())
              gid -= mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1) *
                mdMap->getGlobalDim(2);
          }
          if (mdMap->isBndryPad(tuple(i,j,k)))
            gid = -1;
          if (verbose)
            cout << pid << ": mdVector(" << i << "," << j << "," << k
                 << ") = " << view(i,j,k) << " (should be " << gid << ")"
                 << endl;
          TEST_EQUALITY(view(i,j,k), gid);
        }
      }
    }
  }

  // Check the upper padding along axis 1
  if (numDims > 1)
  {
    if (verbose)
      cout << pid << ": checking upper padding along axis 1" << endl;
    view = mdVector.getUpperPadData(1);
  }
  if (numDims == 2)
  {
    int joff = mdMap->getLocalBounds(1).stop();
    for (int j = 0; j < view.dimension(1); ++j)
    {
      for (int i = 0; i < view.dimension(0); ++i)
      {
        int lid = mdMap->getLocalID(tuple(i,j+joff));
        Sca gid = (Sca) mdMap->getGlobalID(lid);
        if (mdMap->isPeriodic(0))
        {
          Teuchos::Array< dim_type > globalIndex =
            mdMap->getGlobalIndex((size_type)gid);
          if (globalIndex[0] < mdMap->getGlobalBounds(0).start())
            gid += mdMap->getGlobalDim(0);
          if (globalIndex[0] >= mdMap->getGlobalBounds(0).stop())
            gid -= mdMap->getGlobalDim(0);
        }
        if (mdMap->isPeriodic(1))
        {
          Teuchos::Array< dim_type > globalIndex =
            mdMap->getGlobalIndex((size_type)gid);
          if (globalIndex[1] >= mdMap->getGlobalBounds(1).stop())
            gid -= mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1);
        }
        if (mdMap->isBndryPad(tuple(i,j+joff)))
          gid = -1;
        if (verbose)
          cout << pid << ": mdVector(" << i << "," << j+joff << ") = "
               << view(i,j) << " (should be " << gid << ")" << endl;
        TEST_EQUALITY(view(i,j), gid);
      }
    }
  }
  else if (numDims == 3)
  {
    int joff = mdMap->getLocalBounds(1).stop();
    for (int k = 0; k < view.dimension(2); ++k)
    {
      for (int j = 0; j < view.dimension(1); ++j)
      {
        for (int i = 0; i < view.dimension(0); ++i)
        {
          int lid = mdMap->getLocalID(tuple(i,j+joff,k));
          Sca gid = (Sca) mdMap->getGlobalID(lid);
          if (mdMap->isPeriodic(0))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[0] < mdMap->getGlobalBounds(0).start())
              gid += mdMap->getGlobalDim(0);
            if (globalIndex[0] >= mdMap->getGlobalBounds(0).stop())
              gid -= mdMap->getGlobalDim(0);
          }
          if (mdMap->isPeriodic(1))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[1] >= mdMap->getGlobalBounds(1).stop())
              gid -= mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1);
          }
          if (mdMap->isPeriodic(2))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[2] < mdMap->getGlobalBounds(2).start())
              gid += mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1) *
                mdMap->getGlobalDim(2);
            if (globalIndex[2] >= mdMap->getGlobalBounds(2).stop())
              gid -= mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1) *
                mdMap->getGlobalDim(2);
          }
          if (mdMap->isBndryPad(tuple(i,j+joff,k)))
            gid = -1;
          if (verbose)
            cout << pid << ": mdVector(" << i << "," << j+joff << "," << k
                 << ") = " << view(i,j,k) << " (should be " << gid << ")"
                 << endl;
          TEST_EQUALITY(view(i,j,k), gid);
        }
      }
    }
  }

  // Check the lower padding along axis 2
  if (numDims == 3)
  {
    if (verbose)
      cout << pid << ": checking lower padding along axis 2" << endl;
    view = mdVector.getLowerPadData(2);
    for (int k = 0; k < view.dimension(2); ++k)
    {
      for (int j = 0; j < view.dimension(1); ++j)
      {
        for (int i = 0; i < view.dimension(0); ++i)
        {
          int lid = mdMap->getLocalID(tuple(i,j,k));
          Sca gid = (Sca) mdMap->getGlobalID(lid);
          if (mdMap->isPeriodic(0))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[0] < mdMap->getCommPadSize(0))
              gid += mdMap->getGlobalDim(0);
            if (globalIndex[0] >= mdMap->getGlobalBounds(0).stop())
              gid -= mdMap->getGlobalDim(0);
          }
          if (mdMap->isPeriodic(1))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[1] < mdMap->getCommPadSize(1))
              gid += mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1);
            if (globalIndex[1] >= mdMap->getGlobalBounds(1).stop())
              gid -= mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1);
          }
          if (mdMap->isPeriodic(2))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[2] < mdMap->getCommPadSize(2))
              gid += mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1) *
                mdMap->getGlobalDim(2);
          }
          if (mdMap->isBndryPad(tuple(i,j,k)))
            gid = -1;
          if (verbose)
            cout << pid << ": mdVector(" << i << "," << j << "," << k
                 << ") = " << view(i,j,k) << " (should be " << gid << ")"
                 << endl;
          TEST_EQUALITY(view(i,j,k), gid);
        }
      }
    }
  }

  // Check the upper padding along axis 2
  if (numDims == 3)
  {
    if (verbose)
      cout << pid << ": checking upper padding along axis 2" << endl;
    view = mdVector.getUpperPadData(2);
    int koff = mdMap->getLocalBounds(2).stop();
    for (int k = 0; k < view.dimension(2); ++k)
    {
      for (int j = 0; j < view.dimension(1); ++j)
      {
        for (int i = 0; i < view.dimension(0); ++i)
        {
          int lid = mdMap->getLocalID(tuple(i,j,k+koff));
          Sca gid = (Sca) mdMap->getGlobalID(lid);
          if (mdMap->isPeriodic(0))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[0] < mdMap->getCommPadSize(0))
              gid += mdMap->getGlobalDim(0);
            if (globalIndex[0] >= mdMap->getGlobalBounds(0).stop())
              gid -= mdMap->getGlobalDim(0);
          }
          if (mdMap->isPeriodic(1))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[1] < mdMap->getCommPadSize(1))
              gid += mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1);
            if (globalIndex[1] >= mdMap->getGlobalBounds(1).stop())
              gid -= mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1);
          }
          if (mdMap->isPeriodic(2))
          {
            Teuchos::Array< dim_type > globalIndex =
              mdMap->getGlobalIndex((size_type)gid);
            if (globalIndex[2] >= mdMap->getGlobalBounds(2).stop())
              gid -= mdMap->getGlobalDim(0) * mdMap->getGlobalDim(1) *
                mdMap->getGlobalDim(2);
          }
          if (mdMap->isBndryPad(tuple(i,j,k+koff)))
            gid = -1;
          if (verbose)
            cout << pid << ": mdVector(" << i << "," << j << "," << k+koff
                 << ") = " << view(i,j,k) << " (should be " << gid << ")"
                 << endl;
          TEST_EQUALITY(view(i,j,k), gid);
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////

#define UNIT_TEST_GROUP( Sca ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT( MDVector, mdVectorComm, Sca )

UNIT_TEST_GROUP(int)

#if 0
UNIT_TEST_GROUP(long)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)
#endif

}
