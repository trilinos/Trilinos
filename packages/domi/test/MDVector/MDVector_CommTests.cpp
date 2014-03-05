/*
// @HEADER
// ***********************************************************************
//
//            Domi: Multidimensional Datastructures Package
//                 Copyright (2013) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
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
using Domi::splitStringOfIntsWithCommas;
using Domi::TeuchosCommRCP;

string dims          = "8,8";
string axisCommSizes = "-1";
string commPad       = "1,1";
string bndryPad      = "";
string periodic      = "";
bool   verbose       = false;

TEUCHOS_STATIC_SETUP()
{
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.addOutputSetupOptions(true);
  clp.setOption("dims"         , &dims,
                "Comma-separated global dimensions of Field");
  clp.setOption("axisCommSizes", &axisCommSizes,
                "Comma-separated number of processors along each axis");
  clp.setOption("commPad"      , &commPad,
                "Comma-separated list of commPad sizes along each axis");
  clp.setOption("bndryPad"     , &bndryPad,
                "Comma-separated list of bndryPad points on each axis");
  clp.setOption("periodic"     , &periodic,
                "Comma-separated list of axis periodicity flags (use 0,1)");
  clp.setOption("verbose"      , "quiet"       , &verbose,
                "Verbose or quiet output");
}

//
// Templated Unit Test
//

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( MDVector, mdVectorComm, Sca )
{
  // Construct the communicator
  TeuchosCommRCP comm = Teuchos::DefaultComm< int >::getComm();
  int pid = comm->getRank();

  // Convert the command-line arguments into usable arrays
  Array< dim_type > dimVals;
  Array< int > axisCommSizeVals;
  Array< int > commPadVals;
  Array< int > bndryPadVals;
  Array< int > periodicFlags;
  dimVals          = splitStringOfIntsWithCommas(dims         );
  axisCommSizeVals = splitStringOfIntsWithCommas(axisCommSizes);
  commPadVals      = splitStringOfIntsWithCommas(commPad      );
  bndryPadVals     = splitStringOfIntsWithCommas(bndryPad     );
  periodicFlags    = splitStringOfIntsWithCommas(periodic     );

  // Fill out the axisCommSizeVals, if needed
  int numDims = dimVals.size();
  for (int i = axisCommSizeVals.size(); i < numDims; ++i)
    axisCommSizeVals.push_back(-1);

  // Check for 2D or 3D
  if (numDims == 3 && dimVals[2] == 1)
  {
    numDims = 2;
    dimVals.pop_back();
  }
  TEST_ASSERT(numDims == 2 || numDims == 3);
  
  // Print the arrays that will be passed to the MDMap constructor
  if (verbose && pid == 0)
  {
    cout << endl;
    cout << "dims:          " << dimVals          << endl;
    cout << "axisCommSizes: " << axisCommSizeVals << endl;
    cout << "commPad:       " << commPadVals      << endl;
    cout << "bndryPad:      " << bndryPadVals     << endl;
    cout << "periodic:      " << periodicFlags    << endl;
  }

  // Construct the MDComm
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new Domi::MDComm(comm,
                                  axisCommSizeVals(),
                                  periodicFlags() ));

  // Construct the MDMap
  Teuchos::RCP< const Domi::MDMap<> > mdMap =
    Teuchos::rcp(new Domi::MDMap<>(mdComm,
                                   dimVals(),
                                   commPadVals(),
                                   bndryPadVals() ));

  // Construct the MDVector and extract the MDArrayView
  Domi::MDVector< Sca >    mdVector(mdMap);
  Domi::MDArrayView< Sca > mdArray = mdVector.getDataNonConst();

  // Initilize with -1 everywhere
  mdVector.putScalar(-1.0);

  // Assign each owned element the value of its global ID
  if (numDims == 2)
  {
    Domi::Slice iBounds = mdVector.getLocalBounds(0);
    Domi::Slice jBounds = mdVector.getLocalBounds(1);
    for (int j = jBounds.start(); j < jBounds.stop(); ++j)
    {
      for (int i = iBounds.start(); i < iBounds.stop(); ++i)
      {
        int lid = mdMap->getLocalIndex(tuple(i,j));
        mdArray(i,j) = (Sca) mdMap->getGlobalIndex(lid);
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
          int lid = mdMap->getLocalIndex(tuple(i,j,k));
          mdArray(i,j,k) = (Sca) mdMap->getGlobalIndex(lid);
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
  if (numDims == 2)
  {
    Domi::Slice iBounds = mdVector.getLocalBounds(0);
    Domi::Slice jBounds = mdVector.getLocalBounds(1);
    for (int j = jBounds.start(); j < jBounds.stop(); ++j)
    {
      for (int i = iBounds.start(); i < iBounds.stop(); ++i)
      {
        int lid = mdMap->getLocalIndex(tuple(i,j));
        Sca gid = (Sca) mdMap->getGlobalIndex(lid);
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
          int lid = mdMap->getLocalIndex(tuple(i,j,k));
          Sca gid = (Sca) mdMap->getGlobalIndex(lid);
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
  if (numDims == 2)
  {
    for (int j = 0; j < view.dimension(1); ++j)
    {
      for (int i = 0; i < view.dimension(0); ++i)
      {
        int lid = mdMap->getLocalIndex(tuple(i,j));
        Sca gid = (Sca) mdMap->getGlobalIndex(lid);
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
          int lid = mdMap->getLocalIndex(tuple(i,j,k));
          Sca gid = (Sca) mdMap->getGlobalIndex(lid);
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
  if (numDims == 2)
  {
    int ioff = mdMap->getLocalBounds(0).stop();
    for (int j = 0; j < view.dimension(1); ++j)
    {
      for (int i = 0; i < view.dimension(0); ++i)
      {
        int lid = mdMap->getLocalIndex(tuple(i+ioff,j));
        Sca gid = (Sca) mdMap->getGlobalIndex(lid);
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
          int lid = mdMap->getLocalIndex(tuple(i+ioff,j,k));
          Sca gid = (Sca) mdMap->getGlobalIndex(lid);
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
  if (verbose)
    cout << pid << ": checking lower padding along axis 1" << endl;
  view = mdVector.getLowerPadData(1);
  if (numDims == 2)
  {
    for (int j = 0; j < view.dimension(1); ++j)
    {
      for (int i = 0; i < view.dimension(0); ++i)
      {
        int lid = mdMap->getLocalIndex(tuple(i,j));
        Sca gid = (Sca) mdMap->getGlobalIndex(lid);
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
          int lid = mdMap->getLocalIndex(tuple(i,j,k));
          Sca gid = (Sca) mdMap->getGlobalIndex(lid);
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
  if (verbose)
    cout << pid << ": checking upper padding along axis 1" << endl;
  view = mdVector.getUpperPadData(1);
  if (numDims == 2)
  {
    int joff = mdMap->getLocalBounds(1).stop();
    for (int j = 0; j < view.dimension(1); ++j)
    {
      for (int i = 0; i < view.dimension(0); ++i)
      {
        int lid = mdMap->getLocalIndex(tuple(i,j+joff));
        Sca gid = (Sca) mdMap->getGlobalIndex(lid);
        if (mdMap->isBndryPad(tuple(i,j+joff)))
          gid = -1;
        if (verbose)
          cout << pid << ": mdVector(" << i << "," << j+joff << ") = "
               << view(i,j) << " (should be " << gid << ")" << endl;
        TEST_EQUALITY(view(i,j), gid);
      }
    }
  }
  else  // 3D
  {
    int joff = mdMap->getLocalBounds(1).stop();
    for (int k = 0; k < view.dimension(2); ++k)
    {
      for (int j = 0; j < view.dimension(1); ++j)
      {
        for (int i = 0; i < view.dimension(0); ++i)
        {
          int lid = mdMap->getLocalIndex(tuple(i,j+joff,k));
          Sca gid = (Sca) mdMap->getGlobalIndex(lid);
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
          int lid = mdMap->getLocalIndex(tuple(i,j,k));
          Sca gid = (Sca) mdMap->getGlobalIndex(lid);
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
          int lid = mdMap->getLocalIndex(tuple(i,j,k+koff));
          Sca gid = (Sca) mdMap->getGlobalIndex(lid);
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

UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)

}
