
// System includes
#include <assert.h>
#include <string>
#include <sstream>
using std::string;

// Teuchos includes
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_Tuple.hpp"
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::Tuple;
using Teuchos::tuple;

// Domi includes
#include "Domi_ConfigDefs.hpp"
#include "Domi_Utils.hpp"
#include "Domi_MDVector.hpp"
using Domi::splitStringOfIntsWithCommas;

typedef Domi::MDArrayView< const double >::const_iterator citerator;

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  bool success = true;
  try
  {

    Teuchos::CommandLineProcessor clp;
    clp.setDocString
      ("MDVector test code allows the user to specify the parameters for\n"
       "constructing an MDVector for testing the updating of the communica-\n"
       "tion padding.  The 'dims', 'axisCommSizes', 'commPad', 'bndryPad' and\n"
       "'periodic' options should be one or more comma-separated integers.\n"
       "The number of integers specified in the 'dims' option will determine\n"
       "the dimension of the MDMap (1D, 2D, 3D, etc.).  The 'axisCommSizes',\n"
       "'commPad', 'bndryPad' and 'periodic' options can have fewer integers\n"
       "specified than the 'dims' option, and those unspecified values will\n"
       "then receive default values.  The 'periodic' option is just a list of\n"
       "flags, 0 or 1, specifying whether the axis is periodic.\n");

    string dims          = "8,8";
    string axisCommSizes = "-1";
    string commPad       = "1,1";
    string bndryPad      = "";
    string periodic      = "";
    bool   verbose       = false;

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

    // Parse and check the command-line arguments
    Teuchos::CommandLineProcessor::EParseCommandLineReturn
      parseReturn = clp.parse(argc,argv);
    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
      return 0;
    if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
      return 1;

    // Construct the communicator
    Teuchos::RCP< const Teuchos::Comm< int > > comm =
      Teuchos::DefaultComm< int >::getComm();
    int pid = comm->getRank();
    if (verbose && pid == 0)
      cout << "Command line parsed" << endl;

    // Convert the command-line arguments into usable arrays
    Array< int > dimVals;
    Array< int > axisCommSizeVals;
    Array< int > commPadVals;
    Array< int > bndryPadVals;
    Array< int > periodicFlags;
    splitStringOfIntsWithCommas(dims         , dimVals         );
    splitStringOfIntsWithCommas(axisCommSizes, axisCommSizeVals);
    splitStringOfIntsWithCommas(commPad      , commPadVals     );
    splitStringOfIntsWithCommas(bndryPad     , bndryPadVals    );
    splitStringOfIntsWithCommas(periodic     , periodicFlags   );

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
    if (numDims < 2 || numDims > 3)
    {
      if (pid == 0)
        cout << "TEST 1 must be 2D or 3D." << endl;
#ifdef HAVE_MPI
      MPI_Finalize();
#endif
      return 1;
    }

    // Print the arrays that will be passed to the MDMap constructor
    if (verbose && pid == 0)
    {
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
    Teuchos::RCP< const Domi::MDMap< int > > mdMap =
      Teuchos::rcp(new Domi::MDMap< int >(mdComm,
                                          dimVals(),
                                          commPadVals(),
                                          bndryPadVals() ));

    // Construct the MDVector and extract the MDArrayView
    Domi::MDVector< double, int > mdVector(mdMap);
    Domi::MDArrayView<double> mdArray = mdVector.getDataNonConst();

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
          mdArray(i,j) = (double) mdMap->getGlobalIndex(lid);
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
            mdArray(i,j,k) = (double) mdMap->getGlobalIndex(lid);
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
          int    lid = mdMap->getLocalIndex(tuple(i,j));
          double gid = (double) mdMap->getGlobalIndex(lid);
          if (verbose)
            cout << pid << ": mdVector(" << i << "," << j << ") = "
                 << mdArray(i,j) << " (should be " << gid << ")"
                 << endl;
          assert(mdArray(i,j) == gid);
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
            int    lid = mdMap->getLocalIndex(tuple(i,j,k));
            double gid = (double) mdMap->getGlobalIndex(lid);
            if (verbose)
              cout << pid << ": mdVector(" << i << "," << j << "," << k
                   << ") = " << mdArray(i,j,k) << " (should be " << gid << ")"
                   << endl;
            assert(mdArray(i,j,k) == gid);
          }
        }
      }
    }

    // Check the lower padding along axis 0
    if (verbose)
      cout << pid << ": checking lower padding along axis 0" << endl;
    Domi::MDArrayView< const double > view = mdVector.getLowerPadData(0);
    if (numDims == 2)
    {
      for (int j = 0; j < view.dimension(1); ++j)
      {
        for (int i = 0; i < view.dimension(0); ++i)
        {
          int lid = mdMap->getLocalIndex(tuple(i,j));
          double gid = (double) mdMap->getGlobalIndex(lid);
          if (mdMap->isBndryPad(tuple(i,j)))
            gid = -1;
          if (verbose)
            cout << pid << ": mdVector(" << i << "," << j << ") = "
                 << mdArray(i,j) << " (should be " << gid << ")"
                 << endl;
          assert(view(i,j) == gid);
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
            double gid = (double) mdMap->getGlobalIndex(lid);
            if (mdMap->isBndryPad(tuple(i,j,k)))
                gid = -1;
            if (verbose)
              cout << pid << ": mdVector(" << i << "," << j << "," << k
                   << ") = " << view(i,j,k) << " (should be " << gid << ")"
                   << endl;
            assert(view(i,j,k) == gid);
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
          double gid = (double) mdMap->getGlobalIndex(lid);
          if (mdMap->isBndryPad(tuple(i+ioff,j)))
            gid = -1;
          if (verbose)
            cout << pid << ": mdVector(" << i+ioff << "," << j << ") = "
                 << view(i,j) << " (should be " << gid << ")" << endl;
          assert(view(i,j) == gid);
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
            double gid = (double) mdMap->getGlobalIndex(lid);
            if (mdMap->isBndryPad(tuple(i+ioff,j,k)))
              gid = -1;
            if (verbose)
              cout << pid << ": mdVector(" << i+ioff << "," << j << "," << k
                   << ") = " << view(i,j,k) << " (should be " << gid << ")"
                   << endl;
            assert(view(i,j,k) == gid);
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
          double gid = (double) mdMap->getGlobalIndex(lid);
          if (mdMap->isBndryPad(tuple(i,j)))
            gid = -1;
          if (verbose)
            cout << pid << ": mdVector(" << i << "," << j << ") = "
                 << mdArray(i,j) << " (should be " << gid << ")"
                 << endl;
          assert(view(i,j) == gid);
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
            double gid = (double) mdMap->getGlobalIndex(lid);
            if (mdMap->isBndryPad(tuple(i,j,k)))
                gid = -1;
            if (verbose)
              cout << pid << ": mdVector(" << i << "," << j << "," << k
                   << ") = " << view(i,j,k) << " (should be " << gid << ")"
                   << endl;
            assert(view(i,j,k) == gid);
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
          double gid = (double) mdMap->getGlobalIndex(lid);
          if (mdMap->isBndryPad(tuple(i,j+joff)))
            gid = -1;
          if (verbose)
            cout << pid << ": mdVector(" << i << "," << j+joff << ") = "
                 << view(i,j) << " (should be " << gid << ")" << endl;
          assert(view(i,j) == gid);
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
            double gid = (double) mdMap->getGlobalIndex(lid);
            if (mdMap->isBndryPad(tuple(i,j+joff,k)))
              gid = -1;
            if (verbose)
              cout << pid << ": mdVector(" << i << "," << j+joff << "," << k
                   << ") = " << view(i,j,k) << " (should be " << gid << ")"
                   << endl;
            assert(view(i,j,k) == gid);
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
            double gid = (double) mdMap->getGlobalIndex(lid);
            if (mdMap->isBndryPad(tuple(i,j,k)))
              gid = -1;
            if (verbose)
              cout << pid << ": mdVector(" << i << "," << j << "," << k
                   << ") = " << view(i,j,k) << " (should be " << gid << ")"
                   << endl;
            assert(view(i,j,k) == gid);
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
            double gid = (double) mdMap->getGlobalIndex(lid);
            if (mdMap->isBndryPad(tuple(i,j,k+koff)))
              gid = -1;
            if (verbose)
              cout << pid << ": mdVector(" << i << "," << j << "," << k+koff
                   << ") = " << view(i,j,k) << " (should be " << gid << ")"
                   << endl;
            assert(view(i,j,k) == gid);
          }
        }
      }
    }

    if (pid == 0)
      cout << "End Result: TEST PASSED" << endl;
  }
  catch (Teuchos::CommandLineProcessor::HelpPrinted &e)
  {
    return 0;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

}
