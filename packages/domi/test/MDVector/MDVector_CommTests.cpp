
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
       "constructing an MDVector.  The 'test' option specifies an integer ID\n"
       "of astandard suite of tests.  If 'test' is zero, then the user can\n"
       "specifythe parameters for constructing an MDVector and receive output\n"
       "for what that MDVector looks like.  In this case, the 'dims',\n"
       "'axisCommSizes' and 'commPad' options should be one or more comma-\n"
       "separated integers.  The number of integers specified in the 'dims'\n"
       "option will determine the dimension of the MDMap (1D, 2D, 3D, etc.).\n"
       "The 'axisCommSizes' and 'commPad' options can have fewer integers\n"
       "specified than the 'dims' option, and those unspecified values will\n"
       "then receive default values.\n");

    int    test            = 0;
    string filename        = "";
    string dims            = "8,8";
    string axisCommSizes   = "-1";
    string commPad         = "0";
    string bndryPad        = "";
    string periodic        = "";
    string subOrigin       = "";
    string subDims         = "";
    bool   includeBndryPad = false;
    bool   verbose         = false;

    clp.setOption("test"         , &test,
		  "Predefined test specification");
    clp.setOption("writeto"      , &filename,
                  "Write the MDVector to the given filename.  An empty\n"
                  "                                 "
                  "string suppresses output");
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
    clp.setOption("subOrigin"    , &subOrigin,
                  "Comma-separated coordinates of origin of subfield");
    clp.setOption("subDims"      , &subDims,
                  "Comma-separated dimensions of dims of subfield");
    clp.setOption("writeBndryPad", "doNotWriteBndryPad", &includeBndryPad,
                  "Write (or do not write) bndryPad points to the binary file");
    clp.setOption("verbose"      , "quiet"             , &verbose,
		  "Verbose or quiet output");

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

    // Convert the command-line arguments into usable arrays
    Array< int > dimVals;
    Array< int > axisCommSizeVals;
    Array< int > commPadVals;
    Array< int > bndryPadVals;
    Array< int > periodicFlags;
    Array< int > subOriginVals;
    Array< int > subDimsVals;
    splitStringOfIntsWithCommas(dims         , dimVals         );
    splitStringOfIntsWithCommas(axisCommSizes, axisCommSizeVals);
    splitStringOfIntsWithCommas(commPad      , commPadVals     );
    splitStringOfIntsWithCommas(bndryPad     , bndryPadVals    );
    splitStringOfIntsWithCommas(periodic     , periodicFlags   );
    splitStringOfIntsWithCommas(subOrigin    , subOriginVals   );
    splitStringOfIntsWithCommas(subDims      , subDimsVals     );

    int numDims = dimVals.size();
    for (int i = axisCommSizeVals.size(); i < numDims; ++i)
      axisCommSizeVals.push_back(-1);

    if (pid == 0)
    {
      cout << "Command line parsed" << endl;
      cout << "test:          " << test << endl;
    }

    ////////////////////////////////////////////////////////////////////////////////
    // TEST 0: User-specified construction parameters and a printout
    // of the resulting field.
    if (test == 0)
    {
      // Print the arrays that will be passed to the MDMap constructor
      if (pid == 0)
      {
	cout << "dims:          " << dimVals          << endl;
	cout << "axisCommSizes: " << axisCommSizeVals << endl;
	cout << "commPad:       " << commPadVals      << endl;
	cout << "bndryPad:      " << bndryPadVals     << endl;
	cout << "periodic:      " << periodicFlags    << endl;
        cout << "subOrigin      " << subOriginVals    << endl;
        cout << "subDims        " << subDimsVals      << endl;
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

    //   // Subvector?
    //   if (subOriginVals.size() > 0 && subDimsVals.size() > 0)
    //   {
    //     // Construct the subvector
    //     Domi::MDVector subVector =
    //       Teuchos::rcp(new Domi::MdVector("s", "submdVector", mdVector, subOriginVals(),
    //                                   subDimsVals()));
    //     // Print the submdVector
    //     cout << endl << pid << ": " << submdVector->description() << endl;
    //     // Write the submdVector to a file
    //     if (filename.length() > 0)
    //       submdVector->writeBinary(filename, includeBndryPad);
    //   }
    //   else
    //   {
    //     // Print the MdVector
    //     cout << endl << pid << ": " << mdVector->description() << endl;
    //     // Write the MdVector to a file
    //     if (filename.length() > 0)
    //       mdVector->writeBinary(filename, includeBndryPad);
    //   }
    }

    ////////////////////////////////////////////////////////////////////////////////
    // TEST 1: Build an MDVector with size, processor decomposition,
    //         bndryPad points and commPad specified by the user.
    //         Initialize such that owned points are assigned their
    //         global ID, and padding regions are assigned -1.
    //         Perform an updateCommPad() and check that the owned
    //         points and boundary padding points are unchanged and
    //         that communication padding points are updated with the
    //         corresponding global ID.
    if (test == 1)
    {
      // Default commPad should be at least 1,1 to ensure that there is
      // communication
      if (commPad == "0") commPad = "1,1";
      splitStringOfIntsWithCommas(commPad, commPadVals);

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
      if (pid == 0)
      {
	cout << "dims:          " << dimVals          << endl;
	cout << "axisCommSizes: " << axisCommSizeVals << endl;
	cout << "commPad:       " << commPadVals      << endl;
	cout << "bndryPad:      " << bndryPadVals     << endl;
	cout << "periodic:      " << periodicFlags    << endl;
        cout << "subOrigin      " << subOriginVals    << endl;
        cout << "subDims        " << subDimsVals      << endl;
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

#if 0
      // Print the bndryPad view dims and origins
      for (int i = 0; i < numDims; ++i)
      {
	cout << "Axis " << i << endl;
	cout << "  Lower bndryPad origin: " << mdVector.getLowerBndryPadOrigin(i)
             << endl;
	cout << "  Upper bndryPad origin: " << mdVector.getUpperBndryPadOrigin(i)
             << endl;
	cout << "  BndryPad dims: " << mdVector.getBndryPadDims(i)
             << endl;
      }
#endif

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
                     << ") = " << mdArray(i,j) << " (should be " << gid << ")"
                     << endl;
              assert(mdArray(i,j,k) == gid);
            }
	  }
	}
      }

// #if 0
// 	// Left bndryPad points should still be -1
// 	Domi::MdVector f0 =
// 	  Domi::MdVector("f0",
//                      "left bndryPad",
//                      mdVectorPtr,
//                      mdVector.getLowerBndryPadOrigin(0),
//                      mdVector.getBndryPadDims(0));
// 	for (int j = 0; j < f0.getLocalOwnDims(1); ++j)
// 	  for (int i = 0; i < f0.getLocalOwnDims(0); ++i)
// 	  {
// 	    if (verbose)
// 	      cout << pid << ": f0(" << i << "," << j << ") = "
//                    << f0(i,j) << " (should be -1)"
//                    << endl;
// 	    assert(f0(i,j) == -1.0);
// 	  }

// 	// Bottom bndryPad points should still be -1
// 	Domi::MdVector f1 =
// 	  Domi::MdVector("f1",
//                      "bottom bndryPad",
//                      mdVectorPtr,
//                      mdVector.getLowerBndryPadOrigin(1),
//                      mdVector.getBndryPadDims(1));
// 	for (int j = 0; j < f1.getLocalOwnDims(1); ++j)
// 	  for (int i = 0; i < f1.getLocalOwnDims(0); ++i)
// 	  {
// 	    if (verbose)
// 	      cout << pid << ": f1(" << i << "," << j << ") = "
//                    << f1(i,j) << " (should be -1)"
//                    << endl;
// 	    assert(f1(i,j) == -1.0);
// 	  }

// 	// Right bndryPad points should still be -1
// 	Domi::MdVector f2 =
// 	  Domi::MdVector("f2",
//                      "right bndryPad",
//                      mdVectorPtr,
//                      mdVector.getUpperBndryPadOrigin(0),
//                      mdVector.getBndryPadDims(0));
// 	for (int j = 0; j < f2.getLocalOwnDims(1); ++j)
// 	  for (int i = 0; i < f2.getLocalOwnDims(0); ++i)
// 	  {
// 	    if (verbose)
// 	      cout << pid << ": f2(" << i << "," << j << ") = "
//                    << f2(i,j) << " (should be -1)"
//                    << endl;
// 	    assert(f2(i,j) == -1.0);
// 	  }

// 	// Top bndryPad points should still be -1
// 	Domi::MdVector f3 =
// 	  Domi::MdVector("f3",
//                      "top bndryPad",
//                      mdVectorPtr,
//                      mdVector.getUpperBndryPadOrigin(1),
//                      mdVector.getBndryPadDims(1));
// 	for (int j = 0; j < f3.getLocalOwnDims(1); ++j)
// 	  for (int i = 0; i < f3.getLocalOwnDims(0); ++i)
// 	  {
// 	    if (verbose)
// 	      cout << pid << ": f3(" << i << "," << j << ") = "
//                    << f3(i,j) << " (should be -1)"
//                    << endl;
// 	    assert(f3(i,j) == -1.0);
// 	  }

// 	// Left commPad points should now have global IDs
// 	Domi::MdVector h0 =
// 	  Domi::MdVector("h0",
//                      "left commPad",
//                      mdVectorPtr,
//                      0,
//                      Domi::LowerCommPad);
// 	i0 = h0.getLocalOrigin(0);
// 	j0 = h0.getLocalOrigin(1);
//         cout << "MdMap commPad[0] = " << mdMap->getCommPad(0) << endl;
//         cout << "MdMap commPad[1] = " << mdMap->getCommPad(1) << endl;
//         myCommPad = mdMap->getCommPad(1);
//         cout << "MdMap commPad[0] = " << mdMap->getCommPad(0) << endl;
//         cout << "MdMap commPad[1] = " << mdMap->getCommPad(1) << endl;
//         //cout << pid << ": commPad = " << myCommPad << endl;
// 	//for (int j = myCommPad[0]; j < h0.getLocalOwnDims(1)-myCommPad[1]; ++j)
// 	for (int j = mdMap->getCommPad(1)[0]; j < h0.getLocalOwnDims(1)-mdMap->getCommPad(1)[1]; ++j)
// 	  for (int i = 0; i < h0.getLocalOwnDims(0); ++i)
// 	  {
// 	    Domi::Ordinal lid = mdMap->getLocalElement(tuple(i+i0,j+j0));
// 	    Domi::Ordinal gid = mdMap->getGlobalElement(lid);
// 	    Domi::Ordinal global_j = mdMap->getGlobalAxisIndex(gid)[1];
// 	    if (global_j < 0 || global_j > mdVector.getGlobalDims(1)-1)
// 	      gid = -1;
// 	    if (verbose)
// 	      cout << pid << ": h0(" << i << "," << j << ") = "
//                    << h0(i,j) << " (should be " << gid << ")"
//                    << endl;
// 	    assert(h0(i,j) == gid);
// 	  }

// 	// Bottom commPad points should now have global IDs
// 	Domi::MdVector h1 =
// 	  Domi::MdVector("h1",
//                      "bottom commPad",
//                      mdVectorPtr,
//                      1,
//                      Domi::LowerCommPad);
// 	i0 = h1.getLocalOrigin(0);
// 	j0 = h1.getLocalOrigin(1);
//         myCommPad = mdMap->getCommPad(0);
// 	for (int j = 0; j < h1.getLocalOwnDims(1); ++j)
// 	  for (int i = myCommPad[0]; i < h1.getLocalOwnDims(0)-myCommPad[1]; ++i)
// 	  {
// 	    Domi::Ordinal lid = mdMap->getLocalElement(tuple(i+i0,j+j0));
// 	    Domi::Ordinal gid = mdMap->getGlobalElement(lid);
// 	    Domi::Ordinal global_i = mdMap->getGlobalAxisIndex(gid)[0];
// 	    if (global_i < 0 || global_i > mdVector.getGlobalDims(0)-1)
// 	      gid = -1;
// 	    if (verbose)
// 	      cout << pid << ": h1(" << i << "," << j << ") = "
//                    << h1(i,j) << " (should be " << gid << ")"
//                    << endl;
// 	    assert(h1(i,j) == gid);
// 	  }	    

// 	// Right commPad points should now have global IDs
// 	Domi::MdVector h2 =
// 	  Domi::MdVector("h2",
//                      "right commPad",
//                      mdVectorPtr,
//                      0,
//                      Domi::UpperCommPad);
// 	i0 = h2.getLocalOrigin(0);
// 	j0 = h2.getLocalOrigin(1);
//         //myCommPad = mdMap->getCommPad(1);
// 	//for (int j = myCommPad[0]; j < h2.getLocalOwnDims(1)-myCommPad[1]; ++j)
// 	for (int j = mdMap->getCommPad(1)[0]; j < h2.getLocalOwnDims(1)-mdMap->getCommPad(1)[1]; ++j)
// 	  for (int i = 0; i < h2.getLocalOwnDims(0); ++i)
// 	  {
// 	    Domi::Ordinal lid = mdMap->getLocalElement(tuple(i+i0,j+j0));
// 	    Domi::Ordinal gid = mdMap->getGlobalElement(lid);
// 	    Domi::Ordinal global_j = mdMap->getGlobalAxisIndex(gid)[1];
// 	    if (global_j < 0 || global_j > mdVector.getGlobalDims(1)-1)
// 	      gid = -1;
// 	    if (verbose)
// 	      cout << pid << ": h2(" << i << "," << j << ") = "
//                    << h2(i,j) << " (should be " << gid << ")"
//                    << endl;
// 	    assert(h2(i,j) == gid);
// 	  }
// #endif

// 	// Top commPad points should now have global IDs
// 	Domi::MdVector h3 =
// 	  Domi::MdVector("h3",
//                      "top commPad",
//                      mdVectorPtr,
//                      1,
//                      Domi::UpperCommPad);
// 	i0 = h3.getLocalOrigin(0);
// 	j0 = h3.getLocalOrigin(1);
//         myCommPad = mdMap->getCommPad(0);
// 	for (int j = 0; j < h3.getLocalOwnDims(1); ++j)
// 	  for (int i = myCommPad[0]; i < h3.getLocalOwnDims(0)-myCommPad[1]; ++i)
// 	  {
//             cout << "(i,j) = (" << i << "," << j << ")" << endl;
// 	    Domi::Ordinal lid = mdMap->getLocalElement(tuple(i+i0,j+j0));
// 	    Domi::Ordinal gid = mdMap->getGlobalElement(lid);
// 	    Domi::Ordinal global_i = mdMap->getGlobalAxisIndex(gid)[0];
// 	    if (global_i < 0 || global_i > mdVector.getGlobalDims(0)-1)
// 	      gid = -1;
// 	    if (verbose)
// 	      cout << pid << ": h3(" << i << "," << j << ") = "
//                    << h3(i,j) << " (should be " << gid << ")"
//                    << endl;
// 	    assert(h3(i,j) == gid);
// 	  }

//       }
//       else  // 3D
//       {
// 	Domi::Ordinal i0;
// 	Domi::Ordinal j0;
// 	Domi::Ordinal k0;
//         Teuchos::Tuple< int, 2 > iCommPad;
//         Teuchos::Tuple< int, 2 > jCommPad;
//         Teuchos::Tuple< int, 2 > kCommPad;

// 	// Interior data should be unchanged
// 	for (int k = 0; k < mdVector.getLocalOwnDims(2); ++k)
// 	{
// 	  for (int j = 0; j < mdVector.getLocalOwnDims(1); ++j)
// 	  {
// 	    for (int i = 0; i < mdVector.getLocalOwnDims(0); ++i)
// 	    {
// 	      Domi::Ordinal lid = mdMap->getLocalElement(tuple(i,j,k));
// 	      Domi::Ordinal gid = mdMap->getGlobalElement(lid);
// 	      if (verbose)
// 		cout << pid << ": mdVector(" << i << "," << j << "," << k
//                      << ") = " << mdVector(i,j,k) << " (should be " << gid
//                      << ")" << endl;
// 	      assert(mdVector(i,j,k) == gid);
// 	    }
// 	  }
// 	}

// 	// Left bndryPad points should still be -1
// 	Domi::MdVector f0 =
// 	  Domi::MdVector("f0",
//                      "left bndryPad",
//                      mdVectorPtr,
//                      mdVector.getLowerBndryPadOrigin(0),
//                      mdVector.getBndryPadDims(0));
// 	for (int k = 0; k < f0.getLocalOwnDims(2); ++k)
// 	  for (int j = 0; j < f0.getLocalOwnDims(1); ++j)
// 	    for (int i = 0; i < f0.getLocalOwnDims(0); ++i)
// 	    {
// 	      if (verbose)
// 		cout << pid << ": f0(" << i << "," << j << "," << k
//                      << ") = " << f0(i,j,k) << " (should be -1)"
//                      << endl;
// 	      assert(f0(i,j,k) == -1.0);
// 	    }

// 	// Bottom bndryPad points should still be -1
// 	Domi::MdVector f1 =
// 	  Domi::MdVector("f1",
//                      "bottom bndryPad",
//                      mdVectorPtr,
//                      mdVector.getLowerBndryPadOrigin(1),
//                      mdVector.getBndryPadDims(1));
// 	for (int k = 0; k < f1.getLocalOwnDims(2); ++k)
// 	  for (int j = 0; j < f1.getLocalOwnDims(1); ++j)
// 	    for (int i = 0; i < f1.getLocalOwnDims(0); ++i)
// 	    {
// 	      if (verbose)
// 		cout << pid << ": f1(" << i << "," << j << "," << k
//                      << ") = " << f1(i,j,k) << " (should be -1)"
//                      << endl;
// 	      assert(f1(i,j,k) == -1.0);
// 	    }

// 	// Right bndryPad points should still be -1
// 	Domi::MdVector f2 =
// 	  Domi::MdVector("f2",
//                      "right bndryPad",
//                      mdVectorPtr,
//                      mdVector.getUpperBndryPadOrigin(0),
//                      mdVector.getBndryPadDims(0));
// 	for (int k = 0; k < f2.getLocalOwnDims(2); ++k)
// 	  for (int j = 0; j < f2.getLocalOwnDims(1); ++j)
// 	    for (int i = 0; i < f2.getLocalOwnDims(0); ++i)
// 	    {
// 	      if (verbose)
// 		cout << pid << ": f2(" << i << "," << j << "," << k
//                      << ") = " << f2(i,j,k) << " (should be -1)"
//                      << endl;
// 	      assert(f2(i,j,k) == -1.0);
// 	    }

// 	// Top bndryPad points should still be -1
// 	Domi::MdVector f3 =
// 	  Domi::MdVector("f3",
//                      "top bndryPad",
//                      mdVectorPtr,
//                      mdVector.getUpperBndryPadOrigin(1),
//                      mdVector.getBndryPadDims(1));
// 	for (int k = 0; k < f3.getLocalOwnDims(2); ++k)
// 	  for (int j = 0; j < f3.getLocalOwnDims(1); ++j)
// 	    for (int i = 0; i < f3.getLocalOwnDims(0); ++i)
// 	    {
// 	      if (verbose)
// 		cout << pid << ": f3(" << i << "," << j << "," << k
//                      << ") = " << f3(i,j,k) << " (should be -1)"
//                      << endl;
// 	      assert(f3(i,j,k) == -1.0);
// 	    }

// 	// Front bndryPad points should still be -1
// 	Domi::MdVector f4 =
// 	  Domi::MdVector("f4",
//                      "front bndryPad",
//                      mdVectorPtr,
//                      mdVector.getUpperBndryPadOrigin(2),
//                      mdVector.getBndryPadDims(2));
// 	for (int k = 0; k < f4.getLocalOwnDims(2); ++k)
// 	  for (int j = 0; j < f4.getLocalOwnDims(1); ++j)
// 	    for (int i = 0; i < f4.getLocalOwnDims(0); ++i)
// 	    {
// 	      if (verbose)
// 		cout << pid << ": f4(" << i << "," << j << "," << k
//                      << ") = " << f4(i,j,k) << " (should be -1)"
//                      << endl;
// 	      assert(f4(i,j,k) == -1.0);
// 	    }

// 	// Back bndryPad points should still be -1
// 	Domi::MdVector f5 =
// 	  Domi::MdVector("f5",
//                      "back bndryPad",
//                      mdVectorPtr,
//                      mdVector.getUpperBndryPadOrigin(2),
//                      mdVector.getBndryPadDims(2));
// 	for (int k = 0; k < f5.getLocalOwnDims(2); ++k)
// 	  for (int j = 0; j < f5.getLocalOwnDims(1); ++j)
// 	    for (int i = 0; i < f5.getLocalOwnDims(0); ++i)
// 	    {
// 	      if (verbose)
// 		cout << pid << ": f5(" << i << "," << j << "," << k
//                      << ") = " << f5(i,j,k) << " (should be -1)"
//                      << endl;
// 	      assert(f5(i,j,k) == -1.0);
// 	    }

// 	// Left commPad points should now have global IDs
// 	Domi::MdVector h0 =
// 	  Domi::MdVector("h0",
//                      "left commPad",
//                      mdVectorPtr,
//                      0,
//                      Domi::LowerCommPad);
// 	i0 = h0.getLocalOrigin(0);
// 	j0 = h0.getLocalOrigin(1);
// 	k0 = h0.getLocalOrigin(2);
//         jCommPad = mdMap->getCommPad(1);
//         kCommPad = mdMap->getCommPad(2);
// 	for (int k = kCommPad[0]; k < h0.getLocalOwnDims(2)-kCommPad[1]; ++k)
// 	  for (int j = jCommPad[0]; j < h0.getLocalOwnDims(1)-jCommPad[1]; ++j)
// 	    for (int i = 0; i < h0.getLocalOwnDims(0); ++i)
// 	    {
// 	      Domi::Ordinal lid = mdMap->getLocalElement(tuple(i+i0,j+j0,k+k0));
// 	      Domi::Ordinal gid = mdMap->getGlobalElement(lid);
// 	      Domi::Ordinal global_j = mdMap->getGlobalAxisIndex(gid)[1];
// 	      Domi::Ordinal global_k = mdMap->getGlobalAxisIndex(gid)[2];
// 	      if (global_j < 0 || global_j > mdVector.getGlobalDims(1)-1 ||
// 		  global_k < 0 || global_k > mdVector.getGlobalDims(2)-1)
// 		gid = -1;
// 	      if (verbose)
// 		cout << pid << ": h0(" << i << "," << j << "," << k
//                      << ") = " << h0(i,j,k) << " (should be " << gid
//                      << ")" << endl;
// 	      assert(h0(i,j,k) == gid);
// 	    }

// 	// Bottom commPad points should now have global IDs
// 	Domi::MdVector h1 =
// 	  Domi::MdVector("h1",
//                      "bottom commPad",
//                      mdVectorPtr,
//                      1,
//                      Domi::LowerCommPad);
// 	i0 = h1.getLocalOrigin(0);
// 	j0 = h1.getLocalOrigin(1);
// 	k0 = h1.getLocalOrigin(2);
//         iCommPad = mdMap->getCommPad(0);
//         kCommPad = mdMap->getCommPad(2);
// 	for (int k = kCommPad[0]; k < h1.getLocalOwnDims(2)-kCommPad[1]; ++k)
// 	  for (int j = 0; j < h1.getLocalOwnDims(1); ++j)
// 	    for (int i = iCommPad[0]; i < h1.getLocalOwnDims(0)-iCommPad[1]; ++i)
// 	    {
// 	      Domi::Ordinal lid = mdMap->getLocalElement(tuple(i+i0,j+j0,k+k0));
// 	      Domi::Ordinal gid = mdMap->getGlobalElement(lid);
// 	      Domi::Ordinal global_i = mdMap->getGlobalAxisIndex(gid)[0];
// 	      Domi::Ordinal global_k = mdMap->getGlobalAxisIndex(gid)[2];
// 	      if (global_i < 0 || global_i > mdVector.getGlobalDims(0)-1 ||
// 		  global_k < 0 || global_k > mdVector.getGlobalDims(2)-1)
// 		gid = -1;
// 	      if (verbose)
// 		cout << pid << ": h1(" << i << "," << j << "," << k
//                      << ") = " << h1(i,j,k) << " (should be " << gid
//                      << ")" << endl;
// 	      assert(h1(i,j,k) == gid);
// 	    }

// 	// Right commPad points should now have global IDs
// 	Domi::MdVector h2 =
// 	  Domi::MdVector("h2",
//                      "right commPad",
//                      mdVectorPtr,
//                      0,
//                      Domi::UpperCommPad);
// 	i0 = h2.getLocalOrigin(0);
// 	j0 = h2.getLocalOrigin(1);
// 	k0 = h2.getLocalOrigin(2);
//         jCommPad = mdMap->getCommPad(1);
//         kCommPad = mdMap->getCommPad(2);
// 	for (int k = kCommPad[0]; k < h2.getLocalOwnDims(2)-kCommPad[1]; ++k)
// 	  for (int j = jCommPad[0]; j < h2.getLocalOwnDims(1)-jCommPad[1]; ++j)
// 	    for (int i = 0; i < h2.getLocalOwnDims(0); ++i)
// 	    {
// 	      Domi::Ordinal lid = mdMap->getLocalElement(tuple(i+i0,j+j0,k+k0));
// 	      Domi::Ordinal gid = mdMap->getGlobalElement(lid);
// 	      Domi::Ordinal global_j = mdMap->getGlobalAxisIndex(gid)[1];
// 	      Domi::Ordinal global_k = mdMap->getGlobalAxisIndex(gid)[2];
// 	      if (global_j < 0 || global_j > mdVector.getGlobalDims(1)-1 ||
// 		  global_k < 0 || global_k > mdVector.getGlobalDims(2)-1)
// 		gid = -1;
// 	      if (verbose)
// 		cout << pid << ": h2(" << i << "," << j << "," << k
//                      << ") = " << h2(i,j,k) << " (should be " << gid
//                      << ")" << endl;
// 	      assert(h2(i,j,k) == gid);
// 	    }

// 	// Top commPad points should now have global IDs
// 	Domi::MdVector h3 =
// 	  Domi::MdVector("h3",
//                      "top commPad",
//                      mdVectorPtr,
//                      1,
//                      Domi::UpperCommPad);
// 	i0 = h3.getLocalOrigin(0);
// 	j0 = h3.getLocalOrigin(1);
// 	k0 = h3.getLocalOrigin(2);
//         iCommPad = mdMap->getCommPad(0);
//         kCommPad = mdMap->getCommPad(2);
// 	for (int k = kCommPad[0]; k < h3.getLocalOwnDims(2)-kCommPad[1]; ++k)
// 	  for (int j = 0; j < h3.getLocalOwnDims(1); ++j)
// 	    for (int i = iCommPad[0]; i < h3.getLocalOwnDims(0)-iCommPad[1]; ++i)
// 	    {
// 	      Domi::Ordinal lid = mdMap->getLocalElement(tuple(i+i0,j+j0,k+k0));
// 	      Domi::Ordinal gid = mdMap->getGlobalElement(lid);
// 	      Domi::Ordinal global_i = mdMap->getGlobalAxisIndex(gid)[0];
// 	      Domi::Ordinal global_k = mdMap->getGlobalAxisIndex(gid)[2];
// 	      if (global_i < 0 || global_i > mdVector.getGlobalDims(0)-1 ||
// 		  global_k < 0 || global_k > mdVector.getGlobalDims(2)-1)
// 		gid = -1;
// 	      if (verbose)
// 		cout << pid << ": h3(" << i << "," << j << "," << k
//                      << ") = " << h3(i,j,k) << " (should be " << gid
//                      << ")" << endl;
// 	      assert(h3(i,j,k) == gid);
// 	    }

// 	// Front commPad points should now have global IDs
// 	Domi::MdVector h4 =
// 	  Domi::MdVector("h4",
//                      "front commPad",
//                      mdVectorPtr,
//                      2,
//                      Domi::LowerCommPad);
// 	i0 = h4.getLocalOrigin(0);
// 	j0 = h4.getLocalOrigin(1);
// 	k0 = h4.getLocalOrigin(2);
//         iCommPad = mdMap->getCommPad(0);
//         jCommPad = mdMap->getCommPad(1);
// 	for (int k = 0; k < h4.getLocalOwnDims(2); ++k)
// 	  for (int j = jCommPad[0]; j < h4.getLocalOwnDims(1)-jCommPad[1]; ++j)
// 	    for (int i = iCommPad[0]; i < h4.getLocalOwnDims(0)-iCommPad[1]; ++i)
// 	    {
// 	      Domi::Ordinal lid = mdMap->getLocalElement(tuple(i+i0,j+j0,k+k0));
// 	      Domi::Ordinal gid = mdMap->getGlobalElement(lid);
// 	      Domi::Ordinal global_i = mdMap->getGlobalAxisIndex(gid)[0];
// 	      Domi::Ordinal global_j = mdMap->getGlobalAxisIndex(gid)[1];
// 	      if (global_i < 0 || global_i > mdVector.getGlobalDims(0)-1 ||
// 		  global_j < 0 || global_j > mdVector.getGlobalDims(1)-1)
// 		gid = -1;
// 	      if (verbose)
// 		cout << pid << ": h4(" << i << "," << j << "," << k
//                      << ") = " << h4(i,j,k) << " (should be " << gid
//                      << ")" << endl;
// 	      assert(h4(i,j,k) == gid);
// 	    }

// 	// Back commPad points should now have global IDs
// 	Domi::MdVector h5 =
// 	  Domi::MdVector("h5",
//                      "back commPad",
//                      mdVectorPtr,
//                      2,
//                      Domi::UpperCommPad);
// 	i0 = h5.getLocalOrigin(0);
// 	j0 = h5.getLocalOrigin(1);
// 	k0 = h5.getLocalOrigin(2);
//         iCommPad = mdMap->getCommPad(0);
//         jCommPad = mdMap->getCommPad(1);
// 	for (int k = 0; k < h5.getLocalOwnDims(2); ++k)
// 	  for (int j = jCommPad[0]; j < h5.getLocalOwnDims(1)-jCommPad[1]; ++j)
// 	    for (int i = iCommPad[0]; i < h5.getLocalOwnDims(0)-iCommPad[1]; ++i)
// 	    {
// 	      Domi::Ordinal lid = mdMap->getLocalElement(tuple(i+i0,j+j0,k+k0));
// 	      Domi::Ordinal gid = mdMap->getGlobalElement(lid);
// 	      Domi::Ordinal global_i = mdMap->getGlobalAxisIndex(gid)[0];
// 	      Domi::Ordinal global_j = mdMap->getGlobalAxisIndex(gid)[1];
// 	      if (global_i < 0 || global_i > mdVector.getGlobalDims(0)-1 ||
// 		  global_j < 0 || global_j > mdVector.getGlobalDims(1)-1)
// 		gid = -1;
// 	      if (verbose)
// 		cout << pid << ": h5(" << i << "," << j << "," << k
//                      << ") = " << h5(i,j,k) << " (should be " << gid
//                      << ")" << endl;
// 	      assert(h5(i,j,k) == gid);
// 	    }

//        }
    }
    cout << "Test " << test << ": End Result: TEST PASSED on processor "
         << pid << endl;
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
