
// STD includes
#include <iostream>
#include <fstream>
#include <sstream>

// Teuchos includes
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_Tuple.hpp>
using Teuchos::tuple;

// Domi includes
#include <Domi_MDComm.hpp>
#include <Domi_MDMap.hpp>
#include <Domi_MDVector.hpp>
#include <Domi_Slice.hpp>
using Domi::MDArrayView;
using Domi::MDComm;
using Domi::MDMap;
using Domi::MDVector;
using Domi::Slice;

// Macros
#define SCAL double

int main(int argc, char *argv[])
{
  // Construct an MPI session object. If Trilinos was compiled with
  // MPI, then this will initialize MPI properly, and when this object
  // is destructed, MPI will be finalized properly.
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  // Get the default Teuchos::Comm
  Teuchos::RCP< const Teuchos::Comm< int > > comm =
    Teuchos::DefaultComm< int >::getComm();

  // Default problem parameters. These can be over-ridden by command
  // line arguments
  int px  =  -1;     // Number of processors in x-direction
  int py  =  -1;     // Number of processors in y-direction
  int nx  =  33;     // Number of grid points in x-direction
  int ny  =  33;     // Number of grid points in y-direction
  int nt  = 100;     // Number of time steps
  int sf  =   1;     // Screen output frequency
  int ff  =   1;     // File output frequency
  SCAL Re = 50.0;    // Reynolds number
  bool verbose = false;

  // Set up the command-line processor
  Teuchos::CommandLineProcessor clp;
  clp.throwExceptions(false);
  clp.setOption("px", &px, "Number of processors along x-axis");
  clp.setOption("py", &py, "Number of processors along y-axis");
  clp.setOption("nx", &nx, "Global dimension along x-axis"    );
  clp.setOption("ny", &ny, "Global dimension along y-axis"    );
  clp.setOption("nt", &nt, "Number of time steps to take"     );
  clp.setOption("sf", &sf, "Screen output frequency"          );
  clp.setOption("ff", &ff, "File output frequency"            );
  clp.setOption("Re", &Re, "Reynolds number"                  );
  clp.setOption("verbose", "quiet", &verbose, "Verbose or quiet output");

  // Parse the command line
  Teuchos::CommandLineProcessor::EParseCommandLineReturn
    parseReturn = clp.parse(argc,argv);
  if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
    return 0;
  if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL)
    return 1;

  // Sanity check
  if (verbose && (px * py != comm->getSize()))
  {
    cout << "Communicator size " << comm->getSize() << " not compatible with "
         << "processor decomposition (" << px << "," << py << "). Exiting"
         << endl;
    return 1;
  }

  // Print the problem parameters
  if (verbose && (comm->getRank() == 0))
    cout << "Problem Parameters:" << endl
         << "-------------------" << endl
         << "px = " << px << endl
         << "py = " << py << endl
         << "nx = " << nx << endl
         << "ny = " << ny << endl
         << "nt = " << nt << endl
         << "sf = " << sf << endl
         << "ff = " << ff << endl
         << "Re = " << Re << endl;

  // Output the problem parameters to a file
  std::filebuf fb;
  fb.open("params.py", std::ios::out);
  std::ostream os(&fb);
  os << "nx = " << nx << endl
     << "ny = " << ny << endl
     << "nt = " << nt << endl
     << "sf = " << sf << endl
     << "ff = " << ff << endl;
  fb.close();

  // Construct an MDComm
  Teuchos::RCP< const Domi::MDComm > mdComm =
    Teuchos::rcp(new MDComm(comm,             // Teuchos Comm
                            2,                // # dimensions
                            tuple(px, py)));  // # procs in each direction

  // Construct an MDMap suitable for 2nd order finite
  // differencing. This will create a map with no boundary padding and
  // communication padding with a width of 1 in each direction.
  Teuchos::RCP< Domi::MDMap<> > mdMap =
    Teuchos::rcp(new Domi::MDMap<>(mdComm,          // Domi MDComm
                                   tuple(nx, ny),   // Dimensions
                                   tuple( 1,  1))); // Communication padding

  // Construct the structured, distributed field variables, all
  // collocated at the same grid points
  MDVector< SCAL > u(mdMap);        // x component of velocity at step n
  MDVector< SCAL > v(mdMap);        // y component of velocity at step n
  MDVector< SCAL > u_new(mdMap);    // x component of velocity at step n+1
  MDVector< SCAL > v_new(mdMap);    // y component of velocity at step n+1

  // We will need the underlying MDArrayViews to actually index into
  // these fields
  MDArrayView< SCAL > ua     = u.getDataNonConst();
  MDArrayView< SCAL > va     = v.getDataNonConst();
  MDArrayView< SCAL > ua_new = u_new.getDataNonConst();
  MDArrayView< SCAL > va_new = v_new.getDataNonConst();

  // We will need a slice of u along the top boundary in order to
  // enforce the initial condition of a sliding lid
  MDVector< SCAL > u_top(    u    , 1, ny-1);
  MDVector< SCAL > u_top_new(u_new, 1, ny-1);

  // We will need the underlying MDArrayView to actually index into
  // this slice
  MDArrayView< SCAL > ua_top     = u_top.getDataNonConst();
  MDArrayView< SCAL > ua_top_new = u_top_new.getDataNonConst();

  // Set the initial conditions
  u.putScalar(0.0);
  v.putScalar(0.0);
  u_new.putScalar(0.0);
  v_new.putScalar(0.0);

  // Finish the initial conditions by setting the u velocity along the
  // top lid to a value of one
  if (u_top.onSubcommunicator())
  {
    Slice iBounds = u.getLocalBounds(0);
    for (int i = iBounds.start(); i < iBounds.stop(); ++i)
    {
      ua_top(i)     = 1.0;
      ua_top_new(i) = 1.0;
    }
  }

  // Problem parameter setup
  SCAL delta_x = 1.0 / (nx-1);
  SCAL delta_y = 1.0 / (ny-1);
  SCAL Re_h    = Re * std::max(delta_x, delta_y);
  SCAL delta_t = 2.0 * std::min(delta_x, delta_y);
  SCAL cfl     = std::min(delta_x, delta_y) / delta_t;
  if (verbose && (comm->getRank() == 0))
  {
    cout << "Delta-x = " << delta_x << endl
         << "Delta-y = " << delta_y << endl
         << "Re_h    = " << Re_h    << endl
         << "Delta-t = " << delta_t << endl
         << "CFL     = " << cfl     << endl << endl;
  }
  Slice iBounds = mdMap->getLocalInteriorBounds(0);
  Slice jBounds = mdMap->getLocalInteriorBounds(1);
  int ii, jj;

  // Write binary initial condition output files
  cout << "Time step 0, time = 0" << endl;
  if (verbose && (comm->getRank() == 0)) cout << "    Writing u to u0.bin" << endl;
  u.writeBinary("u0.bin");
  if (verbose && (comm->getRank() == 0)) cout << "    Writing v to v0.bin" << endl;
  v.writeBinary("v0.bin");

  // Perform the explicit (Euler) time stepping
  for (int n = 0; n < nt; ++n)
  {
    // Simple output
    if (verbose && ((n+1) % sf == 0))
      cout << "Time step " << n+1 << ", time = " << (n+1) * delta_t << endl;

    // Update the communication padding
    u.updateCommPad();
    v.updateCommPad();

    // Advance the velocities
    for (int j = jBounds.start(); j < jBounds.stop(); ++j)
    {
      for (int i = iBounds.start(); i < iBounds.stop(); ++i)
      {
        // Handle upwinding indexes
        ii = i+1;
        jj = j+1;
        if (ua(i,j) < 0.0) ii = i;
        if (va(i,j) < 0.0) jj = j;
        // X-component
        ua_new(i,j) = ua(i,j) + delta_t *
          ((ua(i+1,j) - 2 * ua(i,j) + ua(i-1,j)) / (delta_x * delta_x) +
           (ua(i,j+1) - 2 * ua(1,j) + ua(i,j-1)) / (delta_y * delta_y) -
           Re * ua(i,j) * (ua(ii,j) - ua(ii-1,j)) / delta_x -
           Re * va(i,j) * (ua(i,jj) - ua(i,jj-1)) / delta_y);

        // Y-component
        va_new(i,j) = va(i,j) + delta_t *
          ((va(i+1,j) - 2 * va(i,j) + va(i-1,j)) / (delta_x * delta_x) +
           (va(i,j+1) - 2 * va(1,j) + va(i,j-1)) / (delta_y * delta_y) -
           Re * ua(i,j) * (va(ii,j) - va(ii-1,j)) / delta_x -
           Re * va(i,j) * (va(i,jj) - va(i,jj-1)) / delta_y);
      }
    }

    // Write binary output files
    if ((n+1) % ff == 0)
    {
      std::stringstream uname;
      uname << "u" << n+1 << ".bin";
      if (verbose && (comm->getRank() == 0)) cout << "    Writing u to " << uname.str() << endl;
      u.writeBinary(uname.str());
      std::stringstream vname;
      vname << "v" << n+1 << ".bin";
      if (verbose && (comm->getRank() == 0)) cout << "    Writing v to " << vname.str() << endl;
      v.writeBinary(vname.str());
    }

    // Prepare for the next time step by swapping the memory for the
    // local field arrays and their "new" counterparts
    MDArrayView< SCAL > temp = ua;
    ua     = ua_new;
    ua_new = temp;
    temp   = va;
    va     = va_new;
    va_new = temp;
  }

  return 0;
}
