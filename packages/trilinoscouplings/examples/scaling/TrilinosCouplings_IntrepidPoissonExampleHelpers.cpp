// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "TrilinosCouplings_IntrepidPoissonExampleHelpers.hpp"
#include <iostream>
#include <sstream>


namespace TrilinosCouplings {
namespace IntrepidPoissonExample {


// Row-major matrix
class Matrix3{
public:
  Matrix3():v_(9,0){}

  Matrix3(const std::vector<double> & vals):v_(vals){};

  // Computes C = this * B
  Matrix3 operator*(const Matrix3& B) const {
    Matrix3 C;
    for(int i=0; i<3; i++)
      for(int j=0; j<3; j++)
        for(int k=0; k<3; k++)
          C(i,j) += (*this)(i,k)* B(k,j);      
    return C;
  }

  double & operator()(int i, int j) { 
    return v_[i*3+j];
  }

  const double & operator()(int i, int j) const{ 
    return v_[i*3+j];
  }


  void print(std::ostream & os) {
    os<<"[ "<<v_[0]<<" "<<v_[1]<<" "<<v_[2]<<" ]\n"
      <<"[ "<<v_[3]<<" "<<v_[4]<<" "<<v_[5]<<" ]\n"
      <<"[ "<<v_[6]<<" "<<v_[7]<<" "<<v_[8]<<" ]"<<std::endl;
  }

  const std::vector<double> & get() { return v_;};

private:
  std::vector<double> v_;
};


Matrix3 DiagonalMatrix(const std::vector<double> &d) {
  Matrix3 rv;
  for(int i=0; i<3; i++)
    rv(i,i) = d[i];
  return rv;
}

Matrix3 X_Rotation(double theta) {
  std::vector<double> v {1.0, 0.0, 0.0, /**/ 0.0, cos(theta), -sin(theta), /**/ 0.0, sin(theta), cos(theta)};
  return Matrix3(v);
}

Matrix3 Y_Rotation(double theta) {
  std::vector<double> v {cos(theta), 0.0, sin(theta), /**/ 0.0, 1.0, 0.0, /**/-sin(theta), 0.0, cos(theta)};
  return Matrix3(v);
}

Matrix3 Z_Rotation(double theta) {
  std::vector<double> v {cos(theta), -sin(theta), 0.0, /**/ sin(theta), cos(theta), 0.0,/**/ 0.0, 0.0, 1.0};
  return Matrix3(v);
}

Matrix3 Transpose(const Matrix3 & A) {
  Matrix3 B;
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      B(j,i) = A(i,j);
  return B;
}



bool use_diffusion_ = false;
double materialTensorOffDiagonalValue_;
  
Matrix3 rotation_, strength_,diff_total_;

std::vector<double> matrix2D_;

bool useDiffusionMatrix() { 
  return use_diffusion_;
}


double getMaterialTensorOffDiagonalValue () {
  if(use_diffusion_)
    throw std::runtime_error("setMaterialTensorOffDiagonalValue has not been called");
  return materialTensorOffDiagonalValue_;
}
void setMaterialTensorOffDiagonalValue (const double newVal) {
  materialTensorOffDiagonalValue_ = newVal;  
}


const std::vector<double>& getDiffusionMatrix() {
  if(!use_diffusion_)
    throw std::runtime_error("setDiffusionRotationStrength has not been called");
  return diff_total_.get();
}



const std::vector<double>& getDiffusionMatrix2D() {
  // Gets the x/y sub-matrix
  if(!use_diffusion_)
    throw std::runtime_error("setDiffusionRotationStrength has not been called");
  matrix2D_.resize(4);
  matrix2D_[0] = diff_total_(0,0);
  matrix2D_[1] = diff_total_(0,1);
  matrix2D_[2] = diff_total_(1,0);
  matrix2D_[3] = diff_total_(1,1);

  return matrix2D_;
}


void setDiffusionRotationAndStrength(const std::vector<double>& theta,  const std::vector<double>& diagonal) {

  rotation_ = Z_Rotation(theta[2]*M_PI/180.0) * Y_Rotation(theta[1]*M_PI/180.0) * X_Rotation(theta[0]*M_PI/180.0);
  strength_ = DiagonalMatrix(diagonal);    
  diff_total_ = rotation_ * strength_ * Transpose(rotation_);
  use_diffusion_ = true;
}



std::string
makeMeshInput (const int nx, const int ny, const int nz)
{
  using std::endl;
  std::ostringstream os;

  TEUCHOS_TEST_FOR_EXCEPTION( nx <= 0 || ny <= 0 || nz <= 0,
    std::invalid_argument, "nx, ny, and nz must all be positive.");

  os << "mesh" << endl
     << "\trectilinear" << endl
     << "\t\tnx = " << nx << endl
     << "\t\tny = " << ny << endl
     << "\t\tnz = " << nz << endl
     << "\t\tbx = 1" << endl
     << "\t\tby = 1" << endl
     << "\t\tbz = 1" << endl
     << "\t\tgmin = 0 0 0" << endl
     << "\t\tgmax = 1 1 1" << endl
     << "\tend" << endl
     << "\tset assign" << endl
     << "\t\tsideset, ilo, 1" << endl
     << "\t\tsideset, jlo, 2" << endl
     << "\t\tsideset, klo, 3" << endl
     << "\t\tsideset, ihi, 4" << endl
     << "\t\tsideset, jhi, 5" << endl
     << "\t\tsideset, khi, 6" << endl
     << "\tend" << endl
     << "end";
  return os.str ();
}

void
setCommandLineArgumentDefaults (int& nx,
                                int& ny,
                                int& nz,
                                std::string& xmlInputParamsFile,
                                std::string& solverName,
                                bool& verbose,
                                bool& debug)
{
  nx = 20;
  ny = 20;
  nz = 20;
  xmlInputParamsFile = "";
  solverName = "GMRES";
  verbose = false;
  debug = false;
}

void
setUpCommandLineArguments (Teuchos::CommandLineProcessor& cmdp,
                           int& nx,
                           int& ny,
                           int& nz,
                           std::string& xmlInputParamsFile,
                           std::string& solverName,
                           double& tol,
                           int& maxNumIters,
                           bool& verbose,
                           bool& debug)
{
  cmdp.setOption ("nx", &nx, "Number of cells along the x dimension");
  cmdp.setOption ("ny", &ny, "Number of cells along the y dimension");
  cmdp.setOption ("nz", &nz, "Number of cells along the z dimension");
  cmdp.setOption ("inputParams", &xmlInputParamsFile, "XML file of input "
                  "parameters, which we read if specified and not \"\".  "
                  "If it has a \"meshInput\" parameter, we use its "
                  "std::string value as the Pamgen mesh specification.  "
                  "Otherwise, we tell Pamgen to make a cube, using "
                  "nx, ny, and nz.");
  cmdp.setOption ("solverName", &solverName, "Name of iterative linear solver "
                  "to use for solving the linear system.  You may use any name "
                  "that Belos::SolverFactory understands.  Examples include "
                  "\"GMRES\" and \"CG\".");
  cmdp.setOption ("tol", &tol, "Tolerance for the linear solve.  If not "
                  "specified, this is read from the input ParameterList (read "
                  "from the XML file).  If specified, this overrides any value "
                  "in the input ParameterList.");
  cmdp.setOption ("maxNumIters", &maxNumIters, "Maximum number of iterations "
                  "in the linear solve.  If not specified, this is read from "
                  "the input ParameterList (read from the XML file).  If "
                  "specified, this overrides any value in the input "
                  "ParameterList.");
  cmdp.setOption ("verbose", "quiet", &verbose,
                  "Whether to print verbose status output.");
  cmdp.setOption ("debug", "release", &debug,
                  "Whether to print copious debugging output to stderr.");
}

void
parseCommandLineArguments (Teuchos::CommandLineProcessor& cmdp,
                           bool& printedHelp,
                           int argc,
                           char* argv[],
                           int& nx,
                           int& ny,
                           int& nz,
                           std::string& xmlInputParamsFile,
                           std::string& solverName,
                           bool& verbose,
                           bool& debug)
{
  using Teuchos::CommandLineProcessor;

  const CommandLineProcessor::EParseCommandLineReturn parseResult =
    cmdp.parse (argc, argv);
  if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED) {
    printedHelp = true;
  }
  else {
    printedHelp = false;
    TEUCHOS_TEST_FOR_EXCEPTION(
      parseResult != CommandLineProcessor::PARSE_SUCCESSFUL,
      std::invalid_argument, "Failed to parse command-line arguments.");
    TEUCHOS_TEST_FOR_EXCEPTION(
      xmlInputParamsFile == "" && (nx <= 0 || ny <= 0 || nz <= 0),
      std::invalid_argument, "If no XML parameters filename is specified (via "
      "--inputParams), then the number of cells along each dimension of the "
      "mesh (--nx, --ny, and --nz) must be positive.");
  }
}

} // namespace IntrepidPoissonExample
} // namespace TrilinosCouplings
