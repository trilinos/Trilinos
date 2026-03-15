// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file  test_01.cpp
    \brief Test StatusTest input mechanism in OptimizationSolver.
*/

#include "ROL_StdObjective.hpp"
#include "ROL_Problem.hpp"
#include "ROL_Solver.hpp"
#include "ROL_BatchManager.hpp"
#include "ROL_UserInputGenerator.hpp"
#include "ROL_Stream.hpp"

#include "ROL_OED_Factory.hpp"
#include "ROL_OED_StdMomentOperator.hpp"

#include "ROL_GlobalMPISession.hpp"

#include <iostream>

template<typename Real>
class LinearRegressionModel : public ROL::StdObjective<Real> {
private:
  std::vector<std::vector<Real>> data_;
  const unsigned nrows_, ncols_;

public:
  LinearRegressionModel() : nrows_(101), ncols_(4) {
    std::vector<std::vector<Real>> data(nrows_, std::vector<Real>(ncols_, 0.0));
    data[  0][0] = -0.217409; data[  0][1] = 65.981300; data[  0][2] =  0.184674; data[  0][3] = -0.597562;
    data[  1][0] = -0.205971; data[  1][1] = 65.977500; data[  1][2] =  0.173927; data[  1][3] = -0.597455;
    data[  2][0] = -0.194621; data[  2][1] = 65.966400; data[  2][2] =  0.163119; data[  2][3] = -0.597137;
    data[  3][0] = -0.183367; data[  3][1] = 65.948400; data[  3][2] =  0.152253; data[  3][3] = -0.596606;
    data[  4][0] = -0.172217; data[  4][1] = 65.923800; data[  4][2] =  0.141332; data[  4][3] = -0.595864;
    data[  5][0] = -0.161179; data[  5][1] = 65.893200; data[  5][2] =  0.130362; data[  5][3] = -0.594910;
    data[  6][0] = -0.150261; data[  6][1] = 65.856800; data[  6][2] =  0.119345; data[  6][3] = -0.593744;
    data[  7][0] = -0.139470; data[  7][1] = 65.815000; data[  7][2] =  0.108286; data[  7][3] = -0.592368;
    data[  8][0] = -0.128814; data[  8][1] = 65.768300; data[  8][2] =  0.097188; data[  8][3] = -0.590780;
    data[  9][0] = -0.118300; data[  9][1] = 65.717100; data[  9][2] =  0.086056; data[  9][3] = -0.588984;
    data[ 10][0] = -0.107936; data[ 10][1] = 65.661600; data[ 10][2] =  0.074893; data[ 10][3] = -0.586977;
    data[ 11][0] = -0.097728; data[ 11][1] = 65.602300; data[ 11][2] =  0.063703; data[ 11][3] = -0.584763;
    data[ 12][0] = -0.087684; data[ 12][1] = 65.539500; data[ 12][2] =  0.052491; data[ 12][3] = -0.582340;
    data[ 13][0] = -0.077811; data[ 13][1] = 65.473600; data[ 13][2] =  0.041260; data[ 13][3] = -0.579711;
    data[ 14][0] = -0.068115; data[ 14][1] = 65.404900; data[ 14][2] =  0.030015; data[ 14][3] = -0.576876;
    data[ 15][0] = -0.058603; data[ 15][1] = 65.333800; data[ 15][2] =  0.018759; data[ 15][3] = -0.573835;
    data[ 16][0] = -0.049281; data[ 16][1] = 65.260700; data[ 16][2] =  0.007496; data[ 16][3] = -0.570591;
    data[ 17][0] = -0.040157; data[ 17][1] = 65.185700; data[ 17][2] = -0.003769; data[ 17][3] = -0.567144;
    data[ 18][0] = -0.031235; data[ 18][1] = 65.109400; data[ 18][2] = -0.015033; data[ 18][3] = -0.563496;
    data[ 19][0] = -0.022523; data[ 19][1] = 65.031900; data[ 19][2] = -0.026292; data[ 19][3] = -0.559647;
    data[ 20][0] = -0.014025; data[ 20][1] = 64.953600; data[ 20][2] = -0.037541; data[ 20][3] = -0.555600;
    data[ 21][0] = -0.005749; data[ 21][1] = 64.874900; data[ 21][2] = -0.048777; data[ 21][3] = -0.551355;
    data[ 22][0] =  0.002301; data[ 22][1] = 64.795900; data[ 22][2] = -0.059996; data[ 22][3] = -0.546914;
    data[ 23][0] =  0.010120; data[ 23][1] = 64.717000; data[ 23][2] = -0.071193; data[ 23][3] = -0.542279;
    data[ 24][0] =  0.017701; data[ 24][1] = 64.638500; data[ 24][2] = -0.082365; data[ 24][3] = -0.537452;
    data[ 25][0] =  0.025040; data[ 25][1] = 64.560600; data[ 25][2] = -0.093508; data[ 25][3] = -0.532433;
    data[ 26][0] =  0.032132; data[ 26][1] = 64.483500; data[ 26][2] = -0.104618; data[ 26][3] = -0.527225;
    data[ 27][0] =  0.038972; data[ 27][1] = 64.407700; data[ 27][2] = -0.115690; data[ 27][3] = -0.521830;
    data[ 28][0] =  0.045556; data[ 28][1] = 64.333200; data[ 28][2] = -0.126722; data[ 28][3] = -0.516250;
    data[ 29][0] =  0.051879; data[ 29][1] = 64.260300; data[ 29][2] = -0.137708; data[ 29][3] = -0.510486;
    data[ 30][0] =  0.057937; data[ 30][1] = 64.189300; data[ 30][2] = -0.148645; data[ 30][3] = -0.504541;
    data[ 31][0] =  0.063725; data[ 31][1] = 64.120300; data[ 31][2] = -0.159530; data[ 31][3] = -0.498416;
    data[ 32][0] =  0.069241; data[ 32][1] = 64.053600; data[ 32][2] = -0.170358; data[ 32][3] = -0.492114;
    data[ 33][0] =  0.074480; data[ 33][1] = 63.989300; data[ 33][2] = -0.181125; data[ 33][3] = -0.485638;
    data[ 34][0] =  0.079439; data[ 34][1] = 63.927700; data[ 34][2] = -0.191828; data[ 34][3] = -0.478989;
    data[ 35][0] =  0.084114; data[ 35][1] = 63.868900; data[ 35][2] = -0.202463; data[ 35][3] = -0.472170;
    data[ 36][0] =  0.088503; data[ 36][1] = 63.813100; data[ 36][2] = -0.213026; data[ 36][3] = -0.465183;
    data[ 37][0] =  0.092602; data[ 37][1] = 63.760400; data[ 37][2] = -0.223513; data[ 37][3] = -0.458031;
    data[ 38][0] =  0.096410; data[ 38][1] = 63.711100; data[ 38][2] = -0.233921; data[ 38][3] = -0.450716;
    data[ 39][0] =  0.099922; data[ 39][1] = 63.665100; data[ 39][2] = -0.244246; data[ 39][3] = -0.443241;
    data[ 40][0] =  0.103138; data[ 40][1] = 63.622700; data[ 40][2] = -0.254484; data[ 40][3] = -0.435608;
    data[ 41][0] =  0.106054; data[ 41][1] = 63.584000; data[ 41][2] = -0.264631; data[ 41][3] = -0.427821;
    data[ 42][0] =  0.108669; data[ 42][1] = 63.549100; data[ 42][2] = -0.274685; data[ 42][3] = -0.419882;
    data[ 43][0] =  0.110981; data[ 43][1] = 63.517900; data[ 43][2] = -0.284641; data[ 43][3] = -0.411793;
    data[ 44][0] =  0.112989; data[ 44][1] = 63.490800; data[ 44][2] = -0.294495; data[ 44][3] = -0.403558;
    data[ 45][0] =  0.114692; data[ 45][1] = 63.467600; data[ 45][2] = -0.304246; data[ 45][3] = -0.395180;
    data[ 46][0] =  0.116087; data[ 46][1] = 63.448600; data[ 46][2] = -0.313888; data[ 46][3] = -0.386662;
    data[ 47][0] =  0.117175; data[ 47][1] = 63.433600; data[ 47][2] = -0.323418; data[ 47][3] = -0.378006;
    data[ 48][0] =  0.117954; data[ 48][1] = 63.422800; data[ 48][2] = -0.332834; data[ 48][3] = -0.369215;
    data[ 49][0] =  0.118424; data[ 49][1] = 63.416200; data[ 49][2] = -0.342131; data[ 49][3] = -0.360294;
    data[ 50][0] =  0.118585; data[ 50][1] = 63.413800; data[ 50][2] = -0.351307; data[ 50][3] = -0.351245;
    data[ 51][0] =  0.118436; data[ 51][1] = 63.415700; data[ 51][2] = -0.360358; data[ 51][3] = -0.342070;
    data[ 52][0] =  0.117978; data[ 52][1] = 63.421700; data[ 52][2] = -0.369281; data[ 52][3] = -0.332775;
    data[ 53][0] =  0.117211; data[ 53][1] = 63.431900; data[ 53][2] = -0.378073; data[ 53][3] = -0.323361;
    data[ 54][0] =  0.116135; data[ 54][1] = 63.446300; data[ 54][2] = -0.386730; data[ 54][3] = -0.313832;
    data[ 55][0] =  0.114752; data[ 55][1] = 63.464800; data[ 55][2] = -0.395250; data[ 55][3] = -0.304191;
    data[ 56][0] =  0.113061; data[ 56][1] = 63.487400; data[ 56][2] = -0.403630; data[ 56][3] = -0.294443;
    data[ 57][0] =  0.111065; data[ 57][1] = 63.514000; data[ 57][2] = -0.411866; data[ 57][3] = -0.284590;
    data[ 58][0] =  0.108764; data[ 58][1] = 63.544500; data[ 58][2] = -0.419956; data[ 58][3] = -0.274636;
    data[ 59][0] =  0.106160; data[ 59][1] = 63.578900; data[ 59][2] = -0.427897; data[ 59][3] = -0.264584;
    data[ 60][0] =  0.103255; data[ 60][1] = 63.617100; data[ 60][2] = -0.435686; data[ 60][3] = -0.254438;
    data[ 61][0] =  0.100050; data[ 61][1] = 63.658900; data[ 61][2] = -0.443319; data[ 61][3] = -0.244202;
    data[ 62][0] =  0.096548; data[ 62][1] = 63.704300; data[ 62][2] = -0.450796; data[ 62][3] = -0.233879;
    data[ 63][0] =  0.092751; data[ 63][1] = 63.753100; data[ 63][2] = -0.458112; data[ 63][3] = -0.223473;
    data[ 64][0] =  0.088662; data[ 64][1] = 63.805300; data[ 64][2] = -0.465266; data[ 64][3] = -0.212988;
    data[ 65][0] =  0.084283; data[ 65][1] = 63.860500; data[ 65][2] = -0.472254; data[ 65][3] = -0.202427;
    data[ 66][0] =  0.079617; data[ 66][1] = 63.918800; data[ 66][2] = -0.479074; data[ 66][3] = -0.191794;
    data[ 67][0] =  0.074667; data[ 67][1] = 63.979900; data[ 67][2] = -0.485724; data[ 67][3] = -0.181093;
    data[ 68][0] =  0.069436; data[ 68][1] = 64.043600; data[ 68][2] = -0.492202; data[ 68][3] = -0.170327;
    data[ 69][0] =  0.063929; data[ 69][1] = 64.109800; data[ 69][2] = -0.498505; data[ 69][3] = -0.159501;
    data[ 70][0] =  0.058148; data[ 70][1] = 64.178200; data[ 70][2] = -0.504630; data[ 70][3] = -0.148619;
    data[ 71][0] =  0.052098; data[ 71][1] = 64.248800; data[ 71][2] = -0.510577; data[ 71][3] = -0.137683;
    data[ 72][0] =  0.045782; data[ 72][1] = 64.321200; data[ 72][2] = -0.516342; data[ 72][3] = -0.126699;
    data[ 73][0] =  0.039205; data[ 73][1] = 64.395200; data[ 73][2] = -0.521923; data[ 73][3] = -0.115670;
    data[ 74][0] =  0.032371; data[ 74][1] = 64.470600; data[ 74][2] = -0.527319; data[ 74][3] = -0.104599;
    data[ 75][0] =  0.025285; data[ 75][1] = 64.547100; data[ 75][2] = -0.532528; data[ 75][3] = -0.093492;
    data[ 76][0] =  0.017951; data[ 76][1] = 64.624500; data[ 76][2] = -0.537547; data[ 76][3] = -0.082351;
    data[ 77][0] =  0.010375; data[ 77][1] = 64.702600; data[ 77][2] = -0.542376; data[ 77][3] = -0.071181;
    data[ 78][0] =  0.002561; data[ 78][1] = 64.781100; data[ 78][2] = -0.547012; data[ 78][3] = -0.059985;
    data[ 79][0] = -0.005485; data[ 79][1] = 64.859600; data[ 79][2] = -0.551453; data[ 79][3] = -0.048768;
    data[ 80][0] = -0.013758; data[ 80][1] = 64.937900; data[ 80][2] = -0.555699; data[ 80][3] = -0.037534;
    data[ 81][0] = -0.022251; data[ 81][1] = 65.015800; data[ 81][2] = -0.559747; data[ 81][3] = -0.026287;
    data[ 82][0] = -0.030960; data[ 82][1] = 65.092800; data[ 82][2] = -0.563596; data[ 82][3] = -0.015030;
    data[ 83][0] = -0.039879; data[ 83][1] = 65.168800; data[ 83][2] = -0.567245; data[ 83][3] = -0.003768;
    data[ 84][0] = -0.049001; data[ 84][1] = 65.243300; data[ 84][2] = -0.570693; data[ 84][3] =  0.007495;
    data[ 85][0] = -0.058320; data[ 85][1] = 65.316100; data[ 85][2] = -0.573937; data[ 85][3] =  0.018756;
    data[ 86][0] = -0.067829; data[ 86][1] = 65.386800; data[ 86][2] = -0.576978; data[ 86][3] =  0.030010;
    data[ 87][0] = -0.077523; data[ 87][1] = 65.455100; data[ 87][2] = -0.579814; data[ 87][3] =  0.041253;
    data[ 88][0] = -0.087394; data[ 88][1] = 65.520600; data[ 88][2] = -0.582444; data[ 88][3] =  0.052482;
    data[ 89][0] = -0.097435; data[ 89][1] = 65.583100; data[ 89][2] = -0.584867; data[ 89][3] =  0.063692;
    data[ 90][0] = -0.107640; data[ 90][1] = 65.642100; data[ 90][2] = -0.587082; data[ 90][3] =  0.074879;
    data[ 91][0] = -0.118001; data[ 91][1] = 65.697200; data[ 91][2] = -0.589088; data[ 91][3] =  0.086040;
    data[ 92][0] = -0.128512; data[ 92][1] = 65.748200; data[ 92][2] = -0.590885; data[ 92][3] =  0.097171;
    data[ 93][0] = -0.139164; data[ 93][1] = 65.794600; data[ 93][2] = -0.592473; data[ 93][3] =  0.108266;
    data[ 94][0] = -0.149951; data[ 94][1] = 65.836000; data[ 94][2] = -0.593850; data[ 94][3] =  0.119324;
    data[ 95][0] = -0.160864; data[ 95][1] = 65.872100; data[ 95][2] = -0.595015; data[ 95][3] =  0.130339;
    data[ 96][0] = -0.171896; data[ 96][1] = 65.902600; data[ 96][2] = -0.595970; data[ 96][3] =  0.141307;
    data[ 97][0] = -0.183039; data[ 97][1] = 65.926900; data[ 97][2] = -0.596712; data[ 97][3] =  0.152226;
    data[ 98][0] = -0.194285; data[ 98][1] = 65.944600; data[ 98][2] = -0.597243; data[ 98][3] =  0.163090;
    data[ 99][0] = -0.205626; data[ 99][1] = 65.955500; data[ 99][2] = -0.597562; data[ 99][3] =  0.173896;
    data[100][0] = -0.205583; data[100][1] = 65.959200; data[100][2] = -0.597668; data[100][3] =  0.173927;
    data_.clear();
    data_.assign(data.begin(),data.end());
  }

  Real value(const std::vector<Real> &theta, Real &tol) override {
    Real val(0);
    unsigned ind = static_cast<unsigned>(ROL::Objective<Real>::getParameter()[0]);
    for (unsigned i = 0u; i < ncols_; ++i)
      val += theta[i] * (data_[ind])[i];
    return val;
  }

  void gradient(std::vector<Real> &g, const std::vector<Real> &theta, Real &tol) override {
    unsigned ind = static_cast<unsigned>(ROL::Objective<Real>::getParameter()[0]);
    g.assign(data_[ind].begin(),data_[ind].end());
  }
};

template<typename Real>
class PolynomialNoise : public ROL::OED::Noise<Real> {
private:
  const Real alpha_;

public:
  PolynomialNoise(Real alpha = Real(1)) : alpha_(alpha) {}

  Real evaluate(const std::vector<Real> &x) const override {
    return std::exp(alpha_ * std::abs(x[0])) / std::exp(alpha_);
  }
};

template<typename Real>
class RegularizationOperator : public ROL::LinearOperator<Real> {
private:
  const Real alpha_;

public:
  RegularizationOperator(Real alpha = Real(1)) : alpha_(alpha) {}

  void apply(ROL::Vector<Real> &Px, const ROL::Vector<Real> &x, Real &tol) const override {
    Px.set(x);
    Px.scale(alpha_);
  }

  void applyInverse(ROL::Vector<Real> &Px, const ROL::Vector<Real> &x, Real &tol) const override {
    Px.set(x);
    Px.scale(static_cast<Real>(1)/alpha_);
  }

  void applyAdjoint(ROL::Vector<Real> &Px, const ROL::Vector<Real> &x, Real &tol) const override {
    apply(Px,x,tol);
  }

  void applyAdjointInverse(ROL::Vector<Real> &Px, const ROL::Vector<Real> &x, Real &tol) const override {
    applyInverse(Px,x,tol);
  }
};

typedef double RealT;

int main(int argc, char *argv[]) {

  ROL::GlobalMPISession mpiSession(&argc, &argv);

  // This little trick lets us print to std::cout only if a (dummy) command-line argument is provided.
  int iprint     = argc - 1;
  ROL::Ptr<std::ostream> outStream;
  ROL::nullstream bhs; // outputs nothing
  if (iprint > 0)
    outStream = ROL::makePtrFromRef(std::cout);
  else
    outStream = ROL::makePtrFromRef(bhs);

  int errorFlag  = 0;

  // *** Test body.

  try {
    std::string filename = "input_03.xml";
    auto parlist = ROL::getParametersFromXmlFile( filename );

    // Setup parameter vector and polynomial model
    const RealT alpha = parlist->sublist("Problem").get("Noise Decay Rate", 5.0);
    auto theta = ROL::makePtr<ROL::StdVector<RealT>>(4,1);
    auto model = ROL::makePtr<LinearRegressionModel<RealT>>();
    auto noise = ROL::makePtr<PolynomialNoise<RealT>>(alpha);

    // Setup experiment sample generator
    std::ofstream ptfile, wtfile;
    ptfile.open("points3.txt");
    wtfile.open("weights3.txt");
    for (unsigned i = 0u; i < 101u; ++i) {
      //ptfile << std::scientific << std::setprecision(16);
      //wtfile << std::scientific << std::setprecision(16);
      ptfile << static_cast<RealT>(i) << std::endl;
      wtfile << static_cast<RealT>(1u) / static_cast<RealT>(101u) << std::endl;
    }
    ptfile.close();
    wtfile.close();
    auto bman = ROL::makePtr<ROL::BatchManager<RealT>>();
    auto sampler = ROL::makePtr<ROL::UserInputGenerator<RealT>>("points3.txt","weights3.txt",101,1,bman);
    auto isampler = ROL::makePtr<ROL::UserInputGenerator<RealT>>("points3.txt","weights3.txt",101,1,bman);

    // Setup factory
    bool homNoise = true;
    std::string regType = "Least Squares";
    std::string ocType = parlist->sublist("OED").get("Optimality Type","A");
    auto type = ROL::OED::StringToRegressionType(regType);
    bool useSVD = parlist->sublist("Problem").get("Use SVD",true);
    RealT tolSVD = parlist->sublist("Problem").get("SVD Drop Tolerance",1e-12);
    auto M = ROL::makePtr<ROL::OED::StdMomentOperator<RealT>>(type,homNoise,noise,useSVD,tolSVD);
    bool addTik = parlist->sublist("Problem").get("Use Tikhonov",false);
    if (addTik) {
      RealT beta  = parlist->sublist("Problem").get("Tikhonov Parameter",1e-4);
      auto P = ROL::makePtr<RegularizationOperator<RealT>>(beta);
      M->setPerturbation(P);
    }

    auto factory = ROL::makePtr<ROL::OED::Factory<RealT>>(model,sampler,theta,M,*parlist);
    if (parlist->sublist("Problem").get("Use Budget Constraint",false)) {
      auto cost = factory->getDesign()->clone();
      cost->setScalar(static_cast<RealT>(1));
      RealT budget = parlist->sublist("Problem").get("Budget",5.0);
      factory->setBudgetConstraint(cost,budget);
    }

    // Generate optimization problem
    auto problem = factory->get(*parlist,sampler);
    problem->setProjectionAlgorithm(*parlist);
    problem->finalize(false,true,*outStream);
    auto test = factory->getDesign()->clone();
    test->randomize(1,2);
    problem->check(true,*outStream,test,0.1);

    // Setup ROL solver
    std::clock_t timer = std::clock();
    auto solver = ROL::makePtr<ROL::Solver<RealT>>(problem,*parlist);
    solver->solve(*outStream);
    *outStream << "  " << ocType << "-optimal design time:      "
               << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC)
               << " seconds" << std::endl;
    factory->profile(*outStream);
    std::stringstream dname;
    dname << ocType << "_optimal_design_ex3";
    factory->printDesign(dname.str());
  }
  catch (std::logic_error& err) {
    *outStream << err.what() << std::endl;
    errorFlag = -1000;
  }; // end try

  if (errorFlag != 0)
    std::cout << "End Result: TEST FAILED" << std::endl;
  else
    std::cout << "End Result: TEST PASSED" << std::endl;

  return 0;

}
