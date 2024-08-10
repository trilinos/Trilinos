// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PDEOPT_MULTIMAT_BRANCHING_PEBBL_H
#define ROL_PDEOPT_MULTIMAT_BRANCHING_PEBBL_H

#include "ROL_PEBBL_Interface.hpp"
#include "../../TOOLS/pdevector.hpp"
#include "MatrixMarket_Tpetra.hpp"

template<class Real>
class MultiMatBranchSub;

template<class Real>
class MultiMatBranching : public ROL::PEBBL::Branching<Real> {
private:
  ROL::Ptr<ROL::Vector<Real>> z0_;

  using ROL::PEBBL::Branching<Real>::verbosity_;
  using ROL::PEBBL::Branching<Real>::outStream_;
  using ROL::PEBBL::Branching<Real>::parlist_;

public:
  MultiMatBranching(const ROL::Ptr<ElasticityFactory<Real>>        &factory,
                    const ROL::Ptr<ROL::ParameterList>             &parlist,
                    const ROL::Ptr<ROL::PEBBL::BranchHelper<Real>> &bHelper,
                    int                                             verbosity = 0,
                    const ROL::Ptr<std::ostream>                   &outStream = ROL::nullPtr)
    : ROL::PEBBL::Branching<Real>(factory,parlist,bHelper,verbosity,outStream) {
    z0_ = factory->buildSolutionVector();
  }

  pebbl::branchSub* blankSub() {
    return new MultiMatBranchSub<Real>(parlist_,ROL::makePtrFromRef<MultiMatBranching<Real>>(*this),verbosity_,outStream_);
  }

//  pebbl::solution* iniitalGuess() {
//
//  }
}; // MultiMatBranching

template <class Real>
class MultiMatBranchSub : public ROL::PEBBL::BranchSub<Real> {
private:
  int method_;
  Real ctol_;
  int T_;
  std::string methodName_;
  ROL::Ptr<ROL::Vector<Real>> c_, x_;
  ROL::Ptr<ROL::Constraint<Real>> con_;
  int numIH_;

  using ROL::PEBBL::BranchSub<Real>::anyChild;
  using ROL::PEBBL::BranchSub<Real>::index_;
  using ROL::PEBBL::BranchSub<Real>::branching_;
  using ROL::PEBBL::BranchSub<Real>::problem0_;
  using ROL::PEBBL::BranchSub<Real>::solution_;
  using ROL::PEBBL::BranchSub<Real>::rndSolution_;
  using ROL::PEBBL::BranchSub<Real>::verbosity_;
  using ROL::PEBBL::BranchSub<Real>::outStream_;

  void round(ROL::Vector<Real> &rx, const ROL::Vector<Real> &x, Real t) const {
    rx.set(x);
    Teuchos::ArrayView<Real> data = getData(rx);
    Real val(0);
    for (auto it = data.begin(); it != data.end(); ++it) {
      val = *it;
      *it = (val < t ? std::floor(val) : std::ceil(val));
    }
  }

  ROL::Ptr<ROL::Vector<Real>> get(ROL::Vector<Real> &x, int ind) const {
    return dynamic_cast<ROL::PartitionedVector<Real>&>(x).get(ind);
  }

  ROL::Ptr<const ROL::Vector<Real>> get(const ROL::Vector<Real> &x, int ind) const {
    return dynamic_cast<const ROL::PartitionedVector<Real>&>(x).get(ind);
  }

  ROL::Ptr<const Tpetra::MultiVector<>> getMultiVector(const ROL::Vector<Real> &x) const {
    try {
      return ROL::dynamicPtrCast<const ROL::TpetraMultiVector<Real>>(get(x,0))->getVector();
    }
    catch (std::exception &e) {
      return ROL::dynamicPtrCast<const PDE_OptVector<Real>>(get(x,0))->getField()->getVector();
    }
  }

  ROL::Ptr<Tpetra::MultiVector<>> getMultiVector(ROL::Vector<Real> &x) const {
    try {
      return ROL::dynamicPtrCast<ROL::TpetraMultiVector<Real>>(get(x,0))->getVector();
    }
    catch (std::exception &e) {
      return ROL::dynamicPtrCast<PDE_OptVector<Real>>(get(x,0))->getField()->getVector();
    }
  }

  Teuchos::ArrayView<Real> getData(ROL::Vector<Real> &x) const {
    return (getMultiVector(x)->getDataNonConst(0))();
  }

  void zeroSlack(ROL::Vector<Real> &x) const {
    size_t nv = dynamic_cast<ROL::PartitionedVector<Real>&>(x).numVectors();
    for (size_t i = 1; i < nv; ++i) {
      get(x,i)->zero();
    }
  }

  void setSlack(ROL::Vector<Real> &x, const ROL::Vector<Real> &c) const {
    size_t nv = dynamic_cast<ROL::PartitionedVector<Real>&>(x).numVectors();
    for (size_t i = 1; i < nv; ++i) {
      get(x,i)->set(*get(c,i-1));
    }
    problem0_->getBoundConstraint()->project(x);
  }

  Real infeasibility(ROL::Vector<Real> &x) {
    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    zeroSlack(x);
    con_->update(x,ROL::UpdateType::Temp);
    con_->value(*c_,x,tol);
    setSlack(x,*c_);
    Real infeas(0);
    if (problem0_->getConstraint()==ROL::nullPtr) {
      x_->set(x);
      problem0_->getPolyhedralProjection()->project(*x_);
      x_->axpy(static_cast<Real>(-1),x);
      infeas = x_->norm();
    }
    else {
      con_->update(x,ROL::UpdateType::Temp);
      con_->value(*c_,x,tol);
      infeas = c_->norm();
    }
    return infeas;
  }

  void setDensity(ROL::Vector<Real> &out, const ROL::Vector<Real> &in) {
    try {
      ROL::Ptr<ROL::Vector<Real>>
        Fx = ROL::staticPtrCast<Filtered_Compliance_Objective<Real>>(
          ROL::staticPtrCast<ROL::SlacklessObjective<Real>>(
            problem0_->getObjective())->getObjective())->applyFilter(*get(in,0));
      Teuchos::ArrayView<Real>  data = getData(out);
      Teuchos::ArrayView<Real> Fdata = (ROL::staticPtrCast<ROL::TpetraMultiVector<Real>>(Fx)->getVector()->getDataNonConst(0))();
      const int nc = data.size()/T_;
      const int nx = 2*std::sqrt(nc/2); // nx = 2*ny
      const int ny =   std::sqrt(nc/2); // nc = nx*ny = 2*ny*ny
      for (int j = 0; j < ny; ++j) { 
        for (int i = 0; i < nx; ++i) {
          data[i+j*nx]  = static_cast<Real>(0.25)*Fdata[j*(nx+1)+i];
          data[i+j*nx] += static_cast<Real>(0.25)*Fdata[j*(nx+1)+(i+1)];
          data[i+j*nx] += static_cast<Real>(0.25)*Fdata[(j+1)*(nx+1)+(i+1)];
          data[i+j*nx] += static_cast<Real>(0.25)*Fdata[(j+1)*(nx+1)+i];
        }
      }
    }
    catch (std::exception &e) {
      out.set(in);
    }
  }

  void printVector(const ROL::Vector<Real> &x, std::string name) const {
    Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<>> vecWriter;
    std::string vecName = name + std::to_string(numIH_) + "_vec.txt";
    vecWriter.writeDenseFile(vecName, *getMultiVector(x));
    std::string mapName = name + std::to_string(numIH_) + "_map.txt";
    vecWriter.writeMapFile(mapName, *getMultiVector(x)->getMap());
  }

public:
  MultiMatBranchSub(const ROL::Ptr<ROL::ParameterList> &parlist,
                    const ROL::Ptr<ROL::PEBBL::Branching<Real>> &branching,
                    int verbosity = 0,
                    const ROL::Ptr<std::ostream> &outStream = ROL::nullPtr)
    : ROL::PEBBL::BranchSub<Real>(branching, verbosity, outStream), numIH_(0) {
    method_ = parlist->sublist("Problem").get("Incumbent Heuristic",0);
    ctol_   = parlist->sublist("Status Test").get("Constraint Tolerance",1e-8);
    std::vector<Real> ym = ROL::getArrayFromStringParameter<Real>(parlist->sublist("Problem"), "Young's Modulus");
    T_ = ym.size();
    methodName_ = "Default";
    if (method_==1)      methodName_ = "Constraint Violation Bisection";
    else if (method_==2) methodName_ = "Neighborhood Rounding";
    if (problem0_->getConstraint()==ROL::nullPtr) {
      x_ = solution_->clone();
      c_   = problem0_->getPolyhedralProjection()->getResidual()->clone();
      con_ = problem0_->getPolyhedralProjection()->getLinearConstraint();
    }
    else {
      c_   = problem0_->getResidualVector()->clone();
      con_ = problem0_->getConstraint();
    }
  }

  MultiMatBranchSub(const MultiMatBranchSub &rpbs)
    : ROL::PEBBL::BranchSub<Real>(rpbs),
      method_(rpbs.method_), ctol_(rpbs.ctol_), T_(rpbs.T_), methodName_(rpbs.methodName_),
      c_(rpbs.c_->clone()), x_(rpbs.x_==ROL::nullPtr ? ROL::nullPtr : rpbs.x_->clone()), con_(rpbs.con_), numIH_(0) {}

  void incumbentHeuristic() {
    std::ios_base::fmtflags flags(outStream_->flags());
    if (verbosity_ > 0) *outStream_ << std::scientific << std::setprecision(8);

    Real tol = std::sqrt(ROL::ROL_EPSILON<Real>());
    const Real zero(0), one(1);
    Teuchos::ArrayView<Real> data = getData(*rndSolution_);
    const int nc = data.size()/T_;
    int cnt(0);
    Real cnorm(0);
    if (method_==1) {
      // If greater than or equal to treshold, then set to 1 else set to 0
      // If threshold == 0, then design will be all 1
      // If threshold == 1, then design will be all 0 unless the component is 1
      const Real half(0.5);
      Real lo(0), up(1), mid(0.5);
      while (up-lo > std::sqrt(ROL::ROL_EPSILON<Real>())) {
        mid = half * (up + lo);
        round(*rndSolution_,*solution_,mid);
        cnorm = infeasibility(*rndSolution_);
        if (cnorm < ctol_) up = mid;
        else               lo = mid;
        if (verbosity_ > 1) {
          *outStream_ << "  cnt = "           << cnt
                      << "  infeasibility = " << cnorm
                      << "  lower bound = "   << lo
                      << "  upper bound = "   << up << std::endl;
        }
        cnt++;
      }
      if (lo == mid) {
        round(*rndSolution_,*solution_,up);
        cnorm = infeasibility(*rndSolution_);
      }
    }
    else if (method_==2) {
      //  -------------------------
      //  | 12| 13| 14| 15| 16| 17|
      //  -------------------------
      //  | 6 | 7 | 8 | 9 | 10| 11|
      //  -------------------------
      //  | 0 | 1 | 2 | 3 | 4 | 5 |
      //  -------------------------
      //rndSolution_->set(*solution_);
      setDensity(*rndSolution_,*solution_);
      const Real nneighbor(2), two(2);
      Real vol0(0), vol(0), d(0), nd(0), itol(0.35);
      for (int i = 0; i < nc; ++i) vol0 += data[i];
      int nx = 2*std::sqrt(nc/2); // nx = 2*ny
      int ny =   std::sqrt(nc/2); // nc = nx*ny = 2*ny*ny
      std::map<std::pair<int,int>,int> nn;
      for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny/2; ++j) {
          d = data[i+j*nx];
          cnt = 0;
          if (d > itol) { // && d < one-itol) {
            for (int k = -1; k < 2; ++k) {
              for (int l = -1; l < 2; ++l) {
                //if (k==0 || l==0) {
                if (k!=0 || l!=0) {
                  int ik = (i+k), jl = (j+l);
                  if ( ik >=0 && ik < nx && jl >= 0 && jl < ny/2 ) {
                    nd = data[ik + jl*nx];
                    if (nd >= one-itol) cnt++;
                  }
                }
              }
            }
          }
          nn.insert({{i,j},cnt});
          ////if (d >= one-itol || cnt >= nneighbor) {
          //if (cnt >= nneighbor) {
          //  data[i+j*nx] = one;
          //  data[i+(ny-1-j)*nx] = one;
          //  vol += (j!=ny-1-j ? two : one);
          //}
          //else {
          //  data[i+j*nx] = zero;
          //  data[i+(ny-1-j)*nx] = zero;
          //}
        }
      }
      for (auto it = nn.begin(); it != nn.end(); ++it) {
        int i = std::get<0>(it->first);
        int j = std::get<1>(it->first);
        cnt = it->second;
        if (cnt >= nneighbor) {
          data[i+j*nx] = one;
          data[i+(ny-1-j)*nx] = one;
          vol += (j!=ny-1-j ? two : one);
        }
        else {
          data[i+j*nx] = zero;
          data[i+(ny-1-j)*nx] = zero;
        }
      }
      if (verbosity_ > 1)
        *outStream_ << "Volume Difference: " << vol-vol0 << std::endl;
      if (vol <= vol0-two) {
        while (vol <= vol0-two) {
          // Get zero cells with nonzero neighbors
          std::multimap<int,std::pair<int,int>> map;
          for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny/2; ++j) {
              d = data[i+j*nx];
              cnt = 0;
              if (d==zero) {
                for (int k = -1; k < 2; ++k) {
                  for (int l = -1; l < 2; ++l) {
                    //if (k==0 || l==0) {
                    if (k!=0 || l!=0) {
                      int ik = (i+k), jl = (j+l);
                      if ( ik >=0 && ik < nx && jl >= 0 && jl < ny/2 ) {
                        nd = data[ik + jl*nx];
                        if (nd == one) cnt++;
                      }
                    }
                  }
                }
                if (cnt > 0) map.insert({cnt,{i,j}});
              }
            }
          }
          std::map<int,std::pair<int,int>>::iterator it = map.end();
          it--;                               // Grab cell with largest number of neighbors
          auto range0 = map.equal_range(it->first);
          for (auto it0 = range0.first; it0 != range0.second; ++it0) {
            int i = std::get<0>(it0->second); // x-index of cell
            int j = std::get<1>(it0->second); // y-index of cell
            data[i+j*nx] = one;               // Set cell to one
            data[i+(ny-1-j)*nx] = one;        // Set symmetric cell to one
            vol += (j!=ny-1-j ? two : one);   // Increase volume
            if (vol > vol0-two) break;
          }
        }
//        // Get zero cells with nonzero neighbors
//        std::multimap<int,std::pair<int,int>> map;
//        std::map<int,int> map2;
//        for (int i = 0; i < nx; ++i) {
//          for (int j = 0; j < ny/2; ++j) {
//            d = data[i+j*nx];
//            cnt = 0;
//            if (d==zero) {
//              for (int k = -1; k < 2; ++k) {
//                for (int l = -1; l < 2; ++l) {
//                  if (k==0 || l==0) {
//                  //if (k!=0 || l!=0) {
//                    int ik = (i+k), jl = (j+l);
//                    if ( ik >=0 && ik < nx && jl >= 0 && jl < ny/2 ) {
//                      nd = data[ik + jl*nx];
//                      if (nd == one) cnt++;
//                    }
//                  }
//                }
//              }
//              if (cnt > 0) {
//                map.insert({cnt,{i,j}});
//                map2.insert({i+j*nx,cnt});
//              }
//            }
//          }
//        }
//        while (vol <= vol0-two) {
//          std::map<int,std::pair<int,int>>::iterator it = map.end();
//          it--;                            // Grab cell with largest number of neighbors
//          int i = std::get<0>(it->second); // x-index of cell
//          int j = std::get<1>(it->second); // y-index of cell
//          auto range0 = map.equal_range(it->first);
//          Real r = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
//          for (auto it0 = range0.first; it0 != range0.second; ++it0) {
//            if (r <= static_cast<Real>(0.5)) {
//              if (i > std::get<0>(it0->second)) {
//                i = std::get<0>(it0->second);
//                j = std::get<1>(it0->second);
//              }
//            }
//            else {
//              if (i < std::get<0>(it0->second)) {
//                i = std::get<0>(it0->second);
//                j = std::get<1>(it0->second);
//              }
//            }
//          }
//          data[i+j*nx] = one;              // Set cell to one
//          data[i+(ny-1-j)*nx] = one;       // Set symmetric cell to one
//          vol += (j!=ny-1-j ? two : one);  // Increase volume
//          map.erase(it);                   // Delete modified cell
//          map2.erase(i+j*nx);              // Delete modified cell
//          for (int k = -1; k < 2; ++k) {
//            for (int l = -1; l < 2; ++l) {
//              if (k==0 || l==0) {
//                int ik = (i+k), jl = (j+l);
//                if ( ik >=0 && ik < nx && jl >= 0 && jl < ny/2 ) { // Check if neighbor indices violate bounds
//                  int index = ik+jl*nx;
//                  auto it2 = map2.find(index);                     // 0-cells that are neighbors of modified cell (i,j)
//                  if (it2 != map2.end()) {                         // If a 0-cell exists
//                    auto range = map.equal_range(it2->second);     // Range in map with key = neighbor count
//                    auto it3 = range.first;                        // Iterator to front of range
//                    bool foundNeighbor = false;
//                    for (it3 = range.first; it3 != range.second; ++it3) {
//                      if (   std::get<0>(it3->second)==ik
//                          && std::get<1>(it3->second)==jl) {      // Find neighbor
//                        foundNeighbor = true;
//                        break;
//                      }
//                    }
//                    if (foundNeighbor) {
//                      int                key = it3->first;        // Increment neighbor counter
//                      std::pair<int,int> val = it3->second;       // Cell (x,y) location
//                      map.erase(it3);                             // Remove {key, val} from map
//                      map.insert({++key,val});                    // Add {++key, val} back to map
//                      it2->second++;                              // Increase count for map2
//                    }
//                  }
//                  else {                                          // If neighbor does not exist in map
//                    if ( data[index] == zero ) {
//                      map.insert({1,{ik,jl}});                    // Add to map
//                      map2.insert({index,1});                     // Add to map2
//                    }
//                  }
//                }
//              }
//            }
//          }
//        }
        //for (auto it = map.rbegin(); it != map.rend(); ++it) {
        //  int i = std::get<0>(it->second);
        //  int j = std::get<1>(it->second);
        //  data[i+j*nx] = one;
        //  vol += one;
        //  if (j != ny-1-j) {
        //    data[i+(ny-1-j)*nx] = one;
        //    vol += one;
        //  }
        //  if (vol > vol0-two) break;
        //}
      }
      if (vol > vol0) {
        // Get nonzero cells with zero neighbors
        std::multimap<int,std::pair<int,int>> map;
        std::map<int,int> map2;
        for (int i = 0; i < nx; ++i) {
          for (int j = 0; j < ny/2; ++j) {
            d = data[i+j*nx];
            cnt = 0;
            if (d==one) {
              for (int k = -1; k < 2; ++k) {
                for (int l = -1; l < 2; ++l) {
                  if (k==0 || l==0) {
                  //if (k!=0 || l!=0) {
                    int ik = (i+k), jl = (j+l);
                    if ( ik >=0 && ik < nx && jl >= 0 && jl < ny/2 ) {
                      nd = data[ik + jl*nx];
                      if (nd == zero) cnt++;
                    }
                  }
                }
              }
              map.insert({cnt,{i,j}});
              map2.insert({i+j*nx,cnt});
            }
          }
        }
        while (vol > vol0) {
          std::map<int,std::pair<int,int>>::iterator it = map.end();
          it--;                            // Grab cell with largest number of neighbors
          int i = std::get<0>(it->second); // x-index of cell
          int j = std::get<1>(it->second); // y-index of cell
          auto range0 = map.equal_range(it->first);
          Real r = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
          for (auto it0 = range0.first; it0 != range0.second; ++it0) {
            if (r <= static_cast<Real>(0.5)) {
              if (i > std::get<0>(it0->second)) {
                i = std::get<0>(it0->second);
                j = std::get<1>(it0->second);
              }
            }
            else {
              if (i < std::get<0>(it0->second)) {
                i = std::get<0>(it0->second);
                j = std::get<1>(it0->second);
              }
            }
          }
          data[i+j*nx] = zero;             // Set cell to zero
          data[i+(ny-1-j)*nx] = zero;      // Set symmetric cell to zero
          vol -= (j!=ny-1-j ? two : one);  // Decrease volume
          map.erase(it);                   // Delete modified cell
          map2.erase(i+j*nx);              // Delete modified cell
          for (int k = -1; k < 2; ++k) {
            for (int l = -1; l < 2; ++l) {
              if (k==0 || l==0) {
                int ik = (i+k), jl = (j+l);
                if ( ik >=0 && ik < nx && jl >= 0 && jl < ny/2 ) { // Check if neighbor indices violate bounds
                  int index = ik+jl*nx;
                  auto it2 = map2.find(index);                     // 0-cells that are neighbors of modified cell (i,j)
                  if (it2 != map2.end()) {                         // If a 0-cell exists
                    auto range = map.equal_range(it2->second);     // Range in map with key = neighbor count
                    auto it3 = range.first;                        // Iterator to front of range
                    bool foundNeighbor = false;
                    for (it3 = range.first; it3 != range.second; ++it3) {
                      if (   std::get<0>(it3->second)==ik
                          && std::get<1>(it3->second)==jl) {      // Find neighbor
                        foundNeighbor = true;
                        break;
                      }
                    }
                    if (foundNeighbor) { 
                      int                key = it3->first;        // Increment neighbor counter
                      std::pair<int,int> val = it3->second;       // Cell (x,y) location
                      map.erase(it3);                             // Remove {key, val} from map
                      map.insert({++key,val});                    // Add {++key, val} back to map
                      it2->second++;                              // Increase count for map2
                    }
                  }
                  else {                                          // If neighbor does not exist in map
                    if ( data[index]==one ) {
                      map.insert({1,{ik,jl}});                    // Add to map
                      map2.insert({index,1});                     // Add to map2
                    }
                  }
                }
              }
            }
          }
        }
        //for (auto it = map.begin(); it != map.end(); ++it) {
        //  int i = std::get<0>(it->second);
        //  int j = std::get<1>(it->second);
        //  data[i+j*nx] = zero;
        //  vol -= one;
        //  if (j != ny-1-j) {
        //    data[i+(ny-1-j)*nx] = zero;
        //    vol -= one;
        //  }
        //  if (vol < vol0) break;
        //}
      }
      if (verbosity_ > 1)
        *outStream_ << "Volume Difference after Modification: " << vol-vol0 << std::endl;
      cnorm = infeasibility(*rndSolution_);
    }
    else {
      Real r(0), sum(0);
      bool oneSet(false);
      while (true) {
        //rndSolution_->set(*solution_);
        setDensity(*rndSolution_,*solution_);
        for (int i = 0; i < nc; ++i) { 
          r = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
          sum = zero;
          oneSet = false;
          for (int t = 0; t < T_; ++t) {
            sum += data[t+i*T_];
            if (r <= sum && !oneSet) {
              data[t+i*T_] = one;
              oneSet = true;
            }
            else {
              data[t+i*T_] = zero;
            }
          }
        }
        cnorm = infeasibility(*rndSolution_);
        cnt++;
        if (cnorm < ctol_) break;
        if (verbosity_ > 1) {
          *outStream_ << "  cnt = " << cnt << "  infeasibility = " << cnorm << std::endl;
        }
      }
    }
    problem0_->getObjective()->update(*rndSolution_,ROL::UpdateType::Temp);
    Real val = problem0_->getObjective()->value(*rndSolution_,tol);
    branching_->foundSolution(new ROL::PEBBL::IntegerSolution<Real>(*rndSolution_,val));
    // Print incumbent
    numIH_++;
    printVector(*rndSolution_,"incumbent");
    printVector(*solution_,"fractional");
    //
    if (verbosity_ > 0) {
      *outStream_ << "MultiMatBranchSub::incumbentHeuristic: " << methodName_ << std::endl;
      *outStream_ << "  Incumbent Value:       " <<   val << std::endl;
      *outStream_ << "  Incumbent Feasibility: " << cnorm << std::endl;
      if (method_ != 1)
        *outStream_ << "  Number of Samples:     " <<   cnt << std::endl;
      outStream_->flags(flags);
    }
  }

  pebbl::branchSub* makeChild(int whichChild = anyChild) override {
    ROL_TEST_FOR_EXCEPTION(whichChild==anyChild, std::logic_error,
      ">>> MultiMatBranchSub::makeChild: whichChild is equal to anyChild!");
    MultiMatBranchSub<Real>* child = new MultiMatBranchSub<Real>(*this);
    child->updateFixed(index_,
      (whichChild==0 ? static_cast<Real>(1) : static_cast<Real>(0)));
    return child;
  }

}; // class MultiMatBranchSub

#endif
