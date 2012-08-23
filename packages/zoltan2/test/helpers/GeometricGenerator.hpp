// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_FilteredIterator.hpp>
#include <Teuchos_ParameterEntry.hpp>
#include <iostream>
#include <ctime>
#include <limits>
#include <climits>
#include <string>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <Tpetra_MultiVector_decl.hpp>
#include <Zoltan2_XpetraMultiVectorInput.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Teuchos_ArrayViewDecl.hpp>

enum shape {SQUARE, RECTANGLE, CIRCLE, CUBE, RECTANGULAR_PRISM, SPHERE};
const std::string shapes[] = {"SQUARE", "RECTANGLE", "CIRCLE", "CUBE", "RECTANGULAR_PRISM", "SPHERE"};
#define SHAPE_COUNT 6

enum distribution {normal, uniform};
const std::string distribution[] = {"distribution", "uniform"};
#define DISTRIBUTION_COUNT 2

#define HOLE_ALLOC_STEP  10
#define MAX_WEIGHT_DIM  10
#define INVALID(STR) "Invalid argument at " + STR
#define INVALIDSHAPE(STR, DIM) "Invalid shape name " + STR + " for " + 	DIM + ".\nValid shapes are \"SQUARE\", \"RECTANGLE\", \"CIRCLE\" for 2D, and \"CUBE\", \"RECTANGULAR_PRISM\", \"SPHERE\" for 3D"

#define INVALID_SHAPE_ARG(SHAPE, REQUIRED) "Invalid argument count for shape " + SHAPE + ". Requires " + REQUIRED + " argument(s)."
#define MAX_ITER_ALLOWED 500

const std::string weight_distribution = "WeightDistribution";
const std::string weight_distribution_string = weight_distribution + "-";
template <typename T>
struct CoordinatePoint {
  T x;
  T y;
  T z;
  /*
	bool isInArea(int dim, T *minCoords, T *maxCoords){
		dim = min(3, dim);
		for (int i = 0; i < dim; ++i){
			bool maybe = true;
			switch(i){
			case 0:
				maybe = x < maxCoords[i] && x > maxCoords[i];
				break;
			case 1:
				maybe = y < maxCoords[i] && y > maxCoords[i];
				break;
			case 2:
				maybe = z < maxCoords[i] && z > maxCoords[i];
				break;
			}
			if (!maybe){
				return false;
			}
		}
		return true;
	}
   */
  CoordinatePoint(){
    x = 0;y=0;z=0;
  }
};

template <typename T>
class Hole{
public:
  CoordinatePoint<T> center;
  T edge1, edge2, edge3;
  Hole(CoordinatePoint<T> center_, T edge1_, T edge2_, T edge3_){
    this->center.x = center_.x;
    this->center.y = center_.y;
    this->center.z = center_.z;
    this->edge1 = edge1_;
    this->edge2 = edge2_;
    this->edge3 = edge3_;
  }
  virtual bool isInArea(CoordinatePoint<T> dot) = 0;
};

template <typename T>
class SquareHole:public Hole<T>{
public:
  SquareHole(CoordinatePoint<T> center_ , T edge_): Hole<T>(center_, edge_, 0 , 0){}

  virtual bool isInArea(CoordinatePoint<T> dot){
    return fabs(dot.x - this->center.x) < this->edge1 / 2 && fabs(dot.y - this->center.y) < this->edge1 / 2;
  }
};

template <typename T>
class RectangleHole:public Hole<T>{
public:
  RectangleHole(CoordinatePoint<T> center_  , T edge_x_, T edge_y_): Hole<T>(center_, edge_x_,  edge_y_, 0){}
  virtual bool isInArea(CoordinatePoint<T> dot){
    return fabs(dot.x - this->center.x) < this->edge1 / 2 && fabs(dot.y - this->center.y) < this->edge2 / 2;
  }
};

template <typename T>
class CircleHole:public Hole<T>{
public:
  CircleHole(CoordinatePoint<T> center_ , T edge_): Hole<T>(center_, edge_, 0 , 0){}
  virtual bool isInArea(CoordinatePoint<T> dot){
    return (dot.x - this->center.x)*(dot.x - this->center.x) + (dot.y - this->center.y) * (dot.y - this->center.y) < this->edge1 * this->edge1;
  }
};

template <typename T>
class CubeHole:public Hole<T>{
public:
  CubeHole(CoordinatePoint<T> center_ , T edge_): Hole<T>(center_, edge_, 0 , 0){}
  virtual bool isInArea(CoordinatePoint<T> dot){
    return fabs(dot.x - this->center.x) < this->edge1 / 2 && fabs(dot.y - this->center.y) < this->edge1 / 2 && fabs(dot.z - this->center.z) < this->edge1 / 2;
  }
};

template <typename T>
class RectangularPrismHole:public Hole<T>{
public:
  RectangularPrismHole(CoordinatePoint<T> center_  , T edge_x_, T edge_y_, T edge_z_): Hole<T>(center_, edge_x_,  edge_y_, edge_z_){}
  virtual bool isInArea(CoordinatePoint<T> dot){
    return fabs(dot.x - this->center.x) < this->edge1 / 2 && fabs(dot.y - this->center.y) < this->edge2 / 2 && fabs(dot.z - this->center.z) < this->edge3 / 2;
  }
};

template <typename T>
class SphereHole:public Hole<T>{
public:
  SphereHole(CoordinatePoint<T> center_ , T edge_): Hole<T>(center_, edge_, 0 , 0){}
  virtual bool isInArea(CoordinatePoint<T> dot){
    return  fabs((dot.x - this->center.x) * (dot.x - this->center.x) * (dot.x - this->center.x)) +
        fabs((dot.y - this->center.y) * (dot.y - this->center.y) * (dot.y - this->center.y)) +
        fabs((dot.z - this->center.z) * (dot.z - this->center.z) * (dot.z - this->center.z))
        <
        this->edge1 * this->edge1 * this->edge1;
  }
};

template <typename T, typename weighttype>
class WeightDistribution{
public:
  virtual weighttype getWeight(CoordinatePoint<T> P)=0;
  virtual weighttype get1DWeight(T x)=0;
  virtual weighttype get2DWeight(T x, T y)=0;
  virtual weighttype get3DWeight(T x, T y, T z)=0;
  WeightDistribution(){};
  virtual ~WeightDistribution(){};
};

/**
 * Expression is a generic following method.
 * a1 (x - x1)^b1 + a2 (y - y1)^b2 + a3 (z - z1)^b3 + c = expression_result
 * if step values are given expression result is applied to a step function as following.
 * expression_result	< step1		value1
 * 						< step2		value2
 * 						< step3		value3
 * 						< step4		value4
 *
 * 						Default values,
 * 							c=1
 * 							a1=a2=a3=0
 * 							x'=y'=z'=0
 * 							b1=b2=b3=0
 * 							steps = NULL
 * 							vals = NULL
 *
 */
template <typename T, typename weighttype>
class SteppedEquation:public WeightDistribution<T, weighttype>{
  T a1,a2,a3;
  T b1,b2,b3;
  T c;
  T x1,y1,z1;

  T *steps;
  weighttype *values;
  int stepCount;
public:
  SteppedEquation(T a1_, T a2_, T a3_, T b1_, T b2_, T b3_, T c_, T x1_, T y1_, T z1_, T *steps_, T *values_, int stepCount_):WeightDistribution<T,weighttype>(){
    this->a1 = a1_;
    this->a2 = a2_;
    this->a3 = a3_;
    this->b1 = b1_;
    this->b2 = b2_;
    this->b3 = b3_;
    this->c = c_;
    this->x1 = x1_;
    this->y1 = y1_;
    this->z1 = z1_;


    this->stepCount = stepCount_;
    if(this->stepCount > 0){
      this->steps = new T[this->stepCount];
      this->values = new T[this->stepCount + 1];

      for (int i = 0; i < this->stepCount; ++i){
        this->steps[i] = steps_[i];
        this->values[i] = values_[i];
      }
      this->values[this->stepCount] = values_[this->stepCount];
    }

  }

  virtual ~SteppedEquation(){
    if(this->stepCount > 0){
      delete [] this->steps;
      delete [] this->values;
    }
  }


  virtual weighttype get1DWeight(T x){
    T expressionRes = this->a1 * pow( (x - this->x1), b1) + c;
    if(this->stepCount > 0){
      for (int i = 0; i < this->stepCount; ++i){
        if (expressionRes < this->steps[i]) return this->values[i];
      }
      return this->values[this->stepCount];
    }
    else {
      return weighttype(expressionRes);
    }
  }

  virtual weighttype get2DWeight(T x, T y){
    T expressionRes = this->a1 * pow( (x - this->x1), b1) + this->a2 * pow( (y - this->y1), b2) + c;
    if(this->stepCount > 0){
      for (int i = 0; i < this->stepCount; ++i){
        if (expressionRes < this->steps[i]) return this->values[i];
      }
      return this->values[this->stepCount];
    }
    else {
      return weighttype(expressionRes);
    }
  }

  void print (T x, T y, T z){
    cout << this->a1 << "*" <<  "math.pow( (" << x  << "- " <<  this->x1 << " ), " << b1 <<")" << "+" <<  this->a2<< "*"<<  "math.pow( (" << y << "-" <<  this->y1 << "), " <<
        b2 << " ) + " << this->a3 << " * math.pow( (" << z << "-" <<  this->z1 << "), " << b3 << ")+ " << c << " == " <<
        this->a1 * pow( (x - this->x1), b1) + this->a2 * pow( (y - this->y1), b2) + this->a3 * pow( (z - this->z1), b3) + c << endl;

  }

  virtual weighttype get3DWeight(T x, T y, T z){
    T expressionRes = this->a1 * pow( (x - this->x1), b1) + this->a2 * pow( (y - this->y1), b2) + this->a3 * pow( (z - this->z1), b3) + this->c;

    //this->print(x,y,z);
    if(this->stepCount > 0){
      for (int i = 0; i < this->stepCount; ++i){
        if (expressionRes < this->steps[i]) {
          //cout << "0exp:" << expressionRes << " step:" << steps[i] << " value:" << values[i] << endl;
          return this->values[i];
        }
      }

      //cout << "1exp:" << expressionRes << " step:" << steps[stepCount] << " value:" << values[stepCount] << endl;
      return this->values[this->stepCount];
    }
    else {
      return weighttype(expressionRes);
    }
  }

  virtual weighttype getWeight(CoordinatePoint<T> p){
    T expressionRes = this->a1 * pow( (p.x - this->x1), b1) + this->a2 * pow( (p.y - this->y1), b2) + this->a3 * pow( (p.z - this->z1), b3);
    if(this->stepCount > 0){
      for (int i = 0; i < this->stepCount; ++i){
        if (expressionRes < this->steps[i]) return this->values[i];
      }
      return this->values[this->stepCount];
    }
    else {
      return weighttype(expressionRes);
    }
  }
};

template <typename T, typename lno_t, typename gno_t>
class CoordinateDistribution{
public:
  gno_t numPoints;
  int dimension;
  lno_t requested;
  gno_t assignedPrevious;
  int worldSize;
  virtual ~CoordinateDistribution(){}

  CoordinateDistribution(gno_t np_, int dim, int worldSize):numPoints(np_), dimension(dim), requested(0), assignedPrevious(0), worldSize(worldSize){}
  virtual CoordinatePoint<T> getPoint(gno_t point_index, unsigned int &state) = 0;
  virtual T getXCenter() = 0;
  virtual T getXRadius() =0;

  void GetPoints(lno_t requestedPointcount, CoordinatePoint<T> *points /*preallocated sized numPoints*/,
      Hole<T> **holes, lno_t holeCount,
      float *sharedRatios_, int myRank){

    for (int i = 0; i < myRank; ++i){
      //cout << "me:" << myRank << " i:" << i << " s:" << sharedRatios_[i]<< endl;
      this->assignedPrevious += lno_t(sharedRatios_[i] * this->numPoints);
      if (i < this->numPoints % this->worldSize){
        this->assignedPrevious += 1;
      }
    }

    this->requested = requestedPointcount;

    unsigned int slice =  UINT_MAX/(this->worldSize);
    unsigned int stateBegin = myRank * slice;

//#ifdef HAVE_ZOLTAN2_OMP
//#pragma omp parallel for
//#endif
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
#endif
    {
      int me = 0;
      int tsize = 1;
#ifdef HAVE_ZOLTAN2_OMP
      me = omp_get_thread_num();
      tsize = omp_get_num_threads();
#endif
      unsigned int state = stateBegin + me * slice/(tsize);

#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
      for(lno_t cnt = 0; cnt < requestedPointcount; ++cnt){
        lno_t iteration = 0;
        while(1){
          if(++iteration > MAX_ITER_ALLOWED) {
            throw "Max number of Iteration is reached for point creation. Check the area criteria or hole coordinates.";
          }
          CoordinatePoint <T> p = this->getPoint( this->assignedPrevious + cnt, &state);

          bool isInHole = false;
          for(lno_t i = 0; i < holeCount; ++i){
            if(holes[i][0].isInArea(p)){
              isInHole = true;
              break;
            }
          }
          if(isInHole) continue;
          points[cnt].x = p.x;

          points[cnt].y = p.y;
          points[cnt].z = p.z;
          break;
        }
      }
    }
//#pragma omp parallel
      /*
    {

      lno_t cnt = 0;
      lno_t iteration = 0;
      while (cnt < requestedPointcount){
        if(++iteration > MAX_ITER_ALLOWED) {
          throw "Max number of Iteration is reached for point creation. Check the area criteria or hole coordinates.";
        }
        CoordinatePoint <T> p = this->getPoint();

        bool isInHole = false;
        for(lno_t i = 0; i < holeCount; ++i){
          if(holes[i][0].isInArea(p)){
            isInHole = true;
            break;
          }
        }
        if(isInHole) continue;
        iteration = 0;
        points[cnt].x = p.x;
        points[cnt].y = p.y;
        points[cnt].z = p.z;
        ++cnt;
      }
    }
    */
  }

  void GetPoints(lno_t requestedPointcount, T **coords/*preallocated sized numPoints*/, lno_t tindex,
      Hole<T> **holes, lno_t holeCount,
      float *sharedRatios_, int myRank){

    for (int i = 0; i < myRank; ++i){
      //cout << "me:" << myRank << " i:" << i << " s:" << sharedRatios_[i]<< endl;
      this->assignedPrevious += lno_t(sharedRatios_[i] * this->numPoints);
      if (i < this->numPoints % this->worldSize){
        this->assignedPrevious += 1;
      }
    }

    this->requested = requestedPointcount;

    unsigned int slice =  UINT_MAX/(this->worldSize);
    unsigned int stateBegin = myRank * slice;
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel
#endif
    {
      int me = 0;
      int tsize = 1;
#ifdef HAVE_ZOLTAN2_OMP
      me = omp_get_thread_num();
      tsize = omp_get_num_threads();
#endif
      unsigned int state = stateBegin + me * (slice/(tsize));
      /*
#pragma omp critical
      {

        cout << "myRank:" << me << " stateBeg:" << stateBegin << " tsize:" << tsize << " state:" << state <<  " slice: " << slice / tsize <<  endl;
      }
      */
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp for
#endif
      for(lno_t cnt = 0; cnt < requestedPointcount; ++cnt){
        lno_t iteration = 0;
        while(1){
          if(++iteration > MAX_ITER_ALLOWED) {
            throw "Max number of Iteration is reached for point creation. Check the area criteria or hole coordinates.";
          }
          CoordinatePoint <T> p = this->getPoint( this->assignedPrevious + cnt, state);

          bool isInHole = false;
          for(lno_t i = 0; i < holeCount; ++i){
            if(holes[i][0].isInArea(p)){
              isInHole = true;
              break;
            }
          }
          if(isInHole) continue;
          coords[0][cnt + tindex] = p.x;
          if(this->dimension > 1){
            coords[1][cnt + tindex] = p.y;
            if(this->dimension > 2){
              coords[2][cnt + tindex] = p.z;
            }
          }
          break;
        }
      }
    }
  }
};

template <typename T, typename lno_t, typename gno_t>
class CoordinateNormalDistribution:public CoordinateDistribution<T,lno_t,gno_t>{
public:
  CoordinatePoint<T> center;
  T standartDevx;
  T standartDevy;
  T standartDevz;


  virtual T getXCenter(){
    return center.x;
  }
  virtual T getXRadius(){
    return standartDevx;
  }

  CoordinateNormalDistribution(gno_t np_, int dim, CoordinatePoint<T> center_ , T sd_x, T sd_y, T sd_z, int worldSize): CoordinateDistribution<T,lno_t,gno_t>(np_,dim,worldSize),
      standartDevx(sd_x), standartDevy(sd_y), standartDevz(sd_z){
    this->center.x = center_.x;
    this->center.y = center_.y;
    this->center.z = center_.z;
  }

  virtual CoordinatePoint<T> getPoint(gno_t pindex, unsigned int &state){

    pindex = 0; // not used in normal distribution.
    CoordinatePoint <T> p;

    for(int i = 0; i < this->dimension; ++i){
      switch(i){
      case 0:
        p.x = normalDist(this->center.x, this->standartDevx, state);
        break;
      case 1:
        p.y = normalDist(this->center.y, this->standartDevy, state);
        break;
      case 2:
        p.z = normalDist(this->center.z, this->standartDevz, state);
        break;
      default:
        throw "unsupported dimension";
      }
    }
    return p;
  }

  virtual ~CoordinateNormalDistribution(){};
private:
  T normalDist(T center, T sd, unsigned int &state) {
    static bool derived=false;
    static T storedDerivation;
    T polarsqrt, normalsquared, normal1, normal2;
    if (!derived) {
      do {
        normal1=2.0*( T(rand_r(&state))/T(RAND_MAX) ) - 1.0;
        normal2=2.0*( T(rand_r(&state))/T(RAND_MAX) ) - 1.0;
        normalsquared=normal1*normal1+normal2*normal2;
      } while ( normalsquared>=1.0 || normalsquared == 0.0);

      polarsqrt=sqrt(-2.0*log(normalsquared)/normalsquared);
      storedDerivation=normal1*polarsqrt;
      derived=true;
      return normal2*polarsqrt*sd + center;
    }
    else {
      derived=false;
      return storedDerivation*sd + center;
    }
  }
};

template <typename T, typename lno_t, typename gno_t>
class CoordinateUniformDistribution:public CoordinateDistribution<T,lno_t, gno_t>{
public:
  T leftMostx;
  T rightMostx;
  T leftMosty;
  T rightMosty;
  T leftMostz;
  T rightMostz;


  virtual T getXCenter(){
    return (rightMostx - leftMostx)/2  + leftMostx;
  }
  virtual T getXRadius(){
    return (rightMostx - leftMostx)/2;
  }


  CoordinateUniformDistribution(gno_t np_, int dim, T l_x, T r_x, T l_y, T r_y, T l_z, T r_z, int worldSize ): CoordinateDistribution<T,lno_t,gno_t>(np_,dim,worldSize),
      leftMostx(l_x), rightMostx(r_x), leftMosty(l_y), rightMosty(r_y), leftMostz(l_z), rightMostz(r_z){}

  virtual ~CoordinateUniformDistribution(){};
  virtual CoordinatePoint<T> getPoint(gno_t pindex, unsigned int &state){


    pindex = 0; //not used in uniform dist.
    CoordinatePoint <T> p;
    for(int i = 0; i < this->dimension; ++i){
      switch(i){
      case 0:
        p.x = uniformDist(this->leftMostx, this->rightMostx, state);
        break;
      case 1:
        p.y = uniformDist(this->leftMosty, this->rightMosty, state);
        break;
      case 2:
        p.z = uniformDist(this->leftMostz, this->rightMostz, state);
        break;
      default:
        throw "unsupported dimension";
      }
    }
    return p;
  }

private:

  T uniformDist(T a, T b, unsigned int &state)
  {
    return (b-a)*(rand_r(&state) / double(RAND_MAX)) + a;
  }
};

template <typename T, typename lno_t, typename gno_t>
class CoordinateGridDistribution:public CoordinateDistribution<T,lno_t,gno_t>{
public:
  T leftMostx;
  T rightMostx;
  T leftMosty;
  T rightMosty;
  T leftMostz;
  T rightMostz;
  gno_t along_X, along_Y, along_Z;
  //T currentX, currentY, currentZ;
  T processCnt;
  int myRank;
  T xstep, ystep, zstep;
  gno_t xshift, yshift, zshift;

  virtual T getXCenter(){
    return (rightMostx - leftMostx)/2  + leftMostx;
  }
  virtual T getXRadius(){
    return (rightMostx - leftMostx)/2;
  }


  CoordinateGridDistribution(gno_t alongX, gno_t alongY, gno_t alongZ, int dim, T l_x, T r_x, T l_y, T r_y, T l_z, T r_z , int myRank_, int worldSize): CoordinateDistribution<T,lno_t,gno_t>(alongX * alongY * alongZ,dim,worldSize),
      leftMostx(l_x), rightMostx(r_x), leftMosty(l_y), rightMosty(r_y), leftMostz(l_z), rightMostz(r_z), myRank(myRank_){
    //currentX = leftMostx, currentY = leftMosty, currentZ = leftMostz;
    this->processCnt = 0;
    this->along_X = alongX; this->along_Y = alongY; this->along_Z = alongZ;

    if(this->along_X > 1)
      this->xstep = (rightMostx - leftMostx) / (alongX - 1);
    else
      this->xstep = 1;
    if(this->along_Y > 1)
      this->ystep = (rightMosty - leftMosty) / (alongY - 1);
    else
      this->ystep = 1;
    if(this->along_Z > 1)
      this->zstep = (rightMostz - leftMostz) / (alongZ - 1);
    else
      this->zstep = 1;
    xshift = 0; yshift = 0; zshift = 0;
  }

  virtual ~CoordinateGridDistribution(){};
  virtual CoordinatePoint<T> getPoint(gno_t pindex, unsigned int &state){
    //lno_t before = processCnt + this->assignedPrevious;
    //cout << "before:" << processCnt << " " << this->assignedPrevious << endl;
    //lno_t xshift = 0, yshift = 0, zshift = 0;

    //lno_t tmp = before % (this->along_X * this->along_Y);
    //xshift  = tmp % this->along_X;
    //yshift = tmp / this->along_X;
    //zshift = before / (this->along_X * this->along_Y);

    state = 0; //not used here
    this->zshift = pindex / (along_X * along_Y);
    this->yshift = (pindex % (along_X * along_Y)) / along_X;
    this->xshift = (pindex % (along_X * along_Y)) % along_X;

    //this->xshift = pindex / (along_Z * along_Y);
   // this->zshift = (pindex % (along_Z * along_Y)) / along_Y;
    //this->yshift = (pindex % (along_Z * along_Y)) % along_Y;

    CoordinatePoint <T> p;
    p.x = xshift * this->xstep + leftMostx;
    p.y = yshift * this->ystep + leftMosty;
    p.z = zshift * this->zstep + leftMostz;
/*
    ++xshift;
    if(xshift == this->along_X){
      ++yshift;
      xshift = 0;
      if(yshift == this->along_Y){
        ++zshift;
        yshift = 0;
      }
    }
*/
    /*
    if(this->processCnt == 0){
      this->xshift = this->assignedPrevious / (along_Z * along_Y);
      //this->yshift = (this->assignedPrevious % (along_X * along_Y)) / along_X;
      this->zshift = (this->assignedPrevious % (along_Z * along_Y)) / along_Y;
      //this->xshift = (this->assignedPrevious % (along_X * along_Y)) % along_X;
      this->yshift = (this->assignedPrevious % (along_Z * along_Y)) % along_Y;
      ++this->processCnt;
    }

    CoordinatePoint <T> p;
    p.x = xshift * this->xstep + leftMostx;
    p.y = yshift * this->ystep + leftMosty;
    p.z = zshift * this->zstep + leftMostz;

    ++yshift;
    if(yshift == this->along_Y){
      ++zshift;
      yshift = 0;
      if(zshift == this->along_Z){
        ++xshift;
        zshift = 0;
      }

    }
    */
    /*
    if(this->requested - 1 > this->processCnt){
      this->processCnt++;
    } else {
      cout << "req:" << this->requested << " pr:" <<this->processCnt << endl;
      cout << "p:" << p.x << " " << p.y << " " << p.z <<  endl ;
      cout << "s:" << xshift << " " << yshift << " " << zshift <<  endl ;
      cout << "st:" << this->xstep << " " << this->ystep << " " << this->zstep <<  endl ;
    }
    */
    return p;
  }

private:

};

template <typename T, typename lno_t, typename gno_t, typename node_t>
class GeometricGenerator {
private:
  Hole<T> **holes; //to represent if there is any hole in the input
  int holeCount;
  int coordinate_dimension;  //dimension of the geometry
  gno_t numGlobalCoords;	//global number of coordinates requested to be created.
  lno_t numLocalCoords;
  float *loadDistributions; //sized as the number of processors, the load of each processor.
  bool loadDistSet;
  bool distinctCoordSet;
  CoordinateDistribution<T, lno_t,gno_t> **coordinateDistributions;
  int distributionCount;
  //CoordinatePoint<T> *points;
  T **coords;
  T **wghts;
  WeightDistribution<T,T> **wd;
  int weight_dimension;  //dimension of the geometry

  //RCP< Tpetra::MultiVector<T, lno_t, gno_t, node_t> >tmVector;
  //RCP< Tpetra::MultiVector<T, lno_t, gno_t, node_t> >tmwVector;
  int worldSize;
  int myRank;
  T minx;
  T maxx;
  T miny;
  T maxy;
  T minz;
  T maxz;
  std::string outfile;



  template <typename tt>
  tt getParamVal( const Teuchos::ParameterEntry& pe, const std::string &paramname){
    tt returnVal;
    try {
      returnVal = Teuchos::getValue<tt>(pe);
    }
    catch (...){
      throw INVALID(paramname);
    }
    return returnVal;
  }

  int countChar (std::string inStr, char countChar){
    int cnt = 0;
    for (unsigned int i = 0; i < inStr.size(); ++i){
      if (inStr[i] == countChar) {
        cnt++;
      }
    }
    return cnt;
  }

  template <typename tt>
  std::string toString(tt obj){
    std::stringstream ss (std::stringstream::in |std::stringstream::out);
    ss << obj;
    std::string tmp = "";
    ss >> tmp;
    return tmp;
  }

  template <typename tt>
  tt fromString(std::string obj){
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    ss << obj;
    tt tmp;
    ss >> tmp;
    if (ss.fail()){
      throw "Cannot convert string " + obj;
    }
    return tmp;
  }


  void splitString(std::string inStr, char splitChar, std::string *outSplittedStr){
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    ss << inStr;
    int cnt = 0;
    while (!ss.eof()){
      std::string tmp = "";
      std::getline(ss, tmp ,splitChar);
      outSplittedStr[cnt++] = tmp;
    }
  }

  /*
	void getDistinctCoordinateDescription(std::string distinctDescription){

		this->distinctCoordinates = new bool[this->dimension];
		if(distinctDescription != ""){
			int argCnt = this->countChar(distinctDescription, ',') + 1;
			if(argCnt != this->dimension) {
				throw "Invalid parameter count for distinct_coordinates. Size should be equal to dimension.";
			}

			std::string *splittedStr = new std::string[argCnt];
			splitString(holeDescription, ',', splittedStr);

			for(int i = 0; i < argCnt; ++i){
				if(splittedStr[i] == "T"){
					distinctCoordinates[i] = true;
				} else if(splittedStr[i] == "F"){
					distinctCoordinates[i] = false;
				} else {
					throw "Invalid parameter " + splittedStr[i] + " for distinct_coordinates.";
				}
			}
			delete []splittedStr;
		}
		else {
			for(int i = 0; i < this->dimension; ++i){
				distinctCoordinates[i] = false;
			}
		}
	}
   */


  void getCoordinateDistributions(std::string coordinate_distributions){
    if(coordinate_distributions == ""){
      throw "At least one distribution function must be provided to coordinate generator.";
    }

    int argCnt = this->countChar(coordinate_distributions, ',') + 1;
    std::string *splittedStr = new std::string[argCnt];
    splitString(coordinate_distributions, ',', splittedStr);
    coordinateDistributions = (CoordinateDistribution<T, lno_t,gno_t> **) malloc(sizeof (CoordinateDistribution<T, lno_t,gno_t> *) * 1);
    for(int i = 0; i < argCnt; ){
      coordinateDistributions = (CoordinateDistribution<T, lno_t,gno_t> **)realloc((void *)coordinateDistributions, (this->distributionCount + 1)* sizeof(CoordinateDistribution<T, lno_t,gno_t> *));

      std::string distName = splittedStr[i++];
      gno_t np_ = 0;
      if(distName == "NORMAL"){
        int reqArg = 5;
        if (this->coordinate_dimension == 3){
          reqArg = 7;
        }
        if(i + reqArg > argCnt) {
          std::string tmp = toString<int>(reqArg);
          throw INVALID_SHAPE_ARG(distName, tmp);
        }
        np_ = fromString<gno_t>(splittedStr[i++]);
        CoordinatePoint<T> pp;

        pp.x = fromString<T>(splittedStr[i++]);
        pp.y = fromString<T>(splittedStr[i++]);
        pp.z = 0;
        if(this->coordinate_dimension == 3){
          pp.z = fromString<T>(splittedStr[i++]);
        }

        T sd_x = fromString<T>(splittedStr[i++]);
        T sd_y = fromString<T>(splittedStr[i++]);
        T sd_z = 0;
        if(this->coordinate_dimension == 3){
          sd_z = fromString<T>(splittedStr[i++]);
        }
        this->coordinateDistributions[this->distributionCount++] = new CoordinateNormalDistribution<T, lno_t,gno_t>(np_, this->coordinate_dimension, pp , sd_x, sd_y, sd_z, worldSize );

      } else if(distName == "UNIFORM" ){
        int reqArg = 5;
        if (this->coordinate_dimension == 3){
          reqArg = 7;
        }
        if(i + reqArg > argCnt) {
          std::string tmp = toString<int>(reqArg);
          throw INVALID_SHAPE_ARG(distName, tmp);
        }
        np_ = fromString<gno_t>(splittedStr[i++]);
        T l_x = fromString<T>(splittedStr[i++]);
        T r_x = fromString<T>(splittedStr[i++]);
        T l_y = fromString<T>(splittedStr[i++]);
        T r_y = fromString<T>(splittedStr[i++]);

        T l_z = 0, r_z = 0;

        if(this->coordinate_dimension == 3){
          l_z = fromString<T>(splittedStr[i++]);
          r_z = fromString<T>(splittedStr[i++]);
        }

        this->coordinateDistributions[this->distributionCount++] = new CoordinateUniformDistribution<T, lno_t,gno_t>( np_,  this->coordinate_dimension, l_x, r_x, l_y, r_y, l_z, r_z, worldSize );
      } else if (distName == "GRID"){
        int reqArg = 6;
        if(this->coordinate_dimension == 3){
          reqArg = 9;
        }
        if(i + reqArg > argCnt) {
          std::string tmp = toString<int>(reqArg);
          throw INVALID_SHAPE_ARG(distName, tmp);
        }

        gno_t np_x = fromString<gno_t>(splittedStr[i++]);
        gno_t np_y = fromString<gno_t>(splittedStr[i++]);
        gno_t np_z = 1;


        if(this->coordinate_dimension == 3){
          np_z = fromString<gno_t>(splittedStr[i++]);
        }

        np_ = np_x * np_y * np_z;
        T l_x = fromString<T>(splittedStr[i++]);
        T r_x = fromString<T>(splittedStr[i++]);
        T l_y = fromString<T>(splittedStr[i++]);
        T r_y = fromString<T>(splittedStr[i++]);

        T l_z = 0, r_z = 0;

        if(this->coordinate_dimension == 3){
          l_z = fromString<T>(splittedStr[i++]);
          r_z = fromString<T>(splittedStr[i++]);
        }

        if(np_x < 1 || np_z < 1 || np_y < 1 ){
          throw "Provide at least 1 point along each dimension for grid test.";
        }
        //cout << "ly:" << l_y << " ry:" << r_y << endl;
        this->coordinateDistributions[this->distributionCount++] = new CoordinateGridDistribution<T, lno_t,gno_t>
        (np_x, np_y,np_z, this->coordinate_dimension, l_x, r_x,l_y, r_y, l_z, r_z , this->myRank, worldSize);

      }
      else {
        std::string tmp = toString<int>(this->coordinate_dimension);
        throw INVALIDSHAPE(distName, tmp);
      }
      this->numGlobalCoords += (gno_t) np_;
    }
    delete []splittedStr;
  }

  void getProcLoadDistributions(std::string proc_load_distributions){

    this->loadDistributions = new float[this->worldSize];
    if(proc_load_distributions == ""){
      float r = 1.0 / this->worldSize;
      for(int i = 0; i < this->worldSize; ++i){
        this->loadDistributions[i] = r;
      }
    } else{


      int argCnt = this->countChar(proc_load_distributions, ',') + 1;
      if(argCnt != this->worldSize) {
        throw "Invalid parameter count load distributions. Given " + toString<int>(argCnt) + " processor size is " + toString<int>(worldSize);
      }
      std::string *splittedStr = new std::string[argCnt];
      splitString(proc_load_distributions, ',', splittedStr);
      for(int i = 0; i < argCnt; ++i){
        this->loadDistributions[i] = fromString<float>(splittedStr[i]);
      }
      delete []splittedStr;


      float sum = 0;
      for(int i = 0; i < this->worldSize; ++i){
        sum += this->loadDistributions[i];
      }
      if (fabs(sum - 1.0) > 10*std::numeric_limits<float>::epsilon()){
        throw "Processor load ratios do not sum to 1.0.";
      }
    }

  }

  void getHoles(std::string holeDescription){
    if(holeDescription == ""){
      return;
    }
    this->holes =  (Hole<T> **) malloc(sizeof (Hole <T>*));
    int argCnt = this->countChar(holeDescription, ',') + 1;
    std::string *splittedStr = new std::string[argCnt];
    splitString(holeDescription, ',', splittedStr);

    for(int i = 0; i < argCnt; ){
      holes = (Hole<T> **)realloc((void *)holes, (this->holeCount + 1)* sizeof(Hole<T> *));

      std::string shapeName = splittedStr[i++];
      if(shapeName == "SQUARE" && this->coordinate_dimension == 2){
        if(i + 3 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "3");
        }
        CoordinatePoint<T> pp;
        pp.x = fromString<T>(splittedStr[i++]);
        pp.y = fromString<T>(splittedStr[i++]);
        T edge = fromString<T>(splittedStr[i++]);
        this->holes[this->holeCount++] = new SquareHole<T>(pp, edge);
      } else if(shapeName == "RECTANGLE" && this->coordinate_dimension == 2){
        if(i + 4 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "4");
        }
        CoordinatePoint<T> pp;
        pp.x = fromString<T>(splittedStr[i++]);
        pp.y = fromString<T>(splittedStr[i++]);
        T edgex = fromString<T>(splittedStr[i++]);
        T edgey = fromString<T>(splittedStr[i++]);

        this->holes[this->holeCount++] = new RectangleHole<T>(pp, edgex,edgey);
      } else if(shapeName == "CIRCLE" && this->coordinate_dimension == 2){
        if(i + 3 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "3");
        }
        CoordinatePoint<T> pp;
        pp.x = fromString<T>(splittedStr[i++]);
        pp.y = fromString<T>(splittedStr[i++]);
        T r = fromString<T>(splittedStr[i++]);
        this->holes[this->holeCount++] = new CircleHole<T>(pp, r);
      }  else if(shapeName == "CUBE" && this->coordinate_dimension == 3){
        if(i + 4 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "4");
        }
        CoordinatePoint<T> pp;
        pp.x = fromString<T>(splittedStr[i++]);
        pp.y = fromString<T>(splittedStr[i++]);
        pp.z = fromString<T>(splittedStr[i++]);
        T edge = fromString<T>(splittedStr[i++]);
        this->holes[this->holeCount++] = new CubeHole<T>(pp, edge);
      }  else if(shapeName == "RECTANGULAR_PRISM" && this->coordinate_dimension == 3){
        if(i + 6 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "6");
        }
        CoordinatePoint<T> pp;
        pp.x = fromString<T>(splittedStr[i++]);
        pp.y = fromString<T>(splittedStr[i++]);
        pp.z = fromString<T>(splittedStr[i++]);
        T edgex = fromString<T>(splittedStr[i++]);
        T edgey = fromString<T>(splittedStr[i++]);
        T edgez = fromString<T>(splittedStr[i++]);
        this->holes[this->holeCount++] = new RectangularPrismHole<T>(pp, edgex, edgey, edgez);

      }  else if(shapeName == "SPHERE" && this->coordinate_dimension == 3){
        if(i + 4 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "4");
        }
        CoordinatePoint<T> pp;
        pp.x = fromString<T>(splittedStr[i++]);
        pp.y = fromString<T>(splittedStr[i++]);
        pp.z = fromString<T>(splittedStr[i++]);
        T r = fromString<T>(splittedStr[i++]);
        this->holes[this->holeCount++] = new SphereHole<T>(pp, r);
      }  else {
        std::string tmp = toString<int>(this->coordinate_dimension);
        throw INVALIDSHAPE(shapeName, tmp);
      }
    }
    delete [] splittedStr;
  }

  void getWeightDistribution(std::string *weight_distribution_arr, int wdimension){
    int wcount = 0;

    this->wd = new WeightDistribution<T,T> *[wdimension];
    for(int ii = 0; ii < MAX_WEIGHT_DIM; ++ii){
      std::string weight_distribution = weight_distribution_arr[ii];
      if(weight_distribution == "") continue;
      if(wcount == wdimension) {
        throw "Weight Dimension is provided as " + toString<int>(wdimension) + ". More weight distribution is provided.";
      }

      int count = this->countChar(weight_distribution, ' ');
      std::string *splittedStr = new string[count + 1];
      this->splitString(weight_distribution, ' ', splittedStr);
      //cout << count << endl;
      T c=1;
      T a1=0,a2=0,a3=0;
      T x1=0,y1=0,z1=0;
      T b1=1,b2=1,b3=1;
      T *steps = NULL;
      T *values= NULL;
      int stepCount = 0;
      int valueCount = 1;

      for (int i = 1; i < count + 1; ++i){
        int assignmentCount = this->countChar(splittedStr[i], '=');
        if (assignmentCount != 1){
          throw "Error at the argument " + splittedStr[i];
        }
        std::string parameter_vs_val[2];
        this->splitString(splittedStr[i], '=', parameter_vs_val);
        std::string parameter = parameter_vs_val[0];
        std::string value = parameter_vs_val[1];
        //cout << "parameter:" << parameter << " value:" << value << endl;

        if (parameter == "a1"){
          a1 = this->fromString<T>(value);
        }
        else if (parameter == "a2"){
          if(this->coordinate_dimension > 1){
            a2 = this->fromString<T>(value);
          }
          else {
            throw  parameter+ " argument is not valid when dimension is " + toString<int>(this->coordinate_dimension);
          }

        }
        else if (parameter == "a3"){
          if(this->coordinate_dimension > 2){
            a3 = this->fromString<T>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + toString<int>(this->coordinate_dimension);
          }
        }
        else if (parameter == "b1"){
          b1 = this->fromString<T>(value);
        }
        else if (parameter == "b2"){
          if(this->coordinate_dimension > 1){
            b2 = this->fromString<T>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + toString<int>(this->coordinate_dimension);
          }
        }
        else if (parameter == "b3"){

          if(this->coordinate_dimension > 2){
            b3 = this->fromString<T>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + toString<int>(this->coordinate_dimension);
          }
        }
        else if (parameter == "c"){
          c = this->fromString<T>(value);
        }
        else if (parameter == "x1"){
          x1 = this->fromString<T>(value);
        }
        else if (parameter == "y1"){
          if(this->coordinate_dimension > 1){
            y1 = this->fromString<T>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + toString<int>(this->coordinate_dimension);
          }
        }
        else if (parameter == "z1"){
          if(this->coordinate_dimension > 2){
            z1 = this->fromString<T>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + toString<int>(this->coordinate_dimension);
          }
        }
        else if (parameter == "steps"){
          stepCount = this->countChar(value, ',') + 1;
          std::string *stepstr = new std::string[stepCount];
          this->splitString(value, ',', stepstr);
          steps = new T[stepCount];
          for (int i = 0; i < stepCount; ++i){
            steps[i] = this->fromString<T>(stepstr[i]);
          }
          delete [] stepstr;
        }
        else if (parameter == "values"){
          valueCount = this->countChar(value, ',') + 1;
          std::string *stepstr = new std::string[valueCount];
          this->splitString(value, ',', stepstr);
          values = new T[valueCount];
          for (int i = 0; i < valueCount; ++i){
            values[i] = this->fromString<T>(stepstr[i]);
          }
          delete [] stepstr;
        }
        else {
          throw "Invalid parameter name at " + splittedStr[i];
        }
      }

      delete []splittedStr;
      if(stepCount + 1!= valueCount){
        throw "Step count: " + this->toString<int>(stepCount) + " must be 1 less than value count: " + this->toString<int>(valueCount);
      }


      this->wd[wcount++] =  new SteppedEquation<T,T>(a1, a2,  a3,  b1,  b2,  b3,  c,  x1,  y1,  z1, steps, values, stepCount);

      if(stepCount > 0){
        delete [] steps;
        delete [] values;

      }
    }
    if(wcount != this->weight_dimension){
      throw "Weight Dimension is provided as " + toString<int>(wdimension) + ". But " + toString<int>(wcount)+" weight distributions are provided.";
    }
  }

  void parseParams(Teuchos::ParameterList params){
    try {
      std::string holeDescription  = "";
      std::string proc_load_distributions = "";
      std::string distinctDescription = "";
      std::string coordinate_distributions = "";
      std::string outfile = "";
      std::string weight_dimension_parameters[MAX_WEIGHT_DIM];
      for (int i = 0; i < MAX_WEIGHT_DIM; ++i){
        weight_dimension_parameters[i] = "";
      }


      for (Teuchos::ParameterList::ConstIterator pit = params.begin(); pit != params.end(); ++pit ){
        const std::string &paramName = params.name(pit);

        const Teuchos::ParameterEntry &pe = params.entry(pit);

        if(paramName.find("hole-") == 0){
          if(holeDescription == ""){
            holeDescription = getParamVal<std::string>(pe, paramName);
          }
          else {
            holeDescription +=","+getParamVal<std::string>(pe, paramName);
          }
        }
        else if(paramName.find("distribution-") == 0){
          if(coordinate_distributions == "")
            coordinate_distributions = getParamVal<std::string>(pe, paramName);
          else
            coordinate_distributions +=  ","+getParamVal<std::string>(pe, paramName);
          //cout << coordinate_distributions << endl;
          //TODO coordinate distribution description
        }

        else if (paramName.find(weight_distribution_string) == 0){
          std::string weight_dist_param = paramName.substr(weight_distribution_string.size());
          int dash_pos = weight_dist_param.find("-");
          std::string distribution_index_string = weight_dist_param.substr(0, dash_pos);
          int distribution_index = fromString<int>(distribution_index_string);

          if(distribution_index >= MAX_WEIGHT_DIM){
            throw "Given distribution index:" + distribution_index_string + " larger than maximum allowed weight dimension:" + toString<int>(MAX_WEIGHT_DIM);
          }
          weight_dimension_parameters[distribution_index] +=  " " + weight_dist_param.substr(dash_pos + 1)+ "="+ getParamVal<std::string>(pe, paramName);
        }
        else if(paramName == "dim"){
          int dim = fromString<int>(getParamVal<std::string>(pe, paramName));
          if(dim < 2 && dim > 3){
            throw INVALID(paramName);
          } else {
            this->coordinate_dimension = dim;
          }
        }
        else if(paramName == "wdim"){
          int dim = fromString<int>(getParamVal<std::string>(pe, paramName));
          if(dim < 1 && dim > MAX_WEIGHT_DIM){
            throw INVALID(paramName);
          } else {
            this->weight_dimension = dim;
          }
        }

        else if(paramName == "proc_load_distributions"){
          proc_load_distributions = getParamVal<std::string>(pe, paramName);
          this->loadDistSet = true;
        }

        else if(paramName == "distinct_coordinates"){
          distinctDescription = getParamVal<std::string>(pe, paramName);
          if(distinctDescription == "T"){
            this->distinctCoordSet = true;
          } else if(distinctDescription == "F"){
            this->distinctCoordSet = true;
          } else {
            throw "Invalid parameter for distinct_coordinates: " + distinctDescription + ". Candiates: T or F." ;
          }
        }

        else if(paramName == "out_file"){
          this->outfile = getParamVal<std::string>(pe, paramName);
        }

        else {
          throw INVALID(paramName);
        }
      }


      if(this->coordinate_dimension == 0){
        throw "Dimension must be provided to coordinate generator.";
      }
      /*
			if(this->globalNumCoords == 0){
				throw "Number of coordinates must be provided to coordinate generator.";
			}
       */
      /*
			if(maxx <= minx ){
				throw "Error: maxx= "+ toString<T>(maxx)+ " and minx=" + toString<T>(minx);
			}
			if(maxy <= miny ){
				throw "Error: maxy= "+ toString<T>(maxy)+ " and miny=" + toString<T>(miny);

			}
			if(this->dimension == 3 && maxz <= minz ){
				throw "Error: maxz= "+ toString<T>(maxz)+ " and minz=" + toString<T>(minz);
			}
       */
      if (this->loadDistSet && this->distinctCoordSet){
        throw "distinct_coordinates and proc_load_distributions parameters cannot be satisfied together.";
      }
      this->getHoles(holeDescription);
      //this->getDistinctCoordinateDescription(distinctDescription);
      this->getProcLoadDistributions(proc_load_distributions);
      this->getCoordinateDistributions(coordinate_distributions);
      this->getWeightDistribution(weight_dimension_parameters, this->weight_dimension);
      /*
			if(this->numGlobalCoords <= 0){
				throw "Must have at least 1 point";
			}
       */
    }
    catch(std::string s){
      throw s;
    }

    catch(char * s){
      throw s;
    }

    catch(char const* s){
      throw s;
    }

  }
public:

  ~GeometricGenerator(){
    if(this->holes){
      for (int i = 0; i < this->holeCount; ++i){
        delete this->holes[i];
      }
      free (this->holes);
    }


    delete []loadDistributions; //sized as the number of processors, the load of each processor.
    //delete []distinctCoordinates; // if processors have different or same range for coordinates to be created.
    if(coordinateDistributions){

      for (int i = 0; i < this->distributionCount; ++i){
        delete this->coordinateDistributions[i];
      }
      free (this->coordinateDistributions);
    }
    if (this->wd){
      for (int i = 0; i < this->weight_dimension; ++i){
        delete this->wd[i];
      }
      delete []this->wd;
    }

    if(this->weight_dimension){
      for(int i = 0; i < this->weight_dimension; ++i)
      delete [] this->wghts[i];
      delete []this->wghts;
    }
    if(this->coordinate_dimension){
      for(int i = 0; i < this->coordinate_dimension; ++i)
      delete [] this->coords[i];
      delete [] this->coords;
    }
    //delete []this->points;
  }

  GeometricGenerator(Teuchos::ParameterList &params, const RCP<const Teuchos::Comm<int> > & comm_){
    this->wd = NULL;
    this->holes = NULL; //to represent if there is any hole in the input
    this->coordinate_dimension = 0;  //dimension of the geometry
    this->weight_dimension = 0;  //dimension of the geometry
    this->worldSize = comm_->getSize(); //comminication world object.
    this->numGlobalCoords = 0;	//global number of coordinates requested to be created.
    this->loadDistributions = NULL; //sized as the number of processors, the load of each processor.
    //this->distinctCoordinates = NULL; // if processors have different or same range for coordinates to be created.
    this->coordinateDistributions = NULL;
    this->holeCount = 0;
    this->distributionCount = 0;
    this->outfile = "";
    //this->points =  NULL;

    /*
		this->minx = 0;
		this->maxx = 0;
		this->miny = 0;
		this->maxy = 0;
		this->minz = 0;
		this->maxz = 0;
     */
    this->loadDistSet = false;
    this->distinctCoordSet = false;
    this->parseParams(params);


    this->myRank = comm_->getRank();

    lno_t myPointCount = 0;
    this->numGlobalCoords = 0;

    gno_t prefixSum = 0;
    for(int i = 0; i < this->distributionCount; ++i){
      for(int ii = 0; ii < worldSize; ++ii){
        lno_t increment  = lno_t (this->coordinateDistributions[i]->numPoints * this->loadDistributions[ii]);
        if (ii < this->coordinateDistributions[i]->numPoints % worldSize){
          increment += 1;
        }
        this->numGlobalCoords += increment;
        if(ii < myRank){
          prefixSum += increment;
        }
      }
      myPointCount += lno_t(this->coordinateDistributions[i]->numPoints * this->loadDistributions[myRank]);
      if (myRank < this->coordinateDistributions[i]->numPoints % worldSize){
        myPointCount += 1;
      }
    }

    this->coords = new T *[this->coordinate_dimension];
    for(int i = 0; i < this->coordinate_dimension; ++i){
      this->coords[i] = new T[myPointCount];
    }

    for (int ii = 0; ii < this->coordinate_dimension; ++ii){
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
      for(int i = 0; i < myPointCount; ++i){
        this->coords[ii][i] = 0;
      }
    }

    this->numLocalCoords = 0;
    srand (myRank + 1);
    for (int i = 0; i < distributionCount; ++i){

      lno_t requestedPointCount = lno_t(this->coordinateDistributions[i]->numPoints *  this->loadDistributions[myRank]);
      if (myRank < this->coordinateDistributions[i]->numPoints % worldSize){
        requestedPointCount += 1;
      }
      //cout << "req:" << requestedPointCount << endl;
      //this->coordinateDistributions[i]->GetPoints(requestedPointCount,this->points + this->numLocalCoords, this->holes, this->holeCount,  this->loadDistributions, myRank);
      this->coordinateDistributions[i]->GetPoints(requestedPointCount,this->coords, this->numLocalCoords, this->holes, this->holeCount,  this->loadDistributions, myRank);
      this->numLocalCoords += requestedPointCount;
    }



    if (this->distinctCoordSet){
      //TODO: Partition and migration.
    }

    if(this->outfile != ""){

      std::ofstream myfile;
      myfile.open ((this->outfile + toString<int>(myRank)).c_str());
      for(lno_t i = 0; i < this->numLocalCoords; ++i){

        myfile << this->coords[0][i];
        if(this->coordinate_dimension > 1){
          myfile << " " << this->coords[1][i];
        }
        if(this->coordinate_dimension > 2){
          myfile << " " << this->coords[2][i];
        }
        myfile << std::endl;
      }
      myfile.close();
    }



    /*
		Zoltan2::XpetraMultiVectorInput < Tpetra::MultiVector<T, lno_t, gno_t, node_t> > xmv (RCP < Tpetra::MultiVector<T, lno_t, gno_t, node_t> > (tmVector));

		RCP< Tpetra::MultiVector<T, lno_t, gno_t, node_t> >tmVector2;
		Zoltan2::PartitioningSolution< Tpetra::MultiVector<T, lno_t, gno_t, node_t> > solution;
		xmv.applyPartitioningSolution<Tpetra::MultiVector<T, lno_t, gno_t, node_t> >(this->tmVector, &tmVector2, solution);
     */

    this->wghts = new T *[this->weight_dimension];
    for(int i = 0; i < this->weight_dimension; ++i){
      this->wghts[i] = new T[this->numLocalCoords];
    }

    for(int ii = 0; ii < this->weight_dimension; ++ii){
      switch(this->coordinate_dimension){
      case 1:
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for (lno_t i = 0; i < this->numLocalCoords; ++i){
          this->wghts[ii][i] = this->wd[ii]->get1DWeight(this->coords[0][i]);
        }
        break;
      case 2:
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for (lno_t i = 0; i < this->numLocalCoords; ++i){
          this->wghts[ii][i] = this->wd[ii]->get2DWeight(this->coords[0][i], this->coords[1][i]);
        }
        break;
      case 3:
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for (lno_t i = 0; i < this->numLocalCoords; ++i){
          this->wghts[ii][i] = this->wd[ii]->get3DWeight(this->coords[0][i], this->coords[1][i], this->coords[2][i]);
        }
        break;
      }
    }
  }

  int getWeightDimension(){
    return this->weight_dimension;
  }
  int getCoordinateDimension(){
    return this->coordinate_dimension;
  }
  lno_t getNumLocalCoords(){
    return this->numLocalCoords;
  }
  gno_t getNumGlobalCoords(){
    return this->numGlobalCoords;
  }

  T **getLocalCoordinatesView(){
    return this->coords;
  }

  T **getLocalWeightsView(){
    return this->wghts;
  }

  void getLocalCoordinatesCopy( T ** c){
    for(int ii = 0; ii < this->coordinate_dimension; ++ii){
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
      for (lno_t i = 0; i < this->numLocalCoords; ++i){
        c[ii][i] = this->coords[ii][i];
      }
    }
  }

  void getLocalWeightsCopy(T **w){
    for(int ii = 0; ii < this->weight_dimension; ++ii){
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
      for (lno_t i = 0; i < this->numLocalCoords; ++i){
        w[ii][i] = this->wghts[ii][i];
      }
    }
  }
};
