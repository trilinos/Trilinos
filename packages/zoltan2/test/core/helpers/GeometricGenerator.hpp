// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GEOMETRICGENERATOR
#define GEOMETRICGENERATOR

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
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Teuchos_ArrayViewDecl.hpp>
#include <Teuchos_RCP.hpp>
#include <Tpetra_Distributor.hpp>
#include <Zoltan2_PartitioningProblem.hpp>


#include <zoltan.h>

#ifdef _MSC_VER
#define NOMINMAX
#include <windows.h>
#endif

using Teuchos::CommandLineProcessor;
using Teuchos::ArrayView;
using std::string;
using Teuchos::ArrayRCP;

namespace GeometricGen{
#define CATCH_EXCEPTIONS(pp) \
        catch (std::runtime_error &e) { \
            std::cout << "Runtime exception returned from " << pp << ": " \
            << e.what() << " FAIL" << std::endl; \
            return -1; \
        } \
        catch (std::logic_error &e) { \
            std::cout << "Logic exception returned from " << pp << ": " \
            << e.what() << " FAIL" << std::endl; \
            return -1; \
        } \
        catch (std::bad_alloc &e) { \
            std::cout << "Bad_alloc exception returned from " << pp << ": " \
            << e.what() << " FAIL" << std::endl; \
            return -1; \
        } \
        catch (std::exception &e) { \
            std::cout << "Unknown exception returned from " << pp << ": " \
            << e.what() << " FAIL" << std::endl; \
            return -1; \
        }




template <typename tMVector_t>
class DOTS{
public:
  std::vector<std::vector<float> > weights;
  tMVector_t *coordinates;
};

template <typename tMVector_t>
int getNumObj(void *data, int *ierr)
{
  *ierr = 0;
  DOTS<tMVector_t> *dots_ = (DOTS<tMVector_t> *) data;
  return dots_->coordinates->getLocalLength();
}
//////////////////////////
template <typename tMVector_t>
void getCoords(void *data, int numGid, int numLid,
  int numObj, ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
  int dim, double *coords_, int *ierr)
{
  typedef typename tMVector_t::scalar_type scalar_t;

  // I know that Zoltan asks for coordinates in gid order.
  if (dim == 3){
    *ierr = 0;
    DOTS<tMVector_t> *dots_ = (DOTS<tMVector_t> *) data;
    double *val = coords_;
    const scalar_t *x = dots_->coordinates->getData(0).getRawPtr();
    const scalar_t *y = dots_->coordinates->getData(1).getRawPtr();
    const scalar_t *z = dots_->coordinates->getData(2).getRawPtr();
    for (int i=0; i < numObj; i++){
      *val++ = static_cast<double>(x[i]);
      *val++ = static_cast<double>(y[i]);
      *val++ = static_cast<double>(z[i]);
    }
  }
  else {
    *ierr = 0;
    DOTS<tMVector_t> *dots_ = (DOTS<tMVector_t> *) data;
    double *val = coords_;
    const scalar_t *x = dots_->coordinates->getData(0).getRawPtr();
    const scalar_t *y = dots_->coordinates->getData(1).getRawPtr();
    for (int i=0; i < numObj; i++){
      *val++ = static_cast<double>(x[i]);
      *val++ = static_cast<double>(y[i]);
    }
  }
}

template <typename tMVector_t>
int getDim(void *data, int *ierr)
{
  *ierr = 0;
  DOTS<tMVector_t> *dots_ = (DOTS<tMVector_t> *) data;
  int dim =  dots_->coordinates->getNumVectors();

  return dim;
}

//////////////////////////
template <typename tMVector_t>
void getObjList(void *data, int numGid, int numLid,
  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids,
  int num_wgts, float *obj_wgts, int *ierr)
{
  typedef typename tMVector_t::global_ordinal_type gno_t;

  *ierr = 0;
  DOTS<tMVector_t> *dots_ = (DOTS<tMVector_t> *) data;

  size_t localLen = dots_->coordinates->getLocalLength();
  const gno_t *ids =
               dots_->coordinates->getMap()->getLocalElementList().getRawPtr();

  if (sizeof(ZOLTAN_ID_TYPE) == sizeof(gno_t))
    memcpy(gids, ids, sizeof(ZOLTAN_ID_TYPE) * localLen);
  else
    for (size_t i=0; i < localLen; i++)
      gids[i] = static_cast<ZOLTAN_ID_TYPE>(ids[i]);

  if (num_wgts > 0){
    float *wgts = obj_wgts;
    for (size_t i=0; i < localLen; i++)
      for (int w=0; w < num_wgts; w++)
        *wgts++ = dots_->weights[w][i];
  }
}


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

const std::string weight_distribution_string = "WeightDistribution-";

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
  virtual ~Hole(){}
};

template <typename T>
class SquareHole:public Hole<T>{
public:
  SquareHole(CoordinatePoint<T> center_ , T edge_): Hole<T>(center_, edge_, 0 , 0){}
  virtual ~SquareHole(){}
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
  virtual ~RectangleHole(){}
};

template <typename T>
class CircleHole:public Hole<T>{
public:
  CircleHole(CoordinatePoint<T> center_ , T edge_): Hole<T>(center_, edge_, 0 , 0){}
  virtual ~CircleHole(){}
  virtual bool isInArea(CoordinatePoint<T> dot){
    return (dot.x - this->center.x)*(dot.x - this->center.x) + (dot.y - this->center.y) * (dot.y - this->center.y) < this->edge1 * this->edge1;
  }
};

template <typename T>
class CubeHole:public Hole<T>{
public:
  CubeHole(CoordinatePoint<T> center_ , T edge_): Hole<T>(center_, edge_, 0 , 0){}
  virtual ~CubeHole(){}
  virtual bool isInArea(CoordinatePoint<T> dot){
    return fabs(dot.x - this->center.x) < this->edge1 / 2 && fabs(dot.y - this->center.y) < this->edge1 / 2 && fabs(dot.z - this->center.z) < this->edge1 / 2;
  }
};

template <typename T>
class RectangularPrismHole:public Hole<T>{
public:
  RectangularPrismHole(CoordinatePoint<T> center_  , T edge_x_, T edge_y_, T edge_z_): Hole<T>(center_, edge_x_,  edge_y_, edge_z_){}
  virtual ~RectangularPrismHole(){}
  virtual bool isInArea(CoordinatePoint<T> dot){
    return fabs(dot.x - this->center.x) < this->edge1 / 2 && fabs(dot.y - this->center.y) < this->edge2 / 2 && fabs(dot.z - this->center.z) < this->edge3 / 2;
  }
};

template <typename T>
class SphereHole:public Hole<T>{
public:
  SphereHole(CoordinatePoint<T> center_ , T edge_): Hole<T>(center_, edge_, 0 , 0){}
  virtual ~SphereHole(){}
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
    std::cout << this->a1 << "*" <<  "math.pow( (" << x  << "- " <<  this->x1 << " ), " << b1 <<")" << "+" <<  this->a2<< "*"<<  "math.pow( (" << y << "-" <<  this->y1 << "), " <<
        b2 << " ) + " << this->a3 << " * math.pow( (" << z << "-" <<  this->z1 << "), " << b3 << ")+ " << c << " == " <<
        this->a1 * pow( (x - this->x1), b1) + this->a2 * pow( (y - this->y1), b2) + this->a3 * pow( (z - this->z1), b3) + c << std::endl;

  }

  virtual weighttype get3DWeight(T x, T y, T z){
    T expressionRes = this->a1 * pow( (x - this->x1), b1) + this->a2 * pow( (y - this->y1), b2) + this->a3 * pow( (z - this->z1), b3) + this->c;

    //this->print(x,y,z);
    if(this->stepCount > 0){
      for (int i = 0; i < this->stepCount; ++i){
        if (expressionRes < this->steps[i]) {
          //std::cout << "0exp:" << expressionRes << " step:" << steps[i] << " value:" << values[i] << std::endl;
          return this->values[i];
        }
      }

      //std::cout << "1exp:" << expressionRes << " step:" << steps[stepCount] << " value:" << values[stepCount] << std::endl;
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

  CoordinateDistribution(gno_t np_, int dim, int wSize) :
    numPoints(np_), dimension(dim), requested(0), assignedPrevious(0),
    worldSize(wSize){}

  virtual CoordinatePoint<T> getPoint(gno_t point_index, unsigned int &state) = 0;
  virtual T getXCenter() = 0;
  virtual T getXRadius() =0;

  void GetPoints(lno_t requestedPointcount, CoordinatePoint<T> *points /*preallocated sized numPoints*/,
      Hole<T> **holes, lno_t holeCount,
      float *sharedRatios_, int myRank){

    for (int i = 0; i < myRank; ++i){
      //std::cout << "me:" << myRank << " i:" << i << " s:" << sharedRatios_[i]<< std::endl;
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
          CoordinatePoint <T> p = this->getPoint( this->assignedPrevious + cnt, state);

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
      //std::cout << "me:" << myRank << " i:" << i << " s:" << sharedRatios_[i]<< std::endl;
      this->assignedPrevious += lno_t(sharedRatios_[i] * this->numPoints);
      if (gno_t(i) < this->numPoints % this->worldSize){
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

        std::cout << "myRank:" << me << " stateBeg:" << stateBegin << " tsize:" << tsize << " state:" << state <<  " slice: " << slice / tsize <<  std::endl;
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

  CoordinateNormalDistribution(gno_t np_, int dim, CoordinatePoint<T> center_ ,
                               T sd_x, T sd_y, T sd_z, int wSize) :
    CoordinateDistribution<T,lno_t,gno_t>(np_,dim,wSize),
    standartDevx(sd_x), standartDevy(sd_y), standartDevz(sd_z)
  {
    this->center.x = center_.x;
    this->center.y = center_.y;
    this->center.z = center_.z;
  }

  virtual CoordinatePoint<T> getPoint(gno_t pindex, unsigned int &state){

    //pindex = 0; // not used in normal distribution.
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
  T normalDist(T center_, T sd, unsigned int &state) {
    T polarsqrt, normalsquared, normal1, normal2;
    do {
      normal1=2.0*( T(rand_r(&state))/T(RAND_MAX) ) - 1.0;
      normal2=2.0*( T(rand_r(&state))/T(RAND_MAX) ) - 1.0;
      normalsquared=normal1*normal1+normal2*normal2;
    } while ( normalsquared>=1.0 || normalsquared == 0.0);

    polarsqrt=sqrt(-2.0*log(normalsquared)/normalsquared);
    return normal2*polarsqrt*sd + center_;
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


  CoordinateUniformDistribution(gno_t np_, int dim, T l_x, T r_x, T l_y, T r_y,
                                T l_z, T r_z, int wSize ) :
      CoordinateDistribution<T,lno_t,gno_t>(np_,dim,wSize),
      leftMostx(l_x), rightMostx(r_x), leftMosty(l_y), rightMosty(r_y),
      leftMostz(l_z), rightMostz(r_z){}

  virtual ~CoordinateUniformDistribution(){};
  virtual CoordinatePoint<T> getPoint(gno_t pindex, unsigned int &state){


    //pindex = 0; //not used in uniform dist.
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


  CoordinateGridDistribution(gno_t alongX, gno_t alongY, gno_t alongZ, int dim,
                             T l_x, T r_x, T l_y, T r_y, T l_z, T r_z ,
                             int myRank_, int wSize) :
      CoordinateDistribution<T,lno_t,gno_t>(alongX * alongY * alongZ,dim,wSize),
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
    //std::cout << "before:" << processCnt << " " << this->assignedPrevious << std::endl;
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
      std::cout << "req:" << this->requested << " pr:" <<this->processCnt << std::endl;
      std::cout << "p:" << p.x << " " << p.y << " " << p.z <<  std::endl ;
      std::cout << "s:" << xshift << " " << yshift << " " << zshift <<  std::endl ;
      std::cout << "st:" << this->xstep << " " << this->ystep << " " << this->zstep <<  std::endl ;
    }
    */
    return p;
  }

private:

};

template <typename scalar_t, typename lno_t, typename gno_t, typename node_t>
class GeometricGenerator {
private:
  Hole<scalar_t> **holes; //to represent if there is any hole in the input
  int holeCount;
  int coordinate_dimension;  //dimension of the geometry
  gno_t numGlobalCoords;	//global number of coordinates requested to be created.
  lno_t numLocalCoords;
  float *loadDistributions; //sized as the number of processors, the load of each processor.
  bool loadDistSet;
  bool distinctCoordSet;
  CoordinateDistribution<scalar_t, lno_t,gno_t> **coordinateDistributions;
  int distributionCount;
  //CoordinatePoint<scalar_t> *points;
  scalar_t **coords;
  scalar_t **wghts;
  WeightDistribution<scalar_t,scalar_t> **wd;
  int numWeightsPerCoord;
  int predistribution;
  RCP<const Teuchos::Comm<int> > comm;
  //RCP< Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> >tmVector;
  //RCP< Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> >tmwVector;
  int worldSize;
  int myRank;
  scalar_t minx;
  scalar_t maxx;
  scalar_t miny;
  scalar_t maxy;
  scalar_t minz;
  scalar_t maxz;
  std::string outfile;
  float perturbation_ratio;

  typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;
  typedef Tpetra::Map<lno_t, gno_t, node_t> tMap_t;


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

  int countChar (std::string inStr, char cntChar){
    int cnt = 0;
    for (unsigned int i = 0; i < inStr.size(); ++i){
      if (inStr[i] == cntChar) {
        cnt++;
      }
    }
    return cnt;
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
    coordinateDistributions = (CoordinateDistribution<scalar_t, lno_t,gno_t> **) malloc(sizeof (CoordinateDistribution<scalar_t, lno_t,gno_t> *) * 1);
    for(int i = 0; i < argCnt; ){
      coordinateDistributions = (CoordinateDistribution<scalar_t, lno_t,gno_t> **)realloc((void *)coordinateDistributions, (this->distributionCount + 1)* sizeof(CoordinateDistribution<scalar_t, lno_t,gno_t> *));

      std::string distName = splittedStr[i++];
      gno_t np_ = 0;
      if(distName == "NORMAL"){
        int reqArg = 5;
        if (this->coordinate_dimension == 3){
          reqArg = 7;
        }
        if(i + reqArg > argCnt) {
          std::string tmp = Teuchos::toString<int>(reqArg);
          throw INVALID_SHAPE_ARG(distName, tmp);
        }
        np_ = fromString<gno_t>(splittedStr[i++]);
        CoordinatePoint<scalar_t> pp;

        pp.x = fromString<scalar_t>(splittedStr[i++]);
        pp.y = fromString<scalar_t>(splittedStr[i++]);
        pp.z = 0;
        if(this->coordinate_dimension == 3){
          pp.z = fromString<scalar_t>(splittedStr[i++]);
        }

        scalar_t sd_x = fromString<scalar_t>(splittedStr[i++]);
        scalar_t sd_y = fromString<scalar_t>(splittedStr[i++]);
        scalar_t sd_z = 0;
        if(this->coordinate_dimension == 3){
          sd_z = fromString<scalar_t>(splittedStr[i++]);
        }
        this->coordinateDistributions[this->distributionCount++] = new CoordinateNormalDistribution<scalar_t, lno_t,gno_t>(np_, this->coordinate_dimension, pp , sd_x, sd_y, sd_z, this->worldSize );

      } else if(distName == "UNIFORM" ){
        int reqArg = 5;
        if (this->coordinate_dimension == 3){
          reqArg = 7;
        }
        if(i + reqArg > argCnt) {
          std::string tmp = Teuchos::toString<int>(reqArg);
          throw INVALID_SHAPE_ARG(distName, tmp);
        }
        np_ = fromString<gno_t>(splittedStr[i++]);
        scalar_t l_x = fromString<scalar_t>(splittedStr[i++]);
        scalar_t r_x = fromString<scalar_t>(splittedStr[i++]);
        scalar_t l_y = fromString<scalar_t>(splittedStr[i++]);
        scalar_t r_y = fromString<scalar_t>(splittedStr[i++]);

        scalar_t l_z = 0, r_z = 0;

        if(this->coordinate_dimension == 3){
          l_z = fromString<scalar_t>(splittedStr[i++]);
          r_z = fromString<scalar_t>(splittedStr[i++]);
        }

        this->coordinateDistributions[this->distributionCount++] = new CoordinateUniformDistribution<scalar_t, lno_t,gno_t>( np_,  this->coordinate_dimension, l_x, r_x, l_y, r_y, l_z, r_z, this->worldSize );
      } else if (distName == "GRID"){
        int reqArg = 6;
        if(this->coordinate_dimension == 3){
          reqArg = 9;
        }
        if(i + reqArg > argCnt) {
          std::string tmp = Teuchos::toString<int>(reqArg);
          throw INVALID_SHAPE_ARG(distName, tmp);
        }

        gno_t np_x = fromString<gno_t>(splittedStr[i++]);
        gno_t np_y = fromString<gno_t>(splittedStr[i++]);
        gno_t np_z = 1;


        if(this->coordinate_dimension == 3){
          np_z = fromString<gno_t>(splittedStr[i++]);
        }

        np_ = np_x * np_y * np_z;
        scalar_t l_x = fromString<scalar_t>(splittedStr[i++]);
        scalar_t r_x = fromString<scalar_t>(splittedStr[i++]);
        scalar_t l_y = fromString<scalar_t>(splittedStr[i++]);
        scalar_t r_y = fromString<scalar_t>(splittedStr[i++]);

        scalar_t l_z = 0, r_z = 0;

        if(this->coordinate_dimension == 3){
          l_z = fromString<scalar_t>(splittedStr[i++]);
          r_z = fromString<scalar_t>(splittedStr[i++]);
        }

        if(np_x < 1 || np_z < 1 || np_y < 1 ){
          throw "Provide at least 1 point along each dimension for grid test.";
        }
        //std::cout << "ly:" << l_y << " ry:" << r_y << std::endl;
        this->coordinateDistributions[this->distributionCount++] = new CoordinateGridDistribution<scalar_t, lno_t,gno_t>
        (np_x, np_y,np_z, this->coordinate_dimension, l_x, r_x,l_y, r_y, l_z, r_z , this->myRank, this->worldSize);

      }
      else {
        std::string tmp = Teuchos::toString<int>(this->coordinate_dimension);
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
        throw "Invalid parameter count load distributions. Given " + Teuchos::toString<int>(argCnt) + " processor size is " + Teuchos::toString<int>(this->worldSize);
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
    this->holes =  (Hole<scalar_t> **) malloc(sizeof (Hole <scalar_t>*));
    int argCnt = this->countChar(holeDescription, ',') + 1;
    std::string *splittedStr = new std::string[argCnt];
    splitString(holeDescription, ',', splittedStr);

    for(int i = 0; i < argCnt; ){
      holes = (Hole<scalar_t> **)realloc((void *)holes, (this->holeCount + 1)* sizeof(Hole<scalar_t> *));

      std::string shapeName = splittedStr[i++];
      if(shapeName == "SQUARE" && this->coordinate_dimension == 2){
        if(i + 3 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "3");
        }
        CoordinatePoint<scalar_t> pp;
        pp.x = fromString<scalar_t>(splittedStr[i++]);
        pp.y = fromString<scalar_t>(splittedStr[i++]);
        scalar_t edge = fromString<scalar_t>(splittedStr[i++]);
        this->holes[this->holeCount++] = new SquareHole<scalar_t>(pp, edge);
      } else if(shapeName == "RECTANGLE" && this->coordinate_dimension == 2){
        if(i + 4 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "4");
        }
        CoordinatePoint<scalar_t> pp;
        pp.x = fromString<scalar_t>(splittedStr[i++]);
        pp.y = fromString<scalar_t>(splittedStr[i++]);
        scalar_t edgex = fromString<scalar_t>(splittedStr[i++]);
        scalar_t edgey = fromString<scalar_t>(splittedStr[i++]);

        this->holes[this->holeCount++] = new RectangleHole<scalar_t>(pp, edgex,edgey);
      } else if(shapeName == "CIRCLE" && this->coordinate_dimension == 2){
        if(i + 3 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "3");
        }
        CoordinatePoint<scalar_t> pp;
        pp.x = fromString<scalar_t>(splittedStr[i++]);
        pp.y = fromString<scalar_t>(splittedStr[i++]);
        scalar_t r = fromString<scalar_t>(splittedStr[i++]);
        this->holes[this->holeCount++] = new CircleHole<scalar_t>(pp, r);
      }  else if(shapeName == "CUBE" && this->coordinate_dimension == 3){
        if(i + 4 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "4");
        }
        CoordinatePoint<scalar_t> pp;
        pp.x = fromString<scalar_t>(splittedStr[i++]);
        pp.y = fromString<scalar_t>(splittedStr[i++]);
        pp.z = fromString<scalar_t>(splittedStr[i++]);
        scalar_t edge = fromString<scalar_t>(splittedStr[i++]);
        this->holes[this->holeCount++] = new CubeHole<scalar_t>(pp, edge);
      }  else if(shapeName == "RECTANGULAR_PRISM" && this->coordinate_dimension == 3){
        if(i + 6 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "6");
        }
        CoordinatePoint<scalar_t> pp;
        pp.x = fromString<scalar_t>(splittedStr[i++]);
        pp.y = fromString<scalar_t>(splittedStr[i++]);
        pp.z = fromString<scalar_t>(splittedStr[i++]);
        scalar_t edgex = fromString<scalar_t>(splittedStr[i++]);
        scalar_t edgey = fromString<scalar_t>(splittedStr[i++]);
        scalar_t edgez = fromString<scalar_t>(splittedStr[i++]);
        this->holes[this->holeCount++] = new RectangularPrismHole<scalar_t>(pp, edgex, edgey, edgez);

      }  else if(shapeName == "SPHERE" && this->coordinate_dimension == 3){
        if(i + 4 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "4");
        }
        CoordinatePoint<scalar_t> pp;
        pp.x = fromString<scalar_t>(splittedStr[i++]);
        pp.y = fromString<scalar_t>(splittedStr[i++]);
        pp.z = fromString<scalar_t>(splittedStr[i++]);
        scalar_t r = fromString<scalar_t>(splittedStr[i++]);
        this->holes[this->holeCount++] = new SphereHole<scalar_t>(pp, r);
      }  else {
        std::string tmp = Teuchos::toString<int>(this->coordinate_dimension);
        throw INVALIDSHAPE(shapeName, tmp);
      }
    }
    delete [] splittedStr;
  }

  void getWeightDistribution(std::string *weight_distribution_arr, int wdimension){
    int wcount = 0;

    this->wd = new WeightDistribution<scalar_t,scalar_t> *[wdimension];
    for(int ii = 0; ii < MAX_WEIGHT_DIM; ++ii){
      std::string weight_distribution = weight_distribution_arr[ii];
      if(weight_distribution == "") continue;
      if(wcount == wdimension) {
        throw "Weight Dimension is provided as " + Teuchos::toString<int>(wdimension) + ". More weight distribution is provided.";
      }

      int count = this->countChar(weight_distribution, ' ');
      std::string *splittedStr = new string[count + 1];
      this->splitString(weight_distribution, ' ', splittedStr);
      //std::cout << count << std::endl;
      scalar_t c=1;
      scalar_t a1=0,a2=0,a3=0;
      scalar_t x1=0,y1=0,z1=0;
      scalar_t b1=1,b2=1,b3=1;
      scalar_t *steps = NULL;
      scalar_t *values= NULL;
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
        //std::cout << "parameter:" << parameter << " value:" << value << std::endl;

        if (parameter == "a1"){
          a1 = this->fromString<scalar_t>(value);
        }
        else if (parameter == "a2"){
          if(this->coordinate_dimension > 1){
            a2 = this->fromString<scalar_t>(value);
          }
          else {
            throw  parameter+ " argument is not valid when dimension is " + Teuchos::toString<int>(this->coordinate_dimension);
          }

        }
        else if (parameter == "a3"){
          if(this->coordinate_dimension > 2){
            a3 = this->fromString<scalar_t>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + Teuchos::toString<int>(this->coordinate_dimension);
          }
        }
        else if (parameter == "b1"){
          b1 = this->fromString<scalar_t>(value);
        }
        else if (parameter == "b2"){
          if(this->coordinate_dimension > 1){
            b2 = this->fromString<scalar_t>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + Teuchos::toString<int>(this->coordinate_dimension);
          }
        }
        else if (parameter == "b3"){

          if(this->coordinate_dimension > 2){
            b3 = this->fromString<scalar_t>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + Teuchos::toString<int>(this->coordinate_dimension);
          }
        }
        else if (parameter == "c"){
          c = this->fromString<scalar_t>(value);
        }
        else if (parameter == "x1"){
          x1 = this->fromString<scalar_t>(value);
        }
        else if (parameter == "y1"){
          if(this->coordinate_dimension > 1){
            y1 = this->fromString<scalar_t>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + Teuchos::toString<int>(this->coordinate_dimension);
          }
        }
        else if (parameter == "z1"){
          if(this->coordinate_dimension > 2){
            z1 = this->fromString<scalar_t>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + Teuchos::toString<int>(this->coordinate_dimension);
          }
        }
        else if (parameter == "steps"){
          stepCount = this->countChar(value, ',') + 1;
          std::string *stepstr = new std::string[stepCount];
          this->splitString(value, ',', stepstr);
          steps = new scalar_t[stepCount];
          for (int j = 0; j < stepCount; ++j){
            steps[j] = this->fromString<scalar_t>(stepstr[j]);
          }
          delete [] stepstr;
        }
        else if (parameter == "values"){
          valueCount = this->countChar(value, ',') + 1;
          std::string *stepstr = new std::string[valueCount];
          this->splitString(value, ',', stepstr);
          values = new scalar_t[valueCount];
          for (int j = 0; j < valueCount; ++j){
            values[j] = this->fromString<scalar_t>(stepstr[j]);
          }
          delete [] stepstr;
        }
        else {
          throw "Invalid parameter name at " + splittedStr[i];
        }
      }

      delete []splittedStr;
      if(stepCount + 1!= valueCount){
        throw "Step count: " + Teuchos::toString<int>(stepCount) + " must be 1 less than value count: " + Teuchos::toString<int>(valueCount);
      }


      this->wd[wcount++] =  new SteppedEquation<scalar_t,scalar_t>(a1, a2,  a3,  b1,  b2,  b3,  c,  x1,  y1,  z1, steps, values, stepCount);

      if(stepCount > 0){
        delete [] steps;
        delete [] values;

      }
    }
    if(wcount != this->numWeightsPerCoord){
      throw "Weight Dimension is provided as " + Teuchos::toString<int>(wdimension) + ". But " + Teuchos::toString<int>(wcount)+" weight distributions are provided.";
    }
  }

  void parseParams(const Teuchos::ParameterList & params){
    try {
      std::string holeDescription  = "";
      std::string proc_load_distributions = "";
      std::string distinctDescription = "";
      std::string coordinate_distributions = "";
      std::string numWeightsPerCoord_parameters[MAX_WEIGHT_DIM];
      for (int i = 0; i < MAX_WEIGHT_DIM; ++i){
        numWeightsPerCoord_parameters[i] = "";
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
          //std::cout << coordinate_distributions << std::endl;
          //TODO coordinate distribution description
        }

        else if (paramName.find(weight_distribution_string) == 0){
          std::string weight_dist_param = paramName.substr(weight_distribution_string.size());
          int dash_pos = weight_dist_param.find("-");
          std::string distribution_index_string = weight_dist_param.substr(0, dash_pos);
          int distribution_index = fromString<int>(distribution_index_string);

          if(distribution_index >= MAX_WEIGHT_DIM){
            throw "Given distribution index:" + distribution_index_string + " larger than maximum allowed number of weights:" + Teuchos::toString<int>(MAX_WEIGHT_DIM);
          }
          numWeightsPerCoord_parameters[distribution_index] +=  " " + weight_dist_param.substr(dash_pos + 1)+ "="+ getParamVal<std::string>(pe, paramName);
        }
        else if(paramName == "dim"){
          int dim = fromString<int>(getParamVal<std::string>(pe, paramName));
          if(dim < 2 || dim > 3){
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
            this->numWeightsPerCoord = dim;
          }
        }
        else if(paramName == "predistribution"){
          int pre_distribution = fromString<int>(getParamVal<std::string>(pe, paramName));
          if(pre_distribution < 0 || pre_distribution > 3){
            throw INVALID(paramName);
          } else {
            this->predistribution = pre_distribution;
          }
        }
        else if(paramName == "perturbation_ratio"){
          float perturbation = fromString<float>(getParamVal<std::string>(pe, paramName));
          if(perturbation < 0 && perturbation > 1){
            throw INVALID(paramName);
          } else {
            this->perturbation_ratio = perturbation;
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
            throw "Invalid parameter for distinct_coordinates: " + distinctDescription + ". Candidates: T or F." ;
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
				throw "Error: maxx= "+ Teuchos::toString<scalar_t>(maxx)+ " and minx=" + Teuchos::toString<scalar_t>(minx);
			}
			if(maxy <= miny ){
				throw "Error: maxy= "+ Teuchos::toString<scalar_t>(maxy)+ " and miny=" + Teuchos::toString<scalar_t>(miny);

			}
			if(this->dimension == 3 && maxz <= minz ){
				throw "Error: maxz= "+ Teuchos::toString<scalar_t>(maxz)+ " and minz=" + Teuchos::toString<scalar_t>(minz);
			}
       */
      if (this->loadDistSet && this->distinctCoordSet){
        throw "distinct_coordinates and proc_load_distributions parameters cannot be satisfied together.";
      }
      this->getHoles(holeDescription);
      //this->getDistinctCoordinateDescription(distinctDescription);
      this->getProcLoadDistributions(proc_load_distributions);
      this->getCoordinateDistributions(coordinate_distributions);
      this->getWeightDistribution(numWeightsPerCoord_parameters, this->numWeightsPerCoord);
      /*
			if(this->numGlobalCoords <= 0){
				throw "Must have at least 1 point";
			}
       */
    }
    catch(std::string &s){
      throw std::runtime_error(s);
    }

    catch(const char *s){
      throw std::runtime_error(s);
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
      for (int i = 0; i < this->numWeightsPerCoord; ++i){
        delete this->wd[i];
      }
      delete []this->wd;
    }

    if(this->numWeightsPerCoord){
      for(int i = 0; i < this->numWeightsPerCoord; ++i)
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

  void print_description(){
    std::cout <<"\nGeometric Generator Parameter File Format:" << std::endl;
    std::cout <<"- dim=Coordinate Dimension: 2 or 3" << std::endl;
    std::cout <<"- Available distributions:" << std::endl;
    std::cout <<"\tUNIFORM: -> distribution-1=UNIFORM,NUMCOORDINATES,XMIN,XMAX,YMIN,YMAX{,ZMIN,ZMAX}" << std::endl;
    std::cout <<"\tGRID: -> distribution-2=GRID,XLENGTH,YLENGTH{,ZLENGTH},XMIN,XMAX,YMIN,YMAX{,ZMIN,ZMAX}" << std::endl;
    std::cout <<"\tNORMAL: -> distribution-3=NORMAL,XCENTER,YCENTER{,ZCENTER},XSD,YSD,{,ZSD}" << std::endl;
    std::cout <<"- wdim=numWeightsPerCoord:  There should be as many weight function as number of weights per coord." << std::endl;
    std::cout <<"- Weight Equation: w = (a1 * (x - x1)^b1) + (a2 * (y - y1)^b2) + (a3 * (z - z1)^b3) + c" << std::endl;
    std::cout << "Parameter settings:" << std::endl;
    std::cout << "\tWeightDistribution-1-a1=a1 " << std::endl;
    std::cout << "\tWeightDistribution-1-a2=a2 " << std::endl;
    std::cout << "\tWeightDistribution-1-a3=a3 " << std::endl;
    std::cout << "\tWeightDistribution-1-b1=b1 " << std::endl;
    std::cout << "\tWeightDistribution-1-b2=b2 " << std::endl;
    std::cout << "\tWeightDistribution-1-b3=b3 " << std::endl;
    std::cout << "\tWeightDistribution-1-x1=x1 " << std::endl;
    std::cout << "\tWeightDistribution-1-y1=y1 " << std::endl;
    std::cout << "\tWeightDistribution-1-z1=z1 " << std::endl;
    std::cout << "\tWeightDistribution-1-c=c " << std::endl;
    std::cout << "\tIt is possible to set step function to the result of weight equation." << std::endl;
    std::cout << "\tWeightDistribution-1-steps=step1,step2,step3:increasing order" << std::endl;
    std::cout << "\tWeightDistribution-1-values=val1,val2,val3,val4." << std::endl;
    std::cout << "\t\tIf w < step1 -> w = val1" << std::endl;
    std::cout << "\t\tElse if w < step2 -> w = val2" << std::endl;
    std::cout << "\t\tElse if w < step3 -> w = val3" << std::endl;
    std::cout << "\t\tElse  -> w = val4" << std::endl;
    std::cout <<"- Holes:" << std::endl;
    std::cout << "\thole-1:SPHERE,XCENTER,YCENTER,ZCENTER,RADIUS (only for dim=3)" << std::endl;
    std::cout << "\thole-2:CUBE,XCENTER,YCENTER,ZCENTER,EDGE (only for dim=3)" << std::endl;
    std::cout << "\thole-3:RECTANGULAR_PRISM,XCENTER,YCENTER,ZCENTER,XEDGE,YEDGE,ZEDGE (only for dim=3)" << std::endl;
    std::cout << "\thole-4:SQUARE,XCENTER,YCENTER,EDGE (only for dim=2)" << std::endl;
    std::cout << "\thole-5:RECTANGLE,XCENTER,YCENTER,XEDGE,YEDGE (only for dim=2)" << std::endl;
    std::cout << "\thole-6:CIRCLE,XCENTER,YCENTER,RADIUS (only for dim=2)" << std::endl;
    std::cout << "- out_file:out_file_path : if provided output will be written to files." << std::endl;
    std::cout << "- proc_load_distributions:ratio_0, ratio_1, ratio_2....ratio_n. Loads of each processor, should be as many as MPI ranks and should sum up to 1." << std::endl;
    std::cout << "- predistribution:distribution_option. Predistribution of the coordinates to the processors. 0 for NONE, 1 RCB, 2 MJ, 3 BLOCK." << std::endl;
    std::cout << "- perturbation_ratio:the percent of the local data which will be perturbed in order to simulate the changes in the dynamic partitioning. Float value between 0 and 1." << std::endl;
  }

  GeometricGenerator(Teuchos::ParameterList &params, const RCP<const Teuchos::Comm<int> > & comm_){
    this->wd = NULL;
    this->comm = comm_;
    this->holes = NULL; //to represent if there is any hole in the input
    this->coordinate_dimension = 0;  //dimension of the geometry
    this->numWeightsPerCoord = 0;
    this->worldSize = comm_->getSize(); //comminication world object.
    this->numGlobalCoords = 0;	//global number of coordinates requested to be created.
    this->loadDistributions = NULL; //sized as the number of processors, the load of each processor.
    //this->distinctCoordinates = NULL; // if processors have different or same range for coordinates to be created.
    this->coordinateDistributions = NULL;
    this->holeCount = 0;
    this->distributionCount = 0;
    this->outfile = "";
    this->predistribution = 0;
    this->perturbation_ratio = 0;
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
    this->myRank = comm_->getRank();

    try {
      this->parseParams(params);
    }
    catch(std::string &s){
      if(myRank == 0){
        print_description();
      }
      throw std::runtime_error(s);
    }

    catch(const char *s){
      if(myRank == 0){
        print_description();
      }
      throw std::runtime_error(s);
    }


    lno_t myPointCount = 0;
    this->numGlobalCoords = 0;

    gno_t prefixSum = 0;
    for(int i = 0; i < this->distributionCount; ++i){
      for(int ii = 0; ii < this->worldSize; ++ii){
        lno_t increment  = lno_t (this->coordinateDistributions[i]->numPoints * this->loadDistributions[ii]);
        if (gno_t(ii) < this->coordinateDistributions[i]->numPoints % this->worldSize){
          increment += 1;
        }
        this->numGlobalCoords += increment;
        if(ii < myRank){
          prefixSum += increment;
        }
      }
      myPointCount += lno_t(this->coordinateDistributions[i]->numPoints * this->loadDistributions[myRank]);
      if (gno_t(myRank) < this->coordinateDistributions[i]->numPoints % this->worldSize){
        myPointCount += 1;
      }
    }

    this->coords = new scalar_t *[this->coordinate_dimension];
    for(int i = 0; i < this->coordinate_dimension; ++i){
      this->coords[i] = new scalar_t[myPointCount];
    }

    for (int ii = 0; ii < this->coordinate_dimension; ++ii){
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
      for(lno_t i = 0; i < myPointCount; ++i){
        this->coords[ii][i] = 0;
      }
    }

    this->numLocalCoords = 0;
    srand ((myRank + 1) * this->numLocalCoords);
    for (int i = 0; i < distributionCount; ++i){

      lno_t requestedPointCount = lno_t(this->coordinateDistributions[i]->numPoints *  this->loadDistributions[myRank]);
      if (gno_t(myRank) < this->coordinateDistributions[i]->numPoints % this->worldSize){
        requestedPointCount += 1;
      }
      //std::cout << "req:" << requestedPointCount << std::endl;
      //this->coordinateDistributions[i]->GetPoints(requestedPointCount,this->points + this->numLocalCoords, this->holes, this->holeCount,  this->loadDistributions, myRank);
      this->coordinateDistributions[i]->GetPoints(requestedPointCount,this->coords, this->numLocalCoords, this->holes, this->holeCount,  this->loadDistributions, myRank);
      this->numLocalCoords += requestedPointCount;
    }

    /*

    if (this->myRank >= 0){
        for(lno_t i = 0; i < this->numLocalCoords; ++i){

          std::cout <<"me:" << this->myRank << " "<< this->coords[0][i];
          if(this->coordinate_dimension > 1){
        	  std::cout << " " << this->coords[1][i];
          }
          if(this->coordinate_dimension > 2){
        	  std::cout  << " " << this->coords[2][i];
          }
          std::cout << std::endl;
        }
    }
	*/



    if (this->predistribution){
    	redistribute();
    }



    int scale = 3;
    if (this->perturbation_ratio > 0){
    	this->perturb_data(this->perturbation_ratio, scale);
    }
    /*
    if (this->myRank >= 0){
    	std::cout << "after distribution" << std::endl;
        for(lno_t i = 0; i < this->numLocalCoords; ++i){

          std::cout <<"me:" << this->myRank << " " << this->coords[0][i];
          if(this->coordinate_dimension > 1){
        	  std::cout << " " << this->coords[1][i];
          }
          if(this->coordinate_dimension > 2){
        	  std::cout  << " " << this->coords[2][i];
          }
          std::cout << std::endl;
        }
    }

    */


    if (this->distinctCoordSet){
      //TODO: Partition and migration.
    }

    if(this->outfile != ""){

      std::ofstream myfile;
      myfile.open ((this->outfile + Teuchos::toString<int>(myRank)).c_str());
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

      if (this->myRank == 0){
    	  std::ofstream gnuplotfile("gnu.gnuplot");
    	  for(int i = 0; i < this->worldSize; ++i){
    		  string s = "splot";
    		  if (this->coordinate_dimension == 2){
    			  s = "plot";
    		  }
    		  if (i > 0){
    			  s = "replot";
    		  }
    		  gnuplotfile << s << " \"" << (this->outfile + Teuchos::toString<int>(i)) << "\"" << std::endl;
    	  }
    	  gnuplotfile  << "pause -1" << std::endl;
    	  gnuplotfile.close();
      }
    }



    /*
		Zoltan2::XpetraMultiVectorAdapter < Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> > xmv (RCP < Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> > (tmVector));

		RCP< Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> >tmVector2;
		Zoltan2::PartitioningSolution< Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> > solution;
		xmv.applyPartitioningSolution<Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> >(this->tmVector, &tmVector2, solution);
     */
    if (this->numWeightsPerCoord > 0){
    	this->wghts = new scalar_t *[this->numWeightsPerCoord];
    	for(int i = 0; i < this->numWeightsPerCoord; ++i){
    		this->wghts[i] = new scalar_t[this->numLocalCoords];
    	}
    }

    for(int ii = 0; ii < this->numWeightsPerCoord; ++ii){
      switch(this->coordinate_dimension){
      case 1:
 	{
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for (lno_t i = 0; i < this->numLocalCoords; ++i){
          this->wghts[ii][i] = this->wd[ii]->get1DWeight(this->coords[0][i]);
        }
	}
        break;
      case 2:
	{
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for (lno_t i = 0; i < this->numLocalCoords; ++i){
          this->wghts[ii][i] = this->wd[ii]->get2DWeight(this->coords[0][i], this->coords[1][i]);
        }
	}
        break;
      case 3:
	{
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
        for (lno_t i = 0; i < this->numLocalCoords; ++i){
          this->wghts[ii][i] = this->wd[ii]->get3DWeight(this->coords[0][i], this->coords[1][i], this->coords[2][i]);
        }
	}
        break;
      }
    }
  }

  //############################################################//
  ///######Start perturbation function##########################//
  //############################################################//
  void perturb_data(float used_perturbation_ratio, int scale){
	  scalar_t *maxCoords= new scalar_t [this->coordinate_dimension];
	  scalar_t *minCoords= new scalar_t [this->coordinate_dimension];
	  for (int dim = 0; dim < this->coordinate_dimension; ++dim){
		  minCoords[dim] = maxCoords[dim] = this->coords[dim][0];
		  for (lno_t i = 1; i < this->numLocalCoords; ++i){
			  if (minCoords[dim] > this->coords[dim][i]){
				  minCoords[dim] = this->coords[dim][i];
			  }

			  if (maxCoords[dim] < this->coords[dim][i]){
				  maxCoords[dim] = this->coords[dim][i];
			  }
		  }




		  scalar_t center = (maxCoords[dim] + minCoords[dim]) / 2;

		  minCoords[dim] = center - (center - minCoords[dim]) * scale;
		  maxCoords[dim] = (maxCoords[dim] - center) * scale + center;

	  }

	  gno_t numLocalToPerturb = gno_t(this->numLocalCoords*used_perturbation_ratio);
	  //std::cout << "numLocalToPerturb :" << numLocalToPerturb  << std::endl;
	  for (int dim = 0; dim < this->coordinate_dimension; ++dim){
		  scalar_t range = maxCoords[dim] - minCoords[dim];
		  for (gno_t i = 0; i < numLocalToPerturb; ++i){
			  this->coords[dim][i] = (rand() / double (RAND_MAX)) * (range) +  minCoords[dim];

		  }
	  }
	  delete []maxCoords;
	  delete []minCoords;
  }

  //############################################################//
  ///######Start Predistribution functions######################//
  //############################################################//

  //Returns the partitioning dimension as even as possible
  void getBestSurface (int remaining, int *dimProcs, int dim, int currentDim, int &bestSurface, int *bestDimProcs){

	  if (currentDim < dim - 1){
		  int ipx = 1;
		  while(ipx <= remaining) {
			  if(remaining % ipx == 0) {
				  int nremain = remaining / ipx;
				  dimProcs[currentDim] = ipx;
				  getBestSurface (nremain, dimProcs, dim, currentDim + 1, bestSurface, bestDimProcs);
			  }
			  ipx++;
		  }
	  }
	  else {
		  dimProcs[currentDim] = remaining;
		  int surface = 0;
		  for (int i = 0; i < dim; ++i) surface += dimProcs[i];
		  if (surface < bestSurface){
			  bestSurface = surface;
			  for (int i = 0; i < dim; ++i) bestDimProcs[i] = dimProcs[i];
		  }
	  }

  }

  //returns min and max coordinates along each dimension
  void getMinMaxCoords(scalar_t *globalMaxCoords, scalar_t *globalMinCoords){
	  scalar_t *maxCoords= new scalar_t [this->coordinate_dimension];
	  scalar_t *minCoords= new scalar_t [this->coordinate_dimension];
	  for (int dim = 0; dim < this->coordinate_dimension; ++dim){
		  minCoords[dim] = maxCoords[dim] = this->coords[dim][0];
		  for (lno_t i = 1; i < this->numLocalCoords; ++i){
			  if (minCoords[dim] > this->coords[dim][i]){
				  minCoords[dim] = this->coords[dim][i];
			  }

			  if (maxCoords[dim] < this->coords[dim][i]){
				  maxCoords[dim] = this->coords[dim][i];
			  }
		  }
	  }

	  reduceAll<int, scalar_t>( *(this->comm), Teuchos::REDUCE_MAX,
			  	  	  	  	  	  this->coordinate_dimension,
			  	  	  	  	  	  maxCoords,
			  	  	  	  	  	  globalMaxCoords);


	  reduceAll<int, scalar_t>( *(this->comm), Teuchos::REDUCE_MIN,
	  			  	  	  	  	  	  this->coordinate_dimension,
	  			  	  	  	  	  	  minCoords,
	  			  	  	  	  	  	  globalMinCoords);

	  delete []maxCoords;
	  delete []minCoords;
  }


  //performs a block partitioning.
  //then distributes the points of the overloaded parts to underloaded parts.
  void blockPartition(int *coordinate_grid_parts){

#ifdef _MSC_VER
	  typedef SSIZE_T ssize_t;
#endif

	  //############################################################//
	  ///getting minimum and maximum coordinates for each dimension///
	  //############################################################//

	  scalar_t *maxCoords= new scalar_t [this->coordinate_dimension];
	  scalar_t *minCoords= new scalar_t [this->coordinate_dimension];
	  //global min and max coordinates in each dimension.
	  this->getMinMaxCoords(maxCoords, minCoords);


	  //############################################################//
	  ///getting the best partitioning number along each dimension ///
	  //############################################################//
	  int remaining = this->worldSize;
	  int coord_dim = this->coordinate_dimension;
	  int *dimProcs = new int[coord_dim];
	  int *bestDimProcs = new int[coord_dim];
	  int currentDim = 0;

	  int bestSurface = 0;
	  dimProcs[0] = remaining;
	  for (int i = 1; i < coord_dim; ++i) dimProcs[i] = 1;
	  for (int i = 0; i < coord_dim; ++i) bestSurface += dimProcs[i];
	  for (int i = 0; i < coord_dim; ++i) bestDimProcs[i] = dimProcs[i];
	  //divides the parts into dimensions as even as possible.
	  getBestSurface ( remaining, dimProcs,  coord_dim,  currentDim, bestSurface, bestDimProcs);


	  delete []dimProcs;

	  //############################################################//
	  ///getting the size of a slice along each dimension ///
	  //############################################################//
	  int *shiftProcCount = new int[coord_dim];
	  //how the consecutive parts along a dimension
	  //differs in the index.

	  int remainingProc = this->worldSize;
	  for (int dim = 0; dim < coord_dim; ++dim){
		  remainingProc = remainingProc /  bestDimProcs[dim];
		  shiftProcCount[dim] = remainingProc;
	  }

	  scalar_t *dim_slices = new scalar_t[coord_dim];
	  for (int dim = 0; dim < coord_dim; ++dim){
		  scalar_t dim_range = maxCoords[dim] - minCoords[dim];
		  scalar_t slice = dim_range / bestDimProcs[dim];
		  dim_slices[dim] = slice;
	  }

	  //############################################################//
	  ///##################Initial part assignments ###############///
	  //############################################################//

	  gno_t *numPointsInParts = new gno_t[this->worldSize];
	  gno_t *numGlobalPointsInParts = new gno_t[this->worldSize];
	  gno_t *numPointsInPartsInclusiveUptoMyIndex = new gno_t[this->worldSize];

	  gno_t *partBegins = new gno_t [this->worldSize];
	  gno_t *partNext = new gno_t[this->numLocalCoords];


	  for (lno_t i = 0; i < this->numLocalCoords; ++i){
		  partNext[i] = -1;
	  }
	  for (int i = 0; i < this->worldSize; ++i) {
		  partBegins[i] = 1;
	  }

	  for (int i = 0; i < this->worldSize; ++i)
		  numPointsInParts[i] = 0;

	  for (lno_t i = 0; i < this->numLocalCoords; ++i){
		  int partIndex = 0;
		  for (int dim = 0; dim < coord_dim; ++dim){
			  int shift = int ((this->coords[dim][i] - minCoords[dim]) / dim_slices[dim]);
			  if (shift >= bestDimProcs[dim]){
				  shift = bestDimProcs[dim] - 1;
			  }
			  shift = shift * shiftProcCount[dim];
			  partIndex += shift;
		  }
		  numPointsInParts[partIndex] += 1;
		  coordinate_grid_parts[i] = partIndex;

		  partNext[i] = partBegins[partIndex];
		  partBegins[partIndex] = i;

	  }

	  //############################################################//
	  ///#########Counting the num points in each  part ###########///
	  //############################################################//
	  reduceAll<int, gno_t>( *(this->comm), Teuchos::REDUCE_SUM,
			  	  	  	  	  	  	  this->worldSize,
			  	  	  	  	  	  	  numPointsInParts,
			  	  	  	  	  	  	  numGlobalPointsInParts);


      Teuchos::scan<int,gno_t>(
    		  *(this->comm), Teuchos::REDUCE_SUM,
    		  this->worldSize,
    		  numPointsInParts,
    		  numPointsInPartsInclusiveUptoMyIndex
      );




	  /*
	  gno_t totalSize = 0;
	  for (int i = 0; i < this->worldSize; ++i){
		  totalSize += numPointsInParts[i];
	  }
	  if (totalSize != this->numLocalCoords){
		  std::cout << "me:" << this->myRank << " ts:" << totalSize << " nl:" << this->numLocalCoords << std::endl;
	  }
	  */


      //std::cout << "me:" << this->myRank << " ilk print" << std::endl;

	  gno_t optimal_num = gno_t(this->numGlobalCoords/double(this->worldSize)+0.5);
#ifdef printparts
	  if (this->myRank == 0){
		  gno_t totalSize = 0;
		  for (int i = 0; i < this->worldSize; ++i){
			  std::cout << "me:" << this->myRank << " NumPoints in part:" << i << " is: " << numGlobalPointsInParts[i] << std::endl;
			  totalSize += numGlobalPointsInParts[i];
		  }
		  std::cout << "Total:" << totalSize << " ng:" << this->numGlobalCoords << std::endl;
		  std::cout << "optimal_num:" << optimal_num << std::endl;
	  }
#endif
	  ssize_t *extraInPart = new ssize_t [this->worldSize];

	  ssize_t extraExchange = 0;
	  for (int i = 0; i < this->worldSize; ++i){
		  extraInPart[i] = numGlobalPointsInParts[i] - optimal_num;
		  extraExchange += extraInPart[i];
	  }
	  if (extraExchange != 0){
		  int addition = -1;
		  if (extraExchange < 0) addition = 1;
		  for (ssize_t i = 0; i < extraExchange; ++i){
			  extraInPart[i % this->worldSize] += addition;
		  }
	  }

	  //############################################################//
	  ///######Check the overloaded and underloaded parts #########///
	  //############################################################//

      int overloadedPartCount = 0;
      int *overloadedPartIndices = new int [this->worldSize];


      int underloadedPartCount = 0;
      int *underloadedPartIndices = new int [this->worldSize];

      for (int i = 0; i < this->worldSize; ++i){
    	  if(extraInPart[i] > 0){
    		  overloadedPartIndices[overloadedPartCount++] = i;
    	  }
    	  else if(extraInPart[i] < 0){
    		  underloadedPartIndices[underloadedPartCount++] = i;
    	  }
      }

      int underloadpartindex = underloadedPartCount - 1;


      int numPartsISendFrom = 0;
      int *mySendFromParts = new int[this->worldSize * 2];
      gno_t *mySendFromPartsCounts = new gno_t[this->worldSize * 2];

      int numPartsISendTo = 0;
      int *mySendParts = new int[this->worldSize * 2];
      gno_t *mySendCountsToParts = new gno_t[this->worldSize * 2];


	  //############################################################//
	  ///######Calculating##########################################//
      ///######*which processors ##################################///
      ///######*which overloaded parts elements should be converted///
      ///######*into which underloaded parts elements #############///
	  //############################################################//
      for (int i = overloadedPartCount - 1; i >= 0; --i){
    	  //get the overloaded part
    	  //the overload
    	  int overloadPart = overloadedPartIndices[i];
    	  gno_t overload = extraInPart[overloadPart];
    	  gno_t myload = numPointsInParts[overloadPart];


    	  //the inclusive load of the processors up to me
    	  gno_t inclusiveLoadUpToMe = numPointsInPartsInclusiveUptoMyIndex[overloadPart];

    	  //the exclusive load of the processors up to me
    	  gno_t exclusiveLoadUptoMe = inclusiveLoadUpToMe - myload;


    	  if (exclusiveLoadUptoMe >= overload){
    		  //this processor does not have to convert anything.
    		  //for this overloaded part.
    		  //set the extra for this processor to zero.
    		  overloadedPartIndices[i] = -1;
    		  extraInPart[overloadPart] = 0;
    		  //then consume underloaded parts.
    		  while (overload > 0){
    			  int nextUnderLoadedPart = underloadedPartIndices[underloadpartindex];
    			  ssize_t underload = extraInPart[nextUnderLoadedPart];
    			  ssize_t left = overload + underload;

    			  if(left >= 0){
    				  extraInPart[nextUnderLoadedPart] = 0;
    				  underloadedPartIndices[underloadpartindex--] = -1;
    			  }
    			  else {
    				  extraInPart[nextUnderLoadedPart] = left;
    			  }
    			  overload = left;
    		  }
    	  }
    	  else if (exclusiveLoadUptoMe < overload){
    		  //if the previous processors load is not enough
    		  //then this processor should convert some of its elements.
    		  gno_t mySendCount = overload - exclusiveLoadUptoMe;
    		  //how much more needed.
    		  gno_t sendAfterMe = 0;
    		  //if my load is not enough
    		  //then the next processor will continue to convert.
    		  if (mySendCount > myload){
    			  sendAfterMe = mySendCount - myload;
    			  mySendCount = myload;
    		  }


    		  //this processors will convert from overloaded part
    		  //as many as mySendCount items.
    		  mySendFromParts[numPartsISendFrom] = overloadPart;
    	      mySendFromPartsCounts[numPartsISendFrom++] = mySendCount;

    	      //first consume underloaded parts for the previous processors.
    		  while (exclusiveLoadUptoMe > 0){
    			  int nextUnderLoadedPart = underloadedPartIndices[underloadpartindex];
    			  ssize_t underload = extraInPart[nextUnderLoadedPart];
    			  ssize_t left = exclusiveLoadUptoMe + underload;

    			  if(left >= 0){
    				  extraInPart[nextUnderLoadedPart] = 0;
    				  underloadedPartIndices[underloadpartindex--] = -1;
    			  }
    			  else {
    				  extraInPart[nextUnderLoadedPart] = left;
    			  }
    			  exclusiveLoadUptoMe = left;
    		  }

    		  //consume underloaded parts for my load.
    		  while (mySendCount > 0){
    			  int nextUnderLoadedPart = underloadedPartIndices[underloadpartindex];
    			  ssize_t underload = extraInPart[nextUnderLoadedPart];
    			  ssize_t left = mySendCount + underload;

    			  if(left >= 0){
    				  mySendParts[numPartsISendTo] = nextUnderLoadedPart;
    				  mySendCountsToParts[numPartsISendTo++] = -underload;

    				  extraInPart[nextUnderLoadedPart] = 0;
    				  underloadedPartIndices[underloadpartindex--] = -1;
    			  }
    			  else {
    				  extraInPart[nextUnderLoadedPart] = left;

    				  mySendParts[numPartsISendTo] = nextUnderLoadedPart;
    				  mySendCountsToParts[numPartsISendTo++] = mySendCount;

    			  }
    			  mySendCount = left;
    		  }
    		  //consume underloaded parts for the load of the processors after my index.
    		  while (sendAfterMe > 0){
    			  int nextUnderLoadedPart = underloadedPartIndices[underloadpartindex];
    			  ssize_t underload = extraInPart[nextUnderLoadedPart];
    			  ssize_t left = sendAfterMe + underload;

    			  if(left >= 0){
    				  extraInPart[nextUnderLoadedPart] = 0;
    				  underloadedPartIndices[underloadpartindex--] = -1;
    			  }
    			  else {
    				  extraInPart[nextUnderLoadedPart] = left;
    			  }
    			  sendAfterMe = left;
    		  }
    	  }
      }


	  //############################################################//
	  ///######Perform actual conversion############################//
	  //############################################################//
      for (int i = 0 ; i < numPartsISendFrom; ++i){

    	  //get the part from which the elements will be converted.
    	  int sendFromPart = mySendFromParts[i];
    	  ssize_t sendCount = mySendFromPartsCounts[i];
    	  while(sendCount > 0){
    		  int partToSendIndex = numPartsISendTo - 1;
    		  int partToSend = mySendParts[partToSendIndex];

    		  int sendCountToThisPart = mySendCountsToParts[partToSendIndex];

    		  //determine which part i should convert to
    		  //and how many to this part.
    		  if (sendCountToThisPart <= sendCount){
    			  mySendParts[partToSendIndex] = 0;
    			  mySendCountsToParts[partToSendIndex] = 0;
    			  --numPartsISendTo;
    			  sendCount -= sendCountToThisPart;
    		  }
    		  else {
    			  mySendCountsToParts[partToSendIndex] = sendCountToThisPart - sendCount;
    			  sendCountToThisPart = sendCount;
    			  sendCount = 0;
    		  }


        	  gno_t toChange = partBegins[sendFromPart];
    		  gno_t previous_begin = partBegins[partToSend];

    		  //do the conversion.
    		  for (int k = 0; k < sendCountToThisPart - 1; ++k){
    			  coordinate_grid_parts[toChange] = partToSend;
    			  toChange = partNext[toChange];
    		  }
    		  coordinate_grid_parts[toChange] = partToSend;

    		  gno_t newBegin = partNext[toChange];
    		  partNext[toChange] = previous_begin;
    		  partBegins[partToSend] = partBegins[sendFromPart];
    		  partBegins[sendFromPart] = newBegin;
    	  }
      }

      //if (this->myRank == 0) std::cout << "4" << std::endl;

#ifdef printparts


      for (int i = 0; i < this->worldSize; ++i) numPointsInParts[i] = 0;

      for (int i = 0; i < this->numLocalCoords; ++i){
    	  numPointsInParts[coordinate_grid_parts[i]] += 1;
      }

	  reduceAll<int, gno_t>( *(this->comm), Teuchos::REDUCE_SUM,
			  	  	  	  	  	  	  this->worldSize,
			  	  	  	  	  	  	  numPointsInParts,
			  	  	  	  	  	  	  numGlobalPointsInParts);
	  if (this->myRank == 0){
		  std::cout << "reassigning" << std::endl;
		  gno_t totalSize = 0;
		  for (int i = 0; i < this->worldSize; ++i){
			  std::cout << "me:" << this->myRank << " NumPoints in part:" << i << " is: " << numGlobalPointsInParts[i] << std::endl;
			  totalSize += numGlobalPointsInParts[i];
		  }
		  std::cout << "Total:" << totalSize << " ng:" << this->numGlobalCoords << std::endl;
	  }
#endif
	  delete []mySendCountsToParts;
	  delete []mySendParts;
	  delete []mySendFromPartsCounts;
	  delete []mySendFromParts;
	  delete []underloadedPartIndices;
	  delete []overloadedPartIndices;
	  delete []extraInPart;
	  delete []partNext;
	  delete []partBegins;
	  delete []numPointsInPartsInclusiveUptoMyIndex;
	  delete []numPointsInParts;
	  delete []numGlobalPointsInParts;

	  delete []shiftProcCount;
	  delete []bestDimProcs;
	  delete []dim_slices;
      delete []minCoords;
      delete []maxCoords;
  }

  //given the part numbers for each local coordinate,
  //distributes the coordinates to the corresponding processors.
  void distribute_points(int *coordinate_grid_parts){

	  Tpetra::Distributor distributor(comm);
	  ArrayView<const int> pIds( coordinate_grid_parts, this->numLocalCoords);
	  gno_t numMyNewGnos = distributor.createFromSends(pIds);


          Kokkos::View<scalar_t*, Kokkos::HostSpace> recvBuf2(
            Kokkos::ViewAllocateWithoutInitializing("recvBuf2"),
            numMyNewGnos);

	  for (int i = 0; i < this->coordinate_dimension; ++i){
              Kokkos::View<scalar_t*, Kokkos::HostSpace> s;
              if (this->numLocalCoords > 0) 
                s = Kokkos::View<scalar_t *, Kokkos::HostSpace>(
                            this->coords[i], this->numLocalCoords); //unmanaged

	      distributor.doPostsAndWaits(s, 1, recvBuf2);

	      delete [] this->coords[i];
	      this->coords[i] = new scalar_t[numMyNewGnos];
	      for (lno_t j = 0; j < numMyNewGnos; ++j){
	    	  this->coords[i][j] = recvBuf2[j];
	      }

	  }
	  this->numLocalCoords = numMyNewGnos;
  }

  //calls MJ for p = numProcs
  int predistributeMJ(int *coordinate_grid_parts){
	  int coord_dim = this->coordinate_dimension;

	  lno_t numLocalPoints = this->numLocalCoords;
	  gno_t numGlobalPoints = this->numGlobalCoords;


	  //T **weight = NULL;
	  //typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tMVector_t;
	  RCP<Tpetra::Map<lno_t, gno_t, node_t> > mp = rcp(
			  new Tpetra::Map<lno_t, gno_t, node_t> (numGlobalPoints, numLocalPoints, 0, comm));

	  Teuchos::Array<Teuchos::ArrayView<const scalar_t> > coordView(coord_dim);



	  for (int i=0; i < coord_dim; i++){
		  if(numLocalPoints > 0){
			  Teuchos::ArrayView<const scalar_t> a(coords[i], numLocalPoints);
			  coordView[i] = a;
		  } else{
			  Teuchos::ArrayView<const scalar_t> a;
			  coordView[i] = a;
		  }
	  }

	  RCP< Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> >tmVector = RCP< Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> >(
			  new Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t>( mp, coordView.view(0, coord_dim), coord_dim));


	  RCP<const tMVector_t> coordsConst = Teuchos::rcp_const_cast<const tMVector_t>(tmVector);
	  std::vector<const scalar_t *> weights;
	  std::vector <int> stride;

	  typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> inputAdapter_t;
	  //inputAdapter_t ia(coordsConst);
	  inputAdapter_t ia(coordsConst,weights, stride);

	  Teuchos::RCP <Teuchos::ParameterList> params ;
	  params =RCP <Teuchos::ParameterList> (new Teuchos::ParameterList, true);


	  params->set("algorithm", "multijagged");
	  params->set("num_global_parts", this->worldSize);

	  //TODO we need to fix the setting parts.
	  //Although MJ sets the parts with
	  //currently the part setting is not correct when migration is done.
	  //params->set("migration_check_option", 2);


	  Zoltan2::PartitioningProblem<inputAdapter_t> *problem;


	  try {
#ifdef HAVE_ZOLTAN2_MPI
		  problem = new Zoltan2::PartitioningProblem<inputAdapter_t>(&ia, params.getRawPtr(),
				  MPI_COMM_WORLD);
#else
		  problem = new Zoltan2::PartitioningProblem<inputAdapter_t>(&ia, params.getRawPtr());
#endif
	  }
	  CATCH_EXCEPTIONS("PartitioningProblem()")

	  try {
		  problem->solve();
	  }
	  CATCH_EXCEPTIONS("solve()")

	  const typename inputAdapter_t::part_t *partIds = problem->getSolution().getPartListView();

	  for (lno_t i = 0; i < this->numLocalCoords;++i){
		  coordinate_grid_parts[i] = partIds[i];
		  //std::cout << "me:" << this->myRank << " i:" << i << " goes to:" << partIds[i] << std::endl;
	  }
	  delete problem;
	  return 0;
  }

  //calls RCP for p = numProcs
  int predistributeRCB(int *coordinate_grid_parts){
	  int rank = this->myRank;
	  int nprocs = this->worldSize;
	  DOTS<tMVector_t> dots_;

	  MEMORY_CHECK(rank==0 || rank==nprocs-1, "After initializing MPI");


	  int nWeights = 0;
	  int debugLevel=0;
	  string memoryOn("memoryOn");
	  string memoryOff("memoryOff");
	  bool doMemory=false;
	  int numGlobalParts = nprocs;
	  int dummyTimer=0;
	  bool remap=0;

	  string balanceCount("balance_object_count");
	  string balanceWeight("balance_object_weight");
	  string mcnorm1("multicriteria_minimize_total_weight");
	  string mcnorm2("multicriteria_balance_total_maximum");
	  string mcnorm3("multicriteria_minimize_maximum_weight");
	  string objective(balanceWeight);   // default

	  // Process command line input
	  CommandLineProcessor commandLine(false, true);
	  //commandLine.setOption("size", &numGlobalCoords,
	  //  "Approximate number of global coordinates.");
	  int input_option = 0;
	  commandLine.setOption("input_option", &input_option,
	    "whether to use mesh creation, geometric generator, or file input");
	  string inputFile = "";

	  commandLine.setOption("input_file", &inputFile,
	    "the input file for geometric generator or file input");


	  commandLine.setOption("size", &numGlobalCoords,
	    "Approximate number of global coordinates.");
	  commandLine.setOption("numParts", &numGlobalParts,
	    "Number of parts (default is one per proc).");
	  commandLine.setOption("nWeights", &nWeights,
	    "Number of weights per coordinate, zero implies uniform weights.");
	  commandLine.setOption("debug", &debugLevel, "Zoltan1 debug level");
	  commandLine.setOption("remap", "no-remap", &remap,
	    "Zoltan1 REMAP parameter; disabled by default for scalability testing");
	  commandLine.setOption("timers", &dummyTimer, "ignored");
	  commandLine.setOption(memoryOn.c_str(), memoryOff.c_str(), &doMemory,
	    "do memory profiling");

	  string doc(balanceCount);
	  doc.append(": ignore weights\n");
	  doc.append(balanceWeight);
	  doc.append(": balance on first weight\n");
	  doc.append(mcnorm1);
	  doc.append(": given multiple weights, balance their total.\n");
	  doc.append(mcnorm3);
	  doc.append(": given multiple weights, "
	             "balance the maximum for each coordinate.\n");
	  doc.append(mcnorm2);
	  doc.append(": given multiple weights, balance the L2 norm of the weights.\n");
	  commandLine.setOption("objective", &objective,  doc.c_str());

	  CommandLineProcessor::EParseCommandLineReturn rc =
	    commandLine.parse(0, NULL);



	  if (rc != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
	    if (rc == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
	      if (rank==0) std::cout << "PASS" << std::endl;
	      return 1;
	    }
	    else {
	      if (rank==0) std::cout << "FAIL" << std::endl;
	      return 0;
	    }
	  }

	  //MEMORY_CHECK(doMemory && rank==0, "After processing parameters");

	  // Create the data structure

	  int coord_dim = this->coordinate_dimension;


	  RCP<Tpetra::Map<lno_t, gno_t, node_t> > mp = rcp(
			  new Tpetra::Map<lno_t, gno_t, node_t> (this->numGlobalCoords, this->numLocalCoords, 0, this->comm));

	  Teuchos::Array<Teuchos::ArrayView<const scalar_t> > coordView(coord_dim);
	  for (int i=0; i < coord_dim; i++){
		  if(numLocalCoords > 0){
			  Teuchos::ArrayView<const scalar_t> a(coords[i], numLocalCoords);
			  coordView[i] = a;
		  } else{
			  Teuchos::ArrayView<const scalar_t> a;
			  coordView[i] = a;
		  }
	  }

	  tMVector_t *tmVector = new tMVector_t( mp, coordView.view(0, coord_dim), coord_dim);

	  dots_.coordinates = tmVector;
	  dots_.weights.resize(nWeights);


	  MEMORY_CHECK(doMemory && rank==0, "After creating input");

	  // Now call Zoltan to partition the problem.

	  float ver;
	  int aok = Zoltan_Initialize(0,NULL, &ver);

	  if (aok != 0){
	    printf("Zoltan_Initialize failed\n");
	    exit(0);
	  }

	  struct Zoltan_Struct *zz;
	  zz = Zoltan_Create(MPI_COMM_WORLD);

	  Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
	  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
	  Zoltan_Set_Param(zz, "CHECK_GEOM", "0");
	  Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
	  Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "0");
	  Zoltan_Set_Param(zz, "RETURN_LISTS", "PART");
	  std::ostringstream oss;
	  oss << numGlobalParts;
	  Zoltan_Set_Param(zz, "NUM_GLOBAL_PARTS", oss.str().c_str());
	  oss.str("");
	  oss << debugLevel;
	  Zoltan_Set_Param(zz, "DEBUG_LEVEL", oss.str().c_str());

	  if (remap)
	    Zoltan_Set_Param(zz, "REMAP", "1");
	  else
	    Zoltan_Set_Param(zz, "REMAP", "0");

	  if (objective != balanceCount){
	    oss.str("");
	    oss << nWeights;
	    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", oss.str().c_str());

	    if (objective == mcnorm1)
	      Zoltan_Set_Param(zz, "RCB_MULTICRITERIA_NORM", "1");
	    else if (objective == mcnorm2)
	      Zoltan_Set_Param(zz, "RCB_MULTICRITERIA_NORM", "2");
	    else if (objective == mcnorm3)
	      Zoltan_Set_Param(zz, "RCB_MULTICRITERIA_NORM", "3");
	  }
	  else{
	    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
	  }

	  Zoltan_Set_Num_Obj_Fn(zz, getNumObj<tMVector_t>, &dots_);
	  Zoltan_Set_Obj_List_Fn(zz, getObjList<tMVector_t>, &dots_);
	  Zoltan_Set_Num_Geom_Fn(zz,  getDim<tMVector_t>, &dots_);
	  Zoltan_Set_Geom_Multi_Fn(zz, getCoords<tMVector_t>, &dots_);

	  int changes, numGidEntries, numLidEntries, numImport, numExport;
	  ZOLTAN_ID_PTR importGlobalGids, importLocalGids;
	  ZOLTAN_ID_PTR exportGlobalGids, exportLocalGids;
	  int *importProcs, *importToPart, *exportProcs, *exportToPart;

	  MEMORY_CHECK(doMemory && rank==0, "Before Zoltan_LB_Partition");


	  aok = Zoltan_LB_Partition(zz, &changes, &numGidEntries, &numLidEntries,
	                            &numImport, &importGlobalGids, &importLocalGids,
	                            &importProcs, &importToPart,
	                            &numExport, &exportGlobalGids, &exportLocalGids,
	                            &exportProcs, &exportToPart);


	  MEMORY_CHECK(doMemory && rank==0, "After Zoltan_LB_Partition");

	  for (lno_t i = 0; i < numLocalCoords; i++)
		  coordinate_grid_parts[i] = exportToPart[i];
	  Zoltan_Destroy(&zz);
	  MEMORY_CHECK(doMemory && rank==0, "After Zoltan_Destroy");

	  delete dots_.coordinates;
	  return 0;
}
  void redistribute(){
	  int *coordinate_grid_parts = new int[this->numLocalCoords];
	  switch (this->predistribution){
	  case 1:
		  this->predistributeRCB(coordinate_grid_parts);
		  break;
	  case 2:

		  this->predistributeMJ(coordinate_grid_parts);
		  break;
	  case 3:
		  //block
		  blockPartition(coordinate_grid_parts);
		  break;
	  }
	  this->distribute_points(coordinate_grid_parts);

	  delete []coordinate_grid_parts;


  }

  //############################################################//
  ///########END Predistribution functions######################//
  //############################################################//


  int getNumWeights(){
    return this->numWeightsPerCoord;
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

  scalar_t **getLocalCoordinatesView(){
    return this->coords;
  }

  scalar_t **getLocalWeightsView(){
    return this->wghts;
  }

  void getLocalCoordinatesCopy( scalar_t ** c){
    for(int ii = 0; ii < this->coordinate_dimension; ++ii){
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
      for (lno_t i = 0; i < this->numLocalCoords; ++i){
        c[ii][i] = this->coords[ii][i];
      }
    }
  }

  void getLocalWeightsCopy(scalar_t **w){
    for(int ii = 0; ii < this->numWeightsPerCoord; ++ii){
#ifdef HAVE_ZOLTAN2_OMP
#pragma omp parallel for
#endif
      for (lno_t i = 0; i < this->numLocalCoords; ++i){
        w[ii][i] = this->wghts[ii][i];
      }
    }
  }
};
}

#endif /* GEOMETRICGENERATOR */
