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

#define HOLE_ALLOC_STEP = 10
#define INVALID(STR) "Invalid argument at " + STR
#define INVALIDSHAPE(STR, DIM) "Invalid shape name " + STR + " for " + 	DIM + ".\nValid shapes are \"SQUARE\", \"RECTANGLE\", \"CIRCLE\" for 2D, and \"CUBE\", \"RECTANGULAR_PRISM\", \"SPHERE\" for 3D"

#define INVALID_SHAPE_ARG(SHAPE, REQUIRED) "Invalid argument count for shape " + SHAPE + ". Requires " + REQUIRED + " argument(s)."
#define MAX_ITER_ALLOWED 500

const std::string equation_with_step_function = "STEPPEDEQUATION";
const std::string equation_with_step_function_parameters = equation_with_step_function + "-";
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
    T expressionRes = this->a1 * pow( (x - this->x1), b1);
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
    T expressionRes = this->a1 * pow( (x - this->x1), b1) + this->a2 * pow( (y - this->y1), b2);
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

template <typename T, typename lno_t>
class CoordinateDistribution{
public:
  lno_t numPoints;
  int dimension;
  lno_t requested;
  lno_t assignedPrevious;
  virtual ~CoordinateDistribution(){}

  CoordinateDistribution(gno_t np_, int dim):numPoints(np_), dimension(dim), requested(0), assignedPrevious(0){}
  virtual CoordinatePoint<T> getPoint() = 0;
  virtual T getXCenter() = 0;
  virtual T getXRadius() =0;

  void GetPoints(lno_t requestedPointcount, CoordinatePoint<T> *points /*preallocated sized numPoints*/, Hole<T> **holes, lno_t holeCount, float *sharedRatios_, int myRank){

    for (int i = 0; i < myRank; ++i){
      //cout << "me:" << myRank << " i:" << i << " s:" << sharedRatios_[i]<< endl;
      this->assignedPrevious += sharedRatios_[i] * this->numPoints;
    }

    this->requested = requestedPointcount;

    int cnt = 0;
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

};

template <typename T, typename lno_t>
class CoordinateNormalDistribution:public CoordinateDistribution<T,lno_t>{
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

  CoordinateNormalDistribution(lno_t np_, int dim, CoordinatePoint<T> center_ , T sd_x, T sd_y, T sd_z ): CoordinateDistribution<T,lno_t>(np_,dim),
      standartDevx(sd_x), standartDevy(sd_y), standartDevz(sd_z){
    this->center.x = center_.x;
    this->center.y = center_.y;
    this->center.z = center_.z;
  }

  virtual CoordinatePoint<T> getPoint(){
    CoordinatePoint <T> p;
    for(int i = 0; i < this->dimension; ++i){
      switch(i){
      case 0:
        p.x = normalDist(this->center.x, this->standartDevx);
        break;
      case 1:
        p.y = normalDist(this->center.y, this->standartDevy);
        break;
      case 2:
        p.z = normalDist(this->center.z, this->standartDevz);
        break;
      default:
        throw "unsupported dimension";
      }
    }
    return p;
  }

  virtual ~CoordinateNormalDistribution(){};
private:
  T normalDist(T mu, T sigma) {
    static bool deviateAvailable=false;        //        flag
    static T storedDeviate;                        //        deviate from previous calculation
    T polar, rsquared, var1, var2;

    //        If no deviate has been stored, the polar Box-Muller transformation is
    //        performed, producing two independent normally-distributed random
    //        deviates.  One is stored for the next round, and one is returned.
    if (!deviateAvailable) {

      //        choose pairs of uniformly distributed deviates, discarding those
      //        that don't fall within the unit circle
      do {
        var1=2.0*( T(rand())/T(RAND_MAX) ) - 1.0;
        var2=2.0*( T(rand())/T(RAND_MAX) ) - 1.0;
        rsquared=var1*var1+var2*var2;
      } while ( rsquared>=1.0 || rsquared == 0.0);

      //        calculate polar tranformation for each deviate
      polar=sqrt(-2.0*log(rsquared)/rsquared);

      //        store first deviate and set flag
      storedDeviate=var1*polar;
      deviateAvailable=true;

      //        return second deviate
      return var2*polar*sigma + mu;
    }

    //        If a deviate is available from a previous call to this function, it is
    //        returned, and the flag is set to false.
    else {
      deviateAvailable=false;
      return storedDeviate*sigma + mu;
    }
  }
};

template <typename T, typename lno_t>
class CoordinateUniformDistribution:public CoordinateDistribution<T,lno_t>{
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


  CoordinateUniformDistribution(lno_t np_, int dim, T l_x, T r_x, T l_y, T r_y, T l_z, T r_z ): CoordinateDistribution<T,lno_t>(np_,dim),
      leftMostx(l_x), rightMostx(r_x), leftMosty(l_y), rightMosty(r_y), leftMostz(l_z), rightMostz(r_z){}

  virtual ~CoordinateUniformDistribution(){};
  virtual CoordinatePoint<T> getPoint(){


    CoordinatePoint <T> p;
    for(int i = 0; i < this->dimension; ++i){
      switch(i){
      case 0:
        p.x = unifRand(this->leftMostx, this->rightMostx);
        break;
      case 1:
        p.y = unifRand(this->leftMosty, this->rightMosty);
        break;
      case 2:
        p.z = unifRand(this->leftMostz, this->rightMostz);
        break;
      default:
        throw "unsupported dimension";
      }
    }
    return p;
  }

private:
  //
  // Generate a random number between 0 and 1
  // return a uniform number in [0,1].
  double unifRand()
  {
    return rand() / double(RAND_MAX);
  }
  //
  // Generate a random number in a real interval.
  // param a one end point of the interval
  // param b the other end of the interval
  // return a inform rand numberin [a,b].
  T unifRand(T a, T b)
  {
    return (b-a)*unifRand() + a;
  }
};

template <typename T, typename lno_t>
class CoordinateGridDistribution:public CoordinateDistribution<T,lno_t>{
public:
  T leftMostx;
  T rightMostx;
  T leftMosty;
  T rightMosty;
  T leftMostz;
  T rightMostz;
  lno_t along_X, along_Y, along_Z;
  //T currentX, currentY, currentZ;
  T processCnt;
  int myRank;
  T xstep, ystep, zstep;


  virtual T getXCenter(){
    return (rightMostx - leftMostx)/2  + leftMostx;
  }
  virtual T getXRadius(){
    return (rightMostx - leftMostx)/2;
  }


  CoordinateGridDistribution(lno_t alongX, lno_t alongY, lno_t alongZ, int dim, T l_x, T r_x, T l_y, T r_y, T l_z, T r_z , int myRank_): CoordinateDistribution<T,lno_t>(alongX * alongY * alongZ,dim),
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

    //cout << xstep << " " << ystep << " " << zstep << endl;

  }

  virtual ~CoordinateGridDistribution(){};
  virtual CoordinatePoint<T> getPoint(){
    lno_t before = processCnt + this->assignedPrevious;
    //cout << "before:" << processCnt << " " << this->assignedPrevious << endl;
    lno_t xshift = 0, yshift = 0, zshift = 0;

    lno_t tmp = before % (this->along_X * this->along_Y);
    yshift  = tmp / this->along_X;
    xshift = tmp % this->along_X;
    zshift = before / (this->along_X * this->along_Y);


    CoordinatePoint <T> p;
    p.x = xshift * this->xstep + leftMostx;
    p.y = yshift * this->ystep + leftMosty;
    p.z = zshift * this->zstep + leftMostz;

    if(this->requested - 1 > this->processCnt)
      this->processCnt++;
    return p;
  }

private:
  //
  // Generate a random number between 0 and 1
  // return a uniform number in [0,1].
  double unifRand()
  {
    return rand() / double(RAND_MAX);
  }
  //
  // Generate a random number in a real interval.
  // param a one end point of the interval
  // param b the other end of the interval
  // return a inform rand numberin [a,b].
  T unifRand(T a, T b)
  {
    return (b-a)*unifRand() + a;
  }
};

template <typename T, typename lno_t, typename gno_t, typename node_t>
class GeometricGenerator {
private:
  Hole<T> **holes; //to represent if there is any hole in the input
  int holeCount;
  int dimension;  //dimension of the geometry
  gno_t numGlobalCoords;	//global number of coordinates requested to be created.
  lno_t numLocalCoords;
  float *loadDistributions; //sized as the number of processors, the load of each processor.
  bool loadDistSet;
  bool distinctCoordSet;
  CoordinateDistribution<T, lno_t> **coordinateDistributions;
  int distributionCount;
  CoordinatePoint<T> *points;

  WeightDistribution<T,T> *wd;

  RCP< Tpetra::MultiVector<T, lno_t, gno_t, node_t> >tmVector;
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
      throw "Cannot convert string " + obj + ".";
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
    coordinateDistributions = (CoordinateDistribution<T, lno_t> **) malloc(sizeof (CoordinateDistribution<T, lno_t> *) * 1);
    for(int i = 0; i < argCnt; ){
      coordinateDistributions = (CoordinateDistribution<T, lno_t> **)realloc((void *)coordinateDistributions, (this->distributionCount + 1)* sizeof(CoordinateDistribution<T, lno_t> *));

      std::string distName = splittedStr[i++];
      lno_t np_ = 0;
      if(distName == "NORMAL"){
        int reqArg = 5;
        if (this->dimension == 3){
          reqArg = 7;
        }
        if(i + reqArg > argCnt) {
          std::string tmp = toString<int>(reqArg);
          throw INVALID_SHAPE_ARG(distName, tmp);
        }
        np_ = fromString<lno_t>(splittedStr[i++]);
        CoordinatePoint<T> pp;

        pp.x = fromString<T>(splittedStr[i++]);
        pp.y = fromString<T>(splittedStr[i++]);
        pp.z = 0;
        if(this->dimension == 3){
          pp.z = fromString<T>(splittedStr[i++]);
        }

        T sd_x = fromString<T>(splittedStr[i++]);
        T sd_y = fromString<T>(splittedStr[i++]);
        T sd_z = 0;
        if(this->dimension == 3){
          sd_z = fromString<T>(splittedStr[i++]);
        }
        this->coordinateDistributions[this->distributionCount++] = new CoordinateNormalDistribution<T, lno_t>(np_, this->dimension, pp , sd_x, sd_y, sd_z );

      } else if(distName == "UNIFORM" ){
        int reqArg = 5;
        if (this->dimension == 3){
          reqArg = 7;
        }
        if(i + reqArg > argCnt) {
          std::string tmp = toString<int>(reqArg);
          throw INVALID_SHAPE_ARG(distName, tmp);
        }
        np_ = fromString<lno_t>(splittedStr[i++]);
        T l_x = fromString<T>(splittedStr[i++]);
        T r_x = fromString<T>(splittedStr[i++]);
        T l_y = fromString<T>(splittedStr[i++]);
        T r_y = fromString<T>(splittedStr[i++]);

        T l_z = 0, r_z = 0;

        if(this->dimension == 3){
          l_z = fromString<T>(splittedStr[i++]);
          r_z = fromString<T>(splittedStr[i++]);
        }

        this->coordinateDistributions[this->distributionCount++] = new CoordinateUniformDistribution<T, lno_t>( np_,  this->dimension, l_x, r_x, l_y, r_y, l_z, r_z );
      } else if (distName == "GRID"){
        int reqArg = 6;
        if(this->dimension == 3){
          reqArg = 9;
        }
        if(i + reqArg > argCnt) {
          std::string tmp = toString<int>(reqArg);
          throw INVALID_SHAPE_ARG(distName, tmp);
        }

        lno_t np_x = fromString<lno_t>(splittedStr[i++]);
        lno_t np_y = fromString<lno_t>(splittedStr[i++]);
        lno_t np_z = 1;


        if(this->dimension == 3){
          np_z = fromString<lno_t>(splittedStr[i++]);
        }

        np_ = np_x * np_y * np_z;
        T l_x = fromString<T>(splittedStr[i++]);
        T r_x = fromString<T>(splittedStr[i++]);
        T l_y = fromString<T>(splittedStr[i++]);
        T r_y = fromString<T>(splittedStr[i++]);

        T l_z = 0, r_z = 0;

        if(this->dimension == 3){
          l_z = fromString<T>(splittedStr[i++]);
          r_z = fromString<T>(splittedStr[i++]);
        }

        if(np_x < 1 || np_z < 1 || np_y < 1 ){
          throw "Provide at least 1 point along each dimension for grid test.";
        }
        this->coordinateDistributions[this->distributionCount++] = new CoordinateGridDistribution<T, lno_t>
        (np_x, np_y,np_z, this->dimension, l_x, r_x,l_y, r_y, l_z, r_z , this->myRank);

      }
      else {
        std::string tmp = toString<int>(this->dimension);
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
    }

    float sum = 0;
    for(int i = 0; i < this->worldSize; ++i){
      sum += this->loadDistributions[i];
    }
    if (fabs(sum - 1.0) > 10*std::numeric_limits<float>::epsilon()){
      throw "Processor load ratios do not sum to 1.0.";
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
      if(shapeName == "SQUARE" && this->dimension == 2){
        if(i + 3 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "3");
        }
        CoordinatePoint<T> pp;
        pp.x = fromString<T>(splittedStr[i++]);
        pp.y = fromString<T>(splittedStr[i++]);
        T edge = fromString<T>(splittedStr[i++]);
        this->holes[this->holeCount++] = new SquareHole<T>(pp, edge);
      } else if(shapeName == "RECTANGLE" && this->dimension == 2){
        if(i + 4 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "4");
        }
        CoordinatePoint<T> pp;
        pp.x = fromString<T>(splittedStr[i++]);
        pp.y = fromString<T>(splittedStr[i++]);
        T edgex = fromString<T>(splittedStr[i++]);
        T edgey = fromString<T>(splittedStr[i++]);

        this->holes[this->holeCount++] = new RectangleHole<T>(pp, edgex,edgey);
      } else if(shapeName == "CIRCLE" && this->dimension == 2){
        if(i + 3 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "3");
        }
        CoordinatePoint<T> pp;
        pp.x = fromString<T>(splittedStr[i++]);
        pp.y = fromString<T>(splittedStr[i++]);
        T r = fromString<T>(splittedStr[i++]);
        this->holes[this->holeCount++] = new CircleHole<T>(pp, r);
      }  else if(shapeName == "CUBE" && this->dimension == 3){
        if(i + 4 > argCnt) {
          throw INVALID_SHAPE_ARG(shapeName, "4");
        }
        CoordinatePoint<T> pp;
        pp.x = fromString<T>(splittedStr[i++]);
        pp.y = fromString<T>(splittedStr[i++]);
        pp.z = fromString<T>(splittedStr[i++]);
        T edge = fromString<T>(splittedStr[i++]);
        this->holes[this->holeCount++] = new CubeHole<T>(pp, edge);
      }  else if(shapeName == "RECTANGULAR_PRISM" && this->dimension == 3){
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

      }  else if(shapeName == "SPHERE" && this->dimension == 3){
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
        std::string tmp = toString<int>(this->dimension);
        throw INVALIDSHAPE(shapeName, tmp);
      }
    }
    delete [] splittedStr;
  }

  void getWeightDistribution(std::string weight_distribution){
    if(weight_distribution == ""){
      return;
    }
    int count = this->countChar(weight_distribution, ' ');
    std::string *splittedStr = new string[count + 1];
    this->splitString(weight_distribution, ' ', splittedStr);
    //cout << count << endl;
    if(splittedStr[0] == "STEPPEDEQUATION"){
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
          if(this->dimension > 1){
            a2 = this->fromString<T>(value);
          }
          else {
            throw  parameter+ " argument is not valid when dimension is " + toString<int>(this->dimension);
          }

        }
        else if (parameter == "a3"){
          if(this->dimension > 2){
            a3 = this->fromString<T>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + toString<int>(this->dimension);
          }
        }
        else if (parameter == "b1"){
          b1 = this->fromString<T>(value);
        }
        else if (parameter == "b2"){
          if(this->dimension > 1){
            b2 = this->fromString<T>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + toString<int>(this->dimension);
          }
        }
        else if (parameter == "b3"){

          if(this->dimension > 2){
            b3 = this->fromString<T>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + toString<int>(this->dimension);
          }
        }
        else if (parameter == "c"){
          c = this->fromString<T>(value);
        }
        else if (parameter == "x1"){
          x1 = this->fromString<T>(value);
        }
        else if (parameter == "y1"){
          if(this->dimension > 1){
            y1 = this->fromString<T>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + toString<int>(this->dimension);
          }
        }
        else if (parameter == "z1"){
          if(this->dimension > 2){
            z1 = this->fromString<T>(value);
          }
          else {
            throw parameter+ " argument is not valid when dimension is " + toString<int>(this->dimension);
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


      this->wd =  new SteppedEquation<T,T>(a1, a2,  a3,  b1,  b2,  b3,  c,  x1,  y1,  z1, steps, values, stepCount);

      if(stepCount > 0){
        delete [] steps;
        delete [] values;
      }
    }
    else if (splittedStr[0] == "ANOTHERWEIGHTCLASS"){

    }
    else {
      throw "Unknown weight distribution name " + splittedStr[0];
    }

  }

  void parseParams(Teuchos::ParameterList params){
    try {
      std::string holeDescription  = "";
      std::string proc_load_distributions = "";
      std::string distinctDescription = "";
      std::string coordinate_distributions = "";
      std::string outfile = "";
      std::string weight_distribution = "";
      std::string stepped_equation_parameters = "";

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

        else if(paramName == "WeightDistribution"){
          weight_distribution = getParamVal<std::string>(pe, paramName);
          //TODO coordinate distribution description
        } else if (paramName.find(equation_with_step_function_parameters) == 0){

          stepped_equation_parameters +=  " " + paramName.substr(equation_with_step_function_parameters.size())+ "="+ getParamVal<std::string>(pe, paramName);
        }
        else if(paramName == "dim"){
          int dim = fromString<int>(getParamVal<std::string>(pe, paramName));
          if(dim < 2 && dim > 3){
            throw INVALID(paramName);
          } else {
            this->dimension = dim;
          }
        }
        /*
				else if(paramName == "num_global_coords"){
					gno_t np = getParamVal<gno_t>(pe, paramName);
					if(np < 1){
						throw INVALID(paramName);
					} else {
						this->globalNumCoords = np;
					}
				}
         */

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
        /*
				else if(paramName == "minx"){
					this->minx = getParamVal<T>(pe, paramName);
				}
				else if(paramName == "maxx"){
					this->maxx = getParamVal<T>(pe, paramName);
				}
				else if(paramName == "miny"){
					this->miny  = getParamVal<T>(pe, paramName);
				}
				else if(paramName == "maxy"){
					this->maxy = getParamVal<T>(pe, paramName);
				}
				else if(paramName == "minz"){
					this->minz = getParamVal<T>(pe, paramName);
				}
				else if(paramName == "maxz"){
					this->maxz = getParamVal<T>(pe, paramName);
				}
         */
        else {
          throw INVALID(paramName);
        }
      }


      if(this->dimension == 0){
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
      if(weight_distribution != equation_with_step_function && stepped_equation_parameters != ""){
        throw "STEPPEDEQUATION parameters are provided without weight distribution is selected.";
      }
      this->getWeightDistribution(weight_distribution+stepped_equation_parameters);
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
      delete this->wd;
    }

    delete []this->points;
  }

  GeometricGenerator(Teuchos::ParameterList &params, const RCP<const Teuchos::Comm<int> > & comm_){
    this->wd = NULL;
    this->holes = NULL; //to represent if there is any hole in the input
    this->dimension = 0;  //dimension of the geometry
    this->worldSize = comm_->getSize(); //comminication world object.
    this->numGlobalCoords = 0;	//global number of coordinates requested to be created.
    this->loadDistributions = NULL; //sized as the number of processors, the load of each processor.
    //this->distinctCoordinates = NULL; // if processors have different or same range for coordinates to be created.
    this->coordinateDistributions = NULL;
    this->holeCount = 0;
    this->distributionCount = 0;
    this->outfile = "";
    this->points =  NULL;

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
        this->numGlobalCoords += this->coordinateDistributions[i]->numPoints * this->loadDistributions[ii];
        if(ii < myRank){
          prefixSum += this->coordinateDistributions[i]->numPoints * this->loadDistributions[ii];
        }
      }
      myPointCount += this->coordinateDistributions[i]->numPoints * this->loadDistributions[myRank];
    }





    this->points = new CoordinatePoint<T> [myPointCount];
    this->numLocalCoords = 0;
    srand ( time(NULL) + myRank);
    for (int i = 0; i < distributionCount; ++i){
      lno_t requestedPointCount = this->coordinateDistributions[i]->numPoints *  this->loadDistributions[myRank];
      this->coordinateDistributions[i]->GetPoints(requestedPointCount,this->points + this->numLocalCoords, this->holes, this->holeCount,  this->loadDistributions, myRank);
      this->numLocalCoords += requestedPointCount;
    }



    if (this->distinctCoordSet){
      //TODO: Partition and migration.
    }

    if(this->outfile != ""){

      std::ofstream myfile;
      myfile.open ((this->outfile + toString<int>(myRank)).c_str());
      for(lno_t i = 0; i < this->numLocalCoords; ++i){

        myfile << this->points[i].x;
        if(this->dimension > 1){
          myfile << " " << this->points[i].y;
        }
        if(this->dimension > 2){
          myfile << " " << this->points[i].z;
        }
        myfile << std::endl;
      }
      myfile.close();
    }

    // target map
    Teuchos::ArrayView<const gno_t> eltList;

    if(this->numLocalCoords > 0){
      gno_t *gnos = new gno_t [this->numLocalCoords];
      for (lno_t i = 0; i < this->numLocalCoords; ++i){
        gnos[i] = i + prefixSum;
      }
      eltList = Teuchos::ArrayView<const gno_t> (gnos, this->numLocalCoords);
    }

    RCP<Tpetra::Map<lno_t, gno_t, node_t> > mp = rcp(
        new Tpetra::Map<lno_t, gno_t, node_t> (this->numGlobalCoords, eltList, 0, comm_));


    T **coords = new T *[this->dimension];
    for(int i = 0; i < this->dimension; ++i){
      coords[i] = new T[this->numLocalCoords];
    }

    for(int i = 0; i < this->numLocalCoords; ++i){
      T tmp = this->points[i].x;
      coords[0][i] = tmp;
      if(this->dimension > 1){
        coords[1][i] = this->points[i].y;
        if(this->dimension > 2){
          coords[2][i] = this->points[i].z;
        }
      }
    }

    Teuchos::Array<Teuchos::ArrayView<const T> > coordView(this->dimension);
    for (int i=0; i < this->dimension; i++){
      if(this->numLocalCoords > 0){
        Teuchos::ArrayView<const T> a(coords[i], this->numLocalCoords);
        coordView[i] = a;
      } else{
        Teuchos::ArrayView<const T> a;
        coordView[i] = a;
      }
    }

    tmVector = RCP< Tpetra::MultiVector<T, lno_t, gno_t, node_t> >(new Tpetra::MultiVector< T, lno_t, gno_t, node_t>( mp,
        coordView.view(0, this->dimension), this->dimension
    ));


    /*
		Zoltan2::XpetraMultiVectorInput < Tpetra::MultiVector<T, lno_t, gno_t, node_t> > xmv (RCP < Tpetra::MultiVector<T, lno_t, gno_t, node_t> > (tmVector));

		RCP< Tpetra::MultiVector<T, lno_t, gno_t, node_t> >tmVector2;
		Zoltan2::PartitioningSolution< Tpetra::MultiVector<T, lno_t, gno_t, node_t> > solution;
		xmv.applyPartitioningSolution<Tpetra::MultiVector<T, lno_t, gno_t, node_t> >(this->tmVector, &tmVector2, solution);
     */
    /*
    switch(this->dimension){
    case 1:
      for (lno_t i = 0; i < this->numLocalCoords; ++i){
        cout << this->wd->get1DWeight(coords[0][i]) << " " ;
      }
      break;
    case 2:
      for (lno_t i = 0; i < this->numLocalCoords; ++i){
        cout << this->wd->get2DWeight(coords[0][i], coords[1][i]) << " " ;
      }
      break;
    case 3:
      for (lno_t i = 0; i < this->numLocalCoords; ++i){
        cout << this->wd->get3DWeight(coords[0][i], coords[1][i], coords[2][i]) << " " ;
      }
      break;
    default:
      cerr <<"error in dimension." << endl;
      exit(1);
    }
     */
    for(int i = 0; i < this->dimension; ++i) delete [] coords[i]; delete []coords;
  }

  RCP< Tpetra::MultiVector<T, lno_t, gno_t, node_t> > getCoordinates(){

    return this->tmVector;
  }

};
