/*Paul
26-Jan-2003 Initial writeup
05-Feb-2003
*/

namespace Tpetra {

// Constructor 1: all values zeroed
//=======================================================================
template<typename OrdinalType, typename ScalarType>
Vector<OrdinalType, ScalarType>::Vector(VectorSpace<OrdinalType, ScalarType> const& vectorSpace)
  
// Constructor 2: values from user array
//=======================================================================
template<typename OrdinalType, typename ScalarType>
Vector<OrdinalType, ScalarType>::Vector(ScalarType* vectorEntries, 
																				VectorSpace<OrdinalType, ScalarType> const& vectorSpace)

// Constructor 3: Copy constructor
//=======================================================================
template<typename OrdinalType, typename ScalarType>
Vector<OrdinalType, ScalarType>::Vector(Vector<OrdinalType, ScalarType> const& vector)

// Destructor
//=======================================================================
template<typename OrdinalType, typename ScalarType>
Vector<OrdinalType, ScalarType>::~Vector()

// combine
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::combineEntries(CombineMode cm, OrdinalType numEntries, 
																										 OrdinalType* indices, ScalarType* values)

// set all to scalar
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::setAllToScalar(ScalarType const value)

// Set all to random
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::setAllToRandom()

// extract copy
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::extractCopy(ScalarType* userArray) const

// extract view
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::extractView(ScalarType** userPointerArray) const

// dot product
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::dotProduct(Vector<OrdinalType, ScalarType> const& x, 
																								 Vector<OrdinalType, ScalarType> const& y) const

// absolute value
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::absoluteValue(Vector<OrdinalType, ScalarType> const& x)

// reciprocal
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::reciprocal(Vector<OrdinalType, ScalarType> const& x)

// scale: this = scalarThis * this
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::scale(ScalarType scalarThis)

// scale: this = scalarX * x
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::scale(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x)

// update: this = scalarThis * this + scalarX * x
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::update(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x, 
																						 ScalarType scalarThis)

// update: x and y
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::update(ScalarType scalarX, Vector<OrdinalType, ScalarType> const& x, 
																						 ScalarType scalarY, Vector<OrdinalType, ScalarType> const& y, 
																						 ScalarType scalarThis)

// 1-norm
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::norm1() const

// 2-norm
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::norm2() const

// Inf-norm
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::normInf() const

// Weighted 2-norm (RMS Norm)
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::normWeighted(Vector<OrdinalType, ScalarType> const& weights) const

// minimum
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::minValue() const

// maximum
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::maxValue() const

// mean
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::meanValue() const

// Vector multiplication (elementwise) 
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::elementwiseMultiply(ScalarType scalarXY, 
																													Vector<OrdinalType, ScalarType> const& x, 
																													Vector<OrdinalType, ScalarType> const& y, 
																													ScalarType scalarThis)

// Reciprocal multiply (elementwise)
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::elementwiseReciprocalMultiply(ScalarType scalarXY, 
																																		Vector<OrdinalType, ScalarType> const& x, 
																																		Vector<OrdinalType, ScalarType> const& y, 
																																		ScalarType scalarThis)

// getSeed
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType Vector<OrdinalType, ScalarType>::getSeed() const

// setSeed
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::setSeed(ScalarType seed)

// [] operator, nonconst version
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType& Vector<OrdinalType, ScalarType>::operator[](OrdinalType index)

// [] operator, const version
//=======================================================================
template<typename OrdinalType, typename ScalarType>
ScalarType const& Vector<OrdinalType, ScalarType>::operator[](OrdinalType index) const

// GetNumMyEntries
//=======================================================================
template<typename OrdinalType, typename ScalarType>
OrdinalType Vector<OrdinalType, ScalarType>::getNumMyEntries() const

// getNumGlobalEntries
//=======================================================================
template<typename OrdinalType, typename ScalarType>
OrdinalType Vector<OrdinalType, ScalarType>::getNumEntries() const

// print
//=======================================================================
template<typename OrdinalType, typename ScalarType>
void Vector<OrdinalType, ScalarType>::print(ostream& os) const

} //namespace Tpetra
