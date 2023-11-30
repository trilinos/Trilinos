#ifndef KRINO_KRINO_KRINO_LIB_AKRI_SIGN_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_SIGN_HPP_

namespace krino {

inline bool sign_change( double f1, double f2 ) {
  return ( (f1 < 0.) ? (f2 >= 0.) : (f2 < 0.) ); // GOMA sign convention
  //return ( (f1 > 0.) ? (f2 <= 0.) : (f2 > 0.) ); // Marching cubes sign convention
}

inline int sign( double f ) {
  return ( (f < 0.) ? -1 : 1 ); // GOMA sign convention
  //return ( (f > 0.) ? 1 : -1 ); // Marching cubes sign convention
}

} // namespace krino


#endif /* KRINO_KRINO_KRINO_LIB_AKRI_SIGN_HPP_ */
