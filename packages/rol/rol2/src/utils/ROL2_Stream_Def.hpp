#pragma once
#ifndef ROL2_STREAM_DEF_HPP
#define ROL2_STREAM_DEF_HPP

namespace ROL2 {

Ptr<std::ostream> makeStreamPtr( std::ostream& os, 
                                 bool          noSuppressOutput ) {
  Ptr<std::ostream> retstream;
  if( noSuppressOutput ) retstream = makePtrFromRef<std::ostream>(os);
  else retstream = makePtr<NullStream>();
  return retstream; 
}

Ptr<std::ostream> makeStreamPtr( Ptr<std::ostream> os, 
                                 bool              noSuppressOutput ) {
  Ptr<std::ostream> retstream;
  if( noSuppressOutput ) retstream = os;
  else retstream = makePtr<NullStream>();
  return retstream; 
}

} // namespace ROL2

#endif // ROL2_STREAM_DEF_HPP

