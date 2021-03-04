#pragma once
#ifndef ROL2_VISITOR_HPP
#define ROL2_VISITOR_HPP

namespace ROL2 {

template<typename _Visitable, typename..._Visitables> 
struct Visitor : Visitor<_Visitable>, Visitor<_Visitables...> {
  using Visitor<_Visitable>::visit;
  using Visitor<_Visitables...>::visit;
};

template<typename _Visitable> 
struct Visitor<_Visitable> {
  virtual ~Visitor<_Visitable>() = default;
  virtual void visit( const _Visitable& ) {};
};

template<typename _Acceptor, typename _Visitor, typename..._Visitors>
struct Visitable : Visitable<_Acceptor,_Visitor>,
                   Visitable<_Acceptor,_Visitors...> {
  using Visitable<_Acceptor,_Visitor>::accept;
  using Visitable<_Acceptor,_Visitors...>::accept;
};

template<typename _Acceptor, typename _Visitor>
struct Visitable<_Acceptor,_Visitor> {
  virtual void accept( _Visitor& visitor ) const {
    static_cast<const _Acceptor&>(*this)->accept(visitor);
 }
};

} // namespace ROL2

#endif // ROL2_VISITOR_HPP

