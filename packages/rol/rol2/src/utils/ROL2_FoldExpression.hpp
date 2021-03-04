#pragma once
#ifndef ROL2_FOLDEXPRESSION_HPP
#define ROL2_FOLDEXPRESSION_HPP

namespace ROL2 {

/** \file ROL2_FoldExpression.hpp
 *  \brief Emulate C++17 fold expressions with C++14
 */

template<typename Function, typename Arg>
auto fold( Function&& f, Arg&& arg ) {
  return std::forward<Arg>(arg);
}

template<typename Function, typename Arg, typename...Args>
auto fold( Function&& f, Arg&& arg, Args&&...args ) {
  return f(std::forward<Arg>(arg), fold(f,std::forward<Args>(args)...));
} 

} // namespace ROL2

#endif // ROL2_FOLDEXPRESSION_HPP

