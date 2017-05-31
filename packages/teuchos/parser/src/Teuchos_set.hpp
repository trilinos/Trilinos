#ifndef TEUCHOS_SET_HPP
#define TEUCHOS_SET_HPP

#include <set>

namespace Teuchos {

template <typename T>
void unite_with(std::set<T>& a, std::set<T> const& b) {
  for (auto& x : b) a.insert(x);
}

template <typename T>
std::set<T> unite(std::set<T> const& a, std::set<T> const& b) {
  auto c = a;
  unite_with(c, b);
  return c;
}

template <typename T>
void subtract_from(std::set<T>& a, std::set<T> const& b) {
  for (auto& x : b) {
    auto it = a.find(x);
    if (it == a.end()) continue;
    a.erase(it);
  }
}

template <typename T>
bool intersects(std::set<T>& a, std::set<T> const& b) {
  for (auto& x : b) {
    auto it = a.find(x);
    if (it != a.end()) return true;
  }
  return false;
}

}

#endif
