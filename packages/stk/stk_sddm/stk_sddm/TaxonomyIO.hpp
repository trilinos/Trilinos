#ifndef STK_SDDM_TAXONOMY_IO_H
#define STK_SDDM_TAXONOMY_IO_H

#include <iosfwd>

namespace stk {
namespace sddm {

class Taxonomy;

void load(std::istream &is, Taxonomy &taxonomy);

std::ostream &printGrammar(std::ostream &os, const Taxonomy &taxonomy);

std::ostream &dump(std::ostream &os, const Taxonomy &taxonomy);
std::ostream &xml(std::ostream &os, const Taxonomy &taxonomy);

std::ostream &operator<<(std::ostream &os, const Taxonomy &taxonomy);
std::ostream &operator<<(std::ostream &os, const Taxon &taxon);

} // namespace sddm
} // namespace stk

#endif // STK_SDDM_TAXONOMY_IO_H
