#ifndef TEUCHOS_YAML_HPP
#define TEUCHOS_YAML_HPP

#include <Teuchos_language.hpp>
#include <Teuchos_reader_tables.hpp>

namespace Teuchos {
namespace yaml {

Language build_language();
LanguagePtr ask_language();
ReaderTablesPtr ask_reader_tables();

}  // end namespace yaml
}  // end namespace Teuchos

#endif
