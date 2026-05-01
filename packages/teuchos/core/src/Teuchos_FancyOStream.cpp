#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_FancyOStream_def.hpp"

#define TEUCHOS_FANCYOSTREAM_INSTANT(CHAR)                                     \
  template class Teuchos::basic_FancyOStream_buf<CHAR,                         \
                                                 std::char_traits<CHAR>>;      \
                                                                               \
  template class Teuchos::basic_FancyOStream<CHAR, std::char_traits<CHAR>>;    \
                                                                               \
  template class Teuchos::basic_OSTab<CHAR, std::char_traits<CHAR>>;           \
                                                                               \
  template Teuchos::RCP<                                                       \
      Teuchos::basic_FancyOStream<CHAR, std::char_traits<CHAR>>>               \
  Teuchos::tab(                                                                \
      const Teuchos::RCP<                                                      \
          Teuchos::basic_FancyOStream<CHAR, std::char_traits<CHAR>>> &,        \
      const int, const std::basic_string<CHAR, std::char_traits<CHAR>>);       \
                                                                               \
  template Teuchos::RCP<                                                       \
      Teuchos::basic_FancyOStream<CHAR, std::char_traits<CHAR>>>               \
  Teuchos::tab(                                                                \
      const Teuchos::RCP<std::basic_ostream<CHAR, std::char_traits<CHAR>>> &,  \
      const int, const std::basic_string<CHAR, std::char_traits<CHAR>>);

TEUCHOS_FANCYOSTREAM_INSTANT(char)
