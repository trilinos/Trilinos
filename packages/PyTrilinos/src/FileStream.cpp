#include "FileStream.h"

FileStream::FileStream(std::FILE* file): self_file(file) {}

std::streambuf::int_type FileStream::overflow(std::streambuf::int_type c) {
  return std::fputc(c, self_file) == EOF?
    std::streambuf::traits_type::eof(): c;
}
