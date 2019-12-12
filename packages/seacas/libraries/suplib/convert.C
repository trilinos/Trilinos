#include <fmt/format.h>

extern "C" {
#if defined(ADDC_)
void convert_(char *str, int64_t *value, int64_t *length, int len)
{
#else
void convert(char *str, int64_t *value, int64_t *length, int len)
{
#endif
  std::string s = fmt::format("{:n}", *value);
  std::copy(s.begin(), s.end(), str);
  int i;
  *length = s.size();
  for (i = s.size(); i < len; i++) {
    str[i] = ' ';
  }
}
}
