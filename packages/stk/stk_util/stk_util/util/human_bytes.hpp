#ifndef STK_HUMAN_BYTES_H
#define STK_HUMAN_BYTES_H

#include <string>

namespace stk {

/*
 * Given a number of bytes, this produces a string of the form
 *    123 B or 123 KB or 123 MB or 123 GB
 * depending on the input number of bytes
 */

std::string human_bytes(size_t bytes);
}
#endif
