"""
Python2and3: Simple set of utilies to write Python code that is Python 2 and
Python 3 compatible.
"""

import sys

#
# Byte array / string / unicode support across Python 2 & 3
#
# Note that the str class in Python 2 is an ASCII string (byte) array and in
# Python 3 it is a Unicode object. For Python 3 code that is backward compatible
# with Python 2, we sometimes need version-specific conversion functions to give
# us the data type we desire. These functions are:
#
#     b(x)    Return a byte array of str x, much like b'<string const>' in
#             Python 3.
#     s(x)    Return a version-specific str object equivalent to x
#     u(x)    return a unicode object equivalent to x, much like
#             u'<string const>' in Python 2.
#     stru()  Returns the stirng 'u' if that would be printed in Python 2
#             where in Python 3 return empty ''.  This is used to build string
#             representations so as genrated by the str() function converting
#             a dict to a string.
#
if sys.version_info < (3,):
  # Python 2
  def b(x): return x
  def s(x): return x
  def u(x): return unicode(x)
  def stru(): return 'u'
else:
  # Python 3
  import codecs
  def b(x): return codecs.latin_1_encode(x)[0]
  def s(x):
     try:
        return x.decode("utf-8")
     except AttributeError:
        return x
  def u(x): return x
  def stru(): return ''


def csvReaderNext(csvReader):
  if sys.version_info < (3,):
    return csvReader.next()
  else:
    return next(csvReader)
