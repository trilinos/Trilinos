"""

Note:
Try to by python2 and python3 compliant
"""
import subprocess # spawning nm
import re         # re matching
import os         # line seperator

class NMParser:
  """Simple NM parser that"""

  # the values are
  nm_option_csv_map = {
    'N' : 'symbol_debug',
    'p' : 'symbol_stack_unwind',
    'R' : 'symbol_ro_data_global',
    'r' : 'symbol_ro_data_local',
    'T' : 'symbol_text_global',
    't' : 'symbol_text_local',
    'u' : 'symbol_unique_global',
  }

  nm_option_desc_map = {
    'N' : 'debugging symbol',
    'p' : 'stack unwind section',
    'R' : 'read only global data',
    'r' : 'read only local data',
    'T' : 'global text section',
    't' : 'local text section',
    'u' : 'unique global symbol',
  }

  nm_re_type_expr = ''.join(nm_option_desc_map)
  nm_re_str = r'^[a-zA-Z0-9]+\s+(?P<size_hex>[a-zA-Z0-9]{2,})\s+(?P<type>[' + nm_re_type_expr + '])\s+'
  nm_re = re.compile(nm_re_str)

  @staticmethod
  def parse_object(filename):
    """
      Simple NM parsing of an object file
      Given an object file, we call nm -aS file

      Next, we parse stdout and match symbol lines corresponding to types
      from nm_option_desc_map.

      Data are aggregated into a dict using the keys from nm_option_desc_map

      The keys are obtained from nm_option_desc_map and enforced inside the regex used
      See nm_re_type_expr, nm_re_str, and nm_re in the static fields of this class
    """
    p = subprocess.Popen(['nm', '-aS', filename],
                         stdout=subprocess.PIPE)
    output = p.communicate()[0]

    nm_counts = dict()

    for line in output.split(os.linesep):
      m = NMParser.nm_re.match(line)
      if m:
        nm_counts[m.group('type')] = nm_counts.get(m.group('type'), 0) + 1
    # return what we found
    return nm_counts

  @staticmethod
  def print_counts(nm_counts,
                   cvs_line=False,
                   csv_header=False):
    for k,v in nm_counts.items():
      print("\"{key}\",{value}".format(key=NMParser.nm_option_desc_map[k],
                                       value=v))
  @staticmethod
  def get_csv_map (nm_counts):
    # create a map of the form: csv_header_str : value
    # loop over the csv_map, which will guarantee we always return the same columns.
    # otherwise, looping over nm_counts will only return csv columns found in this specific file
    # , while the wrapper needs consistent output from all files parsed
    csv_map = { v : nm_counts.get(k,0) for k,v in NMParser.nm_option_csv_map.items() }
    return csv_map

