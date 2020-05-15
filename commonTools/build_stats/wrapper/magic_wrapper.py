#!/usr/bin/env python
'''

Note:
Try to by python2 and python3 compliant

Ideally, keep this file minimal to reduce parsing overhead
Most 'work' is done in the Classes imported
'''

# we need to run subcommands
import os
import csv
import sys
from NMParser import NMParser
from WrapperCommandLineParser import WrapperCommandLineParser
from WrapperOpTimer import WrapperOpTimer

# given a dict of key/val pairs, write them as a CSV line
def write_csv_map(filename,csv_map):
  try:
    with open(filename, 'w') as csvfile:
      writer = csv.DictWriter(csvfile,
                              fieldnames=[ k for k in csv_map ])
      writer.writeheader()
      writer.writerow(csv_map)
  except IOError:
    print("I/O error")

def main(cmdline_args):
  # parse the command line args, find wrapper args and organize the
  # the info into fields in this class
  wcp = WrapperCommandLineParser(cmdline_args)
  # you can 'print()' a WrapperCommandLineParser
  # we could add a verbose mode (----verbose)
  #print(wcp)

  # keep a dict of field : value
  # first do the operation
  # this must be first, as it generates the output file
  csv_map = WrapperOpTimer.time_op(op=wcp.op,
                                   op_output_file=wcp.op_output_file,
                                   output_stats_file=wcp.output_stats_file,
                                   op_args=wcp.op_args)
  # test nm
  # we probably need some to handle the case the .o isn't created
  # as-us, the following will return empty dicts (we parse/validate) the output
  # from NM, so we won't return garbage
  nmp = NMParser.parse_object(wcp.op_output_file)
  # NMParser.print_counts(nmp)
  # add NM's output to our csv data
  # we could move this into the NM parser
  csv_map.update(NMParser.get_csv_map(nmp))

  # ultimately, print the csv data to a file
  # make sure to quote csv columns
  write_csv_map(wcp.output_stats_file, csv_map)


if __name__ == '__main__':
  main(sys.argv)

