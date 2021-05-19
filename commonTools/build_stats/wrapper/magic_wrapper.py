#!/usr/bin/env python3
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
def write_csv_map(filename,csv_map,csv_fields):
  try:
    with open(filename, 'w') as csvfile:
      writer = csv.DictWriter(csvfile,
                              fieldnames=csv_fields,
                              # ignore fields in the csv_map that aren't
                              # in fieldnames
                              extrasaction='ignore')
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
  #
  # WARNING: Be very careful with stdout before these commands.  If the wrapped command
  # has shell redirection it can slurp up Python's output... best to require all messages
  # go after the compiler commnand has completed.
  if wcp.generate_stats():
    (csv_map, returncode) = WrapperOpTimer.time_op(wcp)
    #print("======> Gathering stats...", file=sys.stdout)
  else:
    # only run the command and return the return code
    returncode = 0
    for cmd in wcp.commands:
      returncode |= WrapperOpTimer.run_cmd(cmd)
    #print("##======> NO stats {}".format(wcp.op_output_file), file=sys.stdout)
    return returncode

  if returncode == 0:
    if wcp.parse_nm:
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
    write_csv_map(wcp.output_stats_file,
                  csv_map,
                  csv_fields=wcp.get_output_fields(csv_map))

  # NOTE: Above, we don't write the *.timing file if the build failed because
  # the output target file may not exist!  And we don't want a CSV file entry
  # that does not have all of the fields.  ToDo: It might be nice to have an
  # entry for files that don't build and just put in empty values for NM data
  # so we can display that with the other stats

  return returncode


if __name__ == '__main__':
  sys.exit(main(sys.argv))

