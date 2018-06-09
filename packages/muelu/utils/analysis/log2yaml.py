#!/usr/bin/env python3
"""log2yaml

Usage:
  log2yaml.py -i INPUT [-o OUTPUT] -m <mode>
  log2yaml.py (-h | --help)

Options:
  -h --help                     Show this screen.
  -i FILE --input-file=FILE     Input file
  -o FILE --output-file=FILE    Output file
  -m MODE --mode=MODE           Mode [muelu | albany | drekar]
"""
import logging
import yaml
from   docopt      import docopt

# YAML scheme structure
#
# scheme:                   <string>            # Particular YAML scheme ['muelu' | 'albany' | 'drekar']
# Steps:
#     c_step_1:
#         nl_its:           <int>               # number of nonlinear iterations
#         nl_step_1:
#             its:          <int>               # number of linear iterations
#             res_hist:     <double list>       # residual history
#             setup_time:   <double>            # setup time in s
#             solve_time:   <double>            # solve time in s
#         nl_step_2:
#             its:          <int>
#             res_hist:     <double list>
#             setup_time:   <double>            # setup time in s
#             solve_time:   <double>            # solve time in s
# Timer names:              <string list>
# Total times:
#    <timer_name>:          {MaxOverProcs: <double>, MeanOverCallCounts: <double>, MeanOverProcs: <double>, MinOverProcs: <double>}
# Call counts:
#    <timer_name>:          {MaxOverProcs: <double>, MeanOverCallCounts: <double>, MeanOverProcs: <double>, MinOverProcs: <double>}
# Number of processes:      1
# Statistics collected:     [MinOverProcs, MeanOverProcs, MaxOverProcs, MeanOverCallCounts]
# Time unit:                's'

from muelu_log2yaml     import MueLuLog
from albany_log2yaml    import AlbanyLog

if __name__ == '__main__':

    ## Process input
    options = docopt(__doc__)

    input_file = options['--input-file']
    output_file = options['--output-file']
    if output_file == None:
      output_file = input_file + '.yaml'
    mode = options['--mode']

    ## Validate parameters
    assert(mode == 'muelu' or mode == 'albany' or mode == 'drekar-belos')

    ## Run
    if mode == 'muelu':
        log = MueLuLog()
    elif mode == 'drekar-belos':
        log = DrekarBelosLog()
    else:
        raise

    yaml_data = log.run(input_file)

    ## Save output
    f = open(output_file, 'w')
    yaml.dump(yaml_data, stream=f)
