#!/bin/env python3
import re
import yaml
from log import Log

class MueLuLog(Log):
    """Log parser of MueLu logs.

    Needs to be converged to FSM."""

    def states(self):
        raise RuntimeError("Not a state machine")

    def transitions(self):
        raise RuntimeError("Not a state machine")

    def run(self, filename):
        # construct the YAML string
        mode = {}
        mode['residual'] = False
        mode['timing']   = False

        start_residual  = '\*\*\*\*\* Belos Iterative Solver: '
        mid_residual    = 'Iter.*, \[.*\] :\s*.*'
        end_residual    = '(?<=Number of iterations performed for this solve: )\d*'
        start_timers    = '(?<=TimeMonitor results over )\d*'
        mid_timers      = '.* \(.*\)\s*$'
        end_timers      = '==========='

        timer_names  = []
        timer_times  = {}
        timer_calls  = {}
        timer_serial = None

        through_residual = False

        yaml_string = '{"scheme":"muelu","Steps":{"c_step_0":{"nl_step_0":{'
        with open(filename) as f:
            for line in f:
                if re.search(start_residual, line) != None:
                    assert(timer_serial     == None)
                    assert(mode['residual'] == False)
                    assert(mode['timing']   == False)
                    mode['residual'] = True

                    through_residual = True

                    if yaml_string[-1] != '{':
                        yaml_string += ','
                    yaml_string += '"res_hist":['

                elif re.search(mid_residual, line) != None:
                    assert(mode['residual'] == True)
                    assert(mode['timing']   == False)

                    m = re.search('[^\s]*$', line)
                    res = m.group()
                    if yaml_string[-1] != '[':
                        yaml_string += ','
                    yaml_string += res

                elif re.search(end_residual, line) != None:
                    assert(mode['residual'] == True)
                    assert(mode['timing']   == False)
                    mode['residual'] = False
                    mode['nl_step']  = False

                    m = re.search(end_residual, line)
                    its = m.group()
                    yaml_string += '], "its":' + its + '}}'

                elif re.search(start_timers, line) != None:
                    assert(mode['residual'] == False)
                    assert(mode['timing']   == False)
                    mode['timing'] = True

                    if not through_residual:
                        # Somebody may have given a truncated log with
                        # just timers
                        yaml_string += '"res_hist":[], "its":0}}'

                    m = re.search(start_timers, line)
                    nprocs = m.group()

                    if nprocs == "1":
                        timer_serial = True
                    else:
                        timer_serial = False

                    # Finalize stepping
                    yaml_string += '}'

                    yaml_string += ',"Number of processes":' + nprocs
                    yaml_string += ',"Time unit":s'
                    yaml_string += ',"Statistics collected":["MinOverProcs","MeanOverProcs","MaxOverProcs","MeanOverCallCounts"]'

                elif re.search(end_timers, line) != None:
                    if mode['timing'] == False:
                        # there could be other ======== lines
                        continue

                    assert(mode['residual'] == False)
                    mode['timing'] = False

                    # Timer names
                    yaml_string += ',"Timer names":['
                    for name in timer_names:
                        if yaml_string[-1] != '[':
                            yaml_string += ','
                        yaml_string += '"' + name + '"'
                    yaml_string += ']'

                    # Total times
                    yaml_string += ',"Total times":{'
                    for name in timer_names:
                        if yaml_string[-1] != '{':
                            yaml_string += ','
                        yaml_string += '"' + name + '":{'
                        yaml_string += '"MinOverProcs":'        + timer_times[name]['MinOverProcs']
                        yaml_string += ',"MeanOverProcs":'      + timer_times[name]['MeanOverProcs']
                        yaml_string += ',"MaxOverProcs":'       + timer_times[name]['MaxOverProcs']
                        yaml_string += ',"MeanOverCallCounts":' + timer_times[name]['MeanOverCallCounts']
                        yaml_string += '}'
                    yaml_string += '}'

                    # Call counts
                    yaml_string += ',"Call counts":{'
                    for name in timer_calls:
                        if yaml_string[-1] != '{':
                            yaml_string += ','
                        yaml_string += '"' + name + '":{'
                        yaml_string += '"MinOverProcs":'        + timer_times[name]['MinOverProcs']
                        yaml_string += ',"MeanOverProcs":'      + timer_times[name]['MeanOverProcs']
                        yaml_string += ',"MaxOverProcs":'       + timer_times[name]['MaxOverProcs']
                        yaml_string += ',"MeanOverCallCounts":' + timer_times[name]['MeanOverCallCounts']
                        yaml_string += '}'
                    yaml_string += '}'

                    yaml_string += '}'

                elif re.search(mid_timers, line) != None:
                    if mode['timing'] == False:
                        # there could be other matching lines
                        continue
                    assert(mode['residual'] == False)

                    if re.search('Timer Name', line) != None:
                        # Skip header (in serial it matches the pattern)
                        continue

                    splits = line.split()
                    if timer_serial == True:
                        name = ' '.join(splits[0:-2])
                        name = name.replace('"', '\\"')
                        timer_names.append(name)
                        timer_times[name] = {}
                        timer_times[name]['MinOverProcs']       = splits[-2];
                        timer_times[name]['MeanOverProcs']      = splits[-2];
                        timer_times[name]['MaxOverProcs']       = splits[-2];
                        timer_times[name]['MeanOverCallCounts'] = splits[-2];
                        timer_calls[name] = {}
                        timer_calls[name]['MinOverProcs']       = splits[-1][1:-1];
                        timer_calls[name]['MeanOverProcs']      = splits[-1][1:-1];
                        timer_calls[name]['MaxOverProcs']       = splits[-1][1:-1];
                        timer_calls[name]['MeanOverCallCounts'] = splits[-1][1:-1];
                    else:
                        name = ' '.join(splits[0:-8])
                        name = name.replace('"', '\\"')
                        timer_names.append(name)
                        timer_times[name] = {}
                        timer_times[name]['MinOverProcs']       = splits[-8];
                        timer_times[name]['MeanOverProcs']      = splits[-6];
                        timer_times[name]['MaxOverProcs']       = splits[-4];
                        timer_times[name]['MeanOverCallCounts'] = splits[-2];
                        timer_calls[name] = {}
                        timer_calls[name]['MinOverProcs']       = splits[-7][1:-1];
                        timer_calls[name]['MeanOverProcs']      = splits[-5][1:-1];
                        timer_calls[name]['MaxOverProcs']       = splits[-3][1:-1];
                        timer_calls[name]['MeanOverCallCounts'] = splits[-1][1:-1];

        if timer_serial == None:
            # We did not encounter any timers
            yaml_string += '}}'

        assert(mode['residual'] == False)
        assert(mode['timing']   == False)

        return yaml_string
