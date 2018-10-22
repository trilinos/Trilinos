#!/usr/bin/env python3
import yaml
from log import Log
import logging
import re
from transitions import Machine

class MueLuLog(Log):
    """Log parser of MueLu logs."""

    states = ['none', 'residual', 'timers']

    transitions = [
        { 'trigger': 'n2r',     'source': 'none',       'dest': 'residual'  },
        { 'trigger': 'n2t',     'source': 'none',       'dest': 'timers'    },
        { 'trigger': 'r2n',     'source': 'residual',   'dest': 'none'      },
        { 'trigger': 't2n',     'source': 'timers',     'dest': 'none'      }
    ]

    def run(self, filename):
        """We skip all timer outputs until the last one"""
        logging.basicConfig(level=logging.DEBUG)
        logging.getLogger('transitions').setLevel(logging.WARNING)

        machine = self
        M = Machine(model=machine, states=self.states, transitions=self.transitions, initial='none')

        # construct the YAML string
        start_residual  = '\*\*\*\*\* Belos Iterative Solver: '
        mid_residual    = 'Iter.*, \[.*\] :\s*.*'
        end_residual    = '(?<=Number of iterations performed for this solve: )\d*'
        start_timers    = '(?<=TimeMonitor results over )\d*'
        mid_timers      = '.* \(.*\)\s*$'
        end_timers      = '==========='

        nprocs       = -1
        timer_names  = []
        timer_times  = {}
        timer_calls  = {}
        timer_serial = None

        through_residual = False

        yaml_string = '{"scheme":"muelu","Steps":{"c_step_0":{"nl_step_0":{'
        with open(filename) as f:
            try:
                for line in f:
                    if re.search(start_residual, line) != None:
                        if machine.state != 'none':
                            raise RuntimeError('Wrong state: ' + machine.state)
                        machine.n2r()

                        through_residual = True

                        if yaml_string[-1] != '{':
                            yaml_string += ','
                        yaml_string += '"res_hist":['

                    elif re.search(mid_residual, line) != None:
                        if machine.state != 'residual':
                            raise RuntimeError('Wrong state: ' + machine.state)

                        m = re.search('[^\s]*$', line)
                        res = m.group()
                        if yaml_string[-1] != '[':
                            yaml_string += ','
                        yaml_string += res

                    elif re.search(end_residual, line) != None:
                        if machine.state != 'residual':
                            raise RuntimeError('Wrong state: ' + machine.state)
                        machine.r2n()

                        m = re.search(end_residual, line)
                        its = m.group()
                        yaml_string += '], "its":' + its + '}}'

                    elif re.search(start_timers, line) != None:
                        if machine.state != 'none':
                            raise RuntimeError('Wrong state: ' + machine.state)
                        machine.n2t()

                        # Reset timers
                        timer_names  = []
                        timer_times  = {}
                        timer_calls  = {}

                        m = re.search(start_timers, line)
                        nprocs = m.group()

                        if nprocs == "1":
                            timer_serial = True
                        else:
                            timer_serial = False

                    elif re.search(end_timers, line) != None:
                        if machine.state != 'timers':
                            # there could be other ======== lines
                            continue
                        machine.t2n()

                    elif re.search(mid_timers, line) != None:
                        if machine.state != 'timers':
                            # there could be other matching lines
                            continue

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

                if not through_residual:
                    # Somebody may have given a truncated log with
                    # just timers
                    yaml_string += '"res_hist":[], "its":0}}'

                if not timer_serial == None:
                    # We did encounter some timers

                    # Finalize stepping
                    yaml_string += '}'

                    yaml_string += ',"Number of processes":' + nprocs
                    yaml_string += ',"Time unit":s'
                    yaml_string += ',"Statistics collected":["MinOverProcs","MeanOverProcs","MaxOverProcs","MeanOverCallCounts"]'

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

            except RuntimeError as e:
                raise RuntimeError("Caught an error while parsing on line " + str(lineno) + ": " + e.args[0])

        if timer_serial == None:
            # We did not encounter any timers
            yaml_string += '}}'

        try:
            yaml_data = yaml.load(yaml_string)
        except yaml.parser.ParserError:
            print('Badly formatted YAML string:\n', yaml_string)
            raise RuntimeError('Did you select the right mode?')

        return yaml_data
