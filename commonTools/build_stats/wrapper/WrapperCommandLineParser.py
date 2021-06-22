import os
import sys

class WrapperCommandLineParser:
  """
    Commandline parsing find any wrapper args, determine any output names
  """
  def __init__(self, cmdline_args):
    # base build directory for computing correct relative path
    self.base_build_dir = ''
    # if we write anything out it goes here
    self.output_stats_file = ''
    # if op generates an output file (-o ...)
    self.op_output_file = ''
    # if we perform an operation this is it
    self.op = ''
    self.short_op = ''
    # whether to gather and print a csv_banner
    self.print_csv_banner = False
    # whatever the op's args should be
    self.op_args = []
    # a list of lists of commands to evaluate
    self.commands = []
    # whether we have the output arg
    self.have_output_arg = False
    # ENV control variables
    self.parse_nm = True
    self.output_fields = None

    self.time_cmd = 'not_set'
    # we must parse envs first, because they contain required parameters
    self.parse_env_controls()
    # finally parse the args
    self.parse_cmdline_args(cmdline_args)

  def parse_env_controls(self):
    """Parse control variables from the ENV (rather than command line)

       # REQUIRED
       TRILINOS_BUILD_STATS_TIME_CMD points to a valid GNU Time executable

       # REQUIRED
       TRILINOS_BUILD_STATS_INNER_OP: is the command we are wrapping

       # REQUIRED
       TRILINOS_BUILD_STATS_BASE_DIR : We need to know the `root` of the build
                              tree so we annotate paths correctly (see github
                              PR 8638 for an issue with Makefile builds)
       # OPTIONAL
       TRILINOS_BUILD_STATS_OUTPUT_FIELDS : control what gets written to timing
                               files Can enable only some fields
                               e.g.,
                                     FileName,FileSize,op


    """
    # optional, control which fields we write to a file
    # This does not promise we will not parse all possible fields
    # (That is to say, this does not promise any performance benefits)
    self.output_fields = os.environ.get('TRILINOS_BUILD_STATS_OUTPUT_FIELDS')

    err_msg=''
    # required : TRILINOS_BUILD_STATS_TIME_CMD
    #            TRILINOS_BUILD_STATS_INNER_OP
    #            TRILINOS_BUILD_STATS_BASE_DIR
    if 'TRILINOS_BUILD_STATS_TIME_CMD' not in os.environ:
      err_msg+=os.linesep
      err_msg+=('TRILINOS_BUILD_STATS_TIME_CMD (ENV) is required. CMake should '
                '`find` and set this, if using the build tools manually, locate '
                'GNU Time (typically /usr/bin/time) verify it supports `--format` '
                'and `--output`. Then set TRILINOS_BUILD_STATS_TIME_CMD=/path/to/time')

    if 'TRILINOS_BUILD_STATS_INNER_OP' not in os.environ:
      err_msg+=os.linesep
      err_msg+=('TRILINOS_BUILD_STATS_INNER_OP (ENV) is required. CMake should '
                'set this to a specific operations, e.g., ${CMAKE_C_COMPILER}. '
                'If you are using the tools independently, please see the docs '
                'for examples of how to write the wrapper scripts. E.g., '
                'export TRILINOS_BUILD_STATS_INNER_OP=mpicc')
 
    if 'TRILINOS_BUILD_STATS_BASE_DIR' not in os.environ:
      err_msg+=os.linesep
      err_msg+=('TRILINOS_BUILD_STATS_BASE_DIR (ENV) is required. CMake should '
                'set this to the build directory (top level). If using this script '
                'manually, set this to your build directory (full path). E.g., '
                'export TRILINOS_BUILD_STATS_BASE_DIR=/path/to/build')

    if err_msg:
      sys.stderr.write(err_msg)
      exit(-1)

    # set the required parameters - the dict will throw if these aren't defined
    # but we should have exit() if any errors.
    self.time_cmd = os.environ['TRILINOS_BUILD_STATS_TIME_CMD']
    self.op = os.environ['TRILINOS_BUILD_STATS_INNER_OP']
    self.base_build_dir = os.environ['TRILINOS_BUILD_STATS_BASE_DIR']

    # we name the output as: blah.o.op.timing
    # this will result in blah.ar.timing, blah.mpicc.timing blah.ld.timing...
    self.short_op = os.path.basename(self.op)
    self.output_stats_file = self.short_op + '.timing'

    parse_nm = os.environ.get('TRILINOS_BUILD_STATS_PARSE_NM', "True")
    if parse_nm.lower() == 'true':
      self.parse_nm = True
    elif parse_nm.lower() == 'false':
      self.parse_nm = False
    else:
      msg='ERROR: TRILINOS_BUILD_STATS_PARSE_NM is set to [{}]'.format(parse_nm)
      msg+=', but valid values are True or False. Defaulting to True{}'.format(os.linesep)
      sys.stderr.write(msg);
      self.parse_nm = True

  def __repr__(self):
    return self.lcl_print()

  def __str__(self):
    return self.lcl_print()

  def generate_stats(self):
    return self.have_output_arg

  def lcl_print(self):
    fmt_string = [
      'output_stats_file : {output_stats_file}',
      'op : {op}',
      'op_output_file : {op_output_file}',
      'print_csv_banner : {print_csv_banner}',

    ]
    return '\n'.join(fmt_string).format(
                  output_stats_file=self.output_stats_file,
                  op_output_file=self.op_output_file,
                  op=self.op,
                  print_csv_banner=self.print_csv_banner)

  def get_output_fields(self,csv_map):
    if self.output_fields:
      # this assumes it is a string of comma separated labels
      fields = self.output_fields.split(',')
    else:
      # apply sort here, so the output will be deterministic
      fields = sorted([ k for k in csv_map ])

    return fields

  def generate_commandlets(self, cmdline_args):

    # it seems we need to handle compound commands e.g., && (maybe ||)
    cmdlet = []
    for arg in cmdline_args:
      if arg.strip() == "&&":
        # add the command
        self.commands.append(cmdlet)
        # start a new command
        cmdlet = []
      elif arg.strip() != '':
        cmdlet.append(arg)

    if cmdlet:
      self.commands.append(cmdlet)
    # post - should have all commands broken up into lists of lists (of options)
    return

  def parse_cmdline_arg_helper(self, cmdline_args):

    self.have_output_arg = False
    # we want to do something different for ar, ranlib, or ld.*
    # these commands do not necessarily have a good 'output' arg denoted by -o
    # first try to find -o, if that passes then use it.
    # if not, then do something special based on ar/ranlib/ld.*

    # find the output arg (will raise an exception if not found)
    # we use -o blah.o or -o /path/to/blah.o or none at all
    try:
      output_idx = cmdline_args.index('-o')
      self.op_output_file = cmdline_args[output_idx+1]
      self.output_stats_file = self.op_output_file + '.' + self.output_stats_file
 
      self.have_output_arg = True 
      return

    except:
      pass

    # we failed -o, so try op specific stuff
    if self.short_op.endswith('ar') or self.short_op.endswith('ranlib'):
      for arg in cmdline_args:
        if arg.endswith('.a'):
          self.op_output_file = arg
          self.output_stats_file = arg + '.' + self.output_stats_file
          self.have_output_arg = True
          return
      # we hit this if we can't find a .a
      return


  def parse_cmdline_args(self, cmdline_args):
    wrapper_header_arg = '----get_header'
    print_csv_header=False
    # require that any wrapper arg be the first
    try:
      wrapper_arg_idx = 1
      wrapper_arg = cmdline_args[wrapper_arg_idx]
      if wrapper_arg == wrapper_header_arg:
        self.print_csv_banner=True
        # this isn't implemented....
        sys.stderr.write('----get_header was requested, but is not implemented'
                         '. Doing nothing.')
        exit(0)

      self.parse_cmdline_arg_helper(cmdline_args)

      # Remove the script name
      self.op_args = cmdline_args[1:]
      # we could clean this whole thing up some..
      self.generate_commandlets([self.op] + self.op_args)

    except Exception as e:
      print("Got an error parsing the command line in the compiler wrapper python script")
      print(e)
      raise
      # any error and we give up
      help_msg = ["Compiler wrapper:",
                  "  Usage: wrapper [---base-build-dir=<dir>] ----op=<compiler> [args] | ----get_header",
                  "",
                  "   ----base-build-dir=/path/to/base/build/dir",
                  "   Absolute path to the base project build directory",
                  "   ----op=/path/to/compiler",
                  "   Absolute path to the compiler we are wrapping",
                  "   ----get_header",
                  "   May not be combined with ----op or ----base-build-dir, prints the csv_header generated",
                  "",
                  "  Tool depends on finding a -o <output> option in args",
                  "  statistics will be written to <output>.timing",
                  ]
      print('\n'.join(help_msg))
      #raise
      sys.exit(0)

