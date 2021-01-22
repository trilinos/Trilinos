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
    # whether to gather and print a csv_banner
    self.print_csv_banner = False
    # whatever the op's args should be
    self.op_args = []
    self.parse_cmdline_args(cmdline_args)

  def __repr__(self):
    return self.lcl_print()

  def __str__(self):
    return self.lcl_print()

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

  def parse_cmdline_arg_helper(self, cmdline_args):
    
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
  
      return

    except:
      pass

    # we failed -o, so try op specific stuff
    if self.short_op.endswith('ar') or self.short_op.endswith('ranlib'):
      for arg in cmdline_args:
        if arg.endswith('.a'):
          self.output_stats_file = arg + '.' + self.output_stats_file
          return
      # we hit this if we can't find a .a
      return


  def parse_cmdline_args(self, cmdline_args):
    base_build_dir_arg_prefix = '----base-build-dir='
    wrapper_header_arg = '----get_header'
    wrapper_op_arg_prefix = '----op='
    print_csv_header=False
    have_op=False
    # require that any wrapper arg be the first
    try:
      wrapper_arg_idx = 1
      wrapper_arg = cmdline_args[wrapper_arg_idx]
      if wrapper_arg.startswith(base_build_dir_arg_prefix):
        self.base_build_dir = wrapper_arg.split('=', 1)[1]
        wrapper_arg_idx += 1
      wrapper_arg = cmdline_args[wrapper_arg_idx]
      if wrapper_arg == wrapper_header_arg:
        self.print_csv_banner=True
      elif wrapper_arg.startswith(wrapper_op_arg_prefix):
        self.op = wrapper_arg.split('=', 1)[1]

        # we name the output as: blah.o.op.timing
        # this will result in blah.ar.timing, blah.mpicc.timing blah.ld.timing...
        self.short_op = os.path.basename(self.op)
        self.output_stats_file = short_op + '.timing'

        self.parse_cmdline_arg_helper(cmdline_args)

      else:
        raise Exception('unparseable arguments')

      # Remove the first wrapper_arg_idx+1 args (script name + wrapper args)
      self.op_args = cmdline_args[wrapper_arg_idx+1:]

    except:
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
      sys.exit(0)

