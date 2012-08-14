#!/usr/bin/env python
#-------------------------------------------------------------------------------

import sys
import optparse

#-------------------------------------------------------------------------------

version = "pl2tex.py  Version 0.01 2012-07-19"

description = """
NAME
      pl2tex.py - Tuechos::ParameterList to LaTeX conversion utility
      (version """ + version + """)

SYNOPSIS

      pl2tex.py [OPTIONS] <input_pl.txt> <optional output_pl.tex>

      This utility converts a print out of Tuechos::ParameterList
      from the input_pl.txt to a LaTeX file that can be
      `input`ted/`include`d in another LaTeX document, e.g., a
      user's manual.

DESCRIPTION

      This utility will read in a Tuechos::ParameterList print out
      generated via the ParameterList::print function (termed
      input_pl.txt in this description), and create a LaTeX file
      (termed output_pl.tex) where each ParameterLists (and nested
      ParameterList) is separated into different sections.  Indexing
      of ParameterLists and parameters are automatically generated,
      along with hyperlinking between ParameterLists.

EXAMPLE

      One first needs to generate a text listing of a
      Tuechos::ParameterList, e.g., from a C++ code

      ...
      RCP<IntegratorBuilder<double> > ib = integratorBuilder<double>();
      RCP<const ParameterList> pl = ib->getValidParameters();
      ofstream fout("Rythmos_ParameterList.txt");
      pl->print(fout,Teuchos::ParameterList::PrintOptions().showDoc(true)
                                                           .indent(4));
      fout.close();
      ...

      One can generate the latex to be 'input'ted another latex
      document by executing this script

      % pl2tex.py -r "Integrator Base" --singleDescription VerboseObject \\
           --schematic \\
           ../../../../build/packages/rythmos/test/UnitTest/Rythmos_ParameterList.txt \\
           RythmosUsersManualParameterList_Raw.tex


DESIRED NEW FEATURES
      * ???

AUTHOR(S)
      Curtis Ober,   Sandia National Laboratories, ccober@sandia.gov
"""

#===============================================================================

def getList(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

def getLevel (line):
  description = '''Get the level indentation of line of line.'''
  return len(line) - len(line.lstrip())

def getPLName (line):
  description = '''Get the ParameterList name from line.'''
  return (line.replace('->','')).strip()

def laFix (line):
  return line.replace('_', '\\_').replace('&', '\\&')

def laReplace (line):
  return line.replace('_', '\\_').replace('&', 'and')

class Parameter:
  description = "A class for parameters containing name, default value, and description."
  def __init__(self, line):
      self.line = line
      n, d = line.split('=')
      self.name = n.strip()
      self.default = d.strip()
      self.description = []
  def addDescription(self, line):
      leader = line.find('#')+2
      l = line[leader:].rstrip()
      self.description.append(l)
  def write(self):
      print "     Parameter name = '%s'" % (self.name)
      print "     default = '%s'" % (self.default)
      print "     description = "
      for d in self.description:
        print "         '%s'" % (d)


class PListObj:
  description = "Class of each Tuechos::ParameterList."
  def __init__(self, line, label, parent, parent_label):
      self.name = ''
      self.label = label
      self.level = 0         # Indention level of ParameterList
      self.description = []
      self.verbatim = []
      self.parameters = []
      self.children = []
      self.parent = []
      self.parent.append((parent, parent_label))

      if line.find('->') < 0 :
        print 'Did not find ParameterList name on this line.'
        print '"%s"' % (line)
        sys.exit(1)
      self.name = getPLName(line)
      self.level = getLevel(line)
  def addDescription (self, line):
      leader = line.find('#')+2
      l = line[leader:].rstrip()
      self.description.append(l)
  def addVerbatim (self, line):
      leader = line.find('#')+2
      l = line[leader:].rstrip()
      self.verbatim.append(l)
  def addParam (self, param):
      self.parameters.append(param)
  def addChild (self, line, label):
      tail = line.find('->')
      l = line[:tail].strip()
      self.children.append((l,label))
  def addParent (self, parent, label):
      self.parent.append((parent, label))
  def write(self):
      print "ParameterList name = '%s'" % (self.name)
      print "  label = '%s'" % (self.label)
      print "  level = '%s'" % (self.level)
      print "  parent = " 
      for p in self.parent: print "     '%s, %s'" % (p[0], p[1])
      print "  children = "
      for c in self.children: print "     '%s, %s'" % (c[0], c[1])
      print "  description = "
      for d in self.description: print "     '%s'" % (d)
      print "  parameters = "
      for p in self.parameters: p.write()
      print

  def write_latex2(self, fout):
      l = []
      l.append('\n')
      l.append('% ----------------------------------------------------------\n')
      nf = laFix(self.name)
      l.append('\\subsection{%s}\n' % (nf))
      l.append('\\label{sec:%s}\n' % (laReplace(self.label)))
      l.append('\\index{%s}\n' % (nf))
      l.append('\n')


      # Output description
      l.append('\\begin{tabular}{lp{0.8\\textwidth}}\n')
      l.append('  Description: & \n')
      #l.append('\\noindent \hypertarget{hyp:%s}{Description:}\n' %(self.name))
      for d in self.description:
        l.append('  %s\n' %(d))
      if self.verbatim != []:
        l.append('\\begin{verbatim}\n')
        for v in self.verbatim:
          l.append('  %s\n' %(v))
        l.append('\\end{verbatim}\n')
      l.append('\\\\ & \\\\ \n')


      # Output parents
      l.append('  Parent(s): & \n')
      if self.parent[0][0] == 'ROOT':
        l.append(' ROOT \\\\ \n')
      else:
        for p in self.parent:
          p0f = laFix(p[0])
          p1r = laReplace(p[1])
          l.append('  %s (Section~\\ref{sec:%s})\n' % (p0f, p1r))
          l.append('    \\index{%s} \\newline \n' % (p0f))
        l.append('\\\\ \n')
      l.append('\n')


      # Output children
      l.append('  Child(ren): & \n')
      if self.children == []:
        l.append('None. \\newline \\\\ \n')
      else:
        for c in self.children:
          c0f = laFix(c[0])
          c1r = laReplace(c[1])
          l.append('  %s (Section~\\ref{sec:%s})\n' % (c0f, c1r))
          l.append('    \\index{%s} \\newline \n' % (c0f))
        l.append('\\\\ \n')
      l.append('\n')


      # Output parameters
      l.append('  Parameters: & \n')
      if self.parameters == []:
        l.append('None. \\\\ \n')
      else:
        l.append('\\begin{description}\n')
        for p in self.parameters:
          dall = ''

          i = 0
          while i < len(p.description):
            d = p.description[i]
            ds = d.strip().split()
            if ds[0] == 'Valid' and ds[-1] == 'values:':
              dall += '\n' + laFix(d) + '\n\n'
              i = i+2                       # i+1 should be '{'
              dall += '  \\begin{tabular}{lp{0.4\\textwidth}}\n'
              while p.description[i].strip() != '}':
                d = p.description[i]
                dn = p.description[i+1]
                if dn.strip()[0] == '"' or dn.strip()[0] == '}':
                  dall += laFix(d) + ' & \\\\ \n'
                  i += 1
                else:
                  dall += laFix(d) +' & '+ laFix(dn).strip() +' \\\\ \n'
                  i += 2
              dall += '  \\end{tabular}\n'
              i += 1                        # i+1 should be '}'
            else:
              dall += laFix(d) + '\n'
              i += 1

          dall = dall[:-1]
          pf = laFix(p.name)
          l.append('  \\item[%s = %s] \n%s\n' % (pf, p.default, dall))
          l.append('    \\index{%s!%s}\n' % (self.name, pf))
          l.append('    \\index{%s}\n' % (pf))
        l.append('\\end{description}\n')
        l.append('\\\\ \n')
      l.append('\\end{tabular}\n')


  def write_latex(self, fout):
      l = []
      l.append('\n')
      l.append('% ----------------------------------------------------------\n')
      nf = laFix(self.name)
      l.append('\\subsection{%s}\n' % (nf))
      l.append('\\label{sec:%s}\n' % (laReplace(self.label)))
      l.append('\\index{%s}\n' % (nf))
      l.append('\n')

      # Output description
      l.append('\\begin{list}{}\n')
      l.append('  {\setlength{\leftmargin}{1.0in}\n')
      l.append('   \setlength{\labelwidth}{0.75in}\n')
      l.append('   \setlength{\labelsep}{0.125in}}\n')
      l.append('  \\item[Description:]\n')
      for d in self.description:
        l.append('    %s\n' %(d))
      if self.verbatim != []:
        l.append('\\begin{verbatim}\n')
        for v in self.verbatim:
          l.append('  %s\n' %(v))
        l.append('\\end{verbatim}\n')

      # Output parents
      l.append('  \\item[Parent(s):]\n')
      if self.parent[0][0] == 'ROOT':
        l.append('   ROOT \n')
      else:
        for p in self.parent:
          p0f = laFix(p[0])
          p1r = laReplace(p[1])
          l.append('    %s (Section~\\ref{sec:%s})\n' % (p0f, p1r))
          l.append('      \\index{%s} \n' % (p0f))
          l.append('      \\newline \n')
        l = l[:-1]

      # Output children
      l.append('  \\item[Child(ren):]\n')
      if self.children == []:
        l.append('    None. \n')
      else:
        for c in self.children:
          c0f = laFix(c[0])
          c1r = laReplace(c[1])
          l.append('    %s (Section~\\ref{sec:%s})\n' % (c0f, c1r))
          l.append('      \\index{%s} \n' % (c0f))
          l.append('      \\newline \n')
        l = l[:-1]

      # Output parameters
      l.append('  \\item[Parameters:]\n')
      if self.parameters == []:
        l.append('    None. \n')
      else:
        l.append('    \\begin{description}\n')
        for p in self.parameters:
          dall = ''

          i = 0
          while i < len(p.description):
            d = p.description[i]
            ds = d.strip().split()
            if ds[0] == 'Valid' and ds[-1] == 'values:':
              dall += '\n' + laFix(d) + '\n\n'
              i = i+2                       # i+1 should be '{'
              dall += '      \\begin{tabular}{lp{0.4\\textwidth}}\n'
              while p.description[i].strip() != '}':
                d = p.description[i]
                dn = p.description[i+1]
                if dn.strip()[0] == '"' or dn.strip()[0] == '}':
                  dall += laFix(d) + ' & \\\\ \n'
                  i += 1
                else:
                  dall += laFix(d) +' & '+ laFix(dn).strip() +' \\\\ \n'
                  i += 2
              dall += '      \\end{tabular}\n'
              i += 1                        # i+1 should be '}'
            else:
              dall += laFix(d) + '\n'
              i += 1

          dall = dall[:-1]
          pf = laFix(p.name)
          l.append('      \\item[%s = %s] \n%s\n' % (pf, p.default, dall))
          l.append('        \\index{%s!%s}\n' % (self.name, pf))
          l.append('        \\index{%s}\n' % (pf))
        l.append('\\end{description}\n')
        l.append('\n')
      l.append('\\end{list}\n')


      fout.writelines(l)

#===============================================================================

def main():

  # Define the command line options
  p = optparse.OptionParser(description)

  p.add_option("-r", dest="root", default='ROOT', \
                     action="store", type="string", \
                     help='''Specify the root ParameterList name.''')

  p.add_option("--singleDescription", dest="singleDescription", default=[], \
                     action="callback", type="string", callback=getList, \
                     help='''List of ParameterLists that only a single
                             description is desired on output.
                             ParameterLists are deliminated by commas
                             (",") without spaces.''')

  p.add_option("-s", "--schematic", dest="schematic", default=False, \
                     action="store_true", \
                     help='''Output as a figure a schematic of the
                             ParameterList heirarchy.''')

  #-------------------------------
  # Parse the command line options
  opts, args = p.parse_args()

  plists = []
  fin = open(args[0], 'r')

  current_level = -1
  current_plist = None

  # Setup the Root ParameterList
  line = fin.readline()
  root_line = (getLevel(line)-1)*' ' + opts.root + ' ->'
  current_plist = PListObj(root_line, getPLName(root_line), 'ROOT', 'ROOT')
  plists.append(current_plist)

  while line:

    #print line
    # Check for ParameterList
    if line.find('->') >= 0 :
      current_level = getLevel(line)

      # Should we add ParameterList and is it a duplicate ParameterList?
      add_plist = True
      duplicate_plist = False
      match_plist = None
      for plist in plists:
        name = getPLName(line)
        if name == plist.name:
          duplicate_plist = True
          if name in opts.singleDescription: add_plist = False
          match_plist = plist
          break

      # New plist to add.
      if add_plist:
        next_line = fin.readline()
        if next_line.strip() == '[empty list]':
          print "Warning -- '%s' is an empty list!" % (line.strip())
          line = fin.readline()
        else:
          # Determine which plist is parent.
          parent_plist = None
          label = getPLName(line)
          for plist in reversed(plists):
            if current_level-1 == plist.level:
              parent_plist = plist
              if duplicate_plist: label += '-' + parent_plist.name
              parent_plist.addChild(line, label) 
              break
          current_plist = PListObj(line, label,
                                   parent_plist.name, parent_plist.label)
          plists.append(current_plist)
          line = next_line
        #current_plist.write()

      # plist already in plists.
      else:
        if not (getPLName(line) in opts.singleDescription):
          print "Warning -- '%s' is already in the list!" % (line.strip())
        match_plist.addParent(current_plist.name, current_plist.label)
        line = fin.readline()
        while getLevel(line) >= current_level:
          line = fin.readline()
        current_level = getLevel(line)

      # Check for ParameterList description.
      while getLevel(line) == current_level+1 and line[getLevel(line)] == '#':
        current_plist.addDescription(line)
        line = fin.readline()

      # Check for ParameterList verbatim description.
      while getLevel(line) == current_level+1 and line.strip()=='Description =':
        line = fin.readline()
        while getLevel(line) == current_level+2 and line[getLevel(line)] == '#':
          current_plist.addVerbatim(line)
          line = fin.readline()

      # Check for parameters in ParameterList
      while getLevel(line) == current_level+1 and line.find('=') >= 0:
        current_param = Parameter(line)
        line = fin.readline()
        # Check for Parameter description.
        while getLevel(line) == current_level+2 and line[getLevel(line)] == '#':
          current_param.addDescription(line)
          line = fin.readline()
        #current_param.write()
        current_plist.addParam(current_param)
      #current_plist.write()

    elif line.strip() == '':  # Blank line
      line = fin.readline()
    else:
      print 'Error: Could not parse this line:'
      print "'" + line + "'"
      sys.exit(1)

  fin.close()

  fout = open(args[1], 'w')

  # Write out schematic of ParameterList Heirarchy
  if opts.schematic:
    l = []
    l.append('\\begin{figure} \n')
    l.append('\\centering{} \n')
    l.append('  \\begin{tabular}{p{0.9\\textwidth}}\n')
    current_parent = ''
    parent_count = 1
    for plist in plists:
      if current_parent != plist.parent:
        current_parent = plist.parent
        parent_count = 1
      else: parent_count += 1
      indent = '\hspace*{%gin} ' % (int(plist.level-plists[0].level)*0.2)
      if parent_count < 6:
        nf = laFix(plist.name)
        nr = laReplace(plist.label)
        l.append('%s%s (Section~\\ref{sec:%s})\n' % (indent,nf, nr))
        l.append('    \\index{%s} \\\\ \n' % (nf))
      elif parent_count < 7:
        l.append('%s ... \\\\ \n' % (indent))
    l.append('  \\end{tabular}\n')
    l.append('\\caption{Schematic of ParameterList heirarchy.} \n')
    l.append('\\label{fig:ParameterList-schematic} \n')
    l.append('\\end{figure} \n')
    l.append('\\newpage \n')
    fout.writelines(l)

  for plist in plists:
    plist.write_latex(fout)
  fout.close()

  sys.exit(0)

#-------------------------------------------------------------------------------

if __name__ == "__main__":
  main()
