#!/usr/bin/env python

# This script will read a file and replace any Phalanx Evaluator macros it
# finds with their contents.
import re
import sys
import textwrap


class MacroReplacer:
    """Replace Phalanx Evaluator macros."""

    def __init__(self):
        """Initialize the object.

        Create the variables used to pass data around between members, along
        with the various replacement rules that will be used.
        """
        self.filename = ""
        self.name = ""
        self.line = ""
        self.output = []
        header = """template<typename EvalT, typename Traits>
class {0}
  :
  public PHX::EvaluatorWithBaseImpl<Traits>,
  public PHX::EvaluatorDerived<EvalT, Traits>
{{"""
        public = """
  public:

    {0}(
      const Teuchos::ParameterList& p);

    void
    postRegistrationSetup(
      typename Traits::SetupData d,
      PHX::FieldManager<Traits>& fm);

    void
    evaluateFields(
      typename Traits::EvalData d);
"""
        prePost = """
    void
    preEvaluate(
      typename Traits::PreEvalData d);

    void
    postEvaluate(
      typename Traits::PostEvalData d);
"""
        private = """
  private:

    using ScalarT = typename EvalT::ScalarT;"""
        self.replacements = {}
        self.replacements["PHX_EVALUATOR_CLASS\((.*)\)"] = (
            header + public + private)
        self.replacements["PHX_EVALUATOR_CLASS_PP\((.*)\)"] = (
            header + public + prePost + private)
        self.replacements["PHX_EVALUATOR_CLASS_END"] = (
            """}}; // end of class {0}
""")
        self.replacements["PHX_EVALUATOR_CTOR\((.*),(.*)\)"] = (
            """template<typename EvalT, typename Traits>
{0}<EvalT, Traits>::
{0}(
  const Teuchos::ParameterList& {1})""")
        self.replacements["PHX_EVALUATOR_CTOR_NAMESPACE\((.*),(.*),(.*)\)"] = (
            """template<typename EvalT, typename Traits>
{0}::{1}<EvalT, Traits>::
{1}(
  const Teuchos::ParameterList& {2})""")
        self.replacements["PHX_POST_REGISTRATION_SETUP\((.*),(.*),(.*)\)"] = (
            """template<typename EvalT, typename Traits>
void
{0}<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData {1},
  PHX::FieldManager<Traits>& {2})""")
        self.replacements["PHX_EVALUATE_FIELDS\((.*),(.*)\)"] = (
            """template<typename EvalT, typename Traits>
void
{0}<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData {1})""")
        self.replacements["PHX_PRE_EVALUATE_FIELDS\((.*),(.*)\)"] = (
            """template<typename EvalT, typename Traits>
void
{0}<EvalT, Traits>::
preEvaluate(
  typename Traits::PreEvalData {1})""")
        self.replacements["PHX_POST_EVALUATE_FIELDS\((.*),(.*)\)"] = (
            """template<typename EvalT, typename Traits>
void
{0}<EvalT, Traits>::
postEvaluate(
  typename Traits::PostEvalData {1})""")

    def replaceText(self, regex, text):
        """Replace a regular expression with its corresponding text."""
        indent = self.line.find(regex[:10])
        if indent >= 0:
            match = re.search(regex, self.line)
            if match:
                args = list(match.groups())
                if len(args) > 0:
                    self.name = args[0]
                    replacement = text.format(*args)
                else:
                    replacement = text.format(self.name)
                replacement = textwrap.indent(replacement, " " * indent)
                replacement = replacement[indent:]
                replacement = re.sub(regex, replacement, self.line)
                for i in replacement.splitlines(True):
                    self.output.append(i)
                return True
        return False

    def replaceInFile(self, filename):
        """Make all the necessary replacements in the given file."""
        self.filename = filename
        with open(filename) as f:
            for line in f:
                self.line = line
                foundIt = False
                for (regex, text) in self.replacements.items():
                    foundIt = foundIt or self.replaceText(regex, text)
                if not foundIt:
                    self.output.append(line)

    def writeOutput(self):
        """Overwrite the original file with the modified output."""
        with open(self.filename, 'w') as f:
            for line in self.output:
                f.write(line)


if __name__ == "__main__":
    r = MacroReplacer()
    r.replaceInFile(sys.argv[1])
    r.writeOutput()
