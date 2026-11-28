#!/usr/bin/env python3
"""Migrate application code off the deprecated PHX::TemplateManager API.

Configuring Phalanx with Phalanx_HIDE_DEPRECATED_CODE=ON replaces
PHX::TemplateManager with a form that stores each object with its full type.
Both forms provide get(), set(), has(), build(), buildObjects(), disable(),
isEnabled(), reset(), getAsBase(), forEach() and forEachBase(), so code written
against those methods compiles either way.  This script moves code onto them.

It rewrites the three renamed methods:

    getAsObject<T>()  ->  get<T>()
    disableType<T>()  ->  disable<T>()
    deleteType<T>()   ->  reset<T>()

The replacements work with Phalanx as it is configured today, so the result can
be committed before the flag is turned on.

What this script does NOT do, by design, and reports with file and line instead:

  * It does not rewrite iterator loops.  forEach() skips types with no object,
    while the iterators visited them, so a counter kept alongside the loop has
    to become PHX::mpl::index_v<EvalT,EvalTypes> -- and a loop that counted on
    a null object needs a decision.  Loops over `auto` iterators can only be
    guessed at from the variable name; the compiler finds the rest.
  * It does not rewrite iterator.rcp().  Use getAsBase<T>() or forEachBase().
  * It does not rewrite PHX::TemplateManager<Seq, Base, Foo<_>> written out in
    full.  The template parameters differ between the two forms, so code that
    must build in both needs an #ifdef:

        #ifdef PHALANX_HIDE_DEPRECATED_CODE
          PHX::TemplateManager<Seq, Foo, Base>
        #else
          PHX::TemplateManager<Seq, Base, Foo<_> >
        #endif

    ObjectT takes one template parameter, so a template with more needs an
    alias template.  If the objects reach Base only through a cross cast, use
    an alias that ignores its argument, so that each slot holds a Base.

Code inside a conditional on PHALANX_HIDE_DEPRECATED_CODE is left alone and not
reported, since it has already been dealt with.

The renames match on the method name, as there is no way to tell from the text
what an object's type is.  Review the changes if your own classes have methods
with these names.

Treat a clean run as "the easy ones are done", not as "the migration is
finished": build with Phalanx_HIDE_DEPRECATED_CODE=ON and fix whatever the
compiler still rejects.  The typed form also checks a builder's result when it
is built, so run your tests with the flag on as well.
"""

import argparse
import os
import re
import sys

DEFAULT_EXTENSIONS = [".hpp", ".cpp", ".cc", ".h", ".cxx", ".hxx", ".C", ".impl"]

# The deprecated class's own implementation, which defines the old names, and
# the typed class, whose documentation lists them.
DEFAULT_EXCLUDES = ["Phalanx_TemplateManager_Def.hpp", "Phalanx_TemplateIterator.hpp",
                    "Phalanx_TemplateManager_Typed.hpp"]

RENAMES = {
    "getAsObject": "get",
    "disableType": "disable",
    "deleteType": "reset",
}

RENAME_RE = re.compile(r"\b(" + "|".join(RENAMES) + r")(\s*<)")

# One level of nested template arguments.  A closing `>` has to be followed
# directly by `::iterator`, so that the iterator of a
# std::vector<RCP<Foo_TemplateManager<T> > > does not match.
_ARGS = r"<(?:[^<>]|<[^<>]*>)*>"

RESIDUALS = [
    (re.compile(r"(?<!Vector)TemplateManager\s*" + _ARGS + r"\s*::\s*(?:const_)?iterator\b"),
     "iterator: use forEach() or forEachBase(); replace any loop counter with PHX::mpl::index_v"),
    (re.compile(r"\bFieldManager\s*" + _ARGS + r"\s*::\s*iterator\b"),
     "FieldManager iterator: use the methods templated on the evaluation type, in a PHX::mpl::for_each_type loop"),
    (re.compile(r"\b\w+\s*\.\s*rcp\s*\(\s*\)"),
     "iterator.rcp(): use getAsBase<T>() or forEachBase()"),
    (re.compile(r"\bPHX::TemplateManager\s*<"),
     "template parameters differ between the two forms: needs #ifdef PHALANX_HIDE_DEPRECATED_CODE"),
    # Singular names only: plurals like bc_field_managers are containers of them.
    (re.compile(r"\b(?:\w*tm|fm|\w*[Mm]anager_?|\w*[Ff]actory_?)\s*(?:\.|->)\s*begin\s*\(\s*\)"),
     "possible iterator loop over a template manager or field manager: use forEach() or forEachBase()"),
]

CONDITIONAL_OPEN_RE = re.compile(r"^\s*#\s*if")
CONDITIONAL_CLOSE_RE = re.compile(r"^\s*#\s*endif")
FLAG = "PHALANX_HIDE_DEPRECATED_CODE"


def process(text):
    """Returns (rewritten text, {old name: count}, [(line number, line, why)])."""
    out = []
    counts = {}
    residuals = []
    conditionals = []  # for each open #if, whether it tests the flag

    def rename(match):
        counts[match.group(1)] = counts.get(match.group(1), 0) + 1
        return RENAMES[match.group(1)] + match.group(2)

    # keepends, so that CRLF files stay CRLF
    for number, line in enumerate(text.splitlines(keepends=True), 1):
        if CONDITIONAL_OPEN_RE.match(line):
            conditionals.append(FLAG in line)
        elif CONDITIONAL_CLOSE_RE.match(line):
            if conditionals:
                conditionals.pop()

        if any(conditionals):
            out.append(line)
            continue

        line = RENAME_RE.sub(rename, line)
        out.append(line)

        code = line.split("//", 1)[0]
        for pattern, why in RESIDUALS:
            if pattern.search(code):
                residuals.append((number, line.strip(), why))
                break

    return "".join(out), counts, residuals


def iter_files(paths, extensions, excludes):
    for path in paths:
        if os.path.isfile(path):
            yield path
            continue
        for root, dirs, files in os.walk(path):
            dirs[:] = sorted(d for d in dirs if d not in excludes and not d.startswith("."))
            for name in sorted(files):
                if name not in excludes and os.path.splitext(name)[1] in extensions:
                    yield os.path.join(root, name)


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Migrate code off the deprecated PHX::TemplateManager API.")
    parser.add_argument("path", nargs="+", help="files or directories to scan")
    parser.add_argument("--apply", action="store_true",
                        help="write the changes (the default only reports them)")
    parser.add_argument("--check", action="store_true",
                        help="exit 1 if anything is pending or needs a human, for CI")
    parser.add_argument("--quiet", action="store_true",
                        help="summarize the occurrences that need a human by reason")
    parser.add_argument("--exclude", action="append", default=[],
                        help="file or directory name to skip; may be repeated")
    parser.add_argument("--extension", action="append", default=[],
                        help="additional file extension to scan; may be repeated")
    args = parser.parse_args(argv)

    extensions = DEFAULT_EXTENSIONS + args.extension
    excludes = set(DEFAULT_EXCLUDES + args.exclude)

    files_scanned = 0
    files_changed = 0
    totals = {}
    residuals = []

    for path in iter_files(args.path, extensions, excludes):
        try:
            # surrogateescape round-trips bytes that are not valid UTF-8, and
            # newline="" keeps CRLF files from being rewritten with LF.
            with open(path, encoding="utf-8", errors="surrogateescape", newline="") as fh:
                original = fh.read()
        except OSError as exc:
            print(f"warning: cannot read {path}: {exc}", file=sys.stderr)
            continue
        files_scanned += 1

        updated, counts, found = process(original)
        if counts:
            files_changed += 1
            for name, count in counts.items():
                totals[name] = totals.get(name, 0) + count
            if args.apply:
                with open(path, "w", encoding="utf-8", errors="surrogateescape", newline="") as fh:
                    fh.write(updated)
        residuals.extend((path, number, line, why) for number, line, why in found)

    verb = "rewrote" if args.apply else "would rewrite"
    print(f"scanned {files_scanned} files\n")
    if totals:
        print(f"{verb}:")
        for name in sorted(totals, key=lambda n: -totals[n]):
            print(f"  {totals[name]:>6}  {name} -> {RENAMES[name]}")
        print(f"  {sum(totals.values()):>6}  TOTAL in {files_changed} files\n")
    else:
        print("no mechanical rewrites needed\n")

    if residuals:
        print(f"{len(residuals)} occurrence(s) need a human:")
        if args.quiet:
            by_reason = {}
            for _, _, _, why in residuals:
                by_reason[why] = by_reason.get(why, 0) + 1
            for why, count in sorted(by_reason.items(), key=lambda kv: -kv[1]):
                print(f"  {count:>6}  {why}")
        else:
            for path, number, line, why in residuals:
                print(f"  {path}:{number}\n      {line}\n      -> {why}")
        print()

    if not args.apply and totals:
        print("nothing was modified; re-run with --apply to write these changes")

    pending = bool(totals) and not args.apply
    return 1 if (args.check and (residuals or pending)) else 0


if __name__ == "__main__":
    sys.exit(main())
