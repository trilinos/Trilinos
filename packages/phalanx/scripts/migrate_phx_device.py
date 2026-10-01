#!/usr/bin/env python3
"""Migrate application code off PHX::Device used as an execution space.

PHX::Device is becoming Kokkos::Device<ExecutionSpace, MemorySpace> so that
Phalanx can honour a configured memory space (Kokkos::SharedSpace in
particular).  A Kokkos::Device provides only execution_space, memory_space and
device_type, so every use of PHX::Device that wanted an *execution space* has
to be respelled PHX::exec_space.

PHX::ExecutionSpace is defined as PHX::Device::execution_space, so the
replacements this script makes are correct both before and after that change.
Run it now, against current Phalanx, and the result keeps working either way.

The script also unifies the older names for those two types.  Phalanx has
carried six names for them -- exec_space, ExecSpace, mem_space, MemSpace and the
DefaultExecutionSpace / DefaultMemorySpace pair -- and they collapse onto
PHX::ExecutionSpace and PHX::MemorySpace.  Those renames are unconditional:
the old names only ever meant these types.

What this script does NOT do, by design:

  * It does not touch PHX::Device where it names a device -- the third argument
    of Kokkos::View and friends.  That use stays, and its meaning is exactly
    what the change is for.
  * It does not follow aliases.  `using MyDev = PHX::Device;` hides whether the
    later uses are devices or execution spaces, and no textual tool can tell.
    Those are reported.
  * It does not rewrite templates it has not been told about.  Whether a given
    template's first parameter is an execution space is a fact about that
    template, not something that can be inferred from the text, so the list
    lives in KNOWN_EXEC_SPACE_TEMPLATES below and can be extended with
    --exec-template.

Anything it cannot classify is reported with file and line so you can look at
it.  Treat a clean run as "the easy ones are done", not as "the migration is
finished": build with Phalanx_ENABLE_DEPRECATED_DEVICE_AS_EXECUTION_SPACE=OFF
and fix whatever the compiler still rejects.
"""

import argparse
import os
import re
import sys

DEFAULT_EXTENSIONS = [".hpp", ".cpp", ".cc", ".h", ".cxx", ".hxx", ".C", ".impl"]

# Members an execution space has but Kokkos::Device does not.  execution_space
# and memory_space are deliberately absent: Device has both, so those uses are
# already correct and must not be touched.
BROKEN_MEMBERS = [
    "size_type",
    "array_layout",
    "scratch_memory_space",
    "memory_traits",
]

# Templates whose FIRST parameter is an execution space.  Add your own with
# --exec-template; each is a judgement call about that template's interface.
KNOWN_EXEC_SPACE_TEMPLATES = [
    # Kokkos execution policies
    "Kokkos::RangePolicy", "RangePolicy",
    "Kokkos::MDRangePolicy", "MDRangePolicy",
    "Kokkos::TeamPolicy", "TeamPolicy",
    "Kokkos::WorkGraphPolicy", "WorkGraphPolicy",
]

# Deliberately NOT here, having been checked: every Intrepid2 class template
# takes a DEVICE, not an execution space --
#     template<typename Device, typename outputValueType, typename pointValueType>
#     class Basis { using DeviceType = Device;
#                   using ExecutionSpace = typename DeviceType::execution_space; ... }
# and the same for FunctionSpaceTools, CellTools, OrientationTools,
# ProjectionTools, RealSpaceTools and ArrayTools.  So Intrepid2::Foo<PHX::Device>
# is already right and must be left alone; rewriting it would hand Intrepid2 an
# execution space and silently drop a configured shared memory space.
# panzer::createIntrepid2Basis is the same case despite naming its parameter
# ExecutionSpace -- it forwards straight to Intrepid2::Basis.

# Templates whose relevant parameter is a DEVICE, so PHX::Device is already
# right there and must be left alone.  Checked against the sources: every
# Intrepid2 class template is declared
#     template<typename Device, typename outputValueType, typename pointValueType>
#     class Basis { using DeviceType = Device;
#                   using ExecutionSpace = typename DeviceType::execution_space; }
# and the same for FunctionSpaceTools, CellTools, OrientationTools,
# ProjectionTools, RealSpaceTools, ArrayTools and Cubature.  Matched as
# prefixes, so the concrete bases (Basis_HGRAD_*, Basis_HDIV_*, ...) are covered.
DEVICE_POSITION_PREFIXES = [
    "Intrepid2::",
    "PHX::KokkosViewFactory",
    "KokkosSparse::CrsMatrix",
    # forwards its parameter straight into Intrepid2::Basis, despite naming it
    # ExecutionSpace
    "panzer::createIntrepid2Basis",
    "createIntrepid2Basis",
    # Phalanx's own API that takes a device
    "PHX::getAllocationSize",
    "PHX::MDField",
    "PHX::print",
    "is_device",
    # Kokkos random pools take a DeviceType
    "Kokkos::Random_XorShift64_Pool",
    "Kokkos::Random_XorShift1024_Pool",
]

# Templates that take an execution space but where the right answer is a
# judgement call, not a rename: replacing PHX::Device with PHX::exec_space
# compiles, but silently keeps the execution space's default memory space, so
# the data does not follow a configured shared space.  Reported, never rewritten.
REVIEW_TEMPLATES = {
    "Tpetra::KokkosCompat::KokkosDeviceWrapperNode":
        "takes <ExecutionSpace, MemorySpace = ExecutionSpace::memory_space>; "
        "using PHX::ExecutionSpace alone drops a configured shared memory space, "
        "so decide between <PHX::ExecutionSpace> and "
        "<PHX::ExecutionSpace, PHX::MemorySpace>",
    "KokkosCompat::KokkosDeviceWrapperNode": "see Tpetra::KokkosCompat::KokkosDeviceWrapperNode",
    "KokkosDeviceWrapperNode": "see Tpetra::KokkosCompat::KokkosDeviceWrapperNode",
    # Tpetra's Node parameter is the same question in a different spelling
    "Tpetra::Vector": "the Node parameter determines where Tpetra allocates; decide "
                      "whether it should follow PHX::MemorySpace",
    "Tpetra::MultiVector": "see Tpetra::Vector",
    "Tpetra::CrsMatrix": "see Tpetra::Vector",
    "Tpetra::Map": "see Tpetra::Vector",
}

# An alias hides whether later uses are devices or execution spaces.  Catch the
# spellings that occur in practice, including `typename` and the plain
# using-declaration that just drags the name into scope.
# The headers that DEFINE these names.  Rewriting them turns the definitions
# into self-references (using ExecutionSpace = PHX::ExecutionSpace;), so they
# are skipped.  Matched on basename, so it works wherever Phalanx is checked out.
DEFINING_FILES = {
    "Phalanx_KokkosDeviceTypes.hpp",
    "Phalanx_config.hpp",
    "Phalanx_config.hpp.in",
}

# A file is worth opening if it mentions any name this script can act on.
TRIGGER_TOKENS = ("PHX::Device", "PHX::exec_space", "PHX::ExecSpace",
                  "PHX::mem_space", "PHX::MemSpace")

ALIAS_RE = re.compile(
    r"\busing\s+\w+\s*=\s*(?:typename\s+)?PHX::Device\s*;"
    r"|\btypedef\s+(?:typename\s+)?PHX::Device\s+\w+\s*;"
    r"|\busing\s+PHX::Device\s*;")
ANY_RE = re.compile(r"PHX::Device\b")


class Rewriter:
    def __init__(self, exec_templates):
        self.rules = []
        for member in BROKEN_MEMBERS:
            self.rules.append((
                f"PHX::Device::{member}",
                re.compile(r"\bPHX::Device::" + member + r"\b"),
                f"PHX::ExecutionSpace::{member}",
            ))
        # instance construction: PHX::Device() -> PHX::exec_space()
        self.rules.append((
            "PHX::Device() instance",
            re.compile(r"\bPHX::Device\s*\(\s*\)"),
            "PHX::ExecutionSpace()",
        ))
        # A functor telling Kokkos where to run:
        #     typedef PHX::Device execution_space;
        #     using execution_space = PHX::Device;
        # Unambiguous -- the alias names itself an execution space.
        self.rules.append((
            "typedef/using ... execution_space = PHX::Device",
            re.compile(r"\btypedef\s+(typename\s+)?PHX::Device(\s+execution_space\s*;)"),
            r"typedef \1PHX::ExecutionSpace\2",
        ))
        self.rules.append((
            "using execution_space = PHX::Device",
            re.compile(r"(\busing\s+execution_space\s*=\s*)(?:typename\s+)?PHX::Device(\s*;)"),
            r"\1PHX::ExecutionSpace\2",
        ))

        # PHX::Device::execution_space and ::memory_space are safe -- a
        # Kokkos::Device has both -- but they are two more names for the two
        # types being unified, and the short forms mean exactly the same thing.
        self.rules.append((
            "PHX::Device::execution_space -> PHX::ExecutionSpace",
            re.compile(r"\bPHX::Device::execution_space\b"),
            "PHX::ExecutionSpace",
        ))
        self.rules.append((
            "PHX::Device::memory_space -> PHX::MemorySpace",
            re.compile(r"\bPHX::Device::memory_space\b"),
            "PHX::MemorySpace",
        ))

        # The older spellings of the same two types, unified.  Unambiguous
        # renames: these names only ever meant the execution or memory space.
        for old, new in (("PHX::exec_space", "PHX::ExecutionSpace"),
                         ("PHX::ExecSpace",  "PHX::ExecutionSpace"),
                         ("PHX::mem_space",  "PHX::MemorySpace"),
                         ("PHX::MemSpace",   "PHX::MemorySpace")):
            self.rules.append((
                f"{old} -> {new}",
                re.compile(r"\b" + re.escape(old) + r"\b"),
                new,
            ))

        # panzer::HP::inst().teamPolicy<ScalarT[, Tag], PHX::Device> puts the
        # execution space LAST, so the first-parameter rules below cannot see
        # it.  Rewrite just that trailing argument.
        self.rules.append((
            "teamPolicy<..., PHX::Device> (trailing execution space)",
            # allow one level of nesting: SharedFieldMultTag<0> and friends
            re.compile(r"(\bteamPolicy\s*<(?:[^;<>]|<[^;<>]*>)*,\s*)PHX::Device(\s*>)"),
            r"\1PHX::ExecutionSpace\2",
        ))

        for tmpl in exec_templates:
            esc = re.escape(tmpl)
            self.rules.append((
                f"{tmpl}<PHX::Device",
                re.compile(r"(\b" + esc + r"\s*<\s*)PHX::Device\b"),
                r"\1PHX::ExecutionSpace",
            ))

    def apply(self, text):
        counts = {}
        for name, pattern, repl in self.rules:
            text, n = pattern.subn(repl, text)
            if n:
                counts[name] = counts.get(name, 0) + n
        return text, counts


def classify_residual(line):
    if ALIAS_RE.search(line):
        return "alias of PHX::Device -- follow it by hand"
    for tmpl, why in REVIEW_TEMPLATES.items():
        if re.search(r"\b" + re.escape(tmpl) + r"\s*<[^>]*PHX::Device", line):
            return f"{tmpl}: {why}"
    if re.search(r"(?:View|DynRankView)\s*<[^;]*PHX::Device", line):
        return None  # device position in a View: correct as is
    for prefix in DEVICE_POSITION_PREFIXES:
        if re.search(r"\b" + re.escape(prefix) + r"[\w:]*\s*<[^;]*PHX::Device", line):
            return None  # device position: correct as is
    if re.match(r"\s*(?:\*|//|/\*)", line):
        return None  # comment
    if re.search(r"(?:vector|array|list)\s*<\s*PHX::Device\s*>", line):
        return ("a container of PHX::Device -- if these are execution space "
                "instances (Kokkos::Experimental::partition_space returns "
                "those) it must become PHX::ExecutionSpace")
    return "unrecognized position -- check whether it wants a device or an execution space"


def iter_files(root, extensions, excludes):
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [d for d in dirnames
                       if d not in excludes and not d.startswith(".")]
        for fn in filenames:
            if os.path.splitext(fn)[1] in extensions and fn not in DEFINING_FILES:
                yield os.path.join(dirpath, fn)


def main():
    p = argparse.ArgumentParser(
        description="Respell PHX::Device as PHX::exec_space where it means an execution space.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__)
    p.add_argument("path", nargs="?", default=".", help="directory to scan (default: .)")
    p.add_argument("--apply", action="store_true",
                   help="write the changes; without this nothing is modified")
    p.add_argument("--check", action="store_true",
                   help="report only, for CI; exits non-zero if anything needs attention")
    p.add_argument("--ext", action="append", default=None,
                   help="file extension to scan, repeatable (default: %s)" % " ".join(DEFAULT_EXTENSIONS))
    p.add_argument("--exclude", action="append", default=["build", "CMakeFiles", "node_modules"],
                   help="directory name to skip, repeatable")
    p.add_argument("--exec-template", action="append", default=[],
                   help="extra template whose first parameter is an execution space, repeatable")
    p.add_argument("--quiet", action="store_true", help="suppress the per-site residual listing")
    args = p.parse_args()

    if args.apply and args.check:
        p.error("--apply and --check are mutually exclusive")

    extensions = args.ext if args.ext else DEFAULT_EXTENSIONS
    rewriter = Rewriter(KNOWN_EXEC_SPACE_TEMPLATES + args.exec_template)

    totals, residuals, files_changed, files_scanned = {}, [], 0, 0

    for path in iter_files(args.path, extensions, set(args.exclude)):
        try:
            original = open(path, encoding="utf-8", errors="replace").read()
        except OSError as exc:
            print(f"warning: cannot read {path}: {exc}", file=sys.stderr)
            continue
        files_scanned += 1
        # Any name this script can act on, not just PHX::Device -- the space
        # renames apply to files that never mention a device.
        if not any(tok in original for tok in TRIGGER_TOKENS):
            continue

        updated, counts = rewriter.apply(original)
        if counts:
            files_changed += 1
            for k, v in counts.items():
                totals[k] = totals.get(k, 0) + v
            if args.apply:
                with open(path, "w", encoding="utf-8") as fh:
                    fh.write(updated)

        for n, line in enumerate(updated.splitlines(), 1):
            if ANY_RE.search(line):
                why = classify_residual(line)
                if why:
                    residuals.append((path, n, line.strip(), why))

    verb = "rewrote" if args.apply else "would rewrite"
    print(f"scanned {files_scanned} files\n")
    if totals:
        print(f"{verb}, by rule:")
        for k in sorted(totals, key=lambda k: -totals[k]):
            print(f"  {totals[k]:>6}  {k}")
        print(f"  {sum(totals.values()):>6}  TOTAL in {files_changed} files\n")
    else:
        print("no mechanical rewrites needed\n")

    if residuals:
        print(f"{len(residuals)} occurrence(s) need a human:")
        if not args.quiet:
            for path, n, line, why in residuals:
                print(f"  {path}:{n}\n      {line}\n      -> {why}")
        else:
            by_reason = {}
            for _, _, _, why in residuals:
                by_reason[why] = by_reason.get(why, 0) + 1
            for why, count in sorted(by_reason.items(), key=lambda kv: -kv[1]):
                print(f"  {count:>6}  {why}")
        print()

    if not args.apply and totals:
        print("nothing was modified; re-run with --apply to write these changes")

    return 1 if (residuals and args.check) else 0


if __name__ == "__main__":
    sys.exit(main())
