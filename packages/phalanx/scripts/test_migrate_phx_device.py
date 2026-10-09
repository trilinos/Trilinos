#!/usr/bin/env python3
"""Unit tests for migrate_phx_device.py.

Run them with either of:

    python3 test_migrate_phx_device.py
    python3 -m unittest discover -s packages/phalanx/scripts

Python standard library only, like the script itself, so there is nothing to
install and this can run anywhere Phalanx is configured.

Most of these are regression tests: nearly every case below is a mistake the
script actually made while it was being validated against Phalanx, Panzer and
Drekar.  The comments say which, because the point of a rule is easier to see
from the thing it got wrong than from the rule itself.  The hardest part of
this script is not the rewriting, it is knowing whether a template parameter is
a device or an execution space -- that is a fact about each template, and the
tests below pin down the answer for the ones that are curated.
"""

import importlib.util
import os
import subprocess
import sys
import tempfile
import unittest


SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "migrate_phx_device.py")


def _load():
    spec = importlib.util.spec_from_file_location("migrate_phx_device", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mig = _load()


def rewrite(text):
    """Apply every rule once, the way the script does for a file."""
    return mig.Rewriter(mig.KNOWN_EXEC_SPACE_TEMPLATES).apply(text)[0]


class ExecutionSpaceRenames(unittest.TestCase):
    """The unconditional renames: these names only ever meant one type."""

    def test_old_space_names_collapse(self):
        for old, new in (("PHX::exec_space", "PHX::ExecutionSpace"),
                         ("PHX::ExecSpace", "PHX::ExecutionSpace"),
                         ("PHX::mem_space", "PHX::MemorySpace"),
                         ("PHX::MemSpace", "PHX::MemorySpace")):
            with self.subTest(old=old):
                self.assertEqual(rewrite(f"using T = {old};"),
                                 f"using T = {new};")

    def test_qualified_members_of_device(self):
        self.assertEqual(rewrite("PHX::Device::execution_space e;"),
                         "PHX::ExecutionSpace e;")
        self.assertEqual(rewrite("PHX::Device::memory_space m;"),
                         "PHX::MemorySpace m;")

    def test_members_a_kokkos_device_does_not_have(self):
        # A Kokkos::Device provides only execution_space, memory_space and
        # device_type, so these have to come off the execution space.
        for member in mig.BROKEN_MEMBERS:
            with self.subTest(member=member):
                self.assertEqual(rewrite(f"PHX::Device::{member} x;"),
                                 f"PHX::ExecutionSpace::{member} x;")

    def test_instance_construction(self):
        # A Kokkos::Device is not constructible as an execution space instance,
        # so PHX::Device().fence() has to become PHX::ExecutionSpace().
        self.assertEqual(rewrite("PHX::Device().fence();"),
                         "PHX::ExecutionSpace().fence();")

    def test_functor_execution_space_alias(self):
        # The alias names itself an execution space, so this is unambiguous.
        self.assertEqual(rewrite("typedef PHX::Device execution_space;"),
                         "typedef PHX::ExecutionSpace execution_space;")
        self.assertEqual(rewrite("using execution_space = PHX::Device;"),
                         "using execution_space = PHX::ExecutionSpace;")


class ExecutionSpaceSlots(unittest.TestCase):
    """Kokkos policies: the execution space can be in any argument position."""

    def test_policy_first_position(self):
        for policy in ("RangePolicy", "MDRangePolicy", "TeamPolicy"):
            with self.subTest(policy=policy):
                self.assertEqual(
                    rewrite(f"Kokkos::{policy}<PHX::Device> p;"),
                    f"Kokkos::{policy}<PHX::ExecutionSpace> p;")

    def test_policy_later_position(self):
        # REGRESSION: RangePolicy<LocalOrdinal, PHX::Device> is an index type
        # followed by an execution space.  Matching only the first argument
        # missed these, and a Kokkos::Device in a policy is not rejected -- it
        # fails the is_execution_space test, falls through to the work tag
        # (which only requires an empty type) and the policy silently runs on
        # the DEFAULT execution space.
        self.assertEqual(
            rewrite("Kokkos::RangePolicy<LocalOrdinal, PHX::Device> p;"),
            "Kokkos::RangePolicy<LocalOrdinal, PHX::ExecutionSpace> p;")

    def test_team_policy_trailing_argument_with_nested_brackets(self):
        # REGRESSION: panzer::HP::inst().teamPolicy<ScalarT[, Tag], PHX::Device>
        # puts the execution space LAST, and the real calls carry a nested
        # template argument.  Forbidding nested angle brackets matched 3 of 17
        # Panzer sites.
        self.assertEqual(
            rewrite("auto p = HP::inst().teamPolicy<ScalarT,"
                    "SharedFieldMultTag<0>,PHX::Device>(n);"),
            "auto p = HP::inst().teamPolicy<ScalarT,"
            "SharedFieldMultTag<0>,PHX::ExecutionSpace>(n);")

    def test_team_policy_trailing_argument_without_tag(self):
        self.assertEqual(
            rewrite("auto p = HP::inst().teamPolicy<ScalarT,PHX::Device>(n);"),
            "auto p = HP::inst().teamPolicy<ScalarT,PHX::ExecutionSpace>(n);")


class DeviceSlots(unittest.TestCase):
    """The reverse direction: templates whose parameter really is a device."""

    def test_device_positions_are_left_alone(self):
        # The whole point of the change: PHX::Device in a device slot is
        # correct and must survive untouched.
        for line in ("Kokkos::View<double*,PHX::Device> v;",
                     "Intrepid2::Basis<PHX::Device,double,double> b;",
                     "Kokkos::Random_XorShift64_Pool<PHX::Device> pool;"):
            with self.subTest(line=line):
                self.assertEqual(rewrite(line), line)

    def test_execution_space_in_a_device_slot_becomes_a_device(self):
        # REGRESSION: while PHX::Device WAS an execution space, both spellings
        # named one type and an execution space sitting in a device slot was
        # invisible.  Once they became different types, mixing them produced
        # Basis<Kokkos::Serial> and Basis<Kokkos::Device<Serial,HostSpace>> --
        # unrelated types being assigned to each other.
        cases = [
            ("Intrepid2::Basis<PHX::ExecutionSpace,double,double> b;",
             "Intrepid2::Basis<PHX::Device,double,double> b;"),
            ("Intrepid2::FunctionSpaceTools<PHX::ExecutionSpace> fst;",
             "Intrepid2::FunctionSpaceTools<PHX::Device> fst;"),
            ("createIntrepid2Basis<PHX::ExecutionSpace,double,double>(a,b,c);",
             "createIntrepid2Basis<PHX::Device,double,double>(a,b,c);"),
            ("panzer::createIntrepid2Basis<PHX::ExecutionSpace,double,double>(a,b,c);",
             "panzer::createIntrepid2Basis<PHX::Device,double,double>(a,b,c);"),
            ("b_->getIntrepid2Basis<PHX::ExecutionSpace,double,double>();",
             "b_->getIntrepid2Basis<PHX::Device,double,double>();"),
        ]
        for before, after in cases:
            with self.subTest(before=before):
                self.assertEqual(rewrite(before), after)

    def test_cubature_factory_member_call(self):
        # Intrepid2::DefaultCubatureFactory::create takes a DeviceType, and is
        # reached through a local factory object, so there is no type name to
        # anchor on.
        self.assertEqual(
            rewrite("cubature_factory.create<PHX::ExecutionSpace,double,double>(t,d);"),
            "cubature_factory.create<PHX::Device,double,double>(t,d);")

    def test_qualified_name_is_not_truncated(self):
        # REGRESSION, and the dangerous one: matching the PREFIX of
        # PHX::ExecutionSpace::size_type would yield PHX::Device::size_type --
        # exactly the member a Kokkos::Device does not have.  Caught by
        # PHX::print, which takes any type at all.
        for line in ("std::cout << PHX::print<PHX::ExecutionSpace::size_type>();",
                     "Intrepid2::Basis<PHX::ExecutionSpace::size_type> odd;"):
            with self.subTest(line=line):
                self.assertNotIn("PHX::Device::size_type", rewrite(line))

    def test_templates_taking_any_type_are_left_alone(self):
        # PHX::print<T> and is_device<T> accept an arbitrary type, so an
        # execution space is a perfectly good argument and must not be
        # rewritten.  Asking "is PHX::ExecutionSpace a device" is a legitimate
        # question with the answer "no".
        for line in ("std::cout << PHX::print<PHX::ExecutionSpace>();",
                     "static_assert(is_device<PHX::ExecutionSpace>::value);"):
            with self.subTest(line=line):
                self.assertEqual(rewrite(line), line)


class ResidualClassification(unittest.TestCase):
    """What the script refuses to rewrite, and what it says about it."""

    def test_aliases_are_reported_not_followed(self):
        for line in ("using DeviceSpace = PHX::Device;",
                     "typedef PHX::Device DeviceSpace;",
                     "  using PHX::Device;"):
            with self.subTest(line=line):
                self.assertIn("alias", mig.classify_residual(line))
                self.assertEqual(rewrite(line), line)

    def test_containers_are_reported_distinctly(self):
        # REGRESSION: std::vector<PHX::Device> streams was left alone while the
        # adjacent PHX::Device() became PHX::ExecutionSpace(), which compiled
        # only because the two were the same type at the time.
        why = mig.classify_residual("std::vector<PHX::Device> streams;")
        self.assertIn("container", why)

    def test_tpetra_node_is_a_judgement_call(self):
        for line in ("KokkosDeviceWrapperNode<PHX::Device> node;",
                     "Tpetra::Vector<double,LO,GO,PHX::Device> v;"):
            with self.subTest(line=line):
                self.assertIsNotNone(mig.classify_residual(line))
                self.assertEqual(rewrite(line), line)

    def test_device_positions_and_comments_are_not_residuals(self):
        # REGRESSION: reporting these buried the real cases.  Panzer's count
        # fell from 357 to 48 once they were classified properly.
        for line in ("Intrepid2::Basis<PHX::Device,double,double> b;",
                     "  // PHX::Device is a device here",
                     "  * PHX::Device in a doxygen comment"):
            with self.subTest(line=line):
                self.assertIsNone(mig.classify_residual(line))

    def test_unknown_template_is_reported(self):
        why = mig.classify_residual("SomethingUnknown<PHX::Device> x;")
        self.assertIn("unrecognized", why)

    def test_curated_device_slots_are_not_residuals(self):
        # REGRESSION: the rewriter learned these device slots but the
        # classifier did not, so a fully migrated tree still reported 44
        # "unrecognized" lines -- all of them correct as written.  The residual
        # list is only useful if a human can read it.
        for line in ("b_->getIntrepid2Basis<PHX::Device,double,double>();",
                     "ic = cubature_factory.create<PHX::Device,double,double>(t,o);",
                     "using S = Kokkos::UnorderedMap<Kokkos::pair<int,int>,void,PHX::Device>;",
                     "using ft = FieldTraits<double,C,P,PHX::Device>;",
                     "static_assert(RankCount<C,PHX::Device,P,D>::value == 3,\"x\");"):
            with self.subTest(line=line):
                self.assertIsNone(mig.classify_residual(line))

    def test_device_trait_comparisons_are_not_residuals(self):
        # Phalanx' own MDField plumbing compares against, and defaults into,
        # something already named a device.
        for line in ("static_assert(std::is_same<typename ft::device,PHX::Device>::value,\"x\");",
                     "using device = typename std::conditional<"
                     "!std::is_same<typename prop::device, void>::value,"
                     "typename prop::device, PHX::Device>::type;"):
            with self.subTest(line=line):
                self.assertIsNone(mig.classify_residual(line))

    def test_string_literals_are_not_code(self):
        # A diagnostic message that names the type, seen on its own line.
        self.assertIsNone(
            mig.classify_residual('              "panzer: PHX::Device would repeat "'))

    def test_continuation_lines_are_reported_distinctly(self):
        # The scan is line by line, so a Kokkos::View split across lines shows
        # up as a bare "PHX::Device," with no template in sight.  Saying so
        # beats calling it unrecognized.
        why = mig.classify_residual("               PHX::Device,")
        self.assertIsNotNone(why)
        self.assertIn("split across lines", why)


class FileSelection(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = self.tmp.name
        self.addCleanup(self.tmp.cleanup)

    def write(self, relpath, text):
        path = os.path.join(self.root, relpath)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
        return path

    def found(self, **kwargs):
        kwargs.setdefault("extensions", mig.DEFAULT_EXTENSIONS)
        kwargs.setdefault("excludes", {"build"})
        return {os.path.relpath(p, self.root)
                for p in mig.iter_files(self.root, **kwargs)}

    def test_recurses_into_subdirectories(self):
        self.write("a.hpp", "")
        self.write("deep/nested/b.cpp", "")
        self.assertEqual(self.found(), {"a.hpp", os.path.join("deep", "nested", "b.cpp")})

    def test_honours_extensions(self):
        self.write("a.hpp", "")
        self.write("README.md", "")
        self.assertEqual(self.found(), {"a.hpp"})
        self.assertEqual(self.found(extensions=[".md"]), {"README.md"})

    def test_skips_excluded_and_hidden_directories(self):
        self.write("keep.hpp", "")
        self.write("build/skip.hpp", "")
        self.write(".git/skip.hpp", "")
        self.assertEqual(self.found(), {"keep.hpp"})

    def test_skips_the_headers_that_define_the_names(self):
        # REGRESSION: the script rewrote its own defining header into
        # self-references -- `using ExecutionSpace = PHX::ExecutionSpace;`.
        self.assertIn("Phalanx_KokkosDeviceTypes.hpp", mig.DEFINING_FILES)
        self.write("Phalanx_KokkosDeviceTypes.hpp", "")
        self.write("other.hpp", "")
        self.assertEqual(self.found(), {"other.hpp"})


class CommandLine(unittest.TestCase):
    """End-to-end, including the exit codes a CI gate depends on."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = self.tmp.name
        self.addCleanup(self.tmp.cleanup)

    def write(self, name, text):
        path = os.path.join(self.root, name)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
        return path

    def read(self, name):
        with open(os.path.join(self.root, name), encoding="utf-8") as fh:
            return fh.read()

    def run_script(self, *args):
        return subprocess.run([sys.executable, SCRIPT, self.root, *args],
                              capture_output=True, text=True)

    def test_dry_run_reports_without_writing(self):
        before = "auto x = PHX::exec_space();\n"
        self.write("a.hpp", before)
        res = self.run_script()
        self.assertEqual(res.returncode, 0)
        self.assertEqual(self.read("a.hpp"), before, "dry run must not write")
        self.assertIn("nothing was modified", res.stdout)

    def test_apply_writes(self):
        self.write("a.hpp", "auto x = PHX::exec_space();\n")
        res = self.run_script("--apply")
        self.assertEqual(res.returncode, 0)
        self.assertEqual(self.read("a.hpp"), "auto x = PHX::ExecutionSpace();\n")

    def test_apply_is_idempotent(self):
        self.write("a.hpp", "auto x = PHX::exec_space();\n")
        self.run_script("--apply")
        once = self.read("a.hpp")
        self.run_script("--apply")
        self.assertEqual(self.read("a.hpp"), once)

    def test_check_fails_on_pending_mechanical_rewrites(self):
        # REGRESSION: --check returned non-zero only for residuals, so a tree
        # still full of PHX::exec_space exited 0 and a CI gate let unmigrated
        # code through -- while the help text promised otherwise.
        self.write("a.hpp", "auto x = PHX::exec_space();\n")
        self.assertEqual(self.run_script("--check", "--quiet").returncode, 1)

    def test_check_fails_on_residuals(self):
        self.write("a.hpp", "using DeviceSpace = PHX::Device;\n")
        self.assertEqual(self.run_script("--check", "--quiet").returncode, 1)

    def test_check_passes_on_a_migrated_tree(self):
        self.write("a.hpp",
                   "auto x = PHX::ExecutionSpace();\n"
                   "Kokkos::View<double*,PHX::Device> v;\n"
                   "Intrepid2::Basis<PHX::Device,double,double> b;\n")
        res = self.run_script("--check", "--quiet")
        self.assertEqual(res.returncode, 0, res.stdout + res.stderr)

    def test_check_and_apply_are_mutually_exclusive(self):
        self.write("a.hpp", "")
        self.assertNotEqual(self.run_script("--check", "--apply").returncode, 0)

    def test_file_gate_opens_files_without_the_word_device(self):
        # REGRESSION: the gate only opened files containing "PHX::Device", so a
        # file using nothing but PHX::exec_space was never read.  That reported
        # 28 of Drekar's 99 uses and looked like a success.
        self.write("a.hpp", "using T = PHX::exec_space;\n")
        self.run_script("--apply")
        self.assertEqual(self.read("a.hpp"), "using T = PHX::ExecutionSpace;\n")

    def test_file_gate_opens_already_renamed_files(self):
        # The same trap one step later: a file holding only PHX::ExecutionSpace
        # can still have it sitting in a device slot.
        self.write("a.hpp", "Intrepid2::Basis<PHX::ExecutionSpace> b;\n")
        self.run_script("--apply")
        self.assertEqual(self.read("a.hpp"),
                         "Intrepid2::Basis<PHX::Device> b;\n")

    def test_apply_preserves_bytes_that_are_not_utf8(self):
        # REGRESSION: the file was read with errors="replace", so a byte that is
        # not valid UTF-8 -- a Latin-1 name in a comment, say -- was rewritten
        # as U+FFFD on disk.  The script must not corrupt what it did not come
        # to change.
        path = os.path.join(self.root, "a.hpp")
        raw = b"// Fran\xe7ois\nauto x = PHX::exec_space();\n"
        with open(path, "wb") as fh:
            fh.write(raw)
        self.run_script("--apply")
        with open(path, "rb") as fh:
            after = fh.read()
        self.assertIn(b"\xe7", after, "the non-UTF-8 byte was not preserved")
        self.assertNotIn("\ufffd".encode("utf-8"), after)
        self.assertIn(b"PHX::ExecutionSpace", after, "the rewrite still happened")

    def test_apply_preserves_crlf_line_endings(self):
        # REGRESSION: text-mode writing turned CRLF into LF, so every line of a
        # Windows-style file showed up as changed.
        path = os.path.join(self.root, "a.hpp")
        with open(path, "wb") as fh:
            fh.write(b"auto x = PHX::exec_space();\r\nint y = 0;\r\n")
        self.run_script("--apply")
        with open(path, "rb") as fh:
            after = fh.read()
        self.assertEqual(after.count(b"\r\n"), 2, "CRLF endings were not kept")
        self.assertIn(b"PHX::ExecutionSpace", after)

    def test_extra_exec_template_can_be_supplied(self):
        self.write("a.hpp", "MyPolicy<PHX::Device> p;\n")
        self.run_script("--apply", "--exec-template", "MyPolicy")
        self.assertEqual(self.read("a.hpp"), "MyPolicy<PHX::ExecutionSpace> p;\n")


if __name__ == "__main__":
    unittest.main(verbosity=2)
