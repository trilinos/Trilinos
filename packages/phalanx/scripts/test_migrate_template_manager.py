#!/usr/bin/env python3
"""Unit tests for migrate_template_manager.py.

Run them with either of:

    python3 test_migrate_template_manager.py
    python3 -m unittest discover -s packages/phalanx/scripts

Python standard library only, like the script itself.
"""

import contextlib
import importlib.util
import io
import os
import tempfile
import unittest


SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "migrate_template_manager.py")


def _load():
    spec = importlib.util.spec_from_file_location("migrate_template_manager", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mig = _load()


def rewrite(text):
    return mig.process(text)[0]


def reported(text):
    return [why for _, _, why in mig.process(text)[2]]


def run_quietly(*args):
    with contextlib.redirect_stdout(io.StringIO()) as out:
        status = mig.main(list(args))
    return status, out.getvalue()


class Renames(unittest.TestCase):

    def test_each_rename(self):
        self.assertEqual(rewrite("tm.getAsObject<T>();"), "tm.get<T>();")
        self.assertEqual(rewrite("tm.disableType<T>();"), "tm.disable<T>();")
        self.assertEqual(rewrite("tm.deleteType<T>();"), "tm.reset<T>();")

    def test_dependent_call_keeps_template_keyword(self):
        self.assertEqual(rewrite("(*f)->template getAsObject<EvalT>()->go();"),
                         "(*f)->template get<EvalT>()->go();")

    def test_longer_names_are_left_alone(self):
        for text in ("mygetAsObject<T>();", "getAsObjectFoo<T>();", "disableTypes<T>();"):
            self.assertEqual(rewrite(text), text)

    def test_only_template_calls_are_renamed(self):
        # A method of the same name that takes no template argument is some
        # other class's.
        self.assertEqual(rewrite("x.deleteType(3);"), "x.deleteType(3);")

    def test_idempotent(self):
        once = rewrite("tm.getAsObject<T>(); tm.deleteType<T>();")
        self.assertEqual(rewrite(once), once)

    def test_getAsBase_is_unchanged(self):
        self.assertEqual(rewrite("tm.getAsBase<T>();"), "tm.getAsBase<T>();")

    def test_crlf_is_preserved(self):
        self.assertEqual(rewrite("tm.getAsObject<T>();\r\nx;\r\n"), "tm.get<T>();\r\nx;\r\n")


class Reports(unittest.TestCase):

    def test_iterator(self):
        self.assertEqual(len(reported("Foo_TemplateManager<Traits>::iterator it = m.begin();")), 1)
        self.assertEqual(len(reported("typename Foo_TemplateManager<T<X>>::const_iterator it;")), 1)

    def test_iterator_of_a_vector_of_managers_is_not_reported(self):
        # This loop is over a std::vector, and needs no change.
        self.assertEqual(
            reported("std::vector< RCP<Foo_TemplateManager<Traits> > >::const_iterator eq;"), [])

    def test_vector_template_manager_is_another_class(self):
        self.assertEqual(reported("panzer::VectorTemplateManager<S,B,O>::iterator it;"), [])

    def test_field_manager_iterator(self):
        self.assertEqual(len(reported("typename PHX::FieldManager<Traits>::iterator it;")), 1)

    def test_rcp(self):
        self.assertEqual(len(reported("auto p = it.rcp();")), 1)

    def test_spelled_out_template_manager(self):
        why = reported("class A : public PHX::TemplateManager<S,B,F<_> > {};")
        self.assertEqual(len(why), 1)
        self.assertIn("#ifdef", why[0])

    def test_auto_iterators_are_guessed_from_the_name(self):
        for text in ("for (auto i=fm->begin(); ;)", "auto i = bcstm->begin();",
                     "auto i = my_manager.begin();", "auto i = factory.begin();"):
            self.assertEqual(len(reported(text)), 1, text)
        self.assertEqual(reported("auto i = values.begin();"), [])

    def test_plural_names_are_containers(self):
        for text in ("for (i=bc_field_managers_.begin(); ;)", "a.assign(rowManagers.begin(),x);",
                     "auto i = m_factories.begin();"):
            self.assertEqual(reported(text), [], text)

    def test_comments_are_not_reported(self):
        self.assertEqual(reported("// see PHX::TemplateManager<S,B,O>"), [])

    def test_one_report_per_line(self):
        self.assertEqual(
            len(reported("PHX::TemplateManager<S,B,O>::iterator it = tm.begin();")), 1)

    def test_line_number(self):
        self.assertEqual(mig.process("a;\nb;\nauto p = it.rcp();\n")[2][0][0], 3)


class FlagConditionals(unittest.TestCase):
    """Code that tests the flag has already been dealt with."""

    GUARDED = ("#ifdef PHALANX_HIDE_DEPRECATED_CODE\n"
               "public PHX::TemplateManager<S,F,B> {\n"
               "#else\n"
               "public PHX::TemplateManager<S,B,F<_> > {\n"
               "  x.getAsObject<T>();\n"
               "#endif\n")

    def test_left_alone_and_not_reported(self):
        new, counts, residuals = mig.process(self.GUARDED)
        self.assertEqual(new, self.GUARDED)
        self.assertEqual(counts, {})
        self.assertEqual(residuals, [])

    def test_applies_again_after_endif(self):
        self.assertEqual(rewrite(self.GUARDED + "x.getAsObject<T>();\n"),
                         self.GUARDED + "x.get<T>();\n")

    def test_nested_conditional_inside_the_flag(self):
        text = ("#ifndef PHALANX_HIDE_DEPRECATED_CODE\n#ifdef OTHER\nx.getAsObject<T>();\n"
                "#endif\nx.getAsObject<T>();\n#endif\n")
        self.assertEqual(rewrite(text), text)

    def test_other_conditionals_are_still_migrated(self):
        self.assertEqual(rewrite("#ifdef OTHER\nx.getAsObject<T>();\n#else\ny.deleteType<T>();\n#endif\n"),
                         "#ifdef OTHER\nx.get<T>();\n#else\ny.reset<T>();\n#endif\n")


class CommandLine(unittest.TestCase):

    def setUp(self):
        self.dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.dir.cleanup)
        self.file = self.write("a.cpp", "x.getAsObject<T>();\n")

    def write(self, name, text):
        path = os.path.join(self.dir.name, name)
        with open(path, "w", newline="") as fh:
            fh.write(text)
        return path

    def read(self, path):
        with open(path, newline="") as fh:
            return fh.read()

    def test_default_changes_nothing(self):
        status, out = run_quietly(self.dir.name)
        self.assertEqual(status, 0)
        self.assertIn("would rewrite", out)
        self.assertEqual(self.read(self.file), "x.getAsObject<T>();\n")

    def test_apply(self):
        run_quietly(self.dir.name, "--apply")
        self.assertEqual(self.read(self.file), "x.get<T>();\n")

    def test_check_fails_on_pending_rewrite(self):
        self.assertEqual(run_quietly(self.dir.name, "--check")[0], 1)

    def test_check_fails_on_residual_even_after_apply(self):
        self.write("b.cpp", "auto p = it.rcp();\n")
        self.assertEqual(run_quietly(self.dir.name, "--apply", "--check")[0], 1)

    def test_check_passes_when_clean(self):
        run_quietly(self.dir.name, "--apply")
        self.assertEqual(run_quietly(self.dir.name, "--check")[0], 0)

    def test_extensions_and_excludes(self):
        other = self.write("notes.txt", "x.getAsObject<T>();\n")
        skipped = self.write("skip.cpp", "x.getAsObject<T>();\n")
        own = self.write("Phalanx_TemplateManager_Def.hpp", "x.getAsObject<T>();\n")
        docs = self.write("Phalanx_TemplateManager_Typed.hpp", "x.getAsObject<T>();\n")
        run_quietly(self.dir.name, "--apply", "--exclude", "skip.cpp")
        for path in (other, skipped, own, docs):
            self.assertEqual(self.read(path), "x.getAsObject<T>();\n")

    def test_invalid_utf8_round_trips(self):
        path = os.path.join(self.dir.name, "latin.cpp")
        with open(path, "wb") as fh:
            fh.write(b"// caf\xe9\nx.getAsObject<T>();\n")
        run_quietly(self.dir.name, "--apply")
        with open(path, "rb") as fh:
            self.assertEqual(fh.read(), b"// caf\xe9\nx.get<T>();\n")

    def test_quiet_summarizes(self):
        self.write("b.cpp", "auto p = it.rcp();\n")
        _, out = run_quietly(self.dir.name, "--quiet")
        self.assertNotIn("b.cpp", out)
        self.assertIn("iterator.rcp()", out)


if __name__ == "__main__":
    unittest.main()
