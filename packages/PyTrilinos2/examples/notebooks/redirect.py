"""IPython extension that redirects output"""
from PyTrilinos2 import Teuchos


global _extension_enabled
_extension_enabled = False

global red
red = Teuchos.ostream_redirect()


def start_capture():
    global red
    red.__enter__()


def stop_capture():
    global red
    red.__exit__()


def load_ipython_extension(ip):
    global _extension_enabled
    if _extension_enabled:
        return
    ip.events.register('pre_execute', start_capture)
    ip.events.register('post_execute', stop_capture)
    _extension_enabled = True


def unload_ipython_extension(ip):
    global _extension_enabled
    if not _extension_enabled:
        return

    ip.events.unregister('pre_execute', start_capture)
    ip.events.unregister('post_execute', stop_capture)
    # start_capture was called in pre_execute
    # after unregister we need to call it explicitly:
    stop_capture()
    _extension_enabled = False