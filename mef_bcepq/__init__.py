"""Old name of the package, up to version 0.10. Since 1.0 the program is called
MEF_Masonry and the package `mef_masonry`.

This alias keeps scripts written for the old version running
(`import mef_bcepq`): each `mef_bcepq.X` returns the module `mef_masonry.X`
itself. It will be removed in a future version; replace `mef_bcepq` with
`mef_masonry` in your scripts.
"""
import importlib
import importlib.abc
import importlib.util
import sys
import warnings

warnings.warn("The mef_bcepq package is now called mef_masonry; the old name is only an alias and will be "
              "removed in a future version.", DeprecationWarning, stacklevel=2)

import mef_masonry as _new  # noqa: E402
from mef_masonry import *  # noqa: E402,F401,F403

__version__ = _new.__version__


class _Alias(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    """Resolves mef_bcepq.X as the module mef_masonry.X (the same object)."""

    def find_spec(self, name, path=None, target=None):
        if name.startswith("mef_bcepq."):
            return importlib.util.spec_from_loader(name, self)
        return None

    def create_module(self, spec):
        return importlib.import_module("mef_masonry" + spec.name[len("mef_bcepq"):])

    def exec_module(self, module_):
        pass


if not any(isinstance(f, _Alias) for f in sys.meta_path):
    sys.meta_path.insert(0, _Alias())
