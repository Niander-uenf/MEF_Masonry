"""Nome antigo do pacote, até a versão 0.10. Desde a 1.0, o programa se chama
MEF_Masonry e o pacote, `mef_masonry`.

Este apelido mantém funcionando os roteiros escritos para a versão antiga
(`import mef_bcepq`, `from mef_bcepq.parede import ModeloParede` e assim por
diante): cada `mef_bcepq.X` devolve o próprio módulo `mef_masonry.X`. Ele será
retirado numa versão futura; troque `mef_bcepq` por `mef_masonry` nos roteiros.
"""
import importlib
import importlib.abc
import importlib.util
import sys
import warnings

warnings.warn("O pacote mef_bcepq agora se chama mef_masonry; o nome antigo é só um apelido e será "
              "retirado numa versão futura.", DeprecationWarning, stacklevel=2)

import mef_masonry as _novo  # noqa: E402
from mef_masonry import *  # noqa: E402,F401,F403

__version__ = _novo.__version__


class _Apelido(importlib.abc.MetaPathFinder, importlib.abc.Loader):
    """Resolve mef_bcepq.X como o módulo mef_masonry.X (o mesmo objeto)."""

    def find_spec(self, nome, caminho=None, alvo=None):
        if nome.startswith("mef_bcepq."):
            return importlib.util.spec_from_loader(nome, self)
        return None

    def create_module(self, spec):
        return importlib.import_module("mef_masonry" + spec.name[len("mef_bcepq"):])

    def exec_module(self, modulo):
        pass


if not any(isinstance(f, _Apelido) for f in sys.meta_path):
    sys.meta_path.insert(0, _Apelido())
