"""MEF_Masonry: elementos finitos para alvenaria de blocos em temperatura (interface no navegador).

Uso: python app.py            (abre o navegador em http://127.0.0.1:8765/)
     python app.py --porta 9000 --pasta minhas_analises --sem-navegador
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mef_masonry.app.servidor import principal  # noqa: E402

if __name__ == "__main__":
    # a sensibilidade roda a parede em processos separados; no executável do Windows
    # (PyInstaller), os processos filhos passam por aqui e param no freeze_support
    import multiprocessing
    multiprocessing.freeze_support()
    principal()
