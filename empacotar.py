"""Gera o executável do programa com o PyInstaller (rode no sistema de destino).

Windows: python empacotar.py   ->  dist\\MEF_Masonry\\MEF_Masonry.exe
macOS/Linux: python3 empacotar.py  ->  dist/MEF_Masonry/MEF_Masonry

Instala o PyInstaller se faltar. O executável abre o navegador na interface;
a janela de console mostra o endereço e fica aberta enquanto o programa roda.
As análises são gravadas em resultados_app/, ao lado de onde o programa for aberto.
"""
import os
import subprocess
import sys

AQUI = os.path.dirname(os.path.abspath(__file__))


def main():
    try:
        import PyInstaller  # noqa: F401
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyinstaller"])
    import PyInstaller.__main__
    web = os.path.join(AQUI, "mef_masonry", "app", "web")
    extras = []
    try:                       # solver mais rápido, se estiver instalado (recomendado para a W15)
        import pypardiso  # noqa: F401
        extras = ["--collect-all", "pypardiso"]
        print("pypardiso encontrado: será incluído no executável.")
    except ImportError:
        print("pypardiso não instalado: o executável usará o SuperLU do SciPy.")
    PyInstaller.__main__.run(extras + [
        os.path.join(AQUI, "app.py"),
        "--name", "MEF_Masonry",
        "--onedir", "--console", "--noconfirm", "--clean",
        "--add-data", f"{web}{os.pathsep}mef_masonry/app/web",
        "--hidden-import", "scipy.sparse.linalg",
        "--hidden-import", "scipy.ndimage",
        "--hidden-import", "scipy.interpolate",
        "--hidden-import", "matplotlib.backends.backend_agg",
        "--hidden-import", "mef_masonry.elementos3d",
        "--hidden-import", "mef_masonry.naolinear",
        "--hidden-import", "mef_masonry.visual",
        "--hidden-import", "mef_masonry.parede",
        "--hidden-import", "mef_masonry.dados_parede",
        "--hidden-import", "mef_masonry.biblioteca",
        "--hidden-import", "mef_masonry.sensibilidade",
        "--distpath", os.path.join(AQUI, "dist"),
        "--workpath", os.path.join(AQUI, "build"),
        "--specpath", os.path.join(AQUI, "build"),
    ])
    print("\nPronto:", os.path.join(AQUI, "dist", "MEF_Masonry"))


if __name__ == "__main__":
    main()
