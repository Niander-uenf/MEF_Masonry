"""Builds the executable of the program with PyInstaller (run it on the target system).

Windows: python build_exe.py   ->  dist\\MEF_Masonry\\MEF_Masonry.exe
macOS/Linux: python3 build_exe.py  ->  dist/MEF_Masonry/MEF_Masonry

Installs PyInstaller if it is missing. The executable opens the browser on the interface;
the console window shows the address and stays open while the program runs.
The analyses are saved in results_app/, next to where the program is opened.
"""
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    try:
        import PyInstaller  # noqa: F401
    except ImportError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyinstaller"])
    import PyInstaller.__main__
    web = os.path.join(HERE, "mef_masonry", "app", "web")
    extras = []
    try:                       # faster solver, if installed (recommended for W15)
        import pypardiso  # noqa: F401
        extras = ["--collect-all", "pypardiso"]
        print("pypardiso found: it will be included in the executable.")
    except ImportError:
        print("pypardiso not installed: the executable will use SciPy's SuperLU.")
    PyInstaller.__main__.run(extras + [
        os.path.join(HERE, "app.py"),
        "--name", "MEF_Masonry",
        "--onedir", "--console", "--noconfirm", "--clean",
        "--add-data", f"{web}{os.pathsep}mef_masonry/app/web",
        "--hidden-import", "scipy.sparse.linalg",
        "--hidden-import", "scipy.ndimage",
        "--hidden-import", "scipy.interpolate",
        "--hidden-import", "matplotlib.backends.backend_agg",
        "--hidden-import", "mef_masonry.elements3d",
        "--hidden-import", "mef_masonry.nonlinear",
        "--hidden-import", "mef_masonry.visual",
        "--hidden-import", "mef_masonry.wall",
        "--hidden-import", "mef_masonry.wall_data",
        "--hidden-import", "mef_masonry.library",
        "--hidden-import", "mef_masonry.sensitivity",
        "--distpath", os.path.join(HERE, "dist"),
        "--workpath", os.path.join(HERE, "build"),
        "--specpath", os.path.join(HERE, "build"),
    ])
    print("\nDone:", os.path.join(HERE, "dist", "MEF_Masonry"))


if __name__ == "__main__":
    main()
