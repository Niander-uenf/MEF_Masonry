"""MEF_Masonry: finite elements for block masonry in fire (browser interface).

Usage: python app.py            (opens the browser at http://127.0.0.1:8765/)
       python app.py --port 9000 --folder my_analyses --no-browser
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mef_masonry.app.server import main  # noqa: E402

if __name__ == "__main__":
    # the sensitivity study runs the wall in separate processes; in the Windows executable
    # (PyInstaller), the child processes pass through here and stop at freeze_support
    import multiprocessing
    multiprocessing.freeze_support()
    main()
