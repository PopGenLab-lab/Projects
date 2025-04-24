import sys
from .pyrelfit import pyrelfit as cli

"""
pyrelfit.__main__
~~~~~~~~~~~~~~~~~~~~~
The main entry point for the command line interface.
If installed invoke as ``pyrelfit`` or
as ``python -m pyrelfit`` if not installed.
"""

if __name__ == "__main__":
    sys.exit(cli())