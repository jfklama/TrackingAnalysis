import sys
from unittest.mock import MagicMock

# Prevent ImportError when utils.py does `from pyLCIO import ...`
# pyLCIO requires a full ILD software stack; we mock it for unit tests.
for _mod in ['pyLCIO', 'pyLCIO.IOIMPL', 'pyLCIO.IMPL', 'pyLCIO.EVENT', 'pyLCIO.ROOT']:
    sys.modules[_mod] = MagicMock()
