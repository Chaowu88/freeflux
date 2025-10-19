__author__ = 'Chao Wu'
__version__ = '0.3.7'


from importlib.util import find_spec
from .utils.patch import apply_patch
if find_spec('openopt') is not None:
    apply_patch()

from .core.metabolite import Metabolite
from .core.reaction import Reaction
from .core.model import Model
from .core.emu import EMU
from .core.mdv import MDV, get_natural_MDV, get_substrate_MDV, conv
