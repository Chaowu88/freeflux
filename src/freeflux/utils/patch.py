'''Patch OpenOpt.'''


__author__ = 'Chao Wu'


from importlib.util import find_spec
from pathlib import Path
import re


OPENOPT_MODULES = ['ooIter', 'result', 'runProbSolver']


def _patch_module(spec_dir, mod):
    mod_file = Path(spec_dir) / 'kernel' / f'{mod}.py'
    
    with open(mod_file) as f:
        old_code = f.read()
    
    new_code = re.sub(r'\bclock\b', 'perf_counter', old_code)
    
    if new_code != old_code:
        with open(mod_file, 'w') as f:
            f.write(new_code)
        
        print(f'{mod_file} patched.')
    

def apply_patch():
    spec_dir = find_spec('openopt').submodule_search_locations[0]
    
    for mod in OPENOPT_MODULES:
        _patch_module(spec_dir, mod)