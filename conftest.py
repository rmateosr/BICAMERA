# ABOUTME: Puts em/ on sys.path so pytest works from the repo root.
# ABOUTME: BIC_algorithm_ASM.py imports EMfunctions_ASM flat, the way a script does.
import os
import sys

# The EM driver runs as a script, so Python normally puts em/ at sys.path[0] and a
# plain `import EMfunctions_ASM` resolves. Under pytest the rootdir goes on the path
# instead, so the tests could only be collected by cd-ing into em/ first. This makes
# `pytest` work from anywhere in the repo.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "em"))
