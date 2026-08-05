
import subprocess
from pathlib import Path

def run_rscript(r_script, args):
    r_script = Path(r_script).resolve()
    cmd = ["Rscript", str(r_script)] + [str(a) for a in args]
    subprocess.check_call(cmd)
