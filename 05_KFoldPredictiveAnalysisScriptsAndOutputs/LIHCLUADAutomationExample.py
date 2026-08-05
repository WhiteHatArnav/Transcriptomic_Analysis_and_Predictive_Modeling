import subprocess
import sys
import os
from datetime import datetime

scripts = [
    "/Users/arnavjoshi/Desktop/KFoldLassoCox/LUADKFoldLassoCoxRecurrence/"
    "LUAD_AllRaceTenFoldLassoCoxRecurrence.py",

    "/Users/arnavjoshi/Desktop/KFoldLassoCox/LUADKFoldLassoCoxRecurrence/"
    "LUAD_BlackRaceTenFoldLassoCoxRecurrence.py",

    "/Users/arnavjoshi/Desktop/KFoldLassoCox/LIHCKFoldLassoCoxRecurrence/"
    "LIHC_WhiteRaceTenFoldLassoCoxRecurrence.py",

    "/Users/arnavjoshi/Desktop/KFoldLassoCox/LIHCKFoldLassoCoxRecurrence/"
    "LIHC_AllRaceTenFoldLassoCoxRecurrence.py",
]

LOG_FILE = "KFold_Recurrence_Master_Run.log"

def timestamp():
    return datetime.now().strftime("%Y-%m-%d %H:%M:%S")

def run_and_stream(script_path):
    print(f"\n[{timestamp()}] STARTING: {script_path}\n", flush=True)

    with open(LOG_FILE, "a") as log:
        log.write(f"\n[{timestamp()}] STARTING: {script_path}\n")

        process = subprocess.Popen(
            [sys.executable, "-u", script_path],  # 🔥 UNBUFFERED MODE
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1
        )

        for line in process.stdout:
            print(line, end="", flush=True)   # live to VS Code
            log.write(line)                   # also log

        return_code = process.wait()

        log.write(
            f"\n[{timestamp()}] FINISHED: {script_path} "
            f"(exit code {return_code})\n"
        )

    return return_code

if __name__ == "__main__":
    print(f"[{timestamp()}] ===== MASTER RUN STARTED =====", flush=True)
    with open(LOG_FILE, "a") as log:
        log.write(f"\n[{timestamp()}] ===== MASTER RUN STARTED =====\n")

    for script in scripts:
        if not os.path.exists(script):
            print(f"[ERROR] Script not found: {script}", flush=True)
            continue

        rc = run_and_stream(script)

        if rc != 0:
            print(f"[ERROR] Script failed: {script}", flush=True)
        else:
            print(f"[{timestamp()}] COMPLETED: {script}", flush=True)

    print(f"[{timestamp()}] ===== ALL JOBS FINISHED =====", flush=True)
