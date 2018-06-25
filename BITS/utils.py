import os
import sys
import subprocess
from logzero import logger


def run_command(command):
    """
    General-purpose command executer.
    """
    
    try:
        out = subprocess.check_output(command, shell=True)
    except subprocess.CalledProcessError as proc:
        logger.error(proc.output.decode('utf-8'))
        sys.exit(1)
    else:
        return out.decode('utf-8')


# XXX: TODO: use a pipeflow management system
# (or use "-sync y" option (wait until job complete) of SGE)
def qsub_script(in_script_name, out_script_name, args=(), job_name="run_script", out_log="sge.log", n_core=1, run_directly=False):
    """
    Job submitting function with BITS' template scripts.
    <in_script_name> is the original (not SGE version) script file to be run located in the install path.
    <out_script_name> is a script file to be submitted with headers for SGE.
    <args> is a tuple of arguments for the script (must be in order of use).
    If <run_directly> is True, the script will be run without SGE.
    """

    # Check the template script file
    install_script_root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scripts")   # TODO: where should "install" scripts?
    script_path = os.path.join(install_script_root, in_script_name)
    if not os.path.isfile(script_path):
        logger.error("Script file not found")
        sys.exit(1)

    # Replace parameters and arguments in the template
    command = f"sed -e 's/JOB_NAME/{job_name}/' -e 's/OUT_LOG/{out_log}/' -e 's/N_CORE/{n_core}/' {script_path}"
    for i, arg in enumerate(args, 1):
        command += f" | sed -e 's/\${i}/{arg}/g'"
    command += f" > {out_script_name}"
    run_command(command)

    if run_directly:   # Run without SGE
        print(run_command(f"bash {out_script_name}"))   # TODO: add source ~/.bash_profile before?
    else:   # Run with SGE
        run_qsub_script = os.path.join(install_script_root, "run_qsub.sh")
        print(run_command(f"bash {run_qsub_script} {out_script_name}"))
