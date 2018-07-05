import os
import sys
import subprocess
from interval import interval
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


def make_line(x0, y0, x1, y1, col, width):
    """
    For Plotly.
    Create a line-shape object for Plotly.
    """

    return {'type': 'line', 'xref': 'x', 'yref': 'y',
            'x0': x0, 'y0': y0, 'x1': x1, 'y1': y1,
            'line': {'color': col, 'width': width}}


def interval_len(intvls):
    """
    For pyinterval (https://pyinterval.readthedocs.io/en/latest/), a module for interval.
    Return the sum of the lengths of the intervals in the given interval object.
    """
    
    ret = 0
    for intvl in intvls.components:
        ret += intvl[0][1] - intvl[0][0] + 1
    return ret


def subtract_interval(a_interval, b_interval, length_threshold=0):
    """
    For pyinterval.
    Calculate A - B, where A and B are interval objects.
    Afterwards, remove remaining intervals whose length are less than <length_threshold>.
    """
    
    ret_interval = interval()
    intersection = a_interval & b_interval
    a_index = 0
    i_index = 0

    flag_load_new_interval = True
    while a_index < len(a_interval) and i_index < len(intersection):
        if flag_load_new_interval:
            a_intvl = a_interval[a_index]
        else:
            flag_load_new_interval = False

        if a_intvl[1] < intersection[i_index][0]:
            ret_interval |= interval(a_intvl)
            a_index += 1
        elif a_intvl[0] > intersection[i_index][1]:
            i_index += 1
        else:
            a_start, a_end = a_intvl
            i_start, i_end = intersection[i_index]
            start_contained = True if min(a_start, i_start) == i_start else False
            end_contained = True if max(a_end, i_end) == i_end else False
            if start_contained:
                if end_contained:
                    a_index += 1
                    i_index += 1
                else:
                    # the tail interval that did not intersect
                    # with this B-interval (but maybe with the next B-interval)
                    a_intvl = interval[i_end + 1, a_end]
                    flag_load_new_interval = False
                    i_index += 1
            else:
                if end_contained:
                    ret_interval |= interval[a_start, i_start - 1]
                    a_index += 1
                else:
                    ret_interval |= interval[a_start, i_start - 1]
                    a_intvl = interval[i_end + 1, a_end]
                    flag_load_new_interval = False
                    i_index += 1
    if not flag_load_new_interval:   # add a tail interval if it exists
        ret_interval |= interval(a_intvl)
        a_index += 1
    while a_index < len(a_interval):   # add remaining intervals in A
        ret_interval |= interval(a_interval[a_index])
        a_index += 1

    ret_interval_long = interval()   # length threshold
    for intvl in ret_interval.components:
        if interval_len(intvl) >= length_threshold:
            ret_interval_long |= intvl

    return ret_interval_long
