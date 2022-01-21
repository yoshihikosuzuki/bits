from dataclasses import dataclass
from typing import Callable, Sequence, Optional, List, Mapping
from math import log10
from logzero import logger
from ._pickle import load_pickle, save_pickle
from ._proc import run_command


@dataclass(eq=False, frozen=True)
class Scheduler:
    """Utility for submitting scripts using a job scheduler.

    usage:
      > s = Scheduler("sge", "qsub", "all.q")
      > s.submit("sleep 1s", "sleep.sh", "sleep")

    optional variables:
      @ scheduler_name : Name of job scheduler to be used. "sge" or "slurm".
      @ submit_command : Command to submit a job using the scheduler
      @ queue_name     : Queue/partiion name managed by the scheduler.
      @ prefix_command : Commands to be added before the script.
                         Use for environment-related commands:
                         e.g. "source /path/to/venv/bin/activate"
    """
    scheduler_name: str = "sge"
    submit_command: str = "qsub"
    queue_name: Optional[str] = "all.q"
    prefix_command: Optional[str] = None

    def __post_init__(self):
        assert self.scheduler_name in ("sge", "slurm"), "Unsupported scheduler"

    def submit(self,
               script: str,
               out_script_fname: str,
               job_name: str,
               log_fname: str = "log",
               n_core: int = 1,
               max_cpu_hour: Optional[int] = None,
               max_mem_gb: Optional[int] = None,
               depend_job_ids: Optional[List[str]] = None,
               wait: bool = False) -> str:
        """Generate and submit a script file from a script string.

        positional arguments:
          @ script           : Script (not a file).
          @ out_script_fname : Name of the script file generated and submitted.
          @ job_name         : Display name of the job.
          @ log_fname        : Name of the log file.
          @ n_core           : Number of cores used for the job.
          @ max_cpu_hour     : In hours.
          @ max_mem_gb       : In GB.
          @ depend_job_ids   : Wait to start the script until these jobs finish.
          @ wait             : If True, wait until the script finishes.

        return value:
          @ job_id : Of the script submitted.
        """
        logger.info(f"Submit command(s):\n{script}")
        header = (self.gen_sge_header if self.scheduler_name == "sge"
                  else self.gen_slurm_header)(job_name,
                                              log_fname,
                                              n_core,
                                              max_cpu_hour,
                                              max_mem_gb,
                                              depend_job_ids,
                                              wait)
        with open(out_script_fname, "w") as f:
            f.write('\n'.join(filter(None, ["#!/bin/bash",
                                            header + "\n",
                                            self.prefix_command,
                                            script + "\n"])))
        return (run_command(f"{self.submit_command} {out_script_fname}")
                .split()[2 if self.scheduler_name == "sge" else -1])

    def gen_sge_header(self,
                       job_name: str,
                       log_fname: str,
                       n_core: int,
                       max_cpu_hour: Optional[int],
                       max_mem_gb: Optional[int],
                       depend_job_ids: Optional[List[str]],
                       wait: bool) -> str:
        return '\n'.join(filter(None,
                                [f"#$ -N {job_name}",
                                 f"#$ -o {log_fname}",
                                 "#$ -j y",
                                 "#$ -S /bin/bash",
                                 "#$ -cwd",
                                 f"#$ -q {self.queue_name}"
                                 if self.queue_name is not None else "",
                                 f"#$ -pe smp {n_core}",
                                 f"#$ -l h_cpu={max_cpu_hour}"
                                 if max_cpu_hour is not None else "",
                                 f"#$ -l mem_total={max_mem_gb}G"
                                 if max_mem_gb is not None else "",
                                 f"#$ -hold_jid {','.join(depend_job_ids)}"
                                 if depend_job_ids is not None else "",
                                 f"#$ -sync {'y' if wait else 'n'}"]))

    def gen_slurm_header(self,
                         job_name: str,
                         log_fname: str,
                         n_core: int,
                         max_cpu_hour: Optional[int],
                         max_mem_gb: Optional[int],
                         depend_job_ids: Optional[List[str]],
                         wait: bool) -> str:
        return '\n'.join(filter(None,
                                [f"#SBATCH -J {job_name}",
                                 f"#SBATCH -o {log_fname}",
                                 f"#SBATCH -p {self.queue_name}"
                                 if self.queue_name is not None else "",
                                 "#SBATCH -n 1",
                                 "#SBATCH -N 1",
                                 f"#SBATCH -c {n_core}",
                                 f"#SBATCH -t '{max_cpu_hour}:00:00'"
                                 if max_cpu_hour is not None else "",
                                 f"#SBATCH --mem={max_mem_gb}G"
                                 if max_mem_gb is not None else "",
                                 f"#SBATCH -d afterany:{','.join(depend_job_ids)}"
                                 if depend_job_ids is not None else "",
                                 "#SBATCH --wait" if wait else ""]))


def run_distribute(func: Callable,
                   args: Sequence,
                   shared_args: Mapping,
                   scheduler: Optional[Scheduler],
                   n_distribute: int,
                   n_core: int,
                   max_cpu_hour: Optional[int] = None,
                   max_mem_gb: Optional[int] = None,
                   tmp_dname: str = "tmp_distribute",
                   job_name: str = "job",
                   out_fname: str = "out.pkl",
                   log_level: str = "info") -> List:
    """Distribute `[func(arg, **shared_args, n_core=n_core) for arg in args]`
    into `n_distribute` jobs (`n_core` per job) with a job scheduler.

    positional arguments:
      @ func         : Type must be just like above.
      @ args         : List of arguments to be distributed.
      @ shared_args  : Dict of arguments shared among every job.
                       Add '_' before names of globally shared arguments.
      @ scheduler    : A Scheduler object.
      @ n_distribute : Number of jobs.
      @ n_core       : Number of cores per job.

    optional arguments:
      @ tmp_dname : Directory name for intermediate files.
      @ job_name  : Display job name.
      @ out_fname : Output file name.
      @ log_level : Log level. Must be one of {"info", "debug"}.
    """
    assert isinstance(args, Sequence), \
        "`args` must be a Sequence object"
    assert isinstance(shared_args, Mapping), \
        "`shared_args` must be a Mapping object"
    assert log_level in ("info", "debug"), "Invalid name"
    run_command(f"mkdir -p {tmp_dname}; rm -f {tmp_dname}/*")
    # Save shared arguments as a single pickle object
    shared_args_fname = f"{tmp_dname}/shared_args.pkl"
    save_pickle(shared_args, shared_args_fname)
    # Split and save arguments, and sumit jobs
    n_args_per_job = -(-len(args) // n_distribute)
    job_ids = []
    for i in range(n_distribute):
        index = str(i + 1).zfill(int(log10(n_distribute) + 1))
        _args_fname = f"{tmp_dname}/args.pkl.{index}"
        save_pickle(args[i * n_args_per_job:(i + 1) * n_args_per_job],
                    _args_fname)
        _py_fname = f"{tmp_dname}/scatter.py.{index}"
        with open(_py_fname, 'w') as f:
            f.write(f"""\
import logging
import logzero
from bits.util import load_pickle, save_pickle
from {func.__module__} import {func.__name__}
logzero.loglevel(logging.{"INFO" if log_level == "info" else "DEBUG"})
args = load_pickle("{_args_fname}")
shared_args = load_pickle("{shared_args_fname}")
save_pickle({func.__name__}(args, n_core={n_core}, **shared_args),
            "{tmp_dname}/{out_fname}.{index}")
""")
        job_ids.append(scheduler.submit(f"python {_py_fname}",
                                        f"{tmp_dname}/scatter.sh.{index}",
                                        job_name=f"{job_name}_scatter",
                                        log_fname=f"{tmp_dname}/log.{index}",
                                        n_core=n_core,
                                        max_cpu_hour=max_cpu_hour,
                                        max_mem_gb=max_mem_gb))
    # Merge results
    scheduler.submit("sleep 1s",
                     f"{tmp_dname}/gather.sh",
                     job_name=f"{job_name}_gather",
                     log_fname=f"{tmp_dname}/log.gather",
                     depend_job_ids=job_ids,
                     wait=True)
    merged = []
    script = f"find {tmp_dname} -name '{out_fname}.*' | sort"
    for fname in run_command(script).strip().split('\n'):
        merged += load_pickle(fname)
    return merged
