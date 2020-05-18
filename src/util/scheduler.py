from dataclasses import dataclass
from typing import Optional, Union, List
from logzero import logger
from .proc import run_command


@dataclass(eq=False)
class Scheduler:
    """Utility for submitting jobs using a job scheduler.

    usage:
      > s = Scheduler("sge", "qsub", "all.q")
      > s.submit("sleep 1s", "sleep.sh")

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
               out_fname: str,
               job_name: str = "run_script",
               log_fname: str = "log",
               n_core: int = 1,
               time_limit: Optional[Union[int, str]] = None,
               mem_limit: Optional[int] = None,
               depend_job_ids: List[str] = [],
               wait: bool = False) -> str:
        """Submit `script` after writing it into a file `out_fname`.

        positional arguments:
          @ script         : Script (not a file).
          @ out_fname      : Shell script file name to which `script` is written.
          @ job_name       : Display name of the job.
          @ log_fname      : Log file name.
          @ n_core         : Number of cores used for the job.
          @ time_limit     : In hours. Format depends on the scheduler type.
          @ mem_limit      : In MB.
          @ depend_job_ids : Wait to start the script until these jobs finish.
          @ wait           : If True, wait until the script finishes.

        return value:
          @ job_id : Of the script submitted.
        """
        logger.info(f"Submitting a job: {script}")
        if self.prefix_command is not None:
            script = f"{self.prefix_command}\n{script}"
        with open(out_fname, "w") as f:
            f.write((self._sge_nize if self.scheduler_name == "sge"
                     else self._slurm_nize)(script,
                                            job_name,
                                            log_fname,
                                            n_core,
                                            time_limit,
                                            mem_limit,
                                            depend_job_ids,
                                            wait))
        ret = run_command(f"{self.submit_command} {out_fname}")
        return ret.split()[2] if self.scheduler_name == "sge" else ret.split()[-1]

    def _sge_nize(self,
                  script: str,
                  job_name: str,
                  log_fname: str,
                  n_core: int,
                  time_limit: Optional[Union[int, str]],
                  mem_limit: Optional[int],
                  depend_job_ids: List[str],
                  wait: bool) -> str:
        header = '\n'.join(["#!/bin/bash",
                            f"#$ -N {job_name}",
                            f"#$ -o {log_fname}",
                            "#$ -j y",
                            "#$ -S /bin/bash",
                            "#$ -cwd",
                            "#$ -V",
                            f"#$ -pe smp {n_core}",
                            f"#$ -sync {'y' if wait else 'n'}"])
        if self.queue_name is not None:
            header += f"\n#$ -q {self.queue_name}"
        if time_limit is not None:
            header += f"\n#$ -l h_cpu={time_limit}"
        if mem_limit is not None:
            header += f"\n#$ -l mem_total={mem_limit}"
        if len(depend_job_ids) != 0:
            header += f"\n#$ -hold_jid {','.join(depend_job_ids)}"
        return f"{header}\n\n{script}\n"

    def _slurm_nize(self,
                    script: str,
                    job_name: str,
                    log_fname: str,
                    n_core: int,
                    time_limit: Optional[Union[int, str]],
                    mem_limit: Optional[int],
                    depend_job_ids: List[str],
                    wait: bool) -> str:
        header = '\n'.join(["#!/bin/bash",
                            f"#SBATCH -J {job_name}",
                            f"#SBATCH -o {log_fname}",
                            "#SBATCH -n 1",
                            "#SBATCH -N 1",
                            f"#SBATCH -c {n_core}",
                            f"#SBATCH -t {'24:00:00' if time_limit is None else time_limit}",
                            f"#SBATCH --mem={50000 if mem_limit is None else mem_limit}"])
        if self.queue_name is not None:
            header += f"\n#$ -q {self.queue_name}"
        if len(depend_job_ids) != 0:
            header += f"\n#SBATCH -d afterany:{','.join(depend_job_ids)}"
        if wait:
            header += "\n#SBATCH --wait"
        return f"{header}\n\n{script}\n"
