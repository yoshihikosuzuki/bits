from dataclasses import dataclass
from logzero import logger
from .utils import run_command


@dataclass(repr=False, eq=False)
class Scheduler:
    scheduler_name: str
    submit_command: str
    queue_name: str = None
    out_log: str = "log.stdout"
    err_log: str = "log.stderr"
    merge_log: bool = True

    def __post_init__(self):
        assert self.scheduler_name in ("sge", "slurm"), "Not supported scheduler"

    def submit(self,
               script,
               out_fname,
               job_name="run_script",
               n_core=1,
               time_limit=None,
               mem_limit=None,
               depend=[],
               wait=False):
        with open(out_fname, 'w') as f:
            f.write((self.sge_nize if self.scheduler_name == "sge"
                     else self.slurm_nize)(script,
                                           job_name,
                                           n_core,
                                           time_limit,
                                           mem_limit,
                                           depend,
                                           wait))
        ret = run_command(f"{self.submit_command} {out_fname}")
        return ret.split()[2] if self.scheduler_name == "sge" else ret.split()[-1]

    def sge_nize(self,
                 script,
                 job_name,
                 n_core,
                 time_limit,
                 mem_limit,
                 depend,
                 wait):
        header = '\n'.join([f"#!/bin/bash",
                            f"#$ -N {job_name}",
                            f"#$ -o {self.out_log}",
                            f"#$ -S /bin/bash",
                            f"#$ -cwd",
                            f"#$ -V",
                            f"#$ -pe smp {n_core}",
                            f"#$ -sync {'y' if wait else 'n'}"])
        if self.merge_log:
            header += f"\n#$ -j y"
        else:
            header += f"\n#$ -e {self.err_log}",
        if self.queue_name is not None:
            header += f"\n#$ -q {self.queue_name}"
        if time_limit is not None:
            header += f"\n#$ -l h_cpu={time_limit}"
        if mem_limit is not None:
            header += f"\n#$ -l mem_total={mem_limit}"
        if len(depend) != 0:
            header += f"\n#$ -hold_jid {','.join(depend)}"
        return f"{header}\n\n{script}\n"

    def slurm_nize(self,
                 script,
                 job_name,
                 n_core,
                 time_limit,
                 mem_limit,
                 depend,
                 wait):
        header = '\n'.join([f"#!/bin/bash",
                            f"#SBATCH -J {job_name}",
                            f"#SBATCH -o {self.out_log}",
                            f"#SBATCH -n 1",
                            f"#SBATCH -N 1",
                            f"#SBATCH -c {n_core}",
                            f"#SBATCH -t {'24:00:00' if time_limit is None else time_limit}",
                            f"#SBATCH --mem={50000 if mem_limit is None else mem_limit}"])
        if not self.merge_log:
            header += f"\n#$ -e {self.err_log}",
        if self.queue_name is not None:
            header += f"\n#$ -q {self.queue_name}"
        if len(depend) != 0:
            header += f"\n#SBATCH -d afterany:{','.join(depend)}"
        if wait:
            header += f"\n#SBATCH --wait"
        return f"{header}\n\n{script}\n"


def load_args():
    import argparse
    p = argparse.ArgumentParser(description="Utility for job submission with scheduler.")

    p.add_argument("script_fname",
                   type=str,
                   help="Script file name to be submitted.")

    p.add_argument("scheduler_name",
                   type=str,
                   help="Job scheduler's name.")

    p.add_argument("submit_command",
                   type=str,
                   help="Command to submit a job.")

    p.add_argument("-q",
                   "--queue_name",
                   type=str,
                   default=None,
                   help="Quene name to which jobs will be submitted. [None]")

    p.add_argument("-n",
                   "--job_name",
                   type=str,
                   default="run_script",
                   help="Job name. [run_script]")

    p.add_argument("-o",
                   "--out_log",
                   type=str,
                   default="log.stdout",
                   help="Log file name for standard output. [log.stdout]")

    p.add_argument("-e",
                   "--err_log",
                   type=str,
                   default="log.stderr",
                   help="Log file name for standard error. [log.stderr]")

    p.add_argument("-j",
                   "--merge_log",
                   type=bool,
                   default=True,
                   help="Output stderr in stdout log file. [True]")

    p.add_argument("-p",
                   "--n_core",
                   type=int,
                   default=1,
                   help="Number of cores to be used. [1]")

    p.add_argument("-t",
                   "--time_limit",
                   type=str,
                   default=None,
                   help="Maximum time to be used. [None]")

    p.add_argument("-m",
                   "--mem_limit",
                   type=int,
                   default=None,
                   help="Maximum memory to be used. [None]")

    p.add_argument("-w",
                   "--wait",
                   action="store_true",
                   default=False,
                   help="Wait until the submitted job finishes. [False]")

    return p.parse_args()


if __name__ == "__main__":
    args = load_args()
    s = Scheduler(args.scheduler_name,
                  args.submit_command,
                  args.queue_name,
                  args.out_log,
                  args.err_log,
                  args.merge_log)
    s.submit(f"bash {args.script_fname}",
             f"{args.script_fname}.{args.scheduler_name}",
             args.job_name,
             args.n_core,
             args.time_limit,
             args.mem_limit,
             [],   # NOTE: dependency is not supported
             args.wait)
