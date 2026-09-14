# R-BRMS-CMDSTAN

A Docker container for virtual BRMS runs

## Getting Started

On your virtual machine:

1.  Clone this repository, then:
```bash
cd your/repo/path
```

2.  Build the container image:
```bash
docker build -t r-brms .
```

3. Test build:
```bash
docker run -it --rm -v $(pwd):/project r-brms Rscript brms_cars_test.R
```

## Long-running models

1. create a new `tmux` session on your virtual machine:
```{bash}
tmux new -s myjob
```

2. run your model:
```bash
docker run -it --rm -v $(pwd):/project r-brms Rscript long_running_brms_model_script.R
```

3. detatch from session:
`Ctrl+b, d` to detatch. You can now close your ssh connection to the VM, and the model will keep running.

4: reconnect to session to check on model:
Reestablish ssh, then:
```{bash}
tmux attach -t myjob
```

## Running on an HPC cluster (Slurm + Apptainer)

Shared HPC clusters generally don't run Docker directly (it needs root) and
schedule jobs through Slurm rather than a long-lived VM. Resource requests
(CPUs, memory, walltime) are declared as `#SBATCH` directives inside the
submission script itself, not through a web form — so "pick a script to run"
in an Open OnDemand portal is accurate: that script is the Slurm batch
script.

**1. Build and publish the container image.** A GitHub Actions workflow
(`.github/workflows/docker-publish.yml`) builds the Docker image and pushes
it to `ghcr.io/<your-github-username>/r-brms-contained` automatically on
every push to `main` that touches the `Dockerfile`. No Docker install is
needed on the cluster. The first time this runs, go to the package's
settings on GitHub and set its visibility to Public (or `apptainer remote
login` on the cluster with a GitHub personal access token) so it can be
pulled without extra auth.

**2. On the cluster**, convert the published image to Apptainer/Singularity
format (this is what the `hpc/*.sbatch` scripts do automatically, pulling
into `r-brms.sif` if it's not already present):
```bash
apptainer pull r-brms.sif docker://ghcr.io/<your-github-username>/r-brms-contained:latest
```

**3. Fill in the `#SBATCH` placeholders** in `hpc/test_cars.sbatch` and
`hpc/run_size_bd_modeling.sbatch` (`--partition`, `--account`, and check
`--mem`/`--time`). Find your values with:
```bash
sinfo                                    # available partitions
sacctmgr show associations user=$USER    # your account/allocation
```

**4. Test first, then run for real:**
```bash
sbatch hpc/test_cars.sbatch              # quick smoke test (~minutes)
squeue -u $USER                          # watch it queue/run
cat logs/r-brms-test_<jobid>.out         # check the result

# once that succeeds:
sbatch hpc/run_size_bd_modeling.sbatch
```
After the first real run, check `seff <jobid>` (or `sacct -j <jobid>
--format=MaxRSS,Elapsed`) to see actual memory/time used, and tighten the
`--mem`/`--time` requests in the script — accurate requests generally queue
faster than padded ones.