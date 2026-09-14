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

## Running on an HPC cluster (Slurm + Singularity)

Shared HPC clusters generally don't run Docker directly (it needs root) and
schedule jobs through Slurm rather than a long-lived VM. Resource requests
(CPUs, memory, walltime) are declared as `#SBATCH` directives inside the
submission script itself, not through a web form — so "pick a script to run"
in an Open OnDemand portal is accurate: that script is the Slurm batch
script. The `hpc/*.sbatch` scripts here are written for UCSB's GRIT HPC
cluster (partition `grit_nodes`, no `--account` needed, `singularity`
available); adjust for other clusters as needed.

**1. Build and publish the container image.** A GitHub Actions workflow
(`.github/workflows/docker-publish.yml`) builds the Docker image and pushes
it to `ghcr.io/<your-github-username>/r-brms-contained` automatically on
every push to `main` that touches the `Dockerfile`. No Docker install is
needed on the cluster. The first time this runs, go to the package's
settings on GitHub and set its visibility to Public (or `singularity remote
login` on the cluster with a GitHub personal access token) so it can be
pulled without extra auth.

**2. On the cluster**, connect with `ssh <username>@hpc.grit.ucsb.edu`
(campus network/VPN required; off-campus, hop through the bastion first:
`ssh -J <username>@ssh.grit.ucsb.edu <username>@hpc.grit.ucsb.edu`). This
drops you into a Slurm-backed interactive session automatically. The
`hpc/*.sbatch` scripts pull the published image into `r-brms.sif` on first
run if it's not already present:
```bash
singularity pull r-brms.sif docker://ghcr.io/<your-github-username>/r-brms-contained:latest
```

**3. Transfer data to/from the cluster.** `data/` and `*.rds` are gitignored
(too large/data-specific to version), so getting the input `.RData` onto
the cluster and fitted models back off it is on you — the container just
reads/writes plain files under `/project` (the bind-mounted repo directory),
there's no transfer logic in the job itself. Use `rsync` against a direct
HPC hostname (not the auto-forwarding `hpc.grit.ucsb.edu`, which drops you
straight into a Slurm session rather than a plain shell):
```bash
# push input data up, from your local machine
rsync -avr data/bd_model_env_unified_2026-09-11.RData <username>@ssh.grit.ucsb.edu:~/r-brms-contained/data/


# pull results back down, after a run finishes
rsync -avr <username>@ssh.grit.ucsb.edu:~/r-brms-contained/models/ ./models/
```
Clone/keep the repo in your regular GRIT home directory rather than the
`/home/hpc-scratch` BeeGFS scratch space — scratch is faster but **not
backed up and auto-purges files untouched for 3 months**, which isn't worth
the risk given the data/output sizes here (a few MB to a few 10s of MB).

**4. Test first, then run for real:**
```bash
sbatch hpc/test_cars.sbatch              # quick smoke test (~minutes)
squeue -u $USER                          # watch it queue/run
cat logs/r-brms-test_<jobid>.out         # check the result
```

**5. Estimate runtime before committing to the full model.**
`size_bd_modeling.R` reads `N_CHAINS`/`N_ITER`/`N_WARMUP`/`N_CORES`/
`N_THREADS_PER_CHAIN` from env vars (defaulting to the production values: 4
chains, 4000 iter, 2000 warmup, 1 thread/chain). `hpc/timing_probe.sbatch`
runs the exact same data/formula/family but with 1 chain and 100 warmup +
100 sampling iterations (same 50/50 warmup:sampling ratio as the real run),
and the script itself prints fit time at the end:
```bash
sbatch hpc/timing_probe.sbatch
cat logs/size_bd_timing_probe_<jobid>.out    # look for "Fit time (...)"
```
Extrapolate: `estimated full runtime ≈ probe time × (4000 / 200)`, then pad
it (e.g. x1.5) before setting `-t` in `run_size_bd_modeling.sbatch` — a
short probe's warmup may not fully reach the step size/tree depth the real
2000-iteration warmup settles into, so per-iteration cost can be higher in
the full run than the probe suggests.

**6. Run the real model:**
```bash
sbatch hpc/run_size_bd_modeling.sbatch
```
After it completes, check `seff <jobid>` (or `sacct -j <jobid>
--format=MaxRSS,Elapsed`) to see actual memory/time used, and tighten the
`--mem`/`-t` requests for next time — accurate requests generally queue
faster than padded ones.

## Model structure vs. run parameters

`size_bd_modeling.R` is organized into clearly separated sections: **Model
structure** (family, brms `formula`, priors — what the model *is*), **Run
parameters** (MCMC settings, threading — how it's *run*, overridable via env
vars), and **Execution** (load data, fit, save). Edit the structure section
to change the model; edit the run-parameters section or the env vars passed
into the script (locally, as a timing probe, or full-scale on the cluster)
to change how it runs.

## Multithreading

`brm()` supports within-chain multithreading via
[reduce_sum](https://cran.r-project.org/web/packages/brms/vignettes/brms_threading.html),
splitting each chain's likelihood evaluation across threads — useful when
you have spare CPUs and want to speed up each chain rather than just run
more chains in parallel. Set `N_THREADS_PER_CHAIN` (default 1, i.e. no
threading) alongside `N_CORES`; **total CPUs used is `N_CORES *
N_THREADS_PER_CHAIN`**, so the Slurm `-c`/`--cpus-per-task` request must
equal that product — the `hpc/*.sbatch` scripts export both env vars right
next to the `-c` directive as a reminder to keep them in sync.