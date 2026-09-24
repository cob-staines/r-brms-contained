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

**5. Estimating runtime is unreliable here — watch live progress instead.**
`size_bd_modeling.R` reads `N_CHAINS`/`N_ITER`/`N_WARMUP`/`N_CORES`/
`N_THREADS_PER_CHAIN` from env vars (script defaults: 4 chains, 4000 iter,
2000 warmup, 1 thread/chain; `hpc/run_size_bd_modeling.sbatch` currently
overrides this to 2 threads/chain for production runs - keep this note and
that script's `-c`/`N_CORES`/`N_THREADS_PER_CHAIN` in sync if you change it
again), and `hpc/timing_probe.sbatch` runs the exact same data/formula/family
at a smaller scale.

In practice, linear extrapolation from a short probe doesn't hold for this
model: Stan's warmup adaptation isn't linear in warmup length — its default
windowed adaptation schedule barely engages for a very short warmup (e.g.
10-50 iterations), so a tiny probe can finish suspiciously fast without
doing the real adaptation work, while a slightly larger one can suddenly
hit expensive tree-depth blowups once adaptation properly kicks in. A probe
at `warmup=50` took over 30 minutes here while `warmup=10` took under 2 —
not something you can extrapolate `×N` from.

So instead of trusting an extrapolated number, **watch actual progress live**.
Both scripts run `Rscript` under `stdbuf -oL -eL` (line-buffered, so output
flushes to the log immediately instead of sitting in a block buffer until
the process exits) and `brm()` is called with `refresh` set to print
progress regularly. Submit with a generous `-t`, then:
```bash
sbatch hpc/timing_probe.sbatch        # or hpc/run_size_bd_modeling.sbatch
tail -f logs/<job name>_<jobid>.out   # watch iterations tick by in real time
```
If a run does complete, the script also prints a per-chain warmup/sampling
time breakdown (`$fit$time()`) — useful for seeing whether warmup or
sampling is actually the expensive part.

**Memory is the other thing to check before the real run — and unlike
runtime, it IS safe to extrapolate.** Each chain is a separate OS process
(memory isn't shared between them), so memory scales close to linearly with
chains × kept draws, with no equivalent of the warmup-adaptation
nonlinearity above. `hpc/timing_probe.sbatch` runs `N_CHAINS=4` (matching
production, not fewer) for exactly this reason — since chains run in
parallel this costs no extra wall-clock time, only more requested
CPU/memory, and testing at fewer chains would underestimate real memory
need by up to 4x. After a probe run:
```bash
sacct -j <jobid> --format=MaxRSS
```
A 16GB local VM has already crashed at full production settings (chains=4,
iter=4000, warmup=2000), so treat anything under ~32G as known-insufficient
for the real run — use the probe's `MaxRSS` (and, if you want a second data
point, how it changes as you bump `N_WARMUP`/`N_ITER` in the probe) as the
actual basis for `--mem` in `hpc/run_size_bd_modeling.sbatch`, rather than
guessing.

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

The Bd load response variable is also an env var rather than hardcoded, so
several variants can run concurrently as separate jobs against the same
script without editing it:
- `BD_LOAD_VAR` (default `bd_load_co10_density`) — the response column used
  in the formula, e.g. swap in `bd_load_co50_density` for a different Bd
  load density cutoff.
- `DATA_FILE` (default `data/bd_model_env_unified_2026-09-23.RData`) — the
  input `.RData` to fit against.
- `RUN_TAG` — identifies output files (`models/unified_model_results_zln_<RUN_TAG>_<timestamp>.rds`);
  defaults to `<DATA_FILE basename>_<BD_LOAD_VAR>`, so concurrent runs that
  share a data file but use different `BD_LOAD_VAR` values (or vice versa)
  still get distinct output automatically. Override explicitly if you want a
  more readable tag.

See the `BD_LOAD_VAR` example in `hpc/run_size_bd_modeling.sbatch` for
submitting two such variants side by side on the cluster.

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