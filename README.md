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