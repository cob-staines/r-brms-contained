# R-BRMS-CMDSTAN

A Docker container for VM BRMS runs

## Getting Started

1.  Clone this repository to the VM of interest, then:
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