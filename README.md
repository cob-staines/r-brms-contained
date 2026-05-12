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