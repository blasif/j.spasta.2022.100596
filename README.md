# j.spasta.2022.100596

R-code of the manuscript "A selective view of climatological data and likelihood estimation",
https://doi.org/10.1016/j.spasta.2022.100596.

The R-code was developed and ran using the following software:

**R-version:** 4.1.1  
**Linux distro:** Ubuntu 20.04 LTS

steps to build and run the docker container:

1.  `git clone https://git.math.uzh.ch/reinhard.furrer/j.spasta.2022.100596.git`

2.  `docker build . -t climsimstudy`

3.  `docker run --rm -p 8787:8787 -v /pathToTheRepoFolder:/home/rstudio climsimstudy`

4.  open the web browser and go to localhost:8787

You may need root access for steps 2 and 3.