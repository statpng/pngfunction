
### root password
``` sudo passwd root ```


### install ssh
```
sudo apt update
sudo apt install openssh-server

sudo service ssh start

sudo ufw allow ssh
```

### install R & Rstudio-server
https://posit.co/download/rstudio-server/

```
sudo apt-get install r-base
sudo apt-get install gdebi-core
wget https://download2.rstudio.org/server/focal/amd64/rstudio-server-2023.06.1-524-amd64.deb
sudo gdebi rstudio-server-2023.06.1-524-amd64.deb
```

### install others for installing some R packages
```
sudo apt-get install -y htop
sudo apt-get install -y net-tools
sudo apt-get install -y git
sudo apt-get install -y build-essential libcurl4-gnutls-dev libxml2-dev libssl-dev libfontconfig1-dev
sudo apt-get install -y libzmq3-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev build-essential libcurl4-openssl-dev libxml2-dev libssl-dev libfontconfig1-dev
```



### R packages

```
install.packages(c("tidyverse", "purrr", "parallel", "dplyr", "ggplot2"))
install.packages(c("dirmult", "Ternary", "reshape2", "Rtsne", "mnormt", "quadprog", "compositions", "glmnet"))
```


