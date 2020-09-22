## SCEPTRE package and manuscript code

This repository contains the **sceptre** R package and the code required to reproduce the analyses reported in Katsevich et al. 2020. Users can download and use the **sceptre** package independently of the Katsevich et al. 2020 analysis code.

### Downloading and installing **sceptre**

There are two ways to download and install **sceptre**.

**Way 1**. 

Run the following code within R.

```
library(devtools)
install_github(repo="Timothy-Barry/sceptre_paper", subdir="sceptre")
```

**Way 2**.

Clone the github repository into a suitable location on your machine.

```
git clone https://github.com/Timothy-Barry/sceptre_paper.git
```

Run the bash file *bash_scripts/build_and_install_package.bash*.

```
bash bash_scripts/build_and_install_package.bash
```
Be sure to run this script from within the cloned *sceptre_paper* directory.

### Learning to use **sceptre**

Have a look at the[Getting up-and-running with sceptre](https://github.com/Timothy-Barry/sceptre_paper/sceptre/vignettes/sceptre-basics.html) vignette.
