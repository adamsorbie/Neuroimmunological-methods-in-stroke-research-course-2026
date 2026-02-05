### LIBRARIES

### check if Rtools is installed if necessary


if (Sys.info()['sysname'] == "Windows") {
  if (pkgbuild::has_rtools() == FALSE) {
    print(
      "Rtools is not installed and is required on Windows systems to correctly install packages, please install correct version for your version of R, from: https://cran.r-project.org/bin/windows/Rtools/"
    )
  }
}

if (!require("pacman")) {
  install.packages("pacman", repos = "http://cran.us.r-project.org")
}

packages <- c(
  "BiocManager",
  "tidyverse",
  "ggpubr",
  "ggsci",
  "rstatix",
  "phyloseq",
  "vegan",
  "cowplot",
  "Maaslin2",
  "zoo"
)

pacman::p_load(char = packages)
uninstalled <- packages[!(packages %in% installed.packages()[, "Package"])]

if (length(uninstalled) > 0) {
  print(paste0(
    "installation not completed successfully, please install: ",
    uninstalled
  ))
} else {
  print("Installed packages successfully :)")
}
