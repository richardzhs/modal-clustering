library(meanShiftR)
library(ggplot2)

dir <- getwd()
dir.sep <- "/"
echo <- T
source(paste(dir, 'function_packages', "kernel-density-estimation-functions.R", sep = dir.sep))
source(paste(dir, 'function_packages', "gaussian-mixture-model-functions.R", sep = dir.sep))
source(paste(dir, 'function_packages', "density-optimizing-functions.R", sep = dir.sep))

dat <- read.csv('data/s-originals/s1.csv', header = T, sep = "")
colnames(dat) <- c('f1', 'f2', 'group_id')
M <- as.matrix(dat[,c('f1','f2')])

ggplot(dat, aes(x = f1, y = f2, label = group_id)) +
  geom_point(aes(colour = as.character(group_id)))



sphered.M <- sphere(M)
trial.hs <- c(0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0)
trial.hs <- seq(0.001, 2.0, 0.01)
cv.search.out <- cv.search(sphered.M, trial.par = trial.hs)
h.M <- cv.search.out$opt.smopar


ms <- meanShift(sphered.M, bandwidth = rep(h.M, 2), iterations = 1000000)
max(ms$assignment)

ms$assignment
