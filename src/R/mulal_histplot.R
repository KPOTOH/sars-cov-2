setwd("~/COVID19/")
rm(list=ls(all=TRUE))

# Plot 7
diff2ref <- read.csv('data/mutation_number_distribution.csv')
hist(
  subset(diff2ref, difference < 39)$difference, 
  main = '',
  xlab = 'Расстояние до референса',
  ylab = 'Количество'
)

# Plot 8
pch_diff <- read.csv('data/parent_child_diff.csv')
hist(
  pch_diff$parent_child_diff,
  main = '',
  xlab = 'Расстояние меж отцом и сыном',
  ylab = 'Количество'
)


