# Experiments for All Models
# 2017.8.5

setwd("~/Desktop/final project/codes")

source("shared_functions.R")
source("models/markowitz/t_dependent.R")

# ===== Markowitz with Student t Distr. Model =====
all.gamma <- c(0.00001, 0.00005, 0.0001, 0.0005, 0.001, 0.01, 0.1, 1, 5, 10)
all.data <- GetAllData()
# Daily Freq
RunAllGammas(all.data, all.gamma, 100, 5)

# Weekly Freq
RunAllGammas(GetWeeklyData(all.data), all.gamma, 100, 5)

# Monthly Freq
RunAllGammas(GetMonthlyData(all.data), all.gamma, 100, 5)

