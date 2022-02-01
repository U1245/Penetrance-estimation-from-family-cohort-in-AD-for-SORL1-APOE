# simulated dataset to test the function
mydata <- read.table("data_example.txt", header = TRUE, sep = "\t")

source("run_model.R") #  the run.SORL1_APOE_model needs bped3alleles2alleles function to work.

# example with cut off for beta at times 60, 65 and 70, cnesoring time at 85 years old
run.SORL1_APOE_model(mydata, "time", "status", c(0, 60, 65, 70), 85)

