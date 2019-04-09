owlS01ST4 <- read.csv("rawOWLdata-S01ST4.csv", stringsAsFactors = FALSE)
owlS02ST1 <- read.csv("rawOWLdata-S02ST1.csv", stringsAsFactors = FALSE)

# Time To Charge Ult (s/ult)
hist(owlS01ST4$TTChargeUlt, breaks = 20)
abline(v = mean(owlS01ST4$TTChargeUlt), col="blue")
hist(owlS02ST1$TTChargeUlt, breaks = 20)
abline(v = mean(owlS02ST1$TTChargeUlt), col="blue")

# Number of Ults Per 10 Minutes
hist(owlS01ST4$UltsPer10, breaks = 20)
abline(v = mean(owlS01ST4$UltsPer10), col="blue")
hist(owlS02ST1$UltsPer10, breaks = 20)
abline(v = mean(owlS02ST1$UltsPer10), col="blue")

# Ult Hold Time (avg. seconds between ult gaining and ult use)
hist(owlS01ST4$UltHoldTime, breaks = 20)
abline(v = mean(owlS01ST4$UltHoldTime), col="blue")
hist(owlS02ST1$UltHoldTime, breaks = 20)
abline(v = mean(owlS02ST1$UltHoldTime), col="blue")
