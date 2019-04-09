owlS01ST4 <- read.csv("rawOWLdata-S01ST4.csv", stringsAsFactors = FALSE)
owlS02ST1 <- read.csv("rawOWLdata-S02ST1.csv", stringsAsFactors = FALSE)

hist(owlS01ST4$TT.Charge.Ult..s.ult., breaks = 30)
hist(owlS02ST1$TT.Charge.Ult..s.ult., breaks = 30)

