cov_names=names(example.dat)[1:2]
library(survival)
z=tc(type="PTC", dataset=example.dat, cov_names = cov_names, min_exp_events = 10, min_future_events = 50)
test_that("PTC", {
          expect_equal_to_reference(object=z[[1]], tolerance=.00001, scale=1, file="test_PTC.rds")
          })
