# General Set-up
## Step function of pi(t)
change_time <- c("01/23/2020", "02/04/2020", "02/08/2020")
pi0<- c(1.0, 0.9, 0.5, 0.1)

# Example 1 for different input lengths
NI_complete1 <- c( 41, 41, 41, 45, 62, 131, 200, 270, 375)
RI_complete1 <- c(1, 1, 7, 10)
N1 <- 58.5e6

test_that("Check the input data", {
  expect_error(
    tvt.eSIR(
      RI_complete1 / N1,
      NI_complete1 / N1,
      begin_str = "01/13/2020",
      T_fin = 40,
      pi0 = pi0,
      change_time = change_time,
      dic = FALSE,
      casename = "Hubei_step",
      save_files = FALSE,
      save_mcmc = FALSE,
      save_plot_data = FALSE,
      M = 2E2,
      nburnin = 1E2
    )
  )
})

# Example 2 for correct output
NI_complete2 <- c(41, 45, 62, 131)
RI_complete2 <- c(1, 1, 7, 10)
N2 <- 1E3
test_that("Check the output", {
  expect_output(
    str(tvt.eSIR(
      RI_complete2 / N2,
      (NI_complete2 - RI_complete2) / N2,
      begin_str = "01/13/2020",
      T_fin = 40,
      pi0 = pi0,
      change_time = change_time,
      dic = FALSE,
      casename = "Hubei_step",
      save_files = FALSE,
      save_mcmc = FALSE,
      save_plot_data = FALSE,
      M = 1E2,
      nburnin = 5E1
    )),
    "List"
  )
})




