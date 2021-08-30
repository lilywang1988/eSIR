# General Set-up
## Immunization date and proportion
change_time <- c("04/29/2020")
alpha0 <- c(0.2) # 20% of the susceptible population were found immunized

# Example 1 for different input lengths
NI_complete1 <- c( 41, 41, 41, 45, 62, 131, 200, 270, 375)
RI_complete1 <- c(1, 1, 7, 10)
N1 <- 58.5e6

test_that("Check the input data", {
  expect_error(
    eSAIR(
      RI_complete1 / N1,
      NI_complete1 / N1,
      begin_str = "04/10/2020",
      T_fin = 40,
      alpha0 = alpha0,
      change_time = change_time,
      dic = FALSE,
      casename = "New_York_antibody",
      save_files = FALSE,
      save_mcmc = FALSE,
      save_plot_data = FALSE,
      M = 1E2,
      nburnin = 5E1
    )
  )
})

# Example 2 for correct output
NI_complete2 <- c(41, 45, 62, 131)
RI_complete2 <- c(1, 1, 7, 10)
N2 <- 1E3
test_that("Check the output", {
  expect_output(
    str(eSAIR(
      RI_complete2 / N2,
      (NI_complete2 - RI_complete2) / N2,
      begin_str = "04/10/2020",
      T_fin = 40,
      alpha0 = alpha0,
      change_time = change_time,
      dic = FALSE,
      casename = "New_York_antibody",
      save_files = FALSE,
      save_mcmc = FALSE,
      save_plot_data = FALSE,
      M = 1E2,
      nburnin = 1
    )),
    "List"
  )
})


