gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

lognorm.parm <- function(mu0,var0){
  var <- log(var0 / mu0^2 + 1)
  mu <- log(mu0) - var / 2
  list(mu = mu, var = var)
}

