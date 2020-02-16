# This is the source code for R package eSIR: extended-SIR
# Built on Feb 13, 2020, and last edited on Feb 14, 2020
# Correspondence : Peter X.K. Song, Ph.D. (pxsong@umich.edu)
# Creator: Lili Wang, M.S. (lilywang@umich.edu)
# Model 1: An extended SIR model with a time-varying pi_qbar, a multplier for beta, the transmission rate.
# library(rjags)
# library(gtools) #rdirichlet(n, alpha)
# library(scales) #alpha　function
# library(ggplot2)
# library(chron)
#' SIR model with a time-varying transmission rate
#'
#' SIR model with fixed and known changes in the transmission rate, either stepwise or continuous.
#'
#' In this function we introduce a time-dependent multiplier (boundared between 0 and 1) of the transmisison rate, \eqn{pi_{qbar}(t)}. In this way,  we can edow the transmission some time-dependent changes, either following a step-function or a smooth exponetial function (\eqn{\exp(\lambda_0t)}). The parameters of the function and change points, if any, need to be predefined.
#'
#' @param Y the observed infected proportions by time.
#' @param R the observed removed proportions by time, including death and recovered.
#' @param pi_qbar0 the time dependent function \eqn{\pi_bar{q}(t)} between 0 and 1.
#' @param change_time the change time points for step function pi, defalt value is \code{NULL}.
#' @param exponential logical, whether \eqn{\pi_bar{q}(t)} is exponential \eqn{\exp(-\lambda_0t)} or not; the default is \code{FALSE}.
#' @param lambda0 the rate of decline in the exponential function in \eqn{\exp(-\lambda_0t)}.
#' @param begin_str the character of the starting time, the default is "01/23/2020", which is the starting date that the local government blocked Wuhan City.
#' @param T_fin the maximum follow-up date after the beginning \code{begin_str}, the default is 200.
#' @param nchain the number of MCMC chains called in rjags, the default is 4.
#' @param nadapt the number of iterations for adaptation in the MCMC. The default is 1e4, which is what we suggest to use 1e4 instead.
#' @param M number of draws from each chain, without considering thinning. The default is M=5e3 but we suggest using 5e5 instead.
#' @param thn thinning interval for monitors. Thus, the total number of draws would be \code{round(M/thn)*nchain }. The default is 10.
#' @param nburnin the burn-in period. The default is nburnin=2e3 but we suggest using nburin=2e5.
#' @param dic logical, whether compute the DIC or deviance information criterion.
#' @param file_add the string to denote the location to save the output files and tables.
#' @param save_files logical, whether save (\code{TRUE}) the results or not (\code{FALSE}). This will enable saving the summary table, trace plots, and the plot of the posterior mean of the first derivative of the infection proportion \eqn{\theta_t^I}.
#' @param death_in_R numeric value of average ratio between deaths and cumulative removed subjects. The default is 0.4 within Hubei, and 0.02 outsite Hubei according to the reported data by Feb 11, 2020.
#' @param casename string of the job's name. The default is "tvt.eSIR".
#' @param beta0 the hyperparameter of the mean transmission rate, the default is the one estimated from SARS (0.2586) first-month outbreak.
#' @param gamma0 the hyperparameter of the mean recovery rate (including death), the default is estimated from SARS (0.0821) first-month outbreak.
#' @param R0 the hyperparameter of the mean R0 value. The default is \code{beta0/gamma0}, which can be overwritten by discarding the value set in \code{beta0}.
#' @param gamma0_sd the standard deviance for the prior of the recovery/remove rate, the default is 0.1.
#' @param R0_sd the standard deviance for the prior of R0, the default is 1.
#' @return
#' \item{casename}{casename defined before}
#' \item{incidence_mean}{mean incidence}
#' \item{incidence_ci}{2.5\%, 50\%, and 97.5\% quantiles of the incidences}
#' \item{out_table}{summary table with varibles including the posterior mean of the proportions of the 3 states at their last observation date, and their respective credible inctervals (ci) including the median; the mean and ci of the reporduction number (R0), removed/recovery rate (gamma), transmission rate  (beta)}
#' \item{forecast_infection}{plot to forecast the infection with following lines: the vertial blue line denotes the last observation date; vertial purple line denotes the change point indicating a decrease in infection proportion or the date with 0 value of the posterior mean first-derivative infection proportion \eqn{\theta_t^{\prime I}}; the vertial darkgreen line denotes the deacceleration point of the increasing infection proportion or the date with maximum value of the posterior mean first-derivative infection proportion \eqn{\theta_t^{\prime I}}; the darkgray line denotes the posterior mean of the infection proportion; the red line denotes the posterior median of the infection proportion}
#' \item{forecast_removed}{plot to forecast the removed lines described in \code{forecast_infection}. The meaning of the vertical lines were identical, but the horizontal mean and median were corresponding to the posterior mean and median of the removed state. Moreover, we introduce an additional line for the estimated death proportion, which is based on the input \code{death_in_R}}
#' \item{first_stat_mean}{the mean first stationary date, which is the change point that we observe decline in the infection proportion (\eqn{\theta_t^I}), or its stationary point; it is calculated using the average of the 0-points of all the repeats of posterior draws of the first derivative proportion or \eqn{\theta_t{\prime I}}; this value may be slightly different from the one labeled by the "purple" lines in the \code{forecast_infection} and \code{forecast_removed} two plots, as the latter indicate the 0-value point of the first-derivative the posterior mean of \eqn{\theta_t^I}.}
#' \item{first_stat_ci}{following the definition of \code{first_stat_mean}, but is the corresponding credible interval.}
#' \item{second_stat_mean}{the mean second stationary date, which is the change point that we observe the decline in the increasing spead of the infection proportion (\eqn{\theta_t^I}) or its derivative's stationary point;it is calculated using the average of the stationary values of all the repeats of posterior draws of the first-derivative porportion of infection or \eqn{\theta_t^{\prime I}}; this value may be slightly different from the one labeled by the "darkgreen" lines in the \code{forecast_infection} and \code{forecast_removed} two plots, as the latter indicate the stationary point of the first-derivative the posterior mean of \eqn{\theta_t^I}.}
#' \item{dic_val}{the output of \code{dic.sample()} in \code{rjags}, computing deviance information criterion for model comparison.}
#' @examples
#' NI_complete <- c( 41,41,41,45,62,131,200,270,375,444,549, 729,1052,1423,2714,3554,4903,5806,7153,9074,11177,13522,16678,19665,22112,24953,27100,29631,31728,33366)
#' RI_complete <- c(1,1,7,10,14,20,25,31,34,45,55,71,94,121,152,213,252,345,417,561,650,811,1017,1261,1485,1917,2260,2725,3284,3754)
#' N=58.5e6
#' R <- RI_complete/N
#' Y <- NI_complete/N- R #Jan13->Feb 11
#' ### Step function of pi_qbar(t)
#' change_time <- c("01/23/2020","02/04/2020","02/08/2020")
#' pi_qbar0 <- c(1.0,0.9,0.5,0.1)
#' res.step <-tvt.eSIR(Y,R,begin_str="01/13/2020",T_fin=200,pi_qbar0=pi_qbar0,change_time=change_time,casename="Hubei_step",save_files = T)
#' res.step$forecast_infection
#' ### continuous exponential function of pi_qbar(t)
#' res.exp <- pi.eSIR(Y,R,begin_str="01/13/2020",T_fin=200,pi_qbar0=pi_qbar0,change_time=change_time,exponential=TRUE,lambda0=0.01,casename="Hubei_exp")
#' res.exp$forecast_infection
#' ### without pi_qbar(t)
#' res.nopi <- pi.eSIR(Y,R,begin_str="01/13/2020",T_fin=200,casename="Hubei_nopi")
#' res.nopi$forecast_infection
#'
#'
#'
#' @export
tvt.eSIR <- function (Y,R, pi_qbar0=NULL,change_time=NULL,exponential=FALSE,lambda0=NULL,begin_str="01/23/2020",T_fin=200,nchain=4,nadapt=1e4,M=5e3,thn=10,nburnin=2e3,dic=FALSE,file_add=character(0),save_files=FALSE,death_in_R=0.02,casename="pi.eSIR",beta0=0.2586,gamma0=0.0821,R0=beta0/gamma0,gamma0_sd=0.1, R0_sd=1){

  len <- round(M/thn)*nchain #number of MCMC draws in total

  T_prime <-length(Y)
  if(T_prime!=length(R)) stop("Y and R should be matched.")
  begin <- chron(dates.=begin_str)
  chron_ls <- chron(begin:(begin+T_fin))
  end <- chron(begin:(begin+T_fin))[T_fin]
  message(paste0("The follow-up is from ",begin," to ",end," and the last observed date is ", chron_ls [T_prime],".") )# current data up to this date
  gamma_var <- gamma0_sd^2
  lognorm_gamma_parm <- lognorm.parm(gamma0,gamma_var)
  R0_var <- R0_sd^2
  lognorm_R0_parm <- lognorm.parm(R0,R0_var)

  if(exponential==FALSE & !is.null(change_time) & !is.null(pi_qbar0)){

    #pi_qbar0 <- c(1,0.9,0.5,0.1)
    message("Running for step-function pi(t)")
    if(length(change_time)!=length(pi_qbar0)-1){stop("We need the length of vector change_time to be the length of pi_qbar0 minus 1. ")}
    change_time_chorn<-c(begin-1,chron(dates.=change_time),end)
    pi_qbar <-rep(pi_qbar0,diff(change_time_chorn))

  } else if (exponential==TRUE & !is.null(lambda0)){
    message("Running for exponential-function pi(t)")
    pi_qbar <-exp(-lambda0*(0:(T_fin-1)))
  }else {
    message("Running without pi(t)")

    pi_qbar <- rep(1,T_fin)
  }

  ################ MCMC ##########
  model1.string <-　paste0("
             model{
                   for(t in 2:(T_prime+1)){
                   Km[t-1,1] <- -beta*pi_qbar[t-1]*theta[t-1,1]*theta[t-1,2]
                   Km[t-1,9] <- gamma*theta[t-1,2]
                   Km[t-1,5] <- -Km[t-1,1]-Km[t-1,9]

                   Km[t-1,2] <- -beta*pi_qbar[t-1]*(theta[t-1,1]+0.5*Km[t-1,1])*(theta[t-1,2]+0.5*Km[t-1,5])
                   Km[t-1,10] <- gamma*(theta[t-1,2]+0.5*Km[t-1,5])
                   Km[t-1,6] <- -Km[t-1,2]-Km[t-1,10]

                   Km[t-1,3] <- -beta*pi_qbar[t-1]*(theta[t-1,1]+0.5*Km[t-1,2])*(theta[t-1,2]+0.5*Km[t-1,6])
                   Km[t-1,11] <- gamma*(theta[t-1,2]+0.5*Km[t-1,6])
                   Km[t-1,7] <- -Km[t-1,3]-Km[t-1,11]

                   Km[t-1,4] <- -beta*pi_qbar[t-1]*(theta[t-1,1]+Km[t-1,3])*(theta[t-1,2]+Km[t-1,7])
                   Km[t-1,12] <- gamma*(theta[t-1,2]+Km[t-1,7])
                   Km[t-1,8] <- -Km[t-1,4]-Km[t-1,12]

                   alpha[t-1,1] <- theta[t-1,1]+(Km[t-1,1]+2*Km[t-1,2]+2*Km[t-1,3]+Km[t-1,4])/6
                   alpha[t-1,2] <- theta[t-1,2]+(Km[t-1,5]+2*Km[t-1,6]+2*Km[t-1,7]+Km[t-1,8])/6
                   alpha[t-1,3] <- theta[t-1,3]+(Km[t-1,9]+2*Km[t-1,10]+2*Km[t-1,11]+Km[t-1,12])/6

                   theta[t,1:3] ~ ddirch(k*alpha[t-1,1:3])
                   Y[t-1] ~ dbeta(lambdaY*theta[t,2],lambdaY*(1-theta[t,2]))
                   R[t-1] ~ dbeta(lambdaR*theta[t,3],lambdaR*(1-theta[t,3]))
                  }
                  theta[1,1] <-  1- theta[1,2]- theta[1,3]
                  theta[1,2] ~ dbeta(",1,",",1/max(Y[1],1/N),")
                  theta[1,3] ~ dbeta(",1,",",1/max(R[1],1/N),")
                  gamma ~  dlnorm(",lognorm_gamma_parm[1],",",1/lognorm_gamma_parm[2],")
                  R0 ~ dlnorm(",lognorm_R0_parm[1],",",1/lognorm_R0_parm[2],")
                  beta <- R0*gamma
                  k ~  dgamma(2,0.0001)
                  lambdaY ~ dgamma(2,0.0001)
                  lambdaR ~ dgamma(2,0.0001)
               }
            ")

  model.spec <- textConnection(model1.string)

  posterior <- jags.model(model.spec,data=list('Y'=Y,'R'=R,'T_prime'=T_prime,'pi_qbar'=pi_qbar),n.chains =nchain, n.adapt = nadapt)

  update(posterior,nburnin) #burn-in

  jags_sample <-jags.samples(posterior,c('theta','gamma','R0','beta','Y','lambdaY','lambdaR','k'),n.iter=M*nchain,thin=thn)

  if(dic) {
    dic_val <-dic.samples(posterior,n.iter=M*nchain,thin=thn)
    message(paste0("DIC is: ", dic_val))
  } else dic_val=NULL


  if(save_files) {
    png(paste0(file_add,casename,"theta_p.png"), width = 700, height = 900)
    plot(as.mcmc.list(jags_sample$theta)[[1]][,(1:3)*(T_prime+1)]) # posterior true porbabilities
    dev.off()
    png(paste0(file_add,casename,"R0_p.png"), width = 700, height = 350)
    plot(R0_p<-as.mcmc.list(jags_sample$R0)[[1]])
    dev.off()
    png(paste0(file_add,casename,"gamma_p.png"), width = 700, height = 350)
    plot(gamma_p<-as.mcmc.list(jags_sample$gamma)[[1]])
    dev.off()

    png(paste0(file_add,casename,"beta_p.png"), width = 700, height = 350)
    plot(beta_p<-as.mcmc.list(jags_sample$beta)[[1]])
    dev.off()


    png(paste0(file_add,casename,"lambdaY_p.png"), width = 700, height = 350)
    plot(lambdaY_p<-as.mcmc.list(jags_sample$lambdaY)[[1]])
    dev.off()

    png(paste0(file_add,casename,"lambdaR_p.png"), width = 700, height = 350)
    plot(lambdaR_p<-as.mcmc.list(jags_sample$lambdaR)[[1]])
    dev.off()

    png(paste0(file_add,casename,"k_p.png"), width = 700, height = 350)
    plot(k_p<-as.mcmc.list(jags_sample$k)[[1]])
    dev.off()

  }else{

   R0_p<-as.mcmc.list(jags_sample$R0)[[1]]
   gamma_p<-as.mcmc.list(jags_sample$gamma)[[1]]
   beta_p<-as.mcmc.list(jags_sample$beta)[[1]]
   lambdaY_p<-as.mcmc.list(jags_sample$lambdaY)[[1]]
   lambdaR_p<-as.mcmc.list(jags_sample$lambdaR)[[1]]
   k_p<-as.mcmc.list(jags_sample$k)[[1]]
  }


  theta_p<-array(as.mcmc.list(jags_sample$theta)[[1]],dim=c(len,T_prime+1,3))
  theta_p_mean <- apply(theta_p[,T_prime+1,],2,mean)
  theta_p_ci <- as.vector(apply(theta_p[,T_prime+1,],2,quantile,c(0.025,0.5,0.975)))
  R0_p_mean <- mean(R0_p)
  R0_p_ci <- quantile(R0_p,c(0.025,0.5,0.975))

  gamma_p_mean <- mean(gamma_p)
  gamma_p_ci <- quantile(gamma_p,c(0.025,0.5,0.975))

  beta_p_mean <- mean(beta_p)
  beta_p_ci <- quantile(beta_p,c(0.025,0.5,0.975))

  lambdaY_p_mean <- mean(lambdaY_p)
  lambdaY_p_ci <- quantile(lambdaY_p,c(0.025,0.5,0.975))

  lambdaR_p_mean <- mean(lambdaR_p)
  lambdaR_p_ci <- quantile(lambdaR_p,c(0.025,0.5,0.975))

  k_p_mean <- mean(k_p)
  k_p_ci <- quantile(k_p,c(0.025,0.5,0.975))

  #### Forecast ####
  theta_pp <- array(NA,dim=c(len,T_fin-T_prime,3))
  Y_pp <- matrix(NA,nrow=len,ncol=T_fin-T_prime)
  R_pp <- matrix(NA,nrow=len,ncol=T_fin-T_prime)
  for(l in 1:len){
    thetalt1 <- theta_p[l,T_prime+1,1]
    thetalt2 <- theta_p[l,T_prime+1,2]
    thetalt3 <- theta_p[l,T_prime+1,3]
    betal <- c(beta_p)[l]
    gammal <- c(gamma_p)[l]
    kt <- c(k_p)[l]
    lambdaYl <- c(lambdaY_p)[l]
    lambdaRl <- c(lambdaR_p)[l]
    if(betal<0 |gammal<0 |thetalt1<0 |thetalt2<0 |thetalt3<0) next
    for(t in 1:(T_fin-T_prime)){
      Km <- NULL
      alpha_pp <- NULL
      Km[1] <- -betal*pi_qbar[t+T_prime]*thetalt1*thetalt2
      Km[9] <- gammal*thetalt2
      Km[5] <- -Km[1]-Km[9]

      Km[2] <- -betal*pi_qbar[t+T_prime]*(thetalt1+0.5*Km[1])*(thetalt2+0.5*Km[5])
      Km[10] <- gammal*(thetalt2+0.5*Km[5])
      Km[6] <- -Km[2]-Km[10]

      Km[3] <- -betal*pi_qbar[t+T_prime]*(thetalt1+0.5*Km[2])*(thetalt2+0.5*Km[6])
      Km[11] <- gammal*(thetalt2+0.5*Km[6])
      Km[7] <- -Km[3]-Km[11]

      Km[4] <- -betal*pi_qbar[t+T_prime]*(thetalt1+Km[3])*(thetalt2+Km[7])
      Km[12] <- gammal*(thetalt2+Km[7])
      Km[8] <- -Km[4]-Km[12]

      alpha_pp[1] <- thetalt1+(Km[1]+2*Km[2]+2*Km[3]+Km[4])/6
      alpha_pp[2] <- thetalt2+(Km[5]+2*Km[6]+2*Km[7]+Km[8])/6
      alpha_pp[3] <- thetalt3+(Km[9]+2*Km[10]+2*Km[11]+Km[12])/6

      thetalt_tmp <- rdirichlet(1,kt*c(alpha_pp))
      thetalt1<-theta_pp[l,t,1] <- thetalt_tmp[1]
      thetalt2<-theta_pp[l,t,2] <- thetalt_tmp[2]
      thetalt3<-theta_pp[l,t,3] <- thetalt_tmp[3]

      Y_pp[l,t] <- rbeta(1,lambdaYl*thetalt2,lambdaYl*(1-thetalt2))
      R_pp[l,t] <- rbeta(1,lambdaRl*thetalt3,lambdaRl*(1-thetalt3))
    }
  }

  par(mfrow=c(1,1))
  col2 = gg_color_hue(2)
  Y_band <- data.frame(t(apply(Y_pp,2,quantile,probs=c(0.025,0.5,0.975),na.rm=T)))
  thetaI_band <- data.frame(t(apply(theta_p[,-1,2],2,quantile,probs=c(0.025,0.5,0.975),na.rm=T)))
  Y_mean <- c(colMeans(Y_pp,na.rm = T))
  thetaI_mean <- c(colMeans(theta_p[,-1,2],na.rm = T))

  colnames(Y_band)<- c("lower", "median", "upper")
  colnames(thetaI_band)<- c("lower", "median", "upper")
  data_pre <- data.frame(time=1:T_prime,Y)
  data_post <-data.frame(time=1:T_prime,thetaI_band)
  data_fore <- data.frame(time=(T_prime+1):T_fin,Y_band,Y_mean)

  data_comp<-data.frame(time=1:T_fin,rbind(thetaI_band ,Y_band), phase=c(rep('pre',nrow(thetaI_band)),rep('post',nrow(Y_band))),mean=c(thetaI_mean,Y_mean))

  data_poly<-data.frame(y=c(thetaI_band$upper,rev(thetaI_band$lower),Y_band$upper,rev(Y_band$lower)),x=c(1:T_prime,T_prime:1,(T_prime+1):T_fin,T_fin:(T_prime+1)),phase=c(rep('pre',T_prime*2),rep('post',(T_fin-T_prime)*2)),value=c(rep(col2[1],T_prime*2),rep(col2[2],(T_fin-T_prime)*2)))

  ## First-order derivative check
  thetaS_mat <- cbind(theta_p[,-1,1],theta_pp[,,1])
  thetaI_mat <- cbind(theta_p[,-1,2],theta_pp[,,2])
  dthetaI_mat <- (thetaS_mat*thetaI_mat)*((c(beta_p))%o%(pi_qbar))-thetaI_mat*replicate(T_fin,c(gamma_p))
  dthetaI <- colMeans(dthetaI_mat)
  dthetaI_stationary2 <- (1:T_fin)[which.max(dthetaI)]# first second order derivative=0
  dthetaI_stationary1 <- (dthetaI_stationary2:T_fin)[which.min(dthetaI[dthetaI_stationary2:T_fin]>0)] # first order derivative=0
  dthetaI_stationary2_date <- chron_ls[dthetaI_stationary2]
  dthetaI_stationary1_date <- chron_ls[dthetaI_stationary1]


  incidence_vec <-  rowSums((thetaS_mat*thetaI_mat)*((c(beta_p))%o%(pi_qbar)))
  incidence_mean <-  mean(incidence_vec)

  incidence_ci <- quantile(incidence_vec,c(0.025,0.5,0.975))
  second_order_stationary_vec <- (1:T_fin)[apply(dthetaI_mat,1,which.max)]# first second order derivative=0

  first_order_stationary_vec <- sapply(1:len,function(l){
    (second_order_stationary_vec[l]:T_fin)[which.min(dthetaI_mat[l,second_order_stationary_vec[l]:T_fin]>0)]})
 # first order derivative=0
  second_order_stationary_mean <- mean(second_order_stationary_vec)
  first_order_stationary_mean <- mean(first_order_stationary_vec)

  second_order_stationary_ci <- quantile(second_order_stationary_vec, c(0.025,0.5,0.975))
  first_order_stationary_ci <- quantile(first_order_stationary_vec, c(0.025,0.5,0.975))

  second_order_change_date <- chron_ls[second_order_stationary_mean]
  first_order_change_date <- chron_ls[first_order_stationary_mean]
  second_order_change_date_ci <- chron_ls[second_order_stationary_ci]
  first_order_change_date_ci <- chron_ls[first_order_stationary_ci]

  names(second_order_change_date_ci)<-c("2.5%","50%","97.5%")
  names(first_order_change_date_ci)<-c("2.5%","50%","97.5%")

  if(save_files){
    png(paste0(file_add,casename,"deriv.png"), width = 700, height = 350)
    plot(y=dthetaI,x=chron_ls,type='l',ylab="1st order derivative",main="Infection Proportion")
    abline(h=0,col=2)
    abline(v=begin,col="blue")
    if(exponential==FALSE & !is.null(change_time) & !is.null(pi_qbar0)) abline(v=change_time_chorn[-c(1,length(change_time_chorn))],col="gray")
    dev.off()
  }

  y_text_ht <- max(rbind(thetaI_band ,Y_band))/2
  plot <- ggplot(data = data_poly, aes(x = x, y = y)) +
    geom_polygon(alpha = 0.5,aes(fill=value, group=phase)) +
    labs(title=substitute(paste(casename,": infection forecast with prior ",beta[0],"=",v1,",",gamma[0], "=",v2," and ", R[0],"=",v3), list(casename=casename,v1=format(beta0,digits=3),v2=format(gamma0,digits=3),v3=format(R0,digits=3))),subtitle = substitute(paste("Posterior: ", beta[p],"=",v1,",",gamma[p], "=",v2," and ", R[0],"=",v3), list(v1=format(beta_p_mean,digits=3),v2=format(gamma_p_mean,digits=3),v3=format(R0_p_mean,digits=3))),x = "time", y = "P(Removed)")+
    geom_line(data=data_comp,aes(x=time,y=median),color="red")+geom_vline(xintercept = T_prime,color="blue",show.legend = TRUE)+
    geom_vline(xintercept = dthetaI_stationary2,color="darkgreen",show.legend = TRUE)+
    geom_vline(xintercept = dthetaI_stationary1,color="purple",show.legend = TRUE)+
    geom_line(data=data_comp,aes(x=time,y=mean),color="darkgray")+
    geom_point(data=data_pre,aes(x=time,y=Y))+theme_bw()+
    theme( plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust=0.5))+scale_x_continuous(labels= as.character(chron_ls)[seq(1,T_fin,30)],breaks=seq(1,T_fin,30))+scale_fill_discrete(name="Posterior",labels=c(expression(paste(y[t+1:T],' | ',y[1:t],', ',r[1:t])),expression(paste(theta[1:t]^I,' | ',y[1:t],', ',r[1:t]))))+
    annotate(geom="text", label=as.character(chron(chron_ls[T_prime]),format="mon day"), x=T_prime+12, y=y_text_ht,color="blue")+
    annotate(geom="text", label=as.character(chron(dthetaI_stationary2_date,format="mon day")), x=dthetaI_stationary2+12, y=y_text_ht*1.25,color="darkgreen")+
    annotate(geom="text", label=as.character(chron(dthetaI_stationary1_date,format="mon day")), x=dthetaI_stationary1+12, y=y_text_ht*1.5,color="purple")

  plot_list <- list(data_poly=data_poly,data_comp=data_comp,T_prime=T_prime,dthetaI_stationary2=dthetaI_stationary2,dthetaI_stationary1=dthetaI_stationary1,data_pre=data_pre,dthetaI_stationary2_date,dthetaI_stationary1_date,y_text_ht)

  if(save_files) ggsave(paste0(file_add,casename,"_forecast.png"))

  ### Recovery
  R_band <- data.frame(t(apply(R_pp,2,quantile,probs=c(0.025,0.5,0.975),na.rm=T)))
  thetaR_band <- data.frame(t(apply(theta_p[,-1,3],2,quantile,probs=c(0.025,0.5,0.975),na.rm=T)))
  R_mean <- c(colMeans(R_pp,na.rm = T))
  thetaR_mean <- c(colMeans(theta_p[,-1,3],na.rm = T))

  colnames(R_band)<- c("lower", "median", "upper")
  colnames(thetaR_band)<- c("lower", "median", "upper")
  data_pre_R <- data.frame(time=1:T_prime,R) # previous data
  data_post_R <-data.frame(time=1:T_prime,thetaR_band) # posterior of theta^R
  data_fore_R <- data.frame(time=(T_prime+1):T_fin,R_band,R_mean) # The forecast of R after T_prime

  data_comp_R<-data.frame(time=1:T_fin,rbind(thetaR_band ,R_band), phase=c(rep('pre',nrow(thetaR_band)),rep('post',nrow(R_band))),mean=c(thetaR_mean,R_mean),dead=c(thetaR_mean,R_mean)*death_in_R) # the filled area--polygon

  data_poly_R<-data.frame(y=c(thetaR_band$upper,rev(thetaR_band$lower),R_band$upper,rev(R_band$lower)),x=c(1:T_prime,T_prime:1,(T_prime+1):T_fin,T_fin:(T_prime+1)),phase=c(rep('pre',T_prime*2),rep('post',(T_fin-T_prime)*2)),value=c(rep(col2[1],T_prime*2),rep(col2[2],(T_fin-T_prime)*2)))

  r_text_ht <- max(rbind(thetaR_band ,R_band))/2
  plot2 <- ggplot(data = data_poly_R, aes(x = x, y = y)) +geom_polygon(alpha = 0.5,aes(fill=value, group=phase)) +labs(title=substitute(paste(casename,": removed forecast with prior ",beta[0],"=",v1,",",gamma[0], "=",v2," and ", R[0],"=",v3), list(casename=casename,v1=format(beta0,digits=3),v2=format(gamma0,digits=3),v3=format(R0,digits=3))),subtitle = substitute(paste("posterior: ", beta[p],"=",v1,",",gamma[p], "=",v2," and ", R[0],"=",v3), list(v1=format(beta_p_mean,digits=3),v2=format(gamma_p_mean,digits=3),v3=format(R0_p_mean,digits=3))),x = "time", y = "P(Infected)")+
    geom_line(data=data_comp_R,aes(x=time,y=median),color="red",linetype=1)+geom_line(data=data_comp_R,aes(x=time,y=dead),color="black",linetype=1)+
    geom_vline(xintercept = T_prime,color="blue")+
    geom_vline(xintercept = dthetaI_stationary2,color="darkgreen")+
    geom_vline(xintercept = dthetaI_stationary1,color="purple")+
    geom_line(data=data_comp_R,aes(x=time,y=mean),color="darkgray")+
    geom_point(data=data_pre_R,aes(x=time,y=R))+theme_bw()+
    theme( plot.title = element_text(hjust = 0.5),plot.subtitle = element_text(hjust=0.5))+
    scale_x_continuous(labels= as.character(chron_ls)[seq(1,T_fin,30)],breaks=seq(1,T_fin,30))+
    scale_fill_discrete(name="Posterior",labels=c(expression(paste(r[t+1:T],' | ',y[1:t],', ',r[1:t])),expression(paste(theta[1:t]^R,' | ',y[1:t],', ',r[1:t]))))+
    annotate(geom="text", label=as.character(chron(chron_ls[T_prime]),format="mon day"), x=T_prime+12, y=r_text_ht,color="blue")+annotate(geom="text", label=as.character(chron(dthetaI_stationary2_date,format="mon day")), x=dthetaI_stationary2+12, y=r_text_ht*1.25,color="darkgreen")+
    annotate(geom="text", label=as.character(chron(dthetaI_stationary1_date,format="mon day")), x=dthetaI_stationary1+12, y=r_text_ht*1.5,color="purple")

plot2_list <- list(data_poly_R=data_poly_R,data_comp_R=data_comp_R,T_prime=T_prime,dthetaI_stationary2=dthetaI_stationary2,dthetaI_stationary1=dthetaI_stationary1,data_pre_R=data_pre_R,dthetaI_stationary2_date,dthetaI_stationary1_date)

  if(save_files) ggsave(paste0(file_add,casename,"_forecast2.png"))

out_table<-matrix(c(theta_p_mean,theta_p_ci,R0_p_mean,R0_p_ci,gamma_p_mean,gamma_p_ci,beta_p_mean,beta_p_ci),nrow=1)
#out_table<-matrix(c(theta_p_mean,theta_p_ci,R0_p_mean,R0_p_ci,gamma_p_mean,gamma_p_ci,beta_p_mean,beta_p_ci,k_p_mean,k_p_ci,lambdaY_p_mean,lambdaY_p_ci,lambdaR_p_mean,lambdaR_p_ci,as.character(first_order_change_date),as.character(second_order_change_date)),nrow=1)


  colnames(out_table)<-c("thetaS_p_mean","thetaI_p_mean","thetaR_p_mean","thetaS_p_ci_low","thetaS_p_ci_med","thetaS_p_ci_up","thetaI_p_ci_low","thetaI_p_ci_med","thetaI_p_ci_up","thetaR_p_ci_low","thetaR_p_ci_med","thetaR_p_ci_up","R0_p_mean","R0_p_ci_low","R0_p_ci_med","R0_p_ci_up","gamma_p_mean","gamma_p_ci_low","gamma_p_ci_med","gamma_p_ci_up","beta_p_mean","beta_p_ci_low","beta_p_ci_med","beta_p_ci_up")

  #colnames(out_table)<-c("thetaS_p_mean","thetaI_p_mean","thetaR_p_mean","thetaS_p_ci_low","thetaS_p_ci_med","thetaS_p_ci_up","thetaI_p_ci_low","thetaI_p_ci_med","thetaI_p_ci_up","thetaR_p_ci_low","thetaR_p_ci_med","thetaR_p_ci_up","R0_p_mean","R0_p_ci_low","R0_p_ci_med","R0_p_ci_up","gamma_p_mean","gamma_p_ci_low","gamma_p_ci_med","gamma_p_ci_up","beta_p_mean","beta_p_ci_low","beta_p_ci_med","beta_p_ci_up","k_p_mean","k_p_ci_low","k_p_ci_med","k_p_ci_up","lambdaY_p_mean","lambdaY_p_ci_low","lambdaY_p_ci_med","lambdaY_p_ci_up","lambdaR_p_mean","lambdaR_p_ci_low","lambdaR_p_ci_med","lambdaR_p_ci_up","first_order_change_date","second_order_change_date")

  if(save_files) write.csv(out_table,file=paste0(file_add,casename,"_summary.csv"))

  return(list(casename=casename,incidence_mean=incidence_mean,incidence_ci=incidence_ci,out_table=out_table,forecast_infection=plot,forecast_removed=plot2,first_stat_mean=as.character(first_order_change_date),first_stat_ci=as.character(first_order_change_date_ci),second_stat_mean=as.character(second_order_change_date),second_stat_ci=as.character(second_order_change_date_ci),dic_val=dic_val,plot_list=plot_list,plot2_list=plot2_list,theta_p=theta_p,theta_pp=theta_pp,Y_pp=Y_pp,R_pp=R_pp))
}

if(F){
  NI_complete <- c( 41,41,41,45,62,131,200,270,375,444,549, 729,1052,1423,2714,3554,4903,5806,7153,9074,11177,13522,16678,19665,22112,24953,27100,29631,31728,33366)
  RI_complete <- c(1,1,7,10,14,20,25,31,34,45,55,71,94,121,152,213,252,345,417,561,650,811,1017,1261,1485,1917,2260,2725,3284,3754)
  N=58.5e6
  R <- RI_complete/N
  Y <- NI_complete/N- R #Jan13->Feb 11
  ### Step function of pi_qbar(t)
  change_time <- c("01/23/2020","02/04/2020","02/08/2020")
  pi_qbar0 <- c(1.0,0.9,0.5,0.1)
  res.step <-pi.eSIR(Y,R,begin_str="01/13/2020",T_fin=200,pi_qbar0=pi_qbar0,change_time=change_time,casename="Hubei_step",save_files = T)
  res.step$forecast_infection

  ### continuous exponential function of pi_qbar(t)
  res.exp <- pi.eSIR(Y,R,begin_str="01/13/2020",T_fin=200,pi_qbar0=pi_qbar0,change_time=change_time,exponential=TRUE,lambda0=0.01,casename="Hubei_exp")
  res.exp$forecast_infection

  ### without pi_qbar(t)
  res.nopi <- pi.eSIR(Y,R,begin_str="01/13/2020",T_fin=200,casename="Hubei_nopi")
  res.nopi$forecast_infection
}

#' @noRd
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
#' @noRd
lognorm.parm<-function(mu0,var0){
  var <- log(var0/mu0^2+1)
  mu <- log(mu0)-var/2
  return(round(c(mu,var),3))
}

