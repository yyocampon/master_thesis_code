"
title: Additional functions to help the compute the hypothesis test.  
initial date: 28/03/2024
update date: 28/04/2024
control version: 4
Author: Ocampo-Naranjo Yeison Yovany
Objective: implementing auxiliary functions to help the understandig hypothesis test in easy functions for complex algorithm

Implemented functions:
  - simu_pp: function to simulate point pattern depending of the type clustering model used.
  - simu_pp_control: function to simulate point pattern depending of the type clustering model used with dispersion control.
  - build_theta: function to extract the estimated parameters after using kppm().
  - simulate_pp_models: function to simulate (s) realisations under conditions depending of type of cluster model, estimated parameters.
  - T0_function: 
  - Ti_function: 
  - plot_simulation_results: function to plot the statistic funtions the observed, simulated point patterns under null and altenative models and,
                             central functions on null and alternative models.
  -
  
"


# Necessary libraries -----------------------------------------------------

"The used libraries are not enable:
  - RandomFields: 3.0.62
  - RandomFieldsUtils: 1.2.5
  - spatstat: 3.0.6
  - GET: 0.4
"

{
  library(spatstat)
  library(GET)
  library(RandomFields)
  library(RandomFieldsUtils)
  library(parallel)
  library(fda.usc)
}

# Path to save results
#setwd("C:/Users/Yeison/Documents/YEISON/Maestria/Tesis/Avances_tesis_desarrollo/Scripts/Simulaciones")


{
  # simu_pp -----------------------------------------------------------------
  
  simu_pp <- function(i, params, cluster, model = 'exponential', mu=NULL , nu=0.5, params_cov){
    
    " Function to simulate clustering spatial point pattern.
    Args:
    -------------------
    parameters:
    i: silence variable
    params: vector with cluster parameter, kappa, mu, sigma
    cluster: string with type of clustering model to fit
    model: type of covariance model to use for LGCP
    mu: value with mean for LGCP model
    nu: value for models LGCP
    params_cov: list with variance and scale for covariance model in LGCP
    
    Results:
    pp_sim: pp object with simulated point pattern
  
    "
    stopifnot(cluster %in% c("Thomas", "MatClust", "LGCP"))
    
    # Validar el modelo de covarianza
    stopifnot(model %in% c("exp","exponential", "gauss", "matern"))
    
    # Inicializar el objeto de punto de proceso espacial simulado
    pp_sim <- NULL
    
    # Generar el punto de proceso espacial simulado según el tipo
    pp_sim <- switch(cluster,
                     "Thomas" = rThomas(kappa = params[[1]], scale = params[[2]], mu = params[[3]], nsim = 1),
                     "MatClust" = rMatClust(kappa = params[[1]], scale = params[[2]], mu = params[[3]], nsim = 1),
                     "LGCP" = rLGCP(model = model, mu = mu, param = params_cov),
                     stop("Tipo de proceso puntual no válido")
    )
    
    return(pp_sim)
  }
  
  
  # simu_pp_control ---------------------------------------------------------
  
  simu_pp_control <- function(i, params, cluster, model = 'exponential',
                              nu=0.5, mu=NULL, params_cov=NULL,min_points){
    
    " Function to simulate clustering spatial point pattern with control parameter for overdispersion.
    Args:
    -------------------
    parameters:
    i: silence variable
    params: vector with cluster parameter, kappa, mu, sigma
    cluster: string with type of clustering model to fit
    model: type of covariance model to use for LGCP
    mu: value with mean for LGCP model
    nu: value for models LGCP
    params_cov: list with variance and scale for covariance model in LGCP
    min_point:: value with number of point in pp reference 
    
    Results:
    pp_sim: pp object with simulated point pattern
  
    "
    
    np <- 1
    while(np < 0.8 * min_points || np > 1.2 * min_points) {
      pp_sim <- simu_pp(i = 1, params = params, cluster = cluster, model = model,
                        mu = mu, nu = nu, params_cov = params_cov)
      np <- npoints(pp_sim)
    }
    
    
    return(pp_sim)
  }
  
  # build_theta -------------------------------------------------------------
  
  build_theta <- function(X,cluster=NULL,method = "mincon",model ='exp',nu=0.5,penalised = F){
    
    " Function to compute and extract the parameters of Neyman-Scott-Cox process based on kppm function from spatstat.
    Args:
    -------------------
    parameters:
    X: pp object with the point pattern
    cluster: string with type of model to implement
    method: string with type of estimation parameter method
    model: type of covariance model to use for LGCP
    nu: value for models 
    penalised: T
    
    Results:
    theta: vector with the named values of parameters (addtional strenght cluster y prob siblings in cluster models)
  
    "
    
    k <- kppm(X,trend = ~1,clusters = cluster,model=model,nu =nu,method = method,penalised = penalised)
    
    param_est <- k#summary.kppm(k)
    
    if(cluster != "LGCP"){
      cat(paste0("Cluster model: ",cluster,"\n"))
      #theta <- c(param_est$modelpar,param_est$phi,param_est$psib)
      theta <- c(param_est$modelpar,NA,NA)
      names(theta) <- c(names(param_est$modelpar),"phi","psib")
    }else{
      cat(paste0("LGCox model: ",cluster,"\n"))
      theta <- c(param_est$modelpar)
    }
    return(theta)
  }
  
  
  # simulate_pp_models ------------------------------------------------------
  
  simulate_pp_models <- function(Mod_H, nsimsub, model ='exponential', theta, nu = NULL,min_points) {
    
    " Function to simulate multiple point pattern for Monte Carlo method
    Args:
    -------------------
    parameters:
    Mod_H: string with clustering point pattern
    nsimsub: value with number of simulations of Monte Carlo Method 
    model: string with covariance model
    theta: vector with estimated parameters for clustering model using kppm from spatstat
    min_points: number of point to simulate and controlling overdispersion
    
    Results:
    sim.X: list with pp objects, nsimsub point patterns.
  
    "
    
    if (Mod_H == "LGCP") {
      cat(paste0("Entrando a LGCP para simular ", nsimsub ," patrones \n"))
      sim.X <- lapply(X = 1:nsimsub, FUN = simu_pp_control, model = model, mu = theta[[3]],
                      params_cov = list(theta[[1]],theta[[2]]), nu = nu, cluster = Mod_H,
                      min_points =min_points)
    } else {
      if(Mod_H == 'MatClust'){
        cat(paste0("Entrando a MatClust para simular ", nsimsub ," patrones \n"))
      }else{
        cat("Entrando Thomas para simular ", nsimsub ," patrones \n")
      }
      sim.X <- lapply(X = 1:nsimsub, FUN = simu_pp_control, params = list(theta[1], theta[2], theta[3]), cluster = Mod_H,
                      min_points = min_points)
    }
    return(sim.X)
  }
  
  # T0_function -------------------------------------------------------------
  
  T0_function <- function(X, f.dist, corrector, r_max) {
    
    " Function to compute function based on distance
    Args:
    -------------------
    parameters:
    X: pp object with spatial point pattern
    f.dist: string with selected function based on distance (G, F, J)
    corrector: string border correction to use
    r_max, value with máx value of r in intervale I =(0,r_max)
    
    
    Results:
    T.Statistic: matrix with values of statistic based on distance
  
    "
    # Validate f.dist
    stopifnot(!is.null(f.dist), f.dist %in% c('Kest', 'Lest', 'Jest'))
    
    # Validate corrector
    stopifnot(
      (f.dist %in% c('Kest', 'Lest') & corrector %in% c('border', 'iso', 'Ripley', 'trans')) ||
        (f.dist %in% c("Jest") & corrector %in% c("rs", "km", "cs"))
    )
    
    # Validate r_max
    stopifnot(is.numeric(r_max) && r_max > 0)
    
    T.Statistic <- switch(
      f.dist,
      Kest = Kest(X, correction = corrector, rmax = r_max),
      Lest = Lest(X, correction = corrector, rmax = r_max),
      Jest = {
        T.Statistic <- Jest(X, rmax = r_max, correction = corrector)
        T.Statistic[which(is.na(T.Statistic))] <- 0
        T.Statistic
      }
    )
    
    return(T.Statistic)
  }
  
  # Ti_function -------------------------------------------------------------
  
  Ti_function <- function(X, f.dist, corrector, r_max) {
    
    " Function to compute function based on distance
    Args:
    -------------------
    parameters:
    X: pp object with spatial point pattern
    f.dist: string with selected function based on distance (G, F, J)
    corrector: string border correction to use
    r_max, value with máx value of r in intervale I =(0,r_max)
    
    
    Results:
    T.Statistic: matrix with values of statistic based on distance
  
    "
    
    # Validate f.dist
    stopifnot(!is.null(f.dist), f.dist %in% c('Kest', 'Lest', 'Jest'))
    
    # Validate corrector
    # Validate corrector
    stopifnot(
      (f.dist %in% c('Kest', 'Lest') & corrector %in% c('border', 'iso', 'Ripley', 'trans')) ||
        (f.dist %in% c("Jest") & corrector %in% c("rs", "km", "cs"))
    )
    
    # Validate r_max
    stopifnot(is.numeric(r_max) && r_max > 0)
    
    T.Statistic <- switch(
      f.dist,
      Kest = Kest(X, correction = corrector, rmax = r_max)[[corrector]],
      Lest = Lest(X, correction = corrector, rmax = r_max)[[corrector]],
      Jest = {
        T.Statistic <- Jest(X, rmax = r_max, correction = corrector)[[corrector]]
        T.Statistic[which(is.na(T.Statistic))] <- 0
        T.Statistic
      }
    )
    
    return(T.Statistic)
  }
  
  
  
  # plot_simulation_results -------------------------------------------------
  
  plot_simulation_results <- function(Mod_H0, Mod_H1, sim.M1, sim.M2, MT1, MT2, nsimsub, graphs = F,valor_p) {
    
    " Function to plot all calculated statistic function used
    Args:
    -------------------
    parameters:
    Mod_H0: string with clustering model under null hypothesis
    Mod_H1: string with clustering model under alternative hypothesis
    sim.M1: list with spatial point patterns simulated under null hypothesis
    sim.M2: list with spatial point patterns simulated alternative null hypothesis
    MT1: vector with central functions under null hypothesis
    MT2: vector with central functions under alternative hypothesis
    nsimsub: value with number of simulations
    graphs = F
    valor_p: value from  computed p-value by test
    
    
    Results:
    plot with all plotted curves 
  
    "
    # Validating the 'graphs' parameter
    if (!(graphs %in% c(T, F))) {
      stop("Invalid value for 'graphs'. Use 'T' or 'F'.")
    }
    
    # Computing average curve from H0 and H1
    
    if (graphs == T) {
      plot(sim.M1[, 1], sim.M1[, 3], type = "l", col = "red", lwd = 2, xlab = 'r', ylab = 'J(r)',las=1,
           #main = paste0("Valor-p método Montecarlo: ",valor_p)
           )
      
      for (i in 1:(nsimsub - 1)) {
        lines(sim.M1[, 1], sim.M1[, (i + 3)], col = "lightblue")
        lines(sim.M2[, 1], sim.M2[, (i + 3)], col = "lightsalmon")
      }
      
      lines(sim.M1[, 1], MT1, col = "blue", lwd = 3)  # T_prom H0
      lines(sim.M2[, 1], MT2, col = "red", lwd = 3)  # T_prom H1
      lines(sim.M2[, 1], sim.M2[, 3], col = "black", lwd = 3)  # T_obs
      
      legend("topright", legend = c(bquote(H[0] ~ T[ sim ] ~  ": " ~ .(Mod_H0)),
                                    bquote(H[1] ~ T[ sim ] ~  ": " ~ .(Mod_H1)),
                                    substitute(T[obs], list(obs = "obs"))
      ),
      col = c("blue", "red", "black"), lwd = 3, bty = "n", cex = 1)
    }
  }
  
  
  # pvalfunction ------------------------------------------------------------
  
  pvalfunction <- function(X,r_max =NULL,Mod_H0,Mod_H1,f.dist,corrector,method,nsimsub,
                           params = NULL,model = NULL,nu=NULL,penalised=TRUE,graphs =F){
    " Function to plot all calculated statistic function used
    Args:
    -------------------
    parameters:
    X: pp object with spatial point pattern
    Mod_H0: string with clustering model under null hypothesis
    Mod_H1: string with clustering model under alternative hypothesis
    f.dist: string with selected function based on distance (G, F, J)
    corrector: string border correction to use
    method: string with method to estimate parameters
    model: type of covariance model to use for LGCP
    mu: value with mean for LGCP model
    nu: value for models LGCP
    params: vector with cluster parameter, kappa, mu, sigma
    params_cov: list with variance and scale for covariance model in LGCP
    r_max, value with máx value of r in intervale I =(0,r_max)
    nsimsub: value with number of simulations
    graphs = F
    
    
    Results:
    results: list with p-value and estimated parameters
  
    "
    gc()
    
    k_for_min <- kppm(X = X,trend = ~1,clusters = Mod_H0)
    min_point <- npoints(X) #round(k_for_min$lambda*0.8)
    #max_point <- round(k_for_min$lambda*1.20)
    
    cat(paste0("Estimating parameters under H0: ",Mod_H0," \n"))
    theta0 <- build_theta(X = X,cluster = Mod_H0, method = method, penalised = T)
    cat(paste0("Estimating parameters under H1: ",Mod_H1," \n"))
    theta1 <- build_theta(X = X,cluster = Mod_H1, method = method, penalised = T)
    
    # T0(r) = obs
    T0 <- T0_function(X = X,f.dist = f.dist,corrector = corrector,r_max = r_max)
    r <- T0[["r"]]
    
    # simulating nsubsim using theta[0], theta[1] X11,X12,...,X1(s+1) /  X21,X22,...,X2(s+1)
    sim.X1 <- simulate_pp_models(Mod_H = Mod_H0,nsimsub = nsimsub,
                                 theta = theta0,
                                 min_points= min_point)
    sim.X2 <- simulate_pp_models(Mod_H = Mod_H1,nsimsub = nsimsub,
                                 theta = theta1,
                                 min_points= min_point)
    
    # computing Ti(r) functions, i= 1,2,..., s+1
    Tr_M1 <- cbind(as.matrix(T0),mapply(sim.X1, FUN = Ti_function, f.dist, corrector,r_max))
    Tr_M2 <- cbind(as.matrix(T0),mapply(sim.X2, FUN = Ti_function, f.dist, corrector,r_max))
    
    # Computing average curve from H0 and H1 (central curves)
    MT1 <- rowMeans(Tr_M1[,4:ncol(Tr_M1)])
    MT2 <- rowMeans(Tr_M2[,4:ncol(Tr_M2)])
    
    
    
    # Centering every curve at the central curve (computing differences between Ti and central curve)
    centered_H1 <- (sweep(x = Tr_M2[,3:ncol(Tr_M2)],MARGIN =  1,STATS =  MT2,FUN = "-"))^2 # (T(r)-T1(r))^2 (num)
    centered_H0 <- (sweep(x = Tr_M1[,3:ncol(Tr_M1)],MARGIN =  1,STATS =  MT1,FUN = "-"))^2 # (T(r)-T0(r))^2 (den)
    
    # Computing difference (T(r)-T1(r))^2 - (T(r)-T0(r))^2 to integrate 
    statistic_diff <- cbind(r,arg_diff_int_red =(centered_H1-centered_H0)) 
    
    # Computing int[T(r)-T1(r))^2- (T(r)-T0(r))^2]dr, r in [0,r_max]
    int_difference <- apply(X = statistic_diff[,-1], 2,
                            int.simpson2, x = statistic_diff[,1], method = "CSR")
    
    # Ordering 
    R0_dif <- int_difference[1]
    Ri_dif <- int_difference[-1]
    p_value_diff <- (sum( (Ri_dif <= R0_dif) )+1)/ (nsimsub + 1)
    
    
    # Computing p-value for GET based on statistic_diff
    S_i_diff <- statistic_diff[,-1]
    
    #Construct a curve_set object to use in GET
    cset_diff <- create_curve_set(list(r = r, obs = S_i_diff[,1], sim_m = S_i_diff[,-1]))
    GET_diff <- global_envelope_test(cset_diff, type = "rank", alternative = "two.sided")
    p_value_GET_diff <- attr(GET_diff, "p")
    
    cat("Finish p0-value execution \n")
    results_p <- c(p_value_diff,p_value_GET_diff)
    names(results_p) <- c("p_value_test_diff","p_value_GET_diff")
    
    results <- list("theta0" = theta0,
                    "theta1" = theta1,
                    "H0" = Mod_H0,
                    "H1" = Mod_H1,
                    "p_value"=results_p)
    
    plot_simulation_results(Mod_H0, Mod_H1,Tr_M1, Tr_M2, MT1, MT2, nsimsub ,
                            graphs = graphs,valor_p = round(results$p_value[1],3))
    
    return(results)
  }
  
  
  
  # Graficar pp usando ggplot -----------------------------------------------
  
  
  plot_ppp <- function(pp, main_title = F, additional_text = F) {
    # Crear un dataframe a partir del objeto ppp
    df <- as.data.frame(pp)
    
    # Crear el gráfico ggplot
    g1 <- ggplot(df, aes(x, y)) +
      geom_point(size = 1) +
      theme_bw(base_size = 8) +
      labs(x = "", y = "")+#, title = main_title) +
      #ggtitle(paste0("Patrón observado: ", main_title,additional_text)) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        axis.title.x = element_blank(),      # Eliminar título del eje X
        axis.title.y = element_blank(),      # Eliminar título del eje Y
        axis.text.x = element_blank(),       # Eliminar etiquetas del eje X
        axis.text.y = element_blank(),       # Eliminar etiquetas del eje Y
      )
    return(g1)
  }
  
  
  
  
  
  
}
