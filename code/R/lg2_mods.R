#' @name lg2_mods
#' @aliases disc_log
#' @aliases disc_log_E
#' @aliases disc_LV_comp
#' @aliases disc_LV_E
#' @aliases disc_LV_evol
#' @title Discrete-time logistic population growth, with options for competition, environmental covariate, and trait evolution
#'
#' @description
#' disc_log and disc_log_E are functions for logistic growth in a single species, with user-provided growth rate, initial population size, and strength of intraspecific competition. The function disc_log_E includes an environmental covariate, and growth rate now depends on this covariate as a function of the distance from the optimal trait value. XX disc_LV_comp XX disc_LV_E XX disc_LV_evol XX
#'
#' @param r the user-supplied population growth rate (when the species trait x is at the environmental optimum)
#' @param N0 the initial population size
#' @param alpha the intraspecific interaction coefficient
#' @param E the initial value of the environmental covariate
#' @param x the initial value of the species trait
#' @param P xx
#' @param w xx
#' @param Wmax xx
#' @param h2 xx
#' @return xx
#'
#' @examples
#' # Simulate initial species population growth
#' N1.0 <- rpois(1,10)
#' r1.0 <- 1.67
#' alpha.11 <- 0.00125
#' # Simulation of model for t time steps
#' t <- 30
#' N <- rep(NA, t)
#' N[1] <- N1.0
#' for (i in 2:t) {
#'   N[i] <- disc_log(r=r1.0, N0=N[i-1], alpha=alpha.11)
#' }
#'
#' @rdname lg2_mods
#' @export
disc_log <- function(r, N0, alpha) {
  Nt1 <- (r*N0) / (1+alpha*N0)
  return(Nt1)
}

#' @rdname lg2_mods
#' @export
disc_log_E <- function(r, N0, alpha, E, x) {
  Nt1 <- ((r*exp(-(E-x)^2))*N0) / (1+alpha*N0)
  return(Nt1)
}

#' @rdname lg2_mods
#' @export
disc_LV_comp <- function(r,N0,alpha) {
  num <- r*N0
  denom <- rowSums(alpha %*% diag(N0))
  Nt1 <- mapply(function(x,y) x / (1+y), num, denom)
  return(Nt1)
}

#' @rdname lg2_mods
#' @export
disc_LV_E <- function(r,N0,alpha,E,x) {
  num <- r*exp(-(E-x)^2)*N0
  denom <- rowSums(alpha %*% diag(N0))
  Nt1 <- mapply(function(x,y) x / (1+y), num, denom)
  return(Nt1)
}

#' @rdname lg2_mods
#' @export
disc_LV_evol <- function(N0,alpha,E,x,P,w,Wmax,h2) {
  if(length(N0)==1){N0 = as.matrix(N0)}
  What <- Wmax*sqrt(w/(P+w))
  r <- What*exp((-(((w+(1-h2)*P)/(P+w))*(E-x))^2)/(2*(P+w)))
  num <- r*N0
  denom <- rowSums(alpha %*% diag(N0))
  Nt1 <- mapply(function(x,y) x / (1+y), num, denom)
  return(list(Nt1=Nt1,r=r))
}
