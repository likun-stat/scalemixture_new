

###########################################################################
#  A generic Metropolis sampler.  You have to supply the log likelihood   #
#  function, which need not really be a likelihood function at all.       #
#
#
#  Uppercase K is the size of the blocks of iterations used for
#  adapting the proposal.
#  Lowercase K is the offset to get rid of wild swings in adaptation
#  process that otherwise happen the early
#  iterations.
#

# z ............................... 'data' term in the likelihood functions
# starting.theta ........................................... initial values
# likelihood.fn ..................................... likelihood in-between
# prior.fn ....................................... prior function for theta
# hyper.params ........................................ for prior functions
# n.updates .................................. number of Metropolis updates
# prop.Sigma ................................... proposal covariance matrix

# adapt.cov ......................... whether to update proposal covariance
# return.prop.Sigma.trace........ save proposal covariance from each update
# r.opt ........................................... optimal acceptance rate
# c.0, c.1 .. two coefficients for updating sigma.m and proposal covariance
# K .............................................. adapt every K iterations
#                                                                         #


adaptive.metr <- function(z, starting.theta,
                          likelihood.fn,
                          prior.fn,
                          hyper.params,
                          n.updates,
                          prop.Sigma = NULL,
                          adapt.cov = FALSE,
                          return.prop.Sigma.trace = FALSE,
                          r.opt = .234,
                          c.0 = 10,
                          c.1 = .8,
                          K = 10, ...) {

  call <- match.call()

  eps <- .001
  k <- 3  # the iteration offset

  
  p <- length(starting.theta)

  # If the supplied proposal covariance matrix is either not given or invalid,
  # just use the identity.
  if ((is.null(prop.Sigma)) |
      (length(prop.Sigma) != p^2) |
      (class(try(chol(prop.Sigma), silent=TRUE)) == "try-error")) {
    prop.Sigma <- diag(p)
    cat("Invalid or missing proposal covariance matrix.  Using identity.\n")
  }

  # Initialize sigma.m to the rule of thumb
  sigma.m <- (2.4/p)^2
  r.hat <- 0

  # Initialize prop.C
  prop.C <- chol(prop.Sigma)
  
  # Set up and initialize trace objects
  trace <- matrix(0, n.updates, p)
  sigma.m.trace <- rep(0, n.updates)
  r.trace <- rep(0, n.updates)
  jump.trace <- rep(0, n.updates)

  trace[1, ] <- starting.theta
  sigma.m.trace[1] <- sigma.m
    
  if (return.prop.Sigma.trace) {
    prop.Sigma.trace <- array(0, dim=c(n.updates, p, p))
    prop.Sigma.trace[1, , ] <- prop.Sigma
  }

  # Initialize Metropolis
  theta <- starting.theta
  likelihood <- likelihood.fn(z, theta, ...)
  prior <- prior.fn(theta, hyper.params)
  
  
  #########################################################
  # Begin main loop
  for (i in 2:n.updates) {
    theta.star <- theta + sigma.m * drop(rnorm(p) %*% prop.C)
    prior.star <- prior.fn(theta.star, hyper.params)
    if (prior.star != -Inf) {
      likelihood.star <- likelihood.fn(z, theta.star, ...)
      
      if(is.na(likelihood.star)) likelihood.star <- -Inf
      
      metr.ratio <- exp(prior.star + likelihood.star -
                        prior - likelihood)
      # if(is.na(metr.ratio)){
      #   cat(prior.star, likelihood.star,  prior, likelihood, theta, theta.star,"\n")
      # }
      if (metr.ratio > runif(1)) {
        theta <- theta.star
        prior <- prior.star
        likelihood <- likelihood.star
        jump.trace[i] <- 1
      }
    }

    ########################################################
    # Adapt via my method                                  #
    if ((i %% K) == 0) {
      gamma1 <- c.0 / ((i/K) + k)^(c.1)
      gamma2 <- 1 / ((i/K) + k)^(c.1)
      
      r.hat <- mean(jump.trace[(i - K + 1) : i])

      sigma.m <- exp(log(sigma.m) +
                     gamma1*(r.hat - r.opt))

      if (adapt.cov) {
        prop.Sigma <- prop.Sigma +
                      gamma2*(cov(trace[((i - K + 1) : i), ]) - prop.Sigma)

        while(is(try(chol(prop.Sigma), silent=TRUE), "try-error")) {
          prop.Sigma <- prop.Sigma + eps*diag(p)
          cat("Oops. Proposal covariance matrix is now:\n")
          print(prop.Sigma)
        }
        prop.C <- chol(prop.Sigma)
      }
    }
        
    #                                                      #
    ########################################################


    
    # Update the trace objects
    trace[i, ] <- theta
    sigma.m.trace[i] <- sigma.m
    r.trace[i] <- r.hat
    if (return.prop.Sigma.trace) {
      prop.Sigma.trace[i, , ] <- prop.Sigma
    }

    # Echo every 100 iterations
    if ((i %% 100) == 0)
      cat("Finished", i, "out of", n.updates, "iterations.\n")
  }
  # End main loop
  #########################################################

  # Collect trace objects to return
  res <- list(call=call,
              trace=trace,
              sigma.m.trace=sigma.m.trace,
              r.trace=r.trace,
              acc.prob=mean(jump.trace))
  if (return.prop.Sigma.trace) {
    res <- c(res, prop.Sigma.trace = list(prop.Sigma.trace))
  }

  return(res)

}

#                                                                         #
###########################################################################






###########################################################################
#  A generic Metropolis sampler.  You have to supply the log likelihood   #
#  function, which need not really be a likelihood function at all.       #
#

# z ............................... 'data' term in the likelihood functions
# starting.theta ........................................... initial values
# likelihood.fn ..................................... likelihood in-between
# prior.fn ....................................... prior function for theta
# hyper.params ........................................ for prior functions
# n.updates .................................. number of Metropolis updates
# prop.Sigma ................................... proposal covariance matrix
# sigma.m ..................................... scaling factor for proposal
# verbose ..................................... print out the error message
#                                                                         #

static.metr <- function(z, starting.theta,
                        likelihood.fn,
                        prior.fn,
                        hyper.params,
                        n.updates,
                        prop.Sigma = NULL,
                        sigma.m=NULL, 
                        verbose=FALSE, ...) {

  call <- match.call()

  eps <- .001

  p <- length(starting.theta)

  # If the supplied proposal covariance matrix is either not given or invalid,
  # just use the identity.
  if ((is.null(prop.Sigma)) |
      (length(prop.Sigma) != p^2) |
      (class(try(chol(prop.Sigma), silent=TRUE)) == "try-error")) {
    prop.Sigma <- diag(p)
    if (verbose) cat("Invalid or missing proposal covariance matrix.  Using identity.\n")
  }

  # Initialize sigma.m to the rule of thumb if it's not supplied
  if (is.null(sigma.m)) sigma.m <- (2.4/p)^2

  # Initialize prop.C
  prop.C <- chol(prop.Sigma)
  
  # Set up and initialize trace objects
  trace <- matrix(0, n.updates, p)
  jump.trace <- rep(0, n.updates)

  trace[1, ] <- starting.theta
    
  # Initialize Metropolis
  theta <- starting.theta
  likelihood <- likelihood.fn(z, theta, ...)
  prior <- prior.fn(theta, hyper.params)
  
  
  #########################################################
  # Begin main loop
  for (i in 2:n.updates) {
    theta.star <- theta + sigma.m * drop(rnorm(p) %*% prop.C)
    # cat(theta.star,hyper.params,"\n")
    prior.star <- prior.fn(theta.star, hyper.params)
    if (prior.star != -Inf) {
      likelihood.star <- likelihood.fn(z, theta.star, ...)
      
      if(is.na(likelihood.star)) likelihood.star <- -Inf
      
      metr.ratio <- exp(prior.star + likelihood.star -
                        prior - likelihood)
      if(is.na(metr.ratio)){
        metr.ratio <- 0
      }
      if (metr.ratio > runif(1)) {
        theta <- theta.star
        prior <- prior.star
        likelihood <- likelihood.star
        jump.trace[i] <- 1
      }
    }

    # Update the trace objects
    trace[i, ] <- theta

    # Echo every 100 iterations
    if ((i %% 100) == 0)
      if (verbose) cat("Finished", i, "out of", n.updates, "iterations.\n")
  }
  # End main loop
  #########################################################

  # Collect trace objects to return
  res <- list(call=call, acc.prob=mean(jump.trace[2:n.updates]),
              trace=trace)

  return(res)

}

#                                                                         #
###########################################################################




