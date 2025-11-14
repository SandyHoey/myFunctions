#function to bootstrap parameter confidence intervals from glmer models

#works by simulating data from the model, and the refitting the model with the simulated data and recording the parameters

#currently works with
#' lme4 and glmmTMB
#' Poisson
#' negative binomial
#' binomial
#' zero inflated Poisson
#' zero inflated negative binomial 


#' #loading packages
library(dplyr)
library(ggplot2)

#defining function
boot_param_CI <- function(nsim, model, data){
  
  #so I don't have to load the whole package
  `%like%` <- data.table::`%like%`
  
  
  # lme4 --------------------------------------------------------------------
  ## code for lme4 models
  if(inherits(model, "merMod")){
    #bootstrapping parameter values from model simulations
    betas <- matrix(NA, nrow = nsim,
                    ncol = length(fixef(model)))
    
    for(j in 1:nsim){
      sim_data <- data %>% 
        mutate(y = unlist(simulate(model)))  #simulate response variables from the model
      
      sim_model <- update(model, y ~ ., data = sim_data) #rerun the model with the simulated values as the response
      
      if(is.null(summary(sim_model)$optinfo$conv$lme4$messages) == TRUE){ #if there are no warning messages (model fit fine) then record the parameters
        betas[j,] <- fixef(sim_model)
      }
    }
    
    #seeing how many models fit well
    n_fit <- betas %>% 
      na.omit %>% 
      nrow()
    
    ## Poisson ----
    if(model@call$family == "poisson"){  
      beta_bs <- data.frame(FE = names(fixef(model)),
                            coef = exp(fixef(model)),
                            lower = exp(apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T))), #CI lower bound
                            upper = exp(apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T)))) #CI upper bound
      
      #plotting model coefficients
      (beta_plot <- beta_bs %>% 
          filter(FE != "(Intercept)") %>% 
          ggplot + 
          geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
          geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
          geom_vline(xintercept = 0, lty = "dashed") +
          geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                        vjust = -.6, hjust = .3), size = 3) +
          theme(legend.position = "none") +
          labs(title = paste0("fixed effects (iterations = ", n_fit, ")"),
               subtitle = "Poisson",
               y = "",
               x = "exp(\u03b2)"))
    }
    
    ## negative binomial ----
    if(any(model@call$family %like% "negative.binomial")){  
      beta_bs <- data.frame(FE = names(fixef(model)),
                            coef = exp(fixef(model)),
                            lower = exp(apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T))), #CI lower bound
                            upper = exp(apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T)))) #CI upper bound
      
      #plotting model coefficients
      (beta_plot <- beta_bs %>% 
          filter(FE != "(Intercept)") %>% 
          ggplot + 
          geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
          geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
          geom_vline(xintercept = 0, lty = "dashed") +
          geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                        vjust = -.6, hjust = .3), size = 3) +
          theme(legend.position = "none") +
          labs(title = paste0("fixed effects (iterations = ", n_fit, ")"),
               subtitle = "negative binomial",
               y = "",
               x = "exp(\u03b2)"))
    }
    
    ## binomial ----
    if(model@call$family == "binomial"){  
      
      beta_bs <- data.frame(FE = names(fixef(model)),
                            coef = boot::inv.logit(fixef(model)),
                            lower = boot::inv.logit(apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T))), #CI lower bound
                            upper = boot::inv.logit(apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T)))) #CI upper bound
      
      
      #plotting model coefficients
      (beta_plot <- beta_bs %>% 
          filter(FE != "(Intercept)") %>% 
          ggplot + 
          geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
          geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
          geom_vline(xintercept = 0.5, lty = "dashed") +
          geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                        vjust = -.6, hjust = .3), size = 3) +
          theme(legend.position = "none") +
          labs(title = paste0("fixed effects (iterations = ", n_fit, ")"),
               subtitle = "binomial",
               y = "",
               x = "inv.logit(\u03b2)"))
    }
  }
  
  
  # glmmTMB -----------------------------------------------------------------
  ## code for glmmTMB models
  if(inherits(model, "glmmTMB")){
    
    ## zero inflated ----
    if(length(fixef(model)$zi) != 0){
      #bootstrapping parameter values from model simulations
      betas <- matrix(NA, nrow = nsim,
                      ncol = length(fixef(model)$cond) + length(fixef(model)$zi))
      
      for(j in 1:nsim){
        sim_data <- data %>% 
          mutate(y = unlist(simulate(model)))  #simulate response variables from the model
        
        sim_model <- update(model, y ~ ., data = sim_data) #rerun the model with the simulated values as the response
        
        if(sim_model$fit$convergence == 0 &
           sim_model$sdr$pdHess){ #if model convergence looks good then record the parameters
          betas[j,] <- c(fixef(sim_model)$cond, fixef(sim_model)$zi)
        }
      }
      
      #seeing how many models fit well
      n_fit <- betas %>% 
        na.omit %>% 
        nrow()
      
      ### Poisson ----
      if(model$call$family == "poisson"){  
        beta_bs <- data.frame(FE = c(paste0("cond:", names(fixef(sim_model)$cond)), 
                                     paste0("zi:", names(fixef(sim_model)$zi))),
                              coef = exp(c(fixef(sim_model)$cond, fixef(sim_model)$zi)),
                              lower = exp(apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T))), #CI lower bound
                              upper = exp(apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T)))) #CI upper bound
        
        #plotting model coefficients
        (beta_plot <- beta_bs %>% 
            filter(FE != "cond:(Intercept)",
                   FE != "zi:(Intercept)") %>% 
            ggplot + 
            geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
            geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
            geom_vline(xintercept = 0, lty = "dashed") +
            geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                          vjust = -.6, hjust = .3), size = 3) +
            theme(legend.position = "none") +
            labs(title = paste0("fixed effects (iterations = ", n_fit, ")"),
                 subtitle = "Poisson",
                 y = "",
                 x = "exp(\u03b2)"))
      }
      
      ### negative binomial ----
      if(model$call$family == "nbinom2"){  
        beta_bs <- data.frame(FE = c(paste0("cond:", names(fixef(sim_model)$cond)), 
                                     paste0("zi:", names(fixef(sim_model)$zi))),
                              coef = exp(c(fixef(sim_model)$cond, fixef(sim_model)$zi)),
                              lower = exp(apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T))), #CI lower bound
                              upper = exp(apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T)))) #CI upper bound
        
        #plotting model coefficients
        (beta_plot <- beta_bs %>% 
            filter(FE != "cond:(Intercept)",
                   FE != "zi:(Intercept)") %>% 
            ggplot + 
            geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
            geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
            geom_vline(xintercept = 0, lty = "dashed") +
            geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                          vjust = -.6, hjust = .3), size = 3) +
            theme(legend.position = "none") +
            labs(title = paste0("fixed effects (iterations = ", n_fit, ")"),
                 subtitle = "negative binomial",
                 y = "",
                 x = "exp(\u03b2)"))
      }
    }
    
    ## non-zero inflated models ----
    if(length(fixef(model)$zi) == 0){
      #bootstrapping parameter values from model simulations
      betas <- matrix(NA, nrow = nsim,
                      ncol = length(fixef(model)$cond))
      
      for(j in 1:nsim){
        sim_data <- data %>% 
          mutate(y = unlist(simulate(model)))  #simulate response variables from the model
        
        sim_model <- update(model, y ~ ., data = sim_data) #rerun the model with the simulated values as the response
        
        if(sim_model$fit$convergence == 0 &
           sim_model$sdr$pdHess){ #if there are no warning messages (model fit fine) then record the parameters
          betas[j,] <- fixef(sim_model)$cond
        }
      }
      
      #seeing how many models fit well
      n_fit <- betas %>% 
        na.omit %>% 
        nrow()
      
      ### Poisson ----
      if(model$call$family == "poisson"){  
        beta_bs <- data.frame(FE = names(fixef(sim_model)$cond),
                              coef = exp(fixef(sim_model)$cond),
                              lower = exp(apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T))), #CI lower bound
                              upper = exp(apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T)))) #CI upper bound
        
        #plotting model coefficients
        (beta_plot <- beta_bs %>% 
            filter(FE != "(Intercept)") %>% 
            ggplot + 
            geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
            geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
            geom_vline(xintercept = 0, lty = "dashed") +
            geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                          vjust = -.6, hjust = .3), size = 3) +
            theme(legend.position = "none") +
            labs(title = paste0("fixed effects (iterations = ", n_fit, ")"),
                 subtitle = "Poisson",
                 y = "",
                 x = "exp(\u03b2)"))
      }
      
      ### negative binomial ----
      if(model$call$family == "nbinom2"){  
        beta_bs <- data.frame(FE = names(fixef(sim_model)$cond),
                              coef = exp(fixef(sim_model)$cond),
                              lower = exp(apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T))), #CI lower bound
                              upper = exp(apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T)))) #CI upper bound
        
        #plotting model coefficients
        (beta_plot <- beta_bs %>% 
            filter(FE != "(Intercept)") %>% 
            ggplot + 
            geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
            geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
            geom_vline(xintercept = 0, lty = "dashed") +
            geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                          vjust = -.6, hjust = .3), size = 3) +
            theme(legend.position = "none") +
            labs(title = paste0("fixed effects (iterations = ", n_fit, ")"),
                 subtitle = "negative binomial",
                 y = "",
                 x = "exp(\u03b2)"))
      }
      
      ### binomial ----
      if(model$call$family == "binomial"){  
        
        beta_bs <- data.frame(FE = names(fixef(sim_model)$cond),
                              coef = boot::inv.logit(fixef(sim_model)$cond),
                              lower = boot::inv.logit(apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T))), #CI lower bound
                              upper = boot::inv.logit(apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T)))) #CI upper bound
        
        
        #plotting model coefficients
        (beta_plot <- beta_bs %>% 
            filter(FE != "(Intercept)") %>% 
            ggplot + 
            geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
            geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
            geom_vline(xintercept = 0.5, lty = "dashed") +
            geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                          vjust = -.6, hjust = .3), size = 3) +
            theme(legend.position = "none") +
            labs(title = paste0("fixed effects (iterations = ", n_fit, ")"),
                 subtitle = "binomial",
                 y = "",
                 x = "inv.logit(\u03b2)"))
      }
    }
  }
  
  
  
  return(list(betas, beta_bs, beta_plot))
}

