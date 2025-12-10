#' Bootstrapping confidence intervals for lme4 and glmmTMB
#'
#' @param nsim number of iteration for the bootstrap
#' @param model lme4 (merMod) or glmmTMB model object
#' @param data dataframe used in the model
#' @return returns model parameters for all bootstrap iterations and 95% confidence interval plots and dataframe 
#' @export
#' @import stats
#' @import dplyr
#' @import ggplot2
#' @import patchwork
#' @importFrom dplyr filter
#' @importFrom lme4 fixef
#' @importFrom data.table %like%


#defining function
boot_param_CI <- function(nsim, model, data){
  
  # lme4 --------------------------------------------------------------------
  ## code for lme4 models
  if(inherits(model, "merMod")){
     #bootstrapping parameter values from model simulations
    betas <- matrix(NA, nrow = nsim,
                    ncol = length(fixef(model)))
    
    for(j in 1:nsim){
      sim_data <- data %>% 
        mutate(y = unlist(simulate(model)))  # simulate response variables from the model
      
      sim_model <- update(model, y ~ ., data = sim_data) # rerun the model with the simulated values as the response
      
      if(is.null(summary(sim_model)$optinfo$conv$lme4$messages) == TRUE){ # if there are no warning messages (model fit fine) then record the parameters
        betas[j,] <- fixef(sim_model)
      }
    }
    
    # seeing how many models fit well
    n_fit <- betas %>% 
      na.omit %>% 
      nrow()
    
    ## Poisson ----
    if(model@call$family == "poisson"){  
      beta_bs <- data.frame(FE = names(fixef(model)),
                            coef = exp(fixef(model)),
                            lower = exp(apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T))), #CI lower bound
                            upper = exp(apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T)))) #CI upper bound
      
      # plotting model coefficients
      (beta_plot <- beta_bs %>% 
          dplyr::filter(FE != "(Intercept)") %>% 
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
      
      # plotting model coefficients
      (beta_plot <- beta_bs %>% 
          dplyr::filter(FE != "(Intercept)") %>% 
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
                            coef = plogis(fixef(model)),
                            lower = plogis(apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T))), #CI lower bound
                            upper = plogis(apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T)))) #CI upper bound
      
      
      # plotting model coefficients
      (beta_plot <- beta_bs %>% 
          dplyr::filter(FE != "(Intercept)") %>% 
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
      # bootstrapping parameter values from model simulations
      betas <- matrix(NA, nrow = nsim,
                      ncol = length(fixef(model)$cond) + length(fixef(model)$zi))
      
      for(j in 1:nsim){
        sim_data <- data %>% 
          mutate(y = unlist(simulate(model)))  # simulate response variables from the model
        
        sim_model <- update(model, y ~ ., data = sim_data) #rerun the model with the simulated values as the response
        
        if(sim_model$fit$convergence == 0 &
           sim_model$sdr$pdHess){ # if model convergence looks good then record the parameters
          betas[j,] <- c(fixef(sim_model)$cond, fixef(sim_model)$zi)
        }
      }
      
      # seeing how many models fit well
      n_fit <- betas %>% 
        na.omit %>% 
        nrow()
      
      ### Poisson ----
      if(model$call$family == "poisson"){  
        beta_bs <- data.frame(FE = c(names(fixef(sim_model)$cond), 
                                     names(fixef(sim_model)$zi)),
                              model = c(rep("conditional", length(fixef(sim_model)$cond)), # column for the model type
                                        rep("zi", length(fixef(sim_model)$zi))),
                              coef = c(exp(fixef(sim_model)$cond), plogis(fixef(sim_model)$zi)),
                              lower = apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T)),      #CI lower bound
                              upper = apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T))) %>%  #CI upper bound 
          # back-transforming CI bounds based on the model (con, zi)
          # if the model type is conditional (poisson), the link funciton is log() and is back-transformed with exp()
          # otherwise the model type is zero-inflated (binomial), the link function is logit() and is back-transformed with plogis()
          mutate(lower = if_else(model == "conditional", exp(lower), plogis(lower)),
                 upper = if_else(model == "conditional", exp(upper), plogis(upper))) %>% 
          # turning into a list with a dataframe for each model part (cond, zi) so they plot separately
          group_by(model) %>% 
          group_split()
        
        # setting dataframe names in list
        names(beta_bs) <- c("conditional", "zi")
          
        
        # plotting model coefficients
          # conditional model
        beta_plot_cond <- beta_bs[[1]] %>% 
          dplyr::filter(FE != "(Intercept)") %>% 
          ggplot + 
          geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
          geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
          geom_vline(xintercept = 0, lty = "dashed") +
          geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                        vjust = -.6, hjust = .3), size = 3) +
          theme(legend.position = "none") +
          labs(title = "Conditional model fixed effects",
               subtitle = paste0("Poisson (iterations = ", n_fit, ")"),
               y = "",
               x = "exp(\u03b2)")
        
          # zero inflated model
        beta_plot_zi <- beta_bs[[2]] %>% 
          dplyr::filter(FE != "(Intercept)") %>% 
          ggplot + 
          geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
          geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
          geom_vline(xintercept = 0.5, lty = "dashed") +
          geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                        vjust = -.6, hjust = .3), size = 3) +
          theme(legend.position = "none") +
          labs(title = "Zero inflated model fixed effects",
               subtitle = paste0("Binomial (iterations = ", n_fit, ")"),
               y = "",
               x = "inv.logit(\u03b2)")
        
        #patching plots together
        beta_plot <- beta_plot_cond / beta_plot_zi
      }
      
      ### negative binomial ----
      if(model$call$family == "nbinom2"){  
        beta_bs <- data.frame(FE = c(names(fixef(sim_model)$cond), 
                                     names(fixef(sim_model)$zi)),
                              model = c(rep("conditional", length(fixef(sim_model)$cond)), # column for the model type
                                        rep("zi", length(fixef(sim_model)$zi))), 
                              coef = c(exp(fixef(sim_model)$cond), plogis(fixef(sim_model)$zi)),
                              lower = apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T)),      #CI lower bound
                              upper = apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T))) %>%  #CI upper bound 
          # back-transforming CI bounds based on the model (con, zi)
          # if the model type is conditional (negative binomial), the link funciton is log() and is back-transformed with exp()
          # otherwise the model type is zero-inflated (binomial), the link function is logit() and is back-transformed with plogis()
          mutate(lower = if_else(model == "conditional", exp(lower), plogis(lower)),
                 upper = if_else(model == "conditional", exp(upper), plogis(upper))) %>% 
          # turning into a list with a dataframe for each model part (cond, zi) so they plot separately
          group_by(model) %>% 
          group_split()
        
        # setting dataframe names in list
        names(beta_bs) <- c("conditional", "zi")
        
        
        # plotting model coefficients
          # conditional model
        beta_plot_cond <- beta_bs[[1]] %>% 
          dplyr::filter(FE != "(Intercept)") %>% 
          ggplot + 
          geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
          geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
          geom_vline(xintercept = 0, lty = "dashed") +
          geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                        vjust = -.6, hjust = .3), size = 3) +
          theme(legend.position = "none") +
          labs(title = "Conditional model fixed effects",
               subtitle = paste0("Negative binomial (iterations = ", n_fit, ")"),
               y = "",
               x = "exp(\u03b2)")
        
          # zero inflated model
        beta_plot_zi <- beta_bs[[2]] %>% 
          dplyr::filter(FE != "(Intercept)") %>% 
          ggplot + 
          geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
          geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
          geom_vline(xintercept = 0.5, lty = "dashed") +
          geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                        vjust = -.6, hjust = .3), size = 3) +
          theme(legend.position = "none") +
          labs(title = "Zero inflated model fixed effects",
               subtitle = paste0("Binomial  (iterations = ", n_fit, ")"),
               y = "",
               x = "inv.logit(\u03b2)")
        
        #patching plots together
        beta_plot <- beta_plot_cond / beta_plot_zi
      }
    }
    
    ## non-zero inflated models ----
    if(length(fixef(model)$zi) == 0){
      # bootstrapping parameter values from model simulations
      betas <- matrix(NA, nrow = nsim,
                      ncol = length(fixef(model)$cond))
      
      for(j in 1:nsim){
        sim_data <- data %>% 
          mutate(y = unlist(simulate(model)))  # simulate response variables from the model
        
        sim_model <- update(model, y ~ ., data = sim_data) # rerun the model with the simulated values as the response
        
        if(sim_model$fit$convergence == 0 &
           sim_model$sdr$pdHess){ # if there are no warning messages (model fit fine) then record the parameters
          betas[j,] <- fixef(sim_model)$cond
        }
      }
      
      # seeing how many models fit well
      n_fit <- betas %>% 
        na.omit %>% 
        nrow()
      
      ### Poisson ----
      if(model$call$family == "poisson"){  
        beta_bs <- data.frame(FE = names(fixef(sim_model)$cond),
                              coef = exp(fixef(sim_model)$cond),
                              lower = exp(apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T))), #CI lower bound
                              upper = exp(apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T)))) #CI upper bound
        
        # plotting model coefficients
        (beta_plot <- beta_bs %>% 
            dplyr::filter(FE != "(Intercept)") %>% 
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
        
        # plotting model coefficients
        (beta_plot <- beta_bs %>% 
            dplyr::filter(FE != "(Intercept)") %>% 
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
                              coef = plogis(fixef(sim_model)$cond),
                              lower = plogis(apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T))), #CI lower bound
                              upper = plogis(apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T)))) #CI upper bound
        
        
        # plotting model coefficients
        (beta_plot <- beta_bs %>% 
            dplyr::filter(FE != "(Intercept)") %>% 
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

