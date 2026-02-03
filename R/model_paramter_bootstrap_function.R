#' Bootstrapping confidence intervals for lme4 and glmmTMB
#'
#' @param nsim number of iteration for the bootstrap
#' @param model lme4 (merMod) or glmmTMB model object
#' @param data dataframe used in the model
#' @param newData dataframe used for model predictions
#' @return returns model parameters for all bootstrap iterations and 95% confidence interval plots and dataframe 
#' @export
#' @import stats
#' @import dplyr
#' @import lme4
#' @import ggplot2
#' @import patchwork
#' @importFrom dplyr filter
#' @importFrom data.table %like%


#defining function
boot_param_CI <- function(nsim, model, data, newData = NULL){

  
  # lme4 --------------------------------------------------------------------
  ## code for lme4 models
  if(inherits(model, "merMod")){
    #bootstrapping parameter values from model simulations
    betas <- matrix(NA, nrow = nsim,
                    ncol = length(fixef(model)))
    
    # preparing for model predictions
    if(!is.null(newData)){
      # creating a data frame to put predicted values form each iteration
      pv <- matrix(NA, nrow = nsim,
                   ncol = nrow(newData))
      
      # getting the sd for the random effect if there is one
      # if this is not created, there is no random effect
      re <- as.data.frame(VarCorr(model))$sdcor
    }
    
    
    for(j in 1:nsim){
      sim_data <- data %>% 
        mutate(y = unlist(simulate(model)))  # simulate response variables from the model
      
      sim_model <- update(model, y ~ ., data = sim_data) # rerun the model with the simulated values as the response
      
      if(is.null(summary(sim_model)$optinfo$conv$lme4$messages) == TRUE){ # if there are no warning messages (model fit fine) then record the parameters
        betas[j,] <- fixef(sim_model)
        
        # doing model prediction if new data is provided
        if(!is.null(newData)){
          # if there is no random effect
          if(missing(re)){
            pv[j,] <- predict(sim_model, newData, re.form = NA)
            
          } else(pv[j,] <- predict(sim_model, newData, re.form = NA) + rnorm(1, 0, sd = re))
        }
      }
    }
    
    # seeing how many models fit well
    n_fit <- betas %>% 
      na.omit %>% 
      nrow()
    
    ## Poisson ----
    if(model@call$family == "poisson"){  
      beta_bs <- data.frame(FE = names(fixef(model)),
                            coef = fixef(model),
                            lower = apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T)), #CI lower bound
                            upper = apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T))) #CI upper bound
      
      # plotting model coefficients
      (beta_plot <- beta_bs %>% 
          dplyr::filter(FE != "(Intercept)") %>% 
          ggplot + 
          geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
          # confidence intervals
          geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
          # creating dashed line around 0
          geom_vline(xintercept = 0, lty = "dashed") +
          # adding model coefficient value to plot
          geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                        vjust = -.6, hjust = .3), size = 3) +
          theme(legend.position = "none") +
          labs(title = paste0("fixed effects (iterations = ", n_fit, ")"),
               subtitle = "Poisson",
               y = "",
               x = "\u03b2"))
      
      if(!is.null(newData)){
        # calculating and backtransforming prediciton and prediction interval
        newData$mean <- exp(apply(pv, 2, function(x)mean(x,na.rm = TRUE)))
        newData$lower <- exp(apply(pv, 2, function(x)quantile(x,p = 0.025, na.rm = TRUE)))
        newData$upper <- exp(apply(pv, 2, function(x)quantile(x,p = 0.975, na.rm = TRUE)))
      }
    }
    
    ## negative binomial ----
    if(any(model@call$family %like% "negative.binomial")){  
      beta_bs <- data.frame(FE = names(fixef(model)),
                            coef = fixef(model),
                            lower = apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T)), #CI lower bound
                            upper = apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T))) #CI upper bound
      
      # plotting model coefficients
      (beta_plot <- beta_bs %>% 
          dplyr::filter(FE != "(Intercept)") %>% 
          ggplot + 
          geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
          # confidence intervals
          geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
          # creating dashed line around 0
          geom_vline(xintercept = 0, lty = "dashed") +
          # adding model coefficient value to plot
          geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                        vjust = -.6, hjust = .3), size = 3) +
          theme(legend.position = "none") +
          labs(title = paste0("fixed effects (iterations = ", n_fit, ")"),
               subtitle = "negative binomial",
               y = "",
               x = "\u03b2"))
      
      if(!is.null(newData)){
        # calculating and backtransforming prediciton and prediction interval
        newData$mean <- exp(apply(pv, 2, function(x)mean(x,na.rm = TRUE)))
        newData$lower <- exp(apply(pv, 2, function(x)quantile(x,p = 0.025, na.rm = TRUE)))
        newData$upper <- exp(apply(pv, 2, function(x)quantile(x,p = 0.975, na.rm = TRUE)))
      }
    }
    
    ## binomial ----
    if(model@call$family == "binomial"){  
      
      beta_bs <- data.frame(FE = names(fixef(model)),
                            coef = fixef(model),
                            lower = apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T)), #CI lower bound
                            upper = apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T))) #CI upper bound
      
      
      # plotting model coefficients
      (beta_plot <- beta_bs %>% 
          dplyr::filter(FE != "(Intercept)") %>% 
          ggplot + 
          geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
          # confidence intervals
          geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
          # creating dashed line around 0
          geom_vline(xintercept = 0, lty = "dashed") +
          # adding model coefficient value to plot
          geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                        vjust = -.6, hjust = .3), size = 3) +
          theme(legend.position = "none") +
          labs(title = paste0("fixed effects (iterations = ", n_fit, ")"),
               subtitle = "binomial",
               y = "",
               x = "\u03b2"))
      
      if(!is.null(newData)){
        # calculating and backtransforming prediciton and prediction interval
        newData$mean <- plogis(apply(pv, 2, function(x)mean(x,na.rm = TRUE)))
        newData$lower <- plogis(apply(pv, 2, function(x)quantile(x,p = 0.025, na.rm = TRUE)))
        newData$upper <- plogis(apply(pv, 2, function(x)quantile(x,p = 0.975, na.rm = TRUE)))
      }
    }
  }
  
  
  # glmmTMB -----------------------------------------------------------------
  ## code for glmmTMB models
  if(inherits(model, "glmmTMB")){
    
    ## zero inflated ----
    # conditional model is the count 
    # zero-inflated model is the binomial
    if(length(fixef(model)$zi) != 0){
      # bootstrapping parameter values from model simulations
      betas <- matrix(NA, nrow = nsim,
                      ncol = length(fixef(model)$cond) + length(fixef(model)$zi))
      
      for(j in 1:nsim){
        sim_data <- data %>% 
          mutate(y = unlist(simulate(model)))  # simulate response variables from the model
        
        sim_model <- update(model, y ~ ., data = sim_data) #rerun the model with the simulated values as the response
        
        if(sim_model$fit$convergence == 0 & sim_model$sdr$pdHess){ 
          # if model convergence looks good then record the parameters
          betas[j,] <- c(fixef(sim_model)$cond, fixef(sim_model)$zi)
          
          # doing model prediction if new data is provided
          if(!is.null(newData)){
            # if there is no random effect
            if(missing(re)){
              pv[j,] <- predict(sim_model, newData, re.form = NA)
              
            } else(pv[j,] <- predict(sim_model, newData, re.form = NA) + rnorm(1, 0, sd = re))
          }
        }
      }
      
      # seeing how many models fit well
      n_fit <- betas %>% 
        na.omit %>% 
        nrow()
      
      ### Poisson ----
      if(model$call$family == "poisson"){  
        beta_bs <- data.frame(FE = c(names(fixef(model)$cond), 
                                     names(fixef(model)$zi)),
                              model = c(rep("conditional", length(fixef(model)$cond)), # column for the model type
                                        rep("zi", length(fixef(model)$zi))),
                              coef = c(fixef(model)$cond, fixef(model)$zi),
                              lower = betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T),      #CI lower bound
                              upper = betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T)) %>%  #CI upper bound 
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
          # confidence intervals
          geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
          # creating dashed line at 0
          geom_vline(xintercept = 0, lty = "dashed") +
          # adding model coefficient value to plot
          geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                        vjust = -.6, hjust = .3), size = 3) +
          theme(legend.position = "none") +
          labs(title = "Conditional model fixed effects",
               subtitle = paste0("Poisson (iterations = ", n_fit, ")"),
               y = "",
               x = "\u03b2")
        
        # zero inflated model
        beta_plot_zi <- beta_bs[[2]] %>% 
          dplyr::filter(FE != "(Intercept)") %>% 
          ggplot + 
          geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
          geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
          geom_vline(xintercept = 0, lty = "dashed") +
          geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                        vjust = -.6, hjust = .3), size = 3) +
          theme(legend.position = "none") +
          labs(title = "Zero inflated model fixed effects",
               subtitle = paste0("Binomial (iterations = ", n_fit, ")"),
               y = "",
               x = "\u03b2")
        
        #patching plots together
        beta_plot <- beta_plot_cond / beta_plot_zi
      }
      
      ### negative binomial ----
      if(model$call$family == "nbinom2"){  
        beta_bs <- data.frame(FE = c(names(fixef(model)$cond), 
                                     names(fixef(model)$zi)),
                              model = c(rep("conditional", length(fixef(model)$cond)), # column for the model type
                                        rep("zi", length(fixef(sim_model)$zi))), 
                              coef = c(fixef(model)$cond, fixef(model)$z),
                              lower = apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T)),      #CI lower bound
                              upper = apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T))) %>%  #CI upper bound 
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
          # confidence intervals
          geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
          # creating dashed line at 0
          geom_vline(xintercept = 0, lty = "dashed") +
          # adding model coefficient value to plot
          geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                        vjust = -.6, hjust = .3), size = 3) +
          theme(legend.position = "none") +
          labs(title = "Conditional model fixed effects",
               subtitle = paste0("Negative binomial (iterations = ", n_fit, ")"),
               y = "",
               x = "\u03b2")
        
        # zero inflated model
        beta_plot_zi <- beta_bs[[2]] %>% 
          dplyr::filter(FE != "(Intercept)") %>% 
          ggplot + 
          geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
          # confidence intervals
          geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
          # creating dashed line at 0
          geom_vline(xintercept = 0, lty = "dashed") +
          # adding model coefficient value to plot
          geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                        vjust = -.6, hjust = .3), size = 3) +
          theme(legend.position = "none") +
          labs(title = "Zero inflated model fixed effects",
               subtitle = paste0("Binomial  (iterations = ", n_fit, ")"),
               y = "",
               x = "\u03b2")
        
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
        beta_bs <- data.frame(FE = names(fixef(model)$cond),
                              coef = fixef(model)$cond,
                              lower = apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T)), #CI lower bound
                              upper = apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T))) #CI upper bound
        
        # plotting model coefficients
        (beta_plot <- beta_bs %>% 
            dplyr::filter(FE != "(Intercept)") %>% 
            ggplot + 
            geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
            # confidence intervals
            geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
            # creating dashed line at 0
            geom_vline(xintercept = 0, lty = "dashed") +
            # adding model coefficient value to plot
            geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                          vjust = -.6, hjust = .3), size = 3) +
            theme(legend.position = "none") +
            labs(title = paste0("fixed effects (iterations = ", n_fit, ")"),
                 subtitle = "Poisson",
                 y = "",
                 x = "\u03b2"))
      }
      
      ### negative binomial ----
      if(model$call$family == "nbinom2"){  
        beta_bs <- data.frame(FE = names(fixef(model)$cond),
                              coef = fixef(model)$cond,
                              lower = apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T)), #CI lower bound
                              upper = apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T))) #CI upper bound
        
        # plotting model coefficients
        (beta_plot <- beta_bs %>% 
            dplyr::filter(FE != "(Intercept)") %>% 
            ggplot + 
            geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
            # confidence intervals
            geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
            # creating dashed line at 0
            geom_vline(xintercept = 0, lty = "dashed") +
            # adding model coefficient value to plot
            geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                          vjust = -.6, hjust = .3), size = 3) +
            theme(legend.position = "none") +
            labs(title = paste0("fixed effects (iterations = ", n_fit, ")"),
                 subtitle = "negative binomial",
                 y = "",
                 x = "\u03b2"))
      }
      
      ### binomial ----
      if(model$call$family == "binomial"){  
        
        beta_bs <- data.frame(FE = names(fixef(model)$cond),
                              coef = fixef(model)$cond,
                              lower = apply(betas, 2, function(x) quantile(x, probs = 0.025, na.rm = T)), #CI lower bound
                              upper = apply(betas, 2, function(x) quantile(x, probs = 0.975, na.rm = T))) #CI upper bound
        
        
        # plotting model coefficients
        (beta_plot <- beta_bs %>% 
            dplyr::filter(FE != "(Intercept)") %>% 
            ggplot + 
            geom_point(aes(x = coef, y = FE), colour = "#00BFC4") +
            # confidence intervals
            geom_segment(aes(x = lower, xend = upper, y = FE, yend = FE), colour = "#00BFC4") +
            # creating dashed line at 0
            geom_vline(xintercept = 0, lty = "dashed") +
            # adding model coefficient value to plot
            geom_text(aes(x = coef, y = FE, label = round(coef, 2),
                          vjust = -.6, hjust = .3), size = 3) +
            theme(legend.position = "none") +
            labs(title = paste0("fixed effects (iterations = ", n_fit, ")"),
                 subtitle = "binomial",
                 y = "",
                 x = "\u03b2"))
      }
    }
  }
  
  
  if(!is.null(newData)){
    return(list(betas, beta_bs, beta_plot, newData))
  } else(  return(list(betas, beta_bs, beta_plot))
)
}

