#install.packages('rsconnect')
library(shiny)
library("ggplot2")

# Z test power calculator
power_calculator_proportion_z_test <- function(p1, p0, alpha=0.05, N){
  q0 <- 0.5
  q1 <- 0.5
  p <- p0*q0 + p1*q1
  A <- qnorm(1-alpha/2)*sqrt(p*(1-p)*(1/q0+1/q1))
  B <- sqrt(p0*(1-p0)*(1/q0)+p1*(1-p1)*(1/q1))
  C <- (p0-p1)^2
  beta <- pnorm((sqrt(2*N*C)-A)/B)
  return(round(beta,4))
}

# Simulation calculator - Chi-squared
power_calculator_proportion_simulation_chi <- function(p1, p0, alpha = 0.05, N, runs, seed){
  set.seed(seed)
  pv = replicate(runs, prop.test(c(rbinom(1,N,p0), rbinom(1,N,p1)), c(N,N),cor=F)$p.val)
  beta <- mean(pv <= alpha,  na.rm=TRUE)
  return(beta)
}


# Simulation calculator - Fisher's exact
power_calculator_proportion_simulation_fisher <- function(p1, p0, alpha = 0.05, N, runs, seed){
  set.seed(seed)
  pv = replicate(runs, fisher.test(rbind(c(rbinom(1,N,p0), 460-rbinom(1,N,p0)), c(rbinom(1,N,p1), 460-rbinom(1,N,p1))))$p.val)
  beta <- mean(pv <= alpha,  na.rm=TRUE)
  return(beta)
}

# Test simulation robustness calculator - Chi-squared
seednum <- 100
power_calculator_proportion_robust_chi <- function(p1, p0, alpha = 0.05, N, runs){
  seed_num <- seednum
  random_seed <- matrix(sample(1:10000,seed_num, replace=F))
  beta_matrix <- matrix(0,seed_num,1)
  for(i in 1:seed_num){
    set.seed(random_seed[i,])
    pv = replicate(runs, prop.test(c(rbinom(1,N,p0), rbinom(1,N,p1)), c(N,N),cor=F)$p.val)
    beta <- mean(pv <= alpha,  na.rm=TRUE)
    beta_matrix[i,] <- beta
  }
  return(beta_matrix)
}

# Test simulation robustness calculator - Fisher's exact
power_calculator_proportion_robust_fisher <- function(p1, p0, alpha = 0.05, N, runs){
  seed_num <- seednum
  random_seed <- matrix(sample(1:10000,seed_num, replace=F))
  beta_matrix <- matrix(0,seed_num,1)
  for(i in 1:seed_num){
    set.seed(random_seed[i,])
    pv = replicate(runs, fisher.test(rbind(c(rbinom(1,N,p0), 460-rbinom(1,N,p0)), c(rbinom(1,N,p1), 460-rbinom(1,N,p1))))$p.val)
    beta <- mean(pv <= alpha,  na.rm=TRUE)
    beta_matrix[i,] <- beta
  }
  return(beta_matrix)
}

# Upper confidence interval
U <- function(data){
  ui <- mean(data) + qnorm(0.975)*sd(data)
  return(ui)
}

# Lower confidence interval
L <- function(data){
  li <- mean(data) - qnorm(0.975)*sd(data)
  return(li)
}




# Define UI for application 
ui <- shinyUI(fluidPage(
  headerPanel("Simulation-based validation of sample size calculation for comparing two proportions"),
  fluidRow(
    column(12,
           wellPanel(
             helpText(" This application is specifically designed to assist researchers and practitioners in efficiently estimating the required sample size for comparing proportions between two groups, while considering the desired power and significance level.",  
                      br(), 
                      "In this app, users have the flexibility to choose the number of simulation iterations and define their desired probability of Type 1 error and Power.",
                      br(),
                      br(), 
                      "For analysis, users can select their preferred proportional equality test. For small sample sizes, we recommend using Fisher's exact test due to its accuracy in such scenarios. Alternatively, for large sample sizes, the chi-squared test is recommended for its computational efficiency. To optimize computation time, we have set the simulation iterations for Fisher's exact test relatively smaller.",
                      br(),
                      br(),
                      "The first resulting plot displays the simulation outcomes, providing insights into power analysis and sample size estimation based on predetermined parameters.",
                      br(),
                      "The grey line represents the relationship between power and the number sample size applying the general Z-statistic formula for comparing two proportions of dichotomous variables.",
                      br(),
                      "The light blue line represents the power estimates obtained for each sample size point through simulation (under one random seed). The simulated data follows a binomial distribution for each group, corresponding to the number of subjects.",
                      br(),
                      "The second plot assesses the degree of uncertainty which arises from the randomness of the simulation. Using a range of randomly selected seeds (the number of random seeds labeled as N), the simulation computes the average power and the corresponding number of subjects. A confidence interval is constructed to examine the variability of power. The red line indicates the minimum sample size required to attain the target power level (labeled as P)",
                      br(),
                      "The summary presents the results of power and variance from the uncertainty test conducted through simulations.",
                      br(),
                      br(),
                      "*Notice:",
                      br(),
                      "1. Increasing the number of simulation runs may result in longer computation times.",
                      br(),
                      "2. The uncertainty test may require some time to complete.",
                      br(),
                      "3. After clicking the run button, please allow some time for processing. Kindly refrain from repeatedly pressing the button during this period.",
                      br(),
                      "- The process of generating result plots under the default setting parameters for the Chi-squared test is expected to take around 25 minutes.",
                      br(),
                      "- The process of generating result plots under the default setting parameters for the Fisher's exact test is expected to take around 30 minutes.",
                      br(),
                      "4. Each time you modify the input parameters, please remember to click the reset button."
                      
                      
             )
           ))),
  sidebarLayout(
    sidebarPanel(
      radioButtons("methods",label = "Methods for Testing Equality of Proportions",c("Chi-squared test"="chi", "Fisher's exact test"="fisher")),
      conditionalPanel(condition = "input.methods == 'chi'",
                       numericInput("runs_c", "The Number of Simulation Iterations", value = 1000, min = 10, max = 100000),
                       selectInput("alpha_c", "Significance Level(two-sided)",c("Alpha = 0.01", "Alpha = 0.025", "Alpha = 0.05", "Alpha = 0.10"),selected="Alpha = 0.05"),
                       numericInput("p0_c", "Risk in Group 0 (baseline risk)", value = 0.15, min = 0, max = 1),
                       numericInput("p1_c", "Risk in Group 1 (exposed)", value = 0.09, min = 0, max = 1),
                       sliderInput("target_c", "Power Target", value = .8, min = 0, max = 1),
                       numericInput("maxn_c", "Maximum Number of Subjects(per group)", value = 1000, min = 10, max = 1000000)
      ),
      conditionalPanel(condition = "input.methods == 'fisher'",
                       numericInput("runs_f", "The Number of Simulation Iterations", value = 200, min = 10, max = 100000),
                       selectInput("alpha_f", "Significance Level(two-sided)",c("Alpha = 0.01", "Alpha = 0.025", "Alpha = 0.05", "Alpha = 0.10"),selected="Alpha = 0.05"),
                       numericInput("p0_f", "Risk in Group 0 (baseline risk)", value = 0.15, min = 0, max = 1),
                       numericInput("p1_f", "Risk in Group 1 (exposed)", value = 0.09, min = 0, max = 1),
                       sliderInput("target_f", "Power Target", value = .8, min = 0, max = 1),
                       numericInput("maxn_f", "Maximum Number of Subjects(per group)", value = 1000, min = 10, max = 1000000)
      ),
      actionButton("action","run"),
      actionButton("reset","reset")),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", 
                 plotOutput("powerplot"),
                 column(12, 
                        wellPanel(htmlOutput("nrequired_z_test"),
                                  htmlOutput("nrequired_sim") 
                        ))), 
        tabPanel("Uncertainty", 
                 plotOutput("powerplot_robust"),
                 column(12, 
                        wellPanel(htmlOutput("nrequired_robust")
                        ))),
        tabPanel("Summary",
                 fluidRow(
                   p(withMathJax("$$\\text{The minimum sample size required for each group to achieve the target power:}$$")),
                   wellPanel(htmlOutput("target_N")),
                   p(withMathJax("$$\\text{The closest power value to the target after simulating multiple iterations with different random seeds: }(\\hat{P})$$")),
                   wellPanel(htmlOutput("target_robust")),
                   p(withMathJax("$$\\text{Emperical variance of power from randomized simulation under target: }$$")),
                   wellPanel(htmlOutput("var_ran")),
                   p(withMathJax("$$\\text{Estimated variance of power from sample proportion formula under target: }(\\frac{\\hat{P}*(1-\\hat{P})}{N})$$")),
                   wellPanel(htmlOutput("var"))
                 )))))))



# Define server logic required to draw a histogram
server <- shinyServer(
  function(input, output) {
    
    
    betas <- eventReactive(input$action,{
      
      
      if(input$methods == "chi"){
        p0 <- input$p0_c
        p1 <- input$p1_c
        runs <- input$runs_c
        target <- input$target_c
        maxn <- input$maxn_c
        Ns <- as.matrix(c(seq(1, maxn, 1)))
        seed <- sample(1:10000,1, replace=F)
        alpha <- switch(input$alpha_c,
                        "Alpha = 0.01" = 0.01,
                        "Alpha = 0.025" = 0.025, 
                        "Alpha = 0.05" = 0.05, 
                        "Alpha = 0.10" = 0.10)
        betas_sim <- apply(X=Ns,MARGIN = 1,FUN = power_calculator_proportion_simulation_chi, p1=p1, p0=p0,alpha=alpha, runs=runs, seed=seed)
        betas_z <- apply(X=Ns,MARGIN = 1,FUN = power_calculator_proportion_z_test, p1=p1, p0=p0,alpha=alpha)
        nrequired_z_test <-Ns[which.max(betas_z>=target)]
        Nr <- as.matrix(c(seq(round(nrequired_z_test-0.05*maxn,0),round(nrequired_z_test+0.05*maxn,0),1)))
        betas_robust <- apply(X=Nr,MARGIN = 1,FUN = power_calculator_proportion_robust_chi, p1=p1, p0=p0,alpha=alpha, runs=runs)
        return(list(betas_sim=betas_sim, seed=seed, betas_z=betas_z, Nr=Nr, betas_robust=betas_robust))
      }
      else if(input$methods == "fisher"){
        p0 <- input$p0_f
        p1 <- input$p1_f
        runs <- input$runs_f
        target <- input$target_f
        maxn <- input$maxn_f
        Ns <- as.matrix(c(seq(1, maxn, 1)))
        seed <- sample(1:10000,1, replace=F)
        alpha <- switch(input$alpha_f,
                        "Alpha = 0.01" = 0.01,
                        "Alpha = 0.025" = 0.025, 
                        "Alpha = 0.05" = 0.05, 
                        "Alpha = 0.10" = 0.10)
        betas_sim <- apply(X=Ns,MARGIN = 1,FUN = power_calculator_proportion_simulation_fisher, p1=p1, p0=p0,alpha=alpha, runs=runs, seed=seed)
        betas_z <- apply(X=Ns,MARGIN = 1,FUN = power_calculator_proportion_z_test, p1=p1, p0=p0,alpha=alpha)
        nrequired_z_test <-Ns[which.max(betas_z>=target)]
        Nr <- as.matrix(c(seq(round(nrequired_z_test-0.05*maxn,0),round(nrequired_z_test+0.05*maxn,0),1)))
        betas_robust <- apply(X=Nr,MARGIN = 1,FUN = power_calculator_proportion_robust_fisher, p1=p1, p0=p0,alpha=alpha, runs=runs)
        return(list(betas_sim=betas_sim, seed=seed, betas_z=betas_z, Nr=Nr, betas_robust=betas_robust))
        
      }
    })
    
    
    
    observeEvent(input$action,{
      output$powerplot <- renderPlot({
        if(input$methods == "chi"){
          p0 <- input$p0_c
          p1 <- input$p1_c
          runs <- input$runs_c
          target <- input$target_c
          maxn <- input$maxn_c
          Ns <- as.matrix(c(seq(1, maxn, 1)))
          alpha <- switch(input$alpha_c,
                          "Alpha = 0.01" = 0.01, 
                          "Alpha = 0.025" = 0.025, 
                          "Alpha = 0.05" = 0.05, 
                          "Alpha = 0.10" = 0.10)
          results <- betas()
          plot(NA, ylim=c(0,1), xlim=c(0,maxn*2), main=paste0("Hypothetical Treatment Effect(Difference) = ",abs(round((p1-p0),3))," Percentage Points"),
               ylab="Power (Probability of Statistical Significance)", xlab="Total Number of Subjects")
          lines(Ns*2, results$betas_sim, lwd=4,col="#1B98E059")
          lines(Ns*2, apply(X=Ns,MARGIN = 1,FUN = power_calculator_proportion_z_test, p1=p1, p0=p0,alpha=alpha), lwd=4, col="grey")
          abline(h=target, col="red", lty=2)
          legend("bottomright", inset=.03, title="approaches",
                 c("Chi-squared test Simulation","Z-test"), fill=c("#1B98E059","grey"), horiz=TRUE)
        }
        else if(input$methods == "fisher"){
          p0 <- input$p0_f
          p1 <- input$p1_f
          runs <- input$runs_f
          target <- input$target_f
          maxn <- input$maxn_f
          Ns <- as.matrix(c(seq(1, maxn, 1)))
          alpha <- switch(input$alpha_f,
                          "Alpha = 0.01" = 0.01, 
                          "Alpha = 0.025" = 0.025, 
                          "Alpha = 0.05" = 0.05, 
                          "Alpha = 0.10" = 0.10)
          results <- betas()
          plot(NA, ylim=c(0,1), xlim=c(0,maxn*2), main=paste0("Hypothetical Treatment Effect(Difference) = ",abs(round((p1-p0),3))," Percentage Points"),
               ylab="Power (Probability of Statistical Significance)", xlab="Total Number of Subjects")
          lines(Ns*2, results$betas_sim, lwd=4,col="#1B98E059")
          lines(Ns*2, apply(X=Ns,MARGIN = 1,FUN = power_calculator_proportion_z_test, p1=p1, p0=p0,alpha=alpha), lwd=4, col="grey")
          abline(h=target, col="red", lty=2)
          legend("bottomright", inset=.05, title="approaches",
                 c("Fisher's exact test Simulation","Z-test"), fill=c("#1B98E059","grey"), horiz=TRUE)
        }
      })
      
      output$powerplot_robust <- renderPlot({
        if(input$methods == "chi"){
          target <- input$target_c
          results <- betas()
          # Create DataFrame for Plotting
          DF <- data.frame(Nr=results$Nr,Mean=apply(results$betas_robust,2,mean), U=apply(results$betas_robust,2,U), L=apply(results$betas_robust,2,L))
          # ggplot2 LineGraph with Shading Confidence Interval
          ggplot(DF, aes(Nr*2, Mean)) +                                     
            geom_line(color = "#1B98E059", size = 1) +
            geom_ribbon(aes(ymin=L, ymax=U), alpha=0.1, fill = "#6897bb", 
                        color = "black", linetype = "dotted") +
            geom_hline(yintercept=target, color="red") +
            labs(x="Total Number of Subjects", y="Power (Probability of Statistical Significance)",
                 title = "Confidence Interval of power under different random seeds iterations") +
            theme_bw()
        }
        else if(input$methods == "fisher"){
          target <- input$target_f
          results <- betas()
          # Create DataFrame for Plotting
          DF <- data.frame(Nr=results$Nr,Mean=apply(results$betas_robust,2,mean), U=apply(results$betas_robust,2,U), L=apply(results$betas_robust,2,L))
          # ggplot2 LineGraph with Shading Confidence Interval
          ggplot(DF, aes(Nr*2, Mean)) +                                     
            geom_line(color = "#1B98E059", size = 1) +
            geom_ribbon(aes(ymin=L, ymax=U), alpha=0.1, fill = "#6897bb", 
                        color = "black", linetype = "dotted") +
            geom_hline(yintercept=target, color="red") +
            labs(x="Total Number of Subjects", y="Power (Probability of Statistical Significance)",
                 title = "Confidence Interval of power under different random seeds iterations") +
            theme_bw()
        }
        
      })
      
      
      
      
      output$nrequired_sim <- renderUI({  
        if(input$methods == "chi"){
          p0 <- input$p0_c
          p1 <- input$p1_c
          maxn <- input$maxn_c
          target <- input$target_c
          alpha <- switch(input$alpha_c,
                          "Alpha = 0.01" = 0.01, 
                          "Alpha = 0.025" = 0.025, 
                          "Alpha = 0.05" = 0.05, 
                          "Alpha = 0.10" = 0.10)
          Ns <- as.matrix(c(seq(1, maxn, 1)))
          results<- betas()
          nrequired_sim <-Ns[which.max(results$betas_sim>=target)]
          prequired <- results$betas_sim[which.max(results$betas_sim>=target)]
          seed <- results$seed
          if(prequired<target*0.95){
            str1 <- paste0("The Maximum Number of Subjects may be too small, try to enlarge it.")
          }
          else{
            str1 <- paste0("In order to achieve ", target*100, "% power, you'll need to use a Total sample size of at least ", nrequired_sim*2," under random seed of ", seed, " through simulation.")
          }
          
          HTML(str1)
        }
        else if(input$methods == "fisher"){
          p0 <- input$p0_f
          p1 <- input$p1_f
          maxn <- input$maxn_f
          target <- input$target_f
          alpha <- switch(input$alpha_f,
                          "Alpha = 0.01" = 0.01, 
                          "Alpha = 0.025" = 0.025, 
                          "Alpha = 0.05" = 0.05, 
                          "Alpha = 0.10" = 0.10)
          Ns <- as.matrix(c(seq(1, maxn, 1)))
          results <- betas()
          nrequired_sim <-Ns[which.max(results$betas_sim>=target)]
          prequired <- results$betas_sim[which.max(results$betas_sim>=target)]
          seed <- results$seed
          if(prequired<target*0.95){
            str1 <- paste0("The Maximum Number of Subjects may be too small, try to enlarge it.")
          }
          else{
            str1 <- paste0("In order to achieve ", target*100, "% power, you'll need to use a Total sample size of at least ", nrequired_sim*2," under random seed of ", seed, " through simulation.")
          }
          HTML(str1)
        }
        
        
      })
      output$nrequired_z_test<- renderUI({ 
        if(input$methods == "chi"){
          p0 <- input$p0_c
          p1 <- input$p1_c
          maxn <- input$maxn_c
          target <- input$target_c
          alpha <- switch(input$alpha_c,
                          "Alpha = 0.01" = 0.01, 
                          "Alpha = 0.025" = 0.025, 
                          "Alpha = 0.05" = 0.05, 
                          "Alpha = 0.10" = 0.10)
          Ns <- as.matrix(c(seq(1, maxn, 1)))
          results <- betas()
          nrequired_z_test <-Ns[which.max(results$betas_z>=target)]
          
          str1 <- paste0("In order to achieve ", target*100, "% power, you'll need to use a Total sample size of at least ", nrequired_z_test*2," by Z test general formula.")
          HTML(str1)
        }
        else if(input$methods == "fisher"){
          p0 <- input$p0_f
          p1 <- input$p1_f
          maxn <- input$maxn_f
          target <- input$target_f
          alpha <- switch(input$alpha_f,
                          "Alpha = 0.01" = 0.01, 
                          "Alpha = 0.025" = 0.025, 
                          "Alpha = 0.05" = 0.05, 
                          "Alpha = 0.10" = 0.10)
          Ns <- as.matrix(c(seq(1, maxn, 1)))
          results <- betas()
          nrequired_z_test <-Ns[which.max(results$betas_z>=target)]
          
          str1 <- paste0("In order to achieve ", target*100, "% power, you'll need to use a Total sample size of at least ", nrequired_z_test*2," by Z test general formula.")
          HTML(str1)
        }
      })
      
      output$nrequired_robust <- renderUI({  
        if(input$methods == "chi"){
          target <- input$target_c
          results <- betas()
          Nr <- results$Nr
          Mean <- apply(results$betas_robust,2,mean)
          nrequired_robust <-Nr[which.max(Mean>=target)]
          prequired_robust <- Mean[which.max(Mean>=target)]
          if(prequired_robust<target*0.95){
            str1 <- paste0("The iterations may be too small, try to enlarge the simulation runs.")
          }
          else{
            str1 <- paste0("In order to achieve ", target*100, "% power, you'll need to use a Average Total sample size of at least ", nrequired_robust*2, " through simulation.")
          }
          
          HTML(str1)
        }
        else if(input$methods == "fisher"){
          target <- input$target_f
          results <- betas()
          Nr <- results$Nr
          Mean <- apply(results$betas_robust,2,mean)
          nrequired_robust <-Nr[which.max(Mean>=target)]
          prequired_robust <- Mean[which.max(Mean>=target)]
          if(prequired_robust<target*0.9){
            str1 <- paste0("The iterations may be too small or sample size is large, consider using the chi-squared test.")
            
          }
          else{
            str1 <- paste0("In order to achieve ", target*100, "% power, you'll need to use a Average Total sample size of at least ", nrequired_robust*2, " through simulation.")
          }
          HTML(str1)
        }
      })
      
      output$target_N <- renderUI({  
        if(input$methods == "chi"){      
          target <- input$target_c
          results <- betas()
          Mean <- apply(results$betas_robust,2,mean)
          Nr <- results$Nr
          target_N <- Nr[which.max(Mean>=target)]
          str1 <- paste0(target_N)
          HTML(str1)
        }
        else if(input$methods == "fisher"){      
          target <- input$target_f
          results <- betas()
          Mean <- apply(results$betas_robust,2,mean)
          Nr <- results$Nr
          target_N <- Nr[which.max(Mean>=target)]
          str1 <- paste0(target_N)
          HTML(str1)
        }
      })
      
      output$target_robust <- renderUI({  
        if(input$methods == "chi"){      
          target <- input$target_c
          results <- betas()
          Mean <- apply(results$betas_robust,2,mean)
          target_robust <- Mean[which.max(Mean>=target)]
          str1 <- paste0(target_robust)
          HTML(str1)
        }
        else if(input$methods == "fisher"){      
          target <- input$target_f
          results <- betas()
          Mean <- apply(results$betas_robust,2,mean)
          target_robust <- Mean[which.max(Mean>=target)]
          str1 <- paste0(target_robust)
          HTML(str1)
        }
      })
      
      
      output$var_ran <- renderUI({  
        if(input$methods == "chi"){      
          target <- input$target_c
          results <- betas()
          Mean <- apply(results$betas_robust,2,mean)
          Var <- apply(results$betas_robust,2,var)
          var_ran <- Var[which.max(Mean>=target)]
          str1 <- paste0(round(var_ran,7))
          
          HTML(str1)
        }
        else if(input$methods == "fisher"){      
          target <- input$target_f
          results <- betas()
          Mean <- apply(results$betas_robust,2,mean)
          Var <- apply(results$betas_robust,2,var)
          var_ran <- Var[which.max(Mean>=target)]
          str1 <- paste0(round(var_ran,7))
          
          HTML(str1)
        }
      })
      
      output$var <- renderUI({  
        if(input$methods == "chi"){      
          target <- input$target_c
          results <- betas()
          Nr <- results$Nr
          Mean <- apply(results$betas_robust,2,mean)
          nrequired_robust <-Nr[which.max(Mean>=target)]
          nrequired_beta <- Mean[which.max(Mean>=target)]
          var_ran <-nrequired_beta*(1-nrequired_beta)/seednum
          str1 <- paste0(round(var_ran,7))
          
          HTML(str1)
        }
        else if(input$methods == "fisher"){      
          target <- input$target_f
          results <- betas()
          Nr <- results$Nr
          Mean <- apply(results$betas_robust,2,mean)
          nrequired_robust <-Nr[which.max(Mean>=target)]
          nrequired_beta <- Mean[which.max(Mean>=target)]
          var_ran <-nrequired_beta*(1-nrequired_beta)/seednum
          str1 <- paste0(round(var_ran,7))
          
          HTML(str1)
        }
      })
    })
    
    observeEvent(input$reset,{
      output$powerplot <- renderUI({  
      })
      output$powerplot_robust <- renderUI({  
      })
      output$nrequired_sim <- renderUI({  
      })
      output$nrequired_z_test <- renderUI({  
      })
      output$nrequired_robust <- renderUI({  
      })
      output$target_N <- renderUI({  
      })
      output$target_robust <- renderUI({  
      })
      output$var_ran <- renderUI({  
      })
      output$var <- renderUI({  
      })
      
    })
    
  })


# Run the application 
shinyApp(ui = ui, server = server)
