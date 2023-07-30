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
power_calculator_proportion_robust_chi <- function(p1, p0, alpha = 0.05, N, runs){
  seed_num <- 50
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
  seed_num <- 50
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
                      "Users have the flexibility to select the number of simulation iterations and define their desired probability of Type 1 error and Power.",
                      br(),
                      br(), 
                      "Users can select their preferred proportional equality test for analysis. For small sample sizes, we recommend using Fisher's exact test due to its accuracy in such scenarios. On the other hand, for large sample sizes, the chi-squared test is recommended for its computational efficiency. To balance computation time, we set the simulation iterations for Fisher's exact test relatively smaller.",
                      br(),
                      br(),
                      "The first resulting plot displays the simulation outcomes, providing insights into power analysis and sample size estimation based on predetermined parameters.",
                      br(),
                      "The grey line represents the relationship between power and the number sample size applying the general Z-statistic formula for comparing two proportions of dichotomous variables.",
                      br(),
                      "The light blue line represents the power estimates obtained for each sample size point through simulation (under one random seed). The simulated data follows a binomial distribution for each group, corresponding to the number of subjects.",
                      br(),
                      "The second plot assesses the robustness of the simulation regarding the target power. By using different random seeds, the simulation calculates the average power and corresponding number of subjects. A confidence interval is constructed to examine the variability of power.",
                      br(),
                      "The summary presents the results of power and variance from the robustness test conducted through simulations.",
                      br(),
                      br(),
                      "*Notice:",
                      br(),
                      "Increasing the number of simulation runs may result in longer computation times.",
                      br(),
                      "The robustness test may require some time to complete."
                      ),
             radioButtons("methods",label = "Methods for Testing Equality of Proportions",c("Chi-squared test"="chi", "Fisher's exact test"="fisher"))
           ))),
  sidebarLayout(
  sidebarPanel(
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
        )),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", 
                plotOutput("powerplot"),
                column(12, 
                        wellPanel(htmlOutput("nrequired_z_test"),
                                htmlOutput("nrequired_sim") 
                ))), 
        tabPanel("Robustness", 
               plotOutput("powerplot_robust"),
               column(12, 
                      wellPanel(htmlOutput("nrequired_robust")
               ))),
        tabPanel("Summary",
               fluidRow(
                 p(withMathJax("$$\\text{The minimum sample size required for each group to achieve the target power: (n)}$$")),
                 wellPanel(htmlOutput("target_N")),
                 p(withMathJax("$$\\text{The closest power value to the target after simulating multiple iterations with different random seeds: }(\\hat{p})$$")),
                 wellPanel(htmlOutput("target_robust")),
                 p(withMathJax("$$\\text{Variance of power from randamized simulation under target: }$$")),
                 wellPanel(htmlOutput("var_ran")),
                 p(withMathJax("$$\\text{Variance of power from sample proportion formula under target: }(\\frac{\\hat{p}*(1-\\hat{p})}{n})$$")),
                 wellPanel(htmlOutput("var"))
               )))))))



# Define server logic required to draw a histogram
server <- shinyServer(
  function(input, output) {
    
    # chi_squared
    betas_fun_proportion_sim_chi <- reactive({
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
      betas <- apply(X=Ns,MARGIN = 1,FUN = power_calculator_proportion_simulation_chi, p1=p1, p0=p0,alpha=alpha, runs=runs, seed=seed)
      return(list(betas=betas, seed=seed))
    })
    
    betas_fun_proportion_z_test_chi <- reactive({
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
      betas <- apply(X=Ns,MARGIN = 1,FUN = power_calculator_proportion_z_test, p1=p1, p0=p0,alpha=alpha)
      return(list(betas=betas))
    })
    
    betas_fun_proportion_sim_robust_chi <- reactive({
      p0 <- input$p0_c
      p1 <- input$p1_c
      runs <- input$runs_c
      target <- input$target_c
      maxn <- input$maxn_c
      Ns <- as.matrix(c(seq(1, maxn, 1)))
      results_sim <- betas_fun_proportion_sim_chi()
      nrequired_sim <-Ns[which.max(results_sim$betas>=target)]
      Nr <- as.matrix(c(seq(round(nrequired_sim-0.05*maxn,0),round(nrequired_sim+0.05*maxn,0),1)))
      alpha <- switch(input$alpha_c,
                      "Alpha = 0.01" = 0.01,
                      "Alpha = 0.025" = 0.025, 
                      "Alpha = 0.05" = 0.05, 
                      "Alpha = 0.10" = 0.10)
      betas <- apply(X=Nr,MARGIN = 1,FUN = power_calculator_proportion_robust_chi, p1=p1, p0=p0,alpha=alpha, runs=runs)
      return(list(N=Nr, betas=betas))
    })
    
    # Fisher's exact
    betas_fun_proportion_sim_fisher <- reactive({
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
      betas <- apply(X=Ns,MARGIN = 1,FUN = power_calculator_proportion_simulation_fisher, p1=p1, p0=p0,alpha=alpha, runs=runs, seed=seed)
      return(list(betas=betas, seed=seed))
    })
    
    betas_fun_proportion_z_test_fisher <- reactive({
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
      betas <- apply(X=Ns,MARGIN = 1,FUN = power_calculator_proportion_z_test, p1=p1, p0=p0,alpha=alpha)
      return(list(betas=betas))
    })
    
    betas_fun_proportion_sim_robust_fisher <- reactive({
      p0 <- input$p0_f
      p1 <- input$p1_f
      runs <- input$runs_f
      target <- input$target_f
      maxn <- input$maxn_f
      Ns <- as.matrix(c(seq(1, maxn, 1)))
      results_sim <- betas_fun_proportion_sim_fisher()
      nrequired_sim <-Ns[which.max(results_sim$betas>=target)]
      Nr <- as.matrix(c(seq(round(nrequired_sim-0.05*maxn,0),round(nrequired_sim+0.05*maxn,0),1)))
      alpha <- switch(input$alpha_f,
                      "Alpha = 0.01" = 0.01,
                      "Alpha = 0.025" = 0.025, 
                      "Alpha = 0.05" = 0.05, 
                      "Alpha = 0.10" = 0.10)
      betas <- apply(X=Nr,MARGIN = 1,FUN = power_calculator_proportion_robust_fisher, p1=p1, p0=p0,alpha=alpha, runs=runs)
      return(list(N=Nr, betas=betas))
    })    
    
    
    
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
        results <- betas_fun_proportion_sim_chi()
        plot(NA, ylim=c(0,1), xlim=c(0,maxn*2), main=paste0("Hypothetical Treatment Effect(Difference) = ",abs(round((p1-p0),3))," Percentage Points"),
              ylab="Power (Probability of Statistical Significance)", xlab="Total Number of Subjects")
        lines(Ns*2, results$betas, lwd=4,col="#1B98E059")
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
        results <- betas_fun_proportion_sim_fisher()
        plot(NA, ylim=c(0,1), xlim=c(0,maxn*2), main=paste0("Hypothetical Treatment Effect(Difference) = ",abs(round((p1-p0),3))," Percentage Points"),
             ylab="Power (Probability of Statistical Significance)", xlab="Total Number of Subjects")
        lines(Ns*2, results$betas, lwd=4,col="#1B98E059")
        lines(Ns*2, apply(X=Ns,MARGIN = 1,FUN = power_calculator_proportion_z_test, p1=p1, p0=p0,alpha=alpha), lwd=4, col="grey")
        abline(h=target, col="red", lty=2)
        legend("bottomright", inset=.05, title="approaches",
               c("Fisher's exact test Simulation","Z-test"), fill=c("#1B98E059","grey"), horiz=TRUE)
      }
    })
    
    output$powerplot_robust <- renderPlot({
      if(input$methods == "chi"){
      target <- input$target_c
      results <- betas_fun_proportion_sim_robust_chi()
      # Create DataFrame for Plotting
      DF <- data.frame(N=results$N,Mean=apply(results$betas,2,mean), U=apply(results$betas,2,U), L=apply(results$betas,2,L))
      # ggplot2 LineGraph with Shading Confidence Interval
      ggplot(DF, aes(N*2, Mean)) +                                     
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
        results <- betas_fun_proportion_sim_robust_fisher()
        # Create DataFrame for Plotting
        DF <- data.frame(N=results$N,Mean=apply(results$betas,2,mean), U=apply(results$betas,2,U), L=apply(results$betas,2,L))
        # ggplot2 LineGraph with Shading Confidence Interval
        ggplot(DF, aes(N*2, Mean)) +                                     
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
        results_sim <- betas_fun_proportion_sim_chi()
        nrequired_sim <-Ns[which.max(results_sim$betas>=target)]
        seed <- results_sim$seed
        str1 <- paste0("In order to achieve ", target*100, "% power, you'll need to use a Total sample size of at least ", nrequired_sim*2," under random seed of ", seed, " through simulation.")
      
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
        results_sim <- betas_fun_proportion_sim_fisher()
        nrequired_sim <-Ns[which.max(results_sim$betas>=target)]
        seed <- results_sim$seed
        str1 <- paste0("In order to achieve ", target*100, "% power, you'll need to use a Total sample size of at least ", nrequired_sim*2," under random seed of ", seed, " through simulation.")
        
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
      results_z_test <- betas_fun_proportion_z_test_chi()
      nrequired_z_test <-Ns[which.max(results_z_test$betas>=target)]
      
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
        results_z_test <- betas_fun_proportion_z_test_fisher()
        nrequired_z_test <-Ns[which.max(results_z_test$betas>=target)]
        
        str1 <- paste0("In order to achieve ", target*100, "% power, you'll need to use a Total sample size of at least ", nrequired_z_test*2," by Z test general formula.")
        HTML(str1)
      }
    })
    
    output$nrequired_robust <- renderUI({  
      if(input$methods == "chi"){
      target <- input$target_c
      results <- betas_fun_proportion_sim_robust_chi()
      Nr <- results$N
      Mean <- apply(results$betas,2,mean)
      nrequired_robust <-Nr[which.max(Mean>=target)]
      if(nrequired_robust<target*0.95){
        str1 <- paste0("The iterations are too small, try to enlarge the simulation runs.")
      }
      else{
        str1 <- paste0("In order to achieve ", target*100, "% power, you'll need to use a Average Total sample size of at least ", nrequired_robust*2, " through simulation.")
      }
      
      HTML(str1)
      }
      else if(input$methods == "fisher"){
        target <- input$target_f
        results <- betas_fun_proportion_sim_robust_fisher()
        Nr <- results$N
        Mean <- apply(results$betas,2,mean)
        nrequired_robust <-Nr[which.max(Mean>=target)]
        if(nrequired_robust<target*0.9){
          str1 <- paste0("The iterations are too small or sample size is large, consider using the chi-squared test.")
          
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
      results <- betas_fun_proportion_sim_robust_chi()
      Mean <- apply(results$betas,2,mean)
      Nr <- results$N
      target_N <- Nr[which.max(Mean>=target)]
      str1 <- paste0(target_N)
      HTML(str1)
      }
      else if(input$methods == "fisher"){      
        target <- input$target_f
        results <- betas_fun_proportion_sim_robust_fisher()
        Mean <- apply(results$betas,2,mean)
        Nr <- results$N
        target_N <- Nr[which.max(Mean>=target)]
        str1 <- paste0(target_N)
        HTML(str1)
      }
    })
    
    output$target_robust <- renderUI({  
      if(input$methods == "chi"){      
      target <- input$target_c
      results <- betas_fun_proportion_sim_robust_chi()
      Mean <- apply(results$betas,2,mean)
      target_robust <- Mean[which.max(Mean>=target)]
      str1 <- paste0(target_robust)
      HTML(str1)
      }
      else if(input$methods == "fisher"){      
        target <- input$target_f
        results <- betas_fun_proportion_sim_robust_fisher()
        Mean <- apply(results$betas,2,mean)
        target_robust <- Mean[which.max(Mean>=target)]
        str1 <- paste0(target_robust)
        HTML(str1)
      }
    })
    
    
    output$var_ran <- renderUI({  
      if(input$methods == "chi"){      
      target <- input$target_c
      results <- betas_fun_proportion_sim_robust_chi()
      Mean <- apply(results$betas,2,mean)
      Var <- apply(results$betas,2,var)
      var_ran <- Var[which.max(Mean>=target)]
      str1 <- paste0(round(var_ran,7))
      
      HTML(str1)
      }
      else if(input$methods == "fisher"){      
        target <- input$target_f
        results <- betas_fun_proportion_sim_robust_fisher()
        Mean <- apply(results$betas,2,mean)
        Var <- apply(results$betas,2,var)
        var_ran <- Var[which.max(Mean>=target)]
        str1 <- paste0(round(var_ran,7))
        
        HTML(str1)
      }
    })
    
    output$var <- renderUI({  
      if(input$methods == "chi"){      
      target <- input$target_c
      results <- betas_fun_proportion_sim_robust_chi()
      Nr <- results$N
      Mean <- apply(results$betas,2,mean)
      nrequired_robust <-Nr[which.max(Mean>=target)]
      nrequired_beta <- Mean[which.max(Mean>=target)]
      var_ran <-nrequired_beta*(1-nrequired_beta)/nrequired_robust
      str1 <- paste0(round(var_ran,7))
      
      HTML(str1)
      }
      else if(input$methods == "fisher"){      
        target <- input$target_f
        results <- betas_fun_proportion_sim_robust_fisher()
        Nr <- results$N
        Mean <- apply(results$betas,2,mean)
        nrequired_robust <-Nr[which.max(Mean>=target)]
        nrequired_beta <- Mean[which.max(Mean>=target)]
        var_ran <-nrequired_beta*(1-nrequired_beta)/nrequired_robust
        str1 <- paste0(round(var_ran,7))
        
        HTML(str1)
      }
    })
    
})


# Run the application 
shinyApp(ui = ui, server = server)
