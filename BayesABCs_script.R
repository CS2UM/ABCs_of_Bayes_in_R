#######################################################
###### THE ABCS of Bayesian Estimation in R 
###### Eric Dunford 
###### University of Maryland, College Park
###### Department of Government and Politics
###### edunford@umd.edu
#######################################################


# PART 1: GETTING AN INTUITION FOR HOW THINGS WORK ------------------------------------------------




# Bayesian Analysis: Priors, Likelihoods, and Conjugacy


      # Sample Let's generate a sample of some Binary outcome y
      n = 10
      voted = rbinom(n,1,.4) # True data generating process
      table(voted)
      mean(voted)
  
      # For a prior, we will using a Beta Distribution  (which is conjugate to the binomial distribution)
      a = 3
      b = 2
      prior = rbeta(1000,a,b)
      plot(density(prior))
  
      # Combine distributions: Posterior = prior * likelihood
      # Posterior = Beta(y+a,n-y+b)
      y = sum(voted)
      posterior = rbeta(1000,shape1 = y+a,shape2 = n-y+b)
      mean(posterior)
  
      # Visualization
      curve(dbinom(y,n,x),col="grey",lwd=3,ylim=c(0,.5),ylab="f(x)")
      curve(dbeta(x,a,b)*(1/n),col="orange",lwd=3,add=T)
      curve(dbeta(x,y+a,n-y+b)*(1/n),col="steelblue",lwd=3,add=T)
      legend("topright",c("Data","Prior","Posterior"),lty=1,
             col=c("grey","orange",'steelblue'),box.lwd = .01)
      
      
      
      
# Markov Chain Monte Carlo 
      
      # MARKOV CHAINS
      # A Markov chain wanders around a state space remembering only where it 
      # has been in the last period --- a mathematical "Vegas Rule" --- 
      # Essentially, suppose that there is some target distribution that we’d 
      # like to sample from, but we cannot just draw independent samples from. 
      # Instead we can jump around the space and explore it. Each job is 
      # "proposed" and then some threshold criteria determines whether we
      # actually make the jump or not. 
      
      
  # Monte Carlo INTEGRATION
      
      # Let's examine a simple quadratic function with a [0,1] space
      f <- function(x) x^2 
      curve(f,0,1,lwd=5)
      
      # Now let's sample 10,000 randomly drawn points from a uniform distribution
      set.seed(123)
      n <- 10000
      
      xx <- runif(n) # Random x values
      yy <- runif(n) # Random y values
      
      points(xx, yy, col = "lightgrey", pch = 16) # Plot the values and function
      curve(f, 0, 1, n = 100, col = "red", lwd = 5, add = TRUE)
      
      
      # Now, let's specify a function to compute wheter a point is under the curve
      under_the_curve <- function(x, y) y < x^2 
      z <- under_the_curve(xx, yy)
      
      points(xx[z], yy[z], col = "steelblue", pch = 20)
      
      # Calculate the integral
      sum(z)/n # Monte Carlo integral
      integrate(f,0,1) # Compare to the analytic Intergral
      
      # Not bad! We can use computer to approximate the integrals of some tricky
      # probability statements!
      
      
      
      

      
# METROPOLIS-HASTINGS ALGORITHM (EXAMPLE)
      
      # Say we have some distribution that we want to sample from.
      our_dist <- function(x) .3*dnorm(x, -2, .5) + (1-.3)*dnorm(x, 1.75, 1.5)
      curve(our_dist(x), col="orange", -10, 10, n=301, las=1,lwd=4)
      
      
      # Here we want to define a really simple algorithm to "propose" a random 
      # step. For this, let's draw from a normal distribution with a standard 
      # deviation of 4. Note that the wider the variation, the larger the
      # potential proposal. 
      
      proposal <- function(x) rnorm(1, x, 4)
      proposal(1)
      
      
      # Next, let's craft a function that offers a proposition for a step, and 
      # then only takes that step if its sufficiently close to the previous
      # step (this is where the Markov part comes in)
      
      step <- function(x, target_func, proposal) {
        ## Pick new point
        pick <- proposal(x)
        where_should_I_Jump = target_func(pick)
        where_am_I_currently = target_func(x)
        
        
        ## Acceptance probability:
        should_I_jump <- min(1, where_should_I_Jump/where_am_I_currently)
        
        if(is.nan(should_I_jump)){return(0)} # if NaN, return a 0.
        
        ## Accept new point with some acceptance probability :
        if (runif(1) < should_I_jump){ 
          x <- pick
        }
        ## Returning the point:
        x
      }
      
      # Now, let's iterate through, proposing a step and only taking that step 
      # if it's close to our target distribution. This keeps us where we want to
      # be.
      explore <- function(x, target_func, proposal, nsteps) {
        res <- matrix(NA, nsteps, length(x))
        for (i in seq_len(nsteps)){
          x <- step(x, target_func, proposal) 
          res[i,] <- x # Save the results
        }
        return(res)
      }
      
      
      
      output <- explore(1, # Starting Value 
                  target_func = our_dist, # Target Function (of the distribution we want to explore)
                  proposal, # Our proposal function
                  nsteps = 1000) # Number of steps to take
      
      
      # Plot the time series of each step
      layout(matrix(c(1, 2), 1, 2), widths=c(4, 1))
      par(mar=c(4.1, .5, .5, .5), oma=c(0, 4.1, 0, 0))
      plot(output, type="s", xpd=NA, ylab="Parameter", 
           xlab="Sample", las=1,lwd=2,col="darkgrey",ylim=c(-10,10))
      usr <- par("usr")
      xx <- seq(usr[3], usr[4], length=301)
      plot(our_dist(xx), xx, type="l", yaxs="i", axes=FALSE, 
           xlab="",col="orange",lwd=3)
      
      
      
      # How good did we do in recovering the distribution?
      dev.off()
      hist(output, 100, freq=FALSE, main="", ylim=c(0, .4), las=1,
           xlab="x", ylab="Probability density",col="darkgrey",border="white")
      z <- integrate(our_dist, -Inf, Inf)$value
      curve(our_dist(x) / z, add=TRUE, col="orange", n=200,lwd=3)
      # Not bad considering we only took 1000 steps
      
      
      
      # Now, let's do the same thing but change the variation on the proposal
      # function. 
      out.large_var <- explore(-20, our_dist, function(x) rnorm(1, x,  50), 1000)
      out.small_var <- explore(-20, our_dist, function(x) rnorm(1, x,  .5), 1000)
      
      
      # Look at the time series
      layout(matrix(c(1, 2), 1, 2), widths=c(4, 1))
      par(mar=c(4.1, .5, .5, .5), oma=c(0, 4.1, 0, 0))
      plot(output, type="s", xpd=NA, ylab="Parameter", 
           xlab="Sample", las=1,lwd=2,col="darkgrey",ylim=c(-10,10))
      lines(out.large_var, col="red",lwd=2)
      lines(out.small_var, col="blue",lwd=2)
      plot(our_dist(xx), xx, type="l", yaxs="i", axes=FALSE, 
           xlab="",col="orange",lwd=3)
      
      # Look at the autocorrelation between steps
      par(mfrow=c(2, 3), mar=c(4, 2, 3.5, .5))
      acf(out.small_var, las=1, main="Small steps")
      acf(output, las=1, main="Intermediate steps")
      acf(out.large_var, las=1, main="Large steps")
      pacf(out.small_var, las=1, main="Small steps")
      pacf(output, las=1, main="Intermediate steps")
      pacf(out.large_var, las=1, main="Large steps")
      dev.off()
          
          # If the proposals are too limited, then the algorithm gets stuck; if
          # they're too generous, then it's all over the place.
      
      
      
      # It takes time to explore
      output <- explore(1,target_func = our_dist, proposal, nsteps = 10000) 
      
      par(mfrow=c(2,2))
      plot(output,typ="l",col="grey")
      hist(output, 50, freq=FALSE, main="",ylim=c(0, .3), las=1,xlab="x", ylab="Probability density",col="darkgrey",border="white")
      z <- integrate(our_dist, -Inf, Inf)$value
      curve(our_dist(x) / z, add=TRUE, col="orange", n=200,lwd=3)
      acf(output, las=1,main="")
      pacf(output, las=1,main="")
      
      dev.off()
      
      # Burn-in period
      
          # sometimes it takes a while for a walk to become stationary. Think back
          # to time series --- we only want to draw conclusions from stationary
          # draws.
          
          output2 <- explore(-50, our_dist, function(x) rnorm(1, x,  10), 1000)
          plot(output2[1:100],type="l",lwd=2)
          
          # We want to drop these initial values so that we are only dealing
          # with stationary draws.
          
          par(mfrow=c(1,2))
          plot(output2,type="l",lwd=2)
          plot(output2[-(1:100)],type="l",lwd=2)
          dev.off()
      
          
          

# GIBBS SAMPLER (EXAMPLE)

    # Basic Idea: get a marginal distribution for each parameter by iteratively 
    # conditioning on interim values of the others and continuing this until 
    # samples from this process empirically approximate the desired marginal 
    # distribution. Actually, the Gibbs Sampler is a special case of the MH
    # algorithm.
          
    # In order to use a Gibbs Sampler you have to be able to break up a 
    # multivariate distributions into Univariate distributions (i.e. "marginal" 
    # distributions). This is why conjugacy is soooo important. You can sample 
    # from a KNOWN distribution. Thus, the Gibbs sampler will always win out 
    # when compared to MH. MH has to find where the distribution is, Gibbs
    # samples from it directly. 
          
  
    # Now with a "home made" Gibbs Sampler
    gibbs <- function(start,n, rho){
      
      mat <- matrix(ncol = 2, nrow = n)
      x <- y <- start
      
      mat[1, ] <- c(x, y) #X and Y both == 0
      for (i in 2:n) {
        
        x <- rnorm(1, rho * y, sqrt(1 - rho^2))
        #Drawing from the distribution of x given y
        
        y <- rnorm(1, rho * x, sqrt(1 - rho^2))
        #Drawing from the distribution of y given x
        
        mat[i, ] <- c(x, y) #This matrix serves a memory
        
        #When we iterate back through, we will still have the last value
        #of y which will be involved in the next estimation of x, and so
        #on.
      }
      return(mat)
    }
          
    
    output = gibbs(-10,1000,.5)
    plot(output,xlim=c(-10,10),ylim=c(-10,10),col="grey40")
    x = output[,1];y = output[,2]
    s <-seq(length(x)-1)
    segments(x[s], y[s], x[s+1], y[s+1], col='firebrick')
    
    
    bi.dens<-MASS::kde2d(x,y, n=50)
    plot(output,xlim=c(-10,10),ylim=c(-10,10),col="lightgrey",pch=16)
    contour(bi.dens,xlim=c(-10,10),ylim=c(-10,10),add=T,lwd=3,col="darkorange")
    persp(bi.dens,col="black",border="grey", phi = 20, theta = 50)
      
          
    
    
    
    
    
# PART 2: Using JAGS to Estimate Bayesian Models in R --------------------------------------   
    
    require(R2jags)
    require(mcmcplots)
    require(coda)
    require(hdrcde)
    require(ggmcmc)
    require(ggplot2)
    require(dplyr)
    
# The Basics
    
    
    # Simulate Data
    N <-  150
    x <- rnorm(N,0,1)
    y <- 1 + 3*x + rnorm(N,0,1)
    D <-  data.frame(y,x)
    
    # View Relationship
    qplot(y,x,data=D) + theme_minimal() + geom_smooth()
    
    
    summary(lm(y~x,data=D)) # Frequentist model

    
    # Setting up the data (as a list)
    jags.data <- as.list(D)
    jags.data$N <- N
    
    
    
    # Setting up the Jags Model
    model <- function(){
      
      # Model
      for(i in 1:N){
        y[i] ~ dnorm(mu[i],tau)
        mu[i] <- alpha + beta * x[i]
      }
      
      # Priors
      beta ~ dnorm(0,.001) # variance expressed in terms of "precision" (1/sqrt(.01^2) == 100)
      alpha ~ dnorm(0,.001)
      tau ~ dgamma(1,.1) # hyper-prior
    }
    
    
    
    # Specify the parameters we want to be returned
    params <- c("beta","alpha","tau")
    
    
    fit1 <- jags(data=jags.data,n.iter = 20000,parameters.to.save = params,n.chains = 2,
                 n.burnin = 1000,n.thin = 10,model.file = model)
    
    
    # Examining the output
    
    str(fit1) # a lot going on in here
    
    fit1$BUGSoutput$summary
    
    summary(as.mcmc(fit1)) # as.mcmc() from coda package
             
            # Naive SE = post SD/sqrt(N) -- captures "simulation error"
    
    
    head(fit1$BUGSoutput$sims.matrix)
    dim(fit1$BUGSoutput$sims.matrix)
    
    mean(fit1$BUGSoutput$sims.matrix[,'beta'])
    median(fit1$BUGSoutput$sims.matrix[,'beta'])
    
    
    
    # Visualizing the posterior (the "eyeball test" of convergence)
    par(mfrow=c(2,3))
    hist(fit1$BUGSoutput$sims.matrix[,'beta'],breaks=50,
         col="coral",border="white",main="beta")
    plot(fit1$BUGSoutput$sims.matrix[,'beta'],type="l",col="coral")
    acf(fit1$BUGSoutput$sims.matrix[,'beta'])
    hist(fit1$BUGSoutput$sims.matrix[,'alpha'],breaks=50,
         col="steelblue",border="white",main="alpha")
    plot(fit1$BUGSoutput$sims.matrix[,'alpha'],type="l",col="steelblue")
    acf(fit1$BUGSoutput$sims.matrix[,'alpha'])
    
    dev.off()
    
    
    # Examining othe path of each chain
    beta_path <- fit1$BUGSoutput$sims.array[,,'beta']
    
    head(beta_path)
    
    par(mfrow=c(2,1))
    plot(beta_path[,1],type="l",col="red",ylim=c(2.5,3.5)) # chain 1
    lines(beta_path[,2],type="l",col="blue") # chain 2
    legend("top",c('chain 1','chain 2'),lty=1,col=c('red','blue'),box.lwd = .01)
    
    plot(density(beta_path[,1]),lwd=4,col="red",main="")
    lines(density(beta_path[,2]),lwd=4,col="blue")
    
    
    
# Diagnosis of Convergence --------------------------------------------------------
    
    # From Markov chains theory, one can expect that the chains will eventually 
    # converge to the stationary distribution, but there is no guarantee that
    # the chain converged given the number iterations we specified. 
    
    # How do we know whether the chain has actually converged? We can never be 
    # sure, but there are several tests we can do, both visual and statistical,
    # to see if the chain appears to be converged.
    
    
    # EYE BALL TEST --- Actually look at the posterior 
    # (hist/timeseries/autocorrelation ... see above). Do things look like they 
    # converged (i.e. does the distributions across the two posterior chains
    # find the same place relatively speaking)?
    

    # GEWEKE test ---- put simply, a geweke diag test looks if the the mean in 
    # the first part of the timeseries is statistically different from the mean 
    # in the second part of the time series. If the p-value is significant, then
    # you have a problem. It means our model hasn't converged!
    
        ?geweke.diag()
    
        geweke.diag(fit1$BUGSoutput$sims.matrix[,c('alpha','beta')])
        # Not statistically significant! Which is a great sign
        
        
        
        
    # GELMAN and RUBIN Diagnostic -- for this test, we need to one more than 
    # one chain (with different starting values), we discard some n draws 
    # (burn-in) from all our chains, then we calulate the within-chain and 
    # between-chain variance; from there we calculate the est. var. of the 
    # parameter as a weighted sum of the w/in and between chain var, and
    # finally calculate a scal reduction factor. 
       
         
        gelman.diag(mcmc.list(c1 = as.mcmc(beta_path[,1]),c2 =as.mcmc(beta_path[,2])))
        
        # The results give us the median potential scale reduction factor and its 97.5% quantile
        
        
        gelman.plot(mcmc.list(c1 = as.mcmc(beta_path[,1]),c2 =as.mcmc(beta_path[,2])))
        # potential scale reduction factor changes through the iterations
        
        
        
    # HEIDELBERG and WELCH Diagnostic -- 
        
        # Part 1: 
            # 1. Generate a chain of N iterations and define an alpha level.
            # 2. Calculate the test statistic on the whole chain. Accept or
            # reject null hypothesis that the chain is from a stationary
            # distribution.
            # 3. If null hypothesis is rejected, discard the first 10% of the
            # chain. Calculate the test statistic and accept or reject null.
            # 4. If null hypothesis is rejected, discard the next 10% and
            # calculate the test statistic.
            # 5. Repeat until null hypothesis is accepted or 50% of the chain is
            # discarded. If test still rejects null hypothesis, then the chain
            # fails the test and needs to be run longer.
        
        # Part 2: If the chain passes the first part of the diagnostic, then it
        # takes the part of the chain not discarded from the first part to test
        # the second part. The halfwidth test calculates half the width of the
        # (1 − alpha)% credible interval around the mean. If the ratio of the
        # halfwidth and the mean is lower than some epsilon, then the chain passes the
        # test. Otherwise, we need to run it  longer.
        
        heidel.diag(mcmc.list(c1 = as.mcmc(beta_path[,1]),c2 =as.mcmc(beta_path[,2])))
        
        
    
    # PERSONALLY.... I'm an eye ball and Geweke kind of guy...
        
        
        
        
      
        
        
        
        
        
# [APPLIED EXAMPLE 1] Years of Schooling and Protectionist Preferences ----------------------------------------
   
    
  # Data:
    # The dataset is a cleaned up and modified version of the 1996 American
    # National Election studies as used in Hainmueller and Hiscox (2006). This
    # study investigated the determinants of individuals’ support for free trade
    # or protectionist policies. The outcome variable is a dummy 
    # protectionist—coded as 1 if a respondent expressed a preference for more 
    # protectionist policies, and coded as 0 if a respondent favored free 
    # trade.
    
    # The explanatory variables are (see Table A2 in Hainmueller and
    # Hiscox (2006) for more details):
    
          # - age: the incumbent’s age.
          # - female: a binary indicator for female respondents.
          # - TUmember: a binary indicator for trade union members.
          # - partyid: the respondent’s party identification: coded from 0 “strong
              # Democrat” to 6 “strong Republican”.
          # - ideology: the respondent’s ideology: coded 0 if conservative, 1 if
              # moderate, and 2 if liberal.
          # - schooling: years of full-time education completed.
    
    proc_df <- read.csv("hh_2006.csv")
    
    
    # Examine Data
      head(proc_df)
      summary(proc_df)
      
      # Look at the distribution of each variable
      ggplot(data = reshape2::melt(proc_df), aes(x = factor(value))) + 
        geom_histogram(position = "identity",stat = "count") + 
        facet_wrap(~ variable, scales = "free", as.table = TRUE, nrow = 3) + xlab("") + theme_classic()
    
    
    # MY THEORY: Those who have more schooling are less likely to support
    # protectionist policies.
      
      
      #Frequentist model - as reference
      freq.m <- glm(protectionist~age+female+TUmember+partyid+ideology+schooling,
                    data=proc_df,family=binomial(link="logit"))
      summary(freq.m) 
      
      
      
    # Organize the data as a list (which is how R2Jags reads it)
    jags.data <- as.list(proc_df)
    jags.data$N <- nrow(proc_df)
    
  
    # The Bayesian Model
    HHmodel <- function(){
      
      for(i in 1:N){
        
        protectionist[i] ~ dbern(p[i])  #Distributed Bernoulli
        
        logit(p[i]) <- mu[i]    # Logit link function
        
        # Latent Linear Combination
        mu[i] <- alpha + b[1] * age[i] + b[2] * female[i] + b[3] * TUmember[i] + 
          b[4] * partyid[i] + b[5] * ideology[i] + b[6] * schooling[i]
        
        # Estimates for Model fit
        llh[i] <- log(p[i]) * protectionist[i] + log(1 - p[i]) * (1 - protectionist[i]) 
        
      }
      sumllh <- sum(llh[])
      
      
      #Priors (uninformative)
      alpha ~ dnorm(0,.01)
      
      for(j in 1:6){
        b[j] ~ dnorm(0,.01) 
      }
      
    }
    
    # Parameters to retain
    params <- c("b","alpha","sumllh")
    
    HHfit <- jags(data = jags.data,parameters.to.save = params, n.chains = 2, n.iter = 10000,
                n.burnin = 1000, model.file = HHmodel)
    
    
    # Everything we did before but using a package this time
    mcmcplot(HHfit)
    
    geweke.diag(HHfit) # nothing stands out
    
    # Looks like convergence!
    
    
    # Summary Stats
    summary(as.mcmc(HHfit))
    
    
    
    # another way to visualize the posterior is using the hdrcde package, which
    # shows us the highest density region of the posterior
    hdr.boxplot(HHfit$BUGSoutput$sims.matrix[,'b[6]'])
    hdr.den(HHfit$BUGSoutput$sims.matrix[,'b[6]'])
    
    
    
    
    
    # Let's look at the probabiliy of supporting protectionist policies given a 
    # respondents level of education. Here I'm going to source in some old code 
    # that I put together to look at the predicted probability of Bayesian
    # Limited Dependent Variables. We can source this in directly from Github.
    eval(parse(text = RCurl::getURL("https://raw.githubusercontent.com/edunford/bayes.functions/master/limited.dep.vars.funcs.R")))
    
    
    # Predicted Probabilities....
    pprobs(mcmc.output = HHfit,
           model.data = jags.data,
           params = c("alpha","b"),
           manipulation.var = "schooling",dv="protectionist",
           type = "smoother",points = T)
    
    # Looks like Schooling greatly decreases your odds
    
    
    # Now, let's look at the discrete differences between completing high school and completing college. 
    diffFX(mcmc.output = HHfit,model.data = jags.data,params = c("alpha","b"),
           from = 12,to=16,manipulation.var ="schooling",dv="protectionist",type="ci.distribution")
    
    
    
    # As these functions show, there are easy ways of visualizing the posterior
    # distribution to tell us something about what is going. 
    
          

    
        
# [APPLIED EXAMPLE 2] Grievance Sentiments on Twitter and the Propensity of Conflict around Election Periods in Africa ----------------------------------------   

        # We'll now look at a setting up a hierarchical model using JAGS.
    
    # DATA:
        # The following data looks at the twitter activity of five African 
        # countries 90 days around an election (60 days prior & 30 days after). 
        # The twitter activity variable "grievSentL1" captures the level of 
        # grienvance expressions drawn from a semi-supervise machine coding 
        # technique used to identify "relevant" tweets (i.e. tweets about 
        # mobilizing/opposing the incumbent regime). The goal is to examine the 
        # relationship between expressions of grievances "on the ground" and the
        # propensity for an "instability event" (protest, riots, organized
        # violence, etc.) of occurring. 
    
        # Explanatory Variables: There are a number of country-level variables 
        # that theoretically matter and have been shown to influence the 
        # occurrence of conflict (ethnic exclusion from power, repressive 
        # capatiy of the regime, regime type, economic factors, population), but
        # these indicators don't vary by day --- which our twitter and the 
        # sub-national occurrence of violence does. Thus, we use the 
        # country-level to situated (and contextualize) the underlying 
        # propensity of obeserving an event. In addition, we include "twitter 
        # relevant" controls, such as the log of the total number of tweets, the
        # log of the total number of relevant tweets. All twitter variables are
        # lagged by a day.
    
        
        # Load the data
        electDF <- read.csv(file="elec_viol.csv",stringsAsFactors = F)
        
        head(electDF)
        
        # How is the data distributed?
        plot(as.Date(electDF$date),electDF$activity,pch=16,col=as.factor(electDF$country)) # time and DV
        plot(as.Date(electDF$date),electDF$grievSentL1,pch=16,col=as.factor(electDF$country)) # time and key IV
        hist(electDF$grievSentL1,breaks=50,col="aquamarine",border="white") # hist
        qplot(grievSentL1,activity,data=electDF) + theme_classic() + geom_smooth() # Bivariate relationshp
        
        
        # My Theory: When grievance sentiments are high, the more instability events we observe. 
        
        
        
    # Organizing the data for JAGS
        DF = electDF[complete.cases(electDF),] # Ensure Complete Information (we can't have missing entries)
        DF = DF %>% mutate(year = as.numeric(format.Date(as.Date(date),"%Y")),
                           country_year = paste0(country,"_",year)) %>% select(-year) %>% as.data.frame(.)
        head(DF)
        
        relevant_vars = c("activity","grievSentL1","relevantL1","allL1","activityL1","lngdppc",
             "lnpop","polity2","repress","exclpop","exclgrps","groups")

        jags.data = as.list(DF[, colnames(DF) %in% relevant_vars])
        jags.data$N = nrow(DF)
        jags.data$N.cyears =  DF$country_year %>% unique(.) %>% length(.) 
        cyears = DF$country_year 
        jags.data$cyears = as.numeric(as.factor(cyears)) # country-years to iterate through
        str(jags.data)
        
        
    # Construct the hierarchical model 
        model.hier = function(){
          
          # Level one: country-day
          for(i in 1:N){
            activity[i] ~ dpois(lambda[i])
            log(lambda[i]) <- b.cyears[cyears[i]] + b.sent*grievSentL1[i] + b.rel*relevantL1[i] + b.all*allL1[i] + b.dvLag*activityL1[i]
          }
          
          # Level two: country-year (established here as varying intercepts)
          for(j in 1:N.cyears){ 
            b.cyears[j] ~ dnorm(cyear_mu[j],cyear_tau[j])
            cyear_mu[j] <- b.constant + b.polity*polity2[j] + b.gdp*lngdppc[j] + b.pop*lnpop[j] + b.repress*repress[j] +
              b.exclpop*exclpop[j] + b.exclgrps*exclgrps[j] + b.groups*groups[j]
            
            # uninformative hyperpriors on level 2
            cyear_tau[j] ~ dgamma(1,1) 
            cyear.sigma[j] <- (1 / cyear_tau[j]) # Record variance (1/precision)
          }
          
          
          # Priors -- uninformative
          b.constant ~ dnorm(0,.01)
          b.sent ~ dnorm(0,.01)
          b.dvLag ~ dnorm(0,.01)
          b.rel ~ dnorm(0,.01)
          b.all ~ dnorm(0,.01)
          b.polity ~ dnorm(0,.01)
          b.gdp ~ dnorm(0,.01)
          b.pop ~ dnorm(0,.01)
          b.repress ~ dnorm(0,.01)
          b.exclpop ~ dnorm(0,.01)
          b.exclgrps ~ dnorm(0,.01)
          b.groups ~ dnorm(0,.01)
          
          # Pooled Effects
          pool.cyears <- mean(b.cyears[])
        }
        
        
        params = c("b.cyears","b.sent","b.dvLag","b.rel","b.all","b.constant","b.polity",
                   "b.gdp","b.pop","b.repress","b.exclpop","b.exclgrps","b.groups","pool.cyears")
        
        # Run the model (this will take a moment)
        hier.model.fit <- jags(data = jags.data,n.thin = 100,
                               parameters.to.save = params, n.chains = 2, n.iter = 100000, 
                               n.burnin = 10000, model.file = model.hier)
        
    # Useful to save the posterior object (especially when it takes time to run)
    # save(hier.model.fit,file="~/Desktop/hier.model.fit.Rdata") # save output object
    # load("~/Desktop/hier.model.fit.Rdata")
    
    # Assess convergence
    mcmcplot(hier.model.fit) # Initial results look promising (convergence)
    geweke.diag(hier.model.fit)
    
    
    # Shinystan ----- RSTAN offer a whole new and exciting way of approaching 
    # MCMC with its no U-Turn Hamiltonian, right now we are just going to
    # commandeer it's beautiful posterior visualization capabilities
    
        #  install.packages("shinystan")
        shinystan::launch_shinystan(shinystan::as.shinystan(as.mcmc(hier.model.fit)))
    
    
    
    # Summarize the posterior
    msum <- function(x){
      # Quick clean function to read the output.
      f <- summary(coda::as.mcmc(x))
      out <- cbind(f$statistics[,c(1,2)],f$quantiles[,c(1,5)])
      return(round(out,2))
    }
    
    msum(hier.model.fit)
    
    
    # Presenting Posterior Distributions for our Multi-level model
    rel_cols = names(hier.model.fit$BUGSoutput$sims.array[1,1,])
    rel_cols2 = rel_cols[c(1,3:8,9,16,18)]
    
    # Plot Posterior
    post = as.data.frame(hier.model.fit$BUGSoutput$sims.array[,1,rel_cols2])
    colnames(post)[2:7] = c("Botswana_2014","Kenya_2013","Nigeria_2015",
                            "South Africa_2014","Zambia_2014","Zambia_2015")
    post <- post %>% select(Botswana_2014,Kenya_2013:Zambia_2015,everything())
    colMeans(post)
    
    ggplot(data=reshape2::melt(post)) + facet_wrap(facets=~variable,nrow = 5,ncol=2) +
      geom_histogram(aes(value,fill=ifelse(value>=0,"orange","blue"))) + 
      xlim(-10,10) + theme_minimal() + geom_vline(xintercept = 0,lwd=1,col="grey",alpha=.5) +
      theme(legend.position = "none")
    
    # Zoom in on the posterior of the grievance sentiment variable
    hist(post$b.sent,col="steelblue",border="white",breaks=25,xlim=c(-5,10))
    abline(v=0,h=0,col="grey",lwd=3)
