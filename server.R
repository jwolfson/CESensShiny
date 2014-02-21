library(shiny)

load("C:/Users/Julian/Google Drive/code library/CESensShiny/data/fake.Rdata")
data <- fake
#pS1Y1 <- dnorm(x.pts,mean=muS1Y1,sd=sdS1Y1)
#pS1Y0 <- pS1Y1
S1 <- data$S[data$Z==1]
muS1 <- mean(S1)
sdS1 <- sd(S1)
muS1Y1 <- mean(data$S[data$Z==1&data$Y==1])
sdS1Y1 <- sd(data$S[data$Z==1&data$Y==1])

dS1.np <- density(S1)
pS1.np <- approxfun(dS1.np$x,dS1.np$y)

dS1Y1.np <- density(data$S[data$Z==1&data$Y==1])
pS1Y1.np <- approxfun(dS1Y1.np$x,dS1Y1.np$y)

pY1 <- mean(data$Y[data$Z==1])
pY0 <- mean(data$Y[data$Z==0])
x.pts <- seq(min(S1),max(S1),length.out=100)

# Define server logic
shinyServer(function(input, output) {
    
  doBoot <- function(data,R) {
    allCEs <- sapply(1:R,function(i) {
      rows <- sample(1:nrow(data),nrow(data),replace=TRUE)
      dat <- data[rows,]
      S1 <- dat$S[dat$Z==1]
      muS1 <- mean(S1)
      sdS1 <- sd(S1)
      muS1Y1 <- mean(dat$S[dat$Z==1&dat$Y==1])
      sdS1Y1 <- sd(data$S[dat$Z==1&dat$Y==1])

      dS1.np <- density(S1)
      pS1.np <- approxfun(dS1.np$x,dS1.np$y)
      
      dS1Y1.np <- density(data$S[data$Z==1&data$Y==1])
      pS1Y1.np <- approxfun(dS1Y1.np$x,dS1Y1.np$y)
      
      pY1 <- mean(data$Y[dat$Z==1])
      pY0 <- mean(data$Y[dat$Z==0])
      x.pts <- seq(min(S1),max(S1),length.out=100)

      gam.shape <- muS1^2/sdS1^2
      gam.rate <- muS1/sdS1^2
      
      pS1 <- switch(input$dist,
               "norm" = dnorm(x.pts,mean=muS1,sd=sdS1),
               "gamma" = dgamma(x.pts,shape=gam.shape,rate=gam.rate),
                "np" = pS1.np(x.pts) )
      
      gam.shape <- muS1Y1^2/sdS1Y1^2
      gam.rate <- muS1Y1/sdS1Y1^2
        
      pS1Y1 <- switch(input$dist,
               "norm" = dnorm(x.pts,mean=muS1Y1,sd=sdS1Y1),
               "gamma" = dgamma(x.pts,shape=gam.shape,rate=gam.rate),
                      "np" = pS1Y1.np(x.pts))
      
      shift.mu <- muS1+input$meanshift*sdS1
      scale.sd <- sdS1*input$sdshift
      gam.shape <- shift.mu^2/scale.sd^2
      gam.rate <- shift.mu/scale.sd^2
      
      pS1Y0 <- switch(input$dist, 
             "norm" = dnorm(x.pts,mean=shift.mu,sd=scale.sd),
             "gamma" = dgamma(x.pts,shape=gam.shape,rate=gam.rate),
              "np" = pS1.np(sdS1/scale.sd*(x.pts - shift.mu)+muS1) )
      
      CE <- pS1Y1*pY1/pS1 - pS1Y0*pY0/pS1  
      return(CE)
    })
    return(allCEs)
  }
  
  pS1 <- reactive( {
    gam.shape <- muS1^2/sdS1^2
    gam.rate <- muS1/sdS1^2

    switch(input$dist,
           "norm" = dnorm(x.pts,mean=muS1,sd=sdS1),
           "gamma" = dgamma(x.pts,shape=gam.shape,rate=gam.rate),
           "np" = pS1.np(x.pts))
  })

  pS1Y1 <- reactive( {
    gam.shape <- muS1Y1^2/sdS1Y1^2
    gam.rate <- muS1Y1/sdS1Y1^2
    
    switch(input$dist,
           "norm" = dnorm(x.pts,mean=muS1Y1,sd=sdS1Y1),
           "gamma" = dgamma(x.pts,shape=gam.shape,rate=gam.rate),
           "np" = pS1Y1.np(x.pts))
  })
  
  update.pS1Y0 <- reactive( {
    
    shift.mu <- muS1+input$meanshift*sdS1
    scale.sd <- sdS1*input$sdshift
    gam.shape <- shift.mu^2/scale.sd^2
    gam.rate <- shift.mu/scale.sd^2
    
    
    switch(input$dist, 
      "norm" = dnorm(x.pts,mean=shift.mu,sd=scale.sd),
      "gamma" = dgamma(x.pts,shape=gam.shape,rate=gam.rate),
      "np" = pS1.np(sdS1/scale.sd*(x.pts - shift.mu)+muS1))
  })
    
  output$dens <- renderPlot( {
    hist(S1,probability=TRUE)
    lines(x.pts,pS1(),lwd=2)
    lines(x.pts,update.pS1Y0(),col="red")
    
  })
  
  output$CE <- renderPlot( {
    CE <- pS1Y1()*pY1/pS1() - update.pS1Y0()*pY0/pS1()
    plot(x.pts,CE,type="l",main="Causal Effect = E(Y(1)-Y(0) | S(1))",xlab="S(1)",ylab="CE",ylim=input$yl)
    abline(h=0)
    if(any(abs(na.omit(CE))>1)) {
      text(quantile(x.pts,0.1),mean(CE),"WARNING: Specified sensitivity parameters result in estimated CE outside [-1,1]",col="red",adj=0,cex=1.3) }
    if(input$boot) { 
      boot.CEs <- doBoot(data,input$n.boot)
      boot.SEs <- apply(boot.CEs,1,sd,na.rm=TRUE)
      boot.L <- CE - 1.96*boot.SEs
      boot.U <- CE + 1.96*boot.SEs
      ##boot.L <- apply(boot.CEs,1,quantile,probs=0.025,na.rm=TRUE)
      ##boot.U <- apply(boot.CEs,1,quantile,probs=0.975,na.rm=TRUE)
      lines(x.pts,boot.L,lty="dashed")
      lines(x.pts,boot.U,lty="dashed")
  }

  })  

})