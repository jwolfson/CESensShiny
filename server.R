library(shiny)

# Define the server logic
shinyServer(function(input, output) {
  
  ## Read in the file
  data <- reactive( {
    inFile <- input$file
    if(!is.null(inFile)) {
      dat <- read.csv(inFile$datapath,header=TRUE)
    }
    else {
      if(is.null(input$dataset)) {
        thedataset <- "fake"
      } else {  thedataset <- input$dataset}
      load(paste0("data/",thedataset,".RData")) 
      #load(paste0(thedataset,".RData")) 
      dat <- get(thedataset)  
    }
    dat })
  
  output$datadescription <- renderText( {
    switch(input$dataset,
           "fake" = "Simulated data for testing purposes",
           "modpr8" = "Modified version of antibody titer response to PR-8 flu vaccine",
          "modweiss" = "Modified version of antibody tier response to Weiss flue vaccine",
          "upload" = "Custom user dataset")      
  })
  
  #pS1Y1 <- dnorm(x.pts,mean=muS1Y1,sd=sdS1Y1)
  #pS1Y0 <- pS1Y1
  S1 <- reactive( data()$S[data()$Z==1] )
  muS1 <- reactive( mean(S1()) )
  sdS1 <- reactive( sd(S1()) )
  muS1Y1 <- reactive( mean(data()$S[data()$Z==1&data()$Y==1]) )
  sdS1Y1 <- reactive( sd(data()$S[data()$Z==1&data()$Y==1]) )
  
    dS1.np <- reactive( {
      switch(input$dist,
             "np" = density(S1(),kernel=input$kernel,adjust=input$adjust),
             "normal" = NULL,
             "gamma" = NULL) })

    pS1.np <- reactive( {
    switch(input$dist,
           "np" = approxfun(x=dS1.np()$x,y=dS1.np()$y),
           "normal" = NULL,
           "gamma" = NULL) })

    dS1Y1.np <- reactive( {
      switch(input$dist,
             "np" = density(data()$S[(data()$Z==1)&(data()$Y==1)],kernel=input$kernel,adjust=input$adjust),
             "normal" = NULL,
             "gamma" = NULL) })
  
    pS1Y1.np <- reactive( {
      switch(input$dist,
             "np" = approxfun(dS1Y1.np()$x,dS1Y1.np()$y),
             "normal" = NULL,
             "gamma" = NULL) })
  
  pY1 <- reactive( mean(data()$Y[data()$Z==1]) )
  pY0 <- reactive( mean(data()$Y[data()$Z==0]) )
  x.pts <- reactive( seq(min(S1()),max(S1()),length.out=100) )
  
  doBoot <- function(data,R) {
    allCEs <- lapply(1:R,function(i) {
      rows <- sample(1:nrow(data),nrow(data),replace=TRUE)
      dat <- data.frame(data[rows,])
      S1 <- dat$S[dat$Z==1]
      muS1 <- mean(S1)
      sdS1 <- sd(S1)
      #obsS1Y1 <- observe({print(sum(dat$Z==1&dat$Y==1))})
      if(sum(dat$Z==1&dat$Y==1)>=2) {
      muS1Y1 <- mean(dat$S[(dat$Z==1)&(dat$Y==1)])
      sdS1Y1 <- sd(dat$S[(dat$Z==1)&(dat$Y==1)])

      if(input$dist=="np") {
        dS1.np <- density(S1,kernel=input$kernel,adjust=input$adjust)
        pS1.np <- approxfun(dS1.np$x,dS1.np$y)
      
        dS1Y1.np <- density(dat$S[(dat$Z==1)&(dat$Y==1)],kernel=input$kernel,adjust=input$adjust)
        pS1Y1.np <- approxfun(dS1Y1.np$x,dS1Y1.np$y)
      }
      
      pY1 <- mean(dat$Y[dat$Z==1])
      pY0 <- mean(dat$Y[dat$Z==0])
      x.pts <- seq(min(S1),max(S1),length.out=100)

      gam.shape <- muS1^2/sdS1^2
      gam.rate <- muS1/sdS1^2
      
      pS1 <- switch(input$dist,
               "norm" = dnorm(x.pts,mean=muS1,sd=sdS1),
               "gamma" = suppressWarnings(dgamma(x.pts,shape=gam.shape,rate=gam.rate)),
                "np" = pS1.np(x.pts) )
      
      gam.shape <- muS1Y1^2/sdS1Y1^2
      gam.rate <- muS1Y1/sdS1Y1^2
        
      pS1Y1 <- switch(input$dist,
               "norm" = dnorm(x.pts,mean=muS1Y1,sd=sdS1Y1),
               "gamma" = suppressWarnings(dgamma(x.pts,shape=gam.shape,rate=gam.rate)),
                      "np" = pS1Y1.np(x.pts))
      
      shift.mu <- muS1+input$meanshift*sdS1
      scale.sd <- sdS1*input$sdshift
      gam.shape <- shift.mu^2/scale.sd^2
      gam.rate <- shift.mu/scale.sd^2
      
      pS1Y0 <- switch(input$dist, 
             "norm" = dnorm(x.pts,mean=shift.mu,sd=scale.sd),
             "gamma" = suppressWarnings(dgamma(x.pts,shape=gam.shape,rate=gam.rate)),
              "np" = pS1.np(sdS1/scale.sd*(x.pts - shift.mu)+muS1) )
      
      CE <- pS1Y1*pY1/pS1 - pS1Y0*pY0/pS1
      } else { CE <- rep(NA,length(x.pts)) }
      return(CE)
    })
    return(allCEs)
  }
  
  pS1 <- reactive( {
    gam.shape <- muS1()^2/sdS1()^2
    gam.rate <- muS1()/sdS1()^2

    switch(input$dist,
           "norm" = dnorm(x.pts(),mean=muS1(),sd=sdS1()),
           "gamma" = suppressWarnings(dgamma(x.pts(),shape=gam.shape,rate=gam.rate)),
           "np" = pS1.np()(x.pts()))
  })

  pS1Y1 <- reactive( {
    gam.shape <- muS1Y1()^2/sdS1Y1()^2
    gam.rate <- muS1Y1()/sdS1Y1()^2
    switch(input$dist,
           "norm" = dnorm(x.pts(),mean=muS1Y1(),sd=sdS1Y1()),
           "gamma" = suppressWarnings(dgamma(x.pts(),shape=gam.shape,rate=gam.rate)),
           "np" = pS1Y1.np()(x.pts()) )
  })
  
  update.pS1Y0 <- reactive( {
    
    shift.mu <- muS1()+input$meanshift*sdS1()
    scale.sd <- sdS1()*input$sdshift
    gam.shape <- shift.mu^2/scale.sd^2
    gam.rate <- shift.mu/scale.sd^2
    
    switch(input$dist, 
      "norm" = dnorm(x.pts(),mean=shift.mu,sd=scale.sd),
      "gamma" = dgamma(x.pts(),shape=gam.shape,rate=gam.rate),
      "np" = pS1.np()(sdS1()/scale.sd*(x.pts() - shift.mu)+muS1()))
  })

  CEval <- reactive( {
    pS1Y1()*pY1()/pS1() - update.pS1Y0()*pY0()/pS1()    
  })
  
  output$OOB <- renderText( {
    cevals <- CEval()
    if(all(is.na(cevals))) {
      a <- "WARNING: Chosen sensitivity parameters result in estimated causal effect outside [-1,1]"
    } else if(min(na.omit(cevals))< (-1) | max(na.omit(cevals)) > 1) {
      a <- "WARNING: Chosen sensitivity parameters result in estimated causal effect outside [-1,1]"
    } else { a <- "" } 
    a
  })
  
  output$dens <- renderPlot( {
    hist(S1(),probability=TRUE,main="Histogram of S(1)",xlab="S(1)")
    lines(x.pts(),pS1(),lwd=2)
    lines(x.pts(),update.pS1Y0(),col="red")
  })
  
  output$CE <- renderPlot( {
    plot(x.pts(),CEval(),type="l",lwd=2,
         main="Causal Effect = E(Y(1)-Y(0) | S(1))",
         xlab="S(1)",ylab="Causal Effect",ylim=input$yl)
    abline(h=0)
    if(input$boot) { 
      boot.CEs <- doBoot(data(),input$nboot)
      #boot.SEs <- apply(do.call(rbind,boot.CEs),2,sd,na.rm=TRUE)
      #boot.L <- CE - 1.96*boot.SEs
      #boot.U <- CE + 1.96*boot.SEs
      boot.L <- apply(do.call(rbind,boot.CEs),2,function(x) { quantile(na.omit(x),0.025)})
      boot.U <- apply(do.call(rbind,boot.CEs),2,function(x) { quantile(na.omit(x),0.975)})
      lines(x.pts(),boot.L,lty="dashed")
      lines(x.pts(),boot.U,lty="dashed")
  }

  })
  
})