library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("CESens"),
  
  sidebarPanel(
    fileInput("file","Select data file"),
    selectInput("dist", "Distribution of S(1)",
                list("Normal" = "norm", 
                     "Gamma" = "gamma",
                      "Nonparametric" = "np")),
    sliderInput("meanshift","Mean shift of S(1)|Y(0)=1, in SD",-5,5,0,step=0.05),
    sliderInput("sdshift","SD scale of S(1)|Y(0)=1",0.1,10,1,step=0.05),
    checkboxInput("boot","Compute Bootstrap CIs",value=FALSE),
    numericInput("n.boot","Number of bootstrap resamples",value=500),
    sliderInput("yl","Display range of CE plot",-1,1,c(-1,1),step=0.05)
  ),
  
  # Show the caption and plot of the requested variable against mpg
  mainPanel(
    plotOutput("dens"),
    plotOutput("CE")
  )
))