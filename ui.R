library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("CESens"),
  
  sidebarPanel(
    fileInput("file","Select data file"),
    selectInput("dist", "Distribution of S(1)",
                list("Normal" = "norm", 
                     "Gamma" = "gamma")),
    
    sliderInput("meanshift","Mean shift",-5,5,0,step=0.05),
    sliderInput("sdshift","SD scale",0.1,10,1,step=0.05),
    checkboxInput("boot","Compute Bootstrap CIs",value=FALSE),
    numericInput("n.boot","Number of bootstrap resamples",value=500)
  ),
  
  # Show the caption and plot of the requested variable against mpg
  mainPanel(
    plotOutput("dens"),
    plotOutput("CE")
  )
))