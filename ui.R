library(shiny)

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("CESens"),
  
  sidebarPanel(
    selectInput("dataset", "Dataset",selected="fake",
                choices=list("Simulated" = "fake", 
                             "Modified PR-8 vax" = "modpr8",
                             "Modified Weiss vax" = "modweiss",
                             "Upload CSV..." = "upload")),
    textOutput("datadescription"),
    conditionalPanel(
      condition = "input.dataset =='upload'",
      fileInput("file","Select data file"),
      tags$p('CSV file containing columns named Z, S, and Y')),
    tags$hr(),
    selectInput("dist", "Distribution of S(1)",selected="norm",
                choices=list("Normal" = "norm", 
                     "Gamma" = "gamma",
                      "Nonparametric" = "np")),
    conditionalPanel(
        condition = "input.dist == 'np'",
        selectInput("kernel","Smoothing kernel:",selected="gaussian",
                    choices=list("Gaussian" = "gaussian",
                                 "Epanechnikov" = "epanechnikov",
                                 "Rectangular" = "rectangular",
                                 "Triangular" = "triangular",
                                 "Cosine" = "cosine",
                                 "OptCosine" = "optcosine")),
        sliderInput("adjust","Bandwidth adjustment:",0.1,5,1,step=0.05)),
    
    tags$hr(),
    sliderInput("meanshift","Mean shift of S(1)|Y(0)=1, in SD",-5,5,0,step=0.05),
    sliderInput("sdshift","SD scale of S(1)|Y(0)=1",0.1,10,1,step=0.05),
    strong(textOutput("OOB")),
    tags$hr(),
    checkboxInput("boot","Compute 95% Bootstrap CIs",value=FALSE),
    numericInput("n.boot","Number of bootstrap resamples (min. 250)",value=500,min=250),
    tags$hr(),
    sliderInput("yl","Display range of CE plot",-1,1,c(-1,1),step=0.05)
  ),
  
  mainPanel(
    plotOutput("dens"),
    plotOutput("CE")
  )
))