#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(RColorBrewer)
library(plotly)

brewerCols = brewerCols = brewer.pal(3,'Set1')


d = read.table("combineResults.py.txt",header=1,sep="\t")

d$found = 0
d[d$simulatedRecoveredCount > 0,]$found = 1

d$found5 = 0
d[d$simulatedRecoveredCount >= 5,]$found5 = 1

d$foundHalf = 0
d[d$simulatedRecoveredCount >= d$simulatedCount/2,]$foundHalf = 1
d[d$simulatedCount == 0,]$foundHalf = 0
d = d[d$pctMut <= 0.4,]

plotThing = function(dd,titleText,plotColor)
{
  
  psize = 50 #number of permutations
  ciWidth = 0.3 #width of confidence interval crossbar
  avs = NULL
  q1 = NULL #1st percentile
  q99 = NULL #99th percentile
  levs = levels(factor(dd$pctMut))
  for (l in levs)
  {
    ddd = dd[dd$pctMut == l,]
    zeroes = nrow(ddd[ddd$found == 0,])
    ones = nrow(ddd[ddd$found == 1,])
    avs = c(avs,ones/(ones + zeroes))
    
    #do permuation
    pp = NULL
    for (i in 1:100)
    {
      ss = sample(1:nrow(ddd),size = psize)
      sss = ddd[ss,]
      pp = c(pp,nrow(sss[sss$found ==1,])/psize)
    }
    qs = quantile(pp,c(0.01,0.99))
    q1 = c(q1,qs[1]*100)
    q99 = c(q99,qs[2]*100)
  }
  
  print(cbind(as.numeric(levels(factor(dd$pctMut))),as.numeric(avs),as.numeric(q1),as.numeric(q99)))
  plot(as.numeric(levs)*100,avs*100,xlab="Percent Reads Mutated",ylab="Percent of Simulated Samples",main=titleText,pch=20,ylim=c(0,100),col=plotColor)
  lines(as.numeric(levs)*100,avs*100,col=plotColor)
  for (i in 1:length(avs))
  {
    if(q99[i] - q1[i] == 0)
    {
      next
    }
    m = avs[i]
    x = as.numeric(levs[i])*100
    segments(x,q1[i],x,q99[i],col=plotColor)
    segments(x-ciWidth,q1[i],x+ciWidth,q1[i],col=plotColor)
    segments(x-ciWidth,q99[i],x+ciWidth,q99[i],col=plotColor)
  }
}


pointsThing = function(dd,dataCol,dataColor)
{
  
  psize = 50 #number of permutations
  ciWidth = 0.5 #width of confidence interval crossbar
  avs = NULL
  q1 = NULL #1st percentile
  q99 = NULL #99th percentile
  levs = levels(factor(dd$pctMut))
  for (l in levs)
  {
    ddd = dd[dd$pctMut == l,]
    zeroes = nrow(ddd[ddd[,dataCol] == 0,])
    ones = nrow(ddd[ddd[,dataCol] == 1,])
    avs = c(avs,ones/(ones + zeroes))
    
    #do permuation
    pp = NULL
    for (i in 1:100)
    {
      ss = sample(1:nrow(ddd),size = psize)
      sss = ddd[ss,]
      pp = c(pp,nrow(sss[sss[,dataCol] ==1,])/psize)
    }
    qs = quantile(pp,c(0.01,0.99))
    q1 = c(q1,qs[1]*100)
    q99 = c(q99,qs[2]*100)
  }
  
  print(paste('data:',dataCol))
  print(cbind(as.numeric(levels(factor(dd$pctMut))),as.numeric(avs),as.numeric(q1),as.numeric(q99)))
  points(as.numeric(levs)*100,avs*100,pch=20,ylim=c(0,100),col=dataColor)
  lines(as.numeric(levs)*100,avs*100,col=dataColor)
  for (i in 1:length(avs))
  {
    if(q99[i] - q1[i] == 0)
    {
      next
    }
    m = avs[i]
    x = as.numeric(levs[i])*100
    segments(x,q1[i],x,q99[i],col=dataColor)
    segments(x-ciWidth,q1[i],x+ciWidth,q1[i],col=dataColor)
    segments(x-ciWidth,q99[i],x+ciWidth,q99[i],col=dataColor)
  }
}


# Define UI for application that draws a histogram
ui <- fluidPage(
  fluidRow(
    fluidRow(
      column(12,align="center",
         h3("Probability of observing an editing event")
      )
    ),
    fluidRow(
      column(6,align="center",
      sliderInput('trans_eff','Transfection efficiency (%)',min=0,max=100,value=50,step=1,round=0),
      sliderInput('seq_depth','Sequencing depth',min=0,max=300,value=20,step=10),
      sliderInput('min_obs','Min number of observed reads',min=0,max=10,value=1,step=1)
      ),
      column(6,align="center",
             plotlyOutput("modelPlot")
      )
    )
  ),
  hr(),
  
  fluidRow(
    column(12,align="center",
   h3("Simulated Detection of Editing Events"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
      sidebarPanel(
         selectInput(
           inputId = "depth",
           label="Sequencing Depth:",
           selected = 10,
           choices = c(10,20,30,50,100)
           
         )
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot")
      )
   )
   )
  )
)


server <- function(input, output) {
   
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R
     print(paste('head of d is',head(d$depth)))
     print(paste('got here input depth is ',input$depth))
     dd = d[d$depth == input$depth,]
     plotThing(dd,paste0(input$depth,'x Coverage'),brewerCols[1])
     pointsThing(dd,'found5',brewerCols[2])
     pointsThing(dd,'foundHalf',brewerCols[3])
     legend('bottomright',fill=brewerCols,legend=c('>=1','>=5','50%'),title='Number/Percent of Simulated Events\nRecovered',bty='n')
     
   })
   output$modelPlot <- renderPlotly({
     trans_eff = input$trans_eff/100
     seq_depth = input$seq_depth
     min_obs = input$min_obs
     mutFreq = seq(0,1,by=0.01)
     pvals = 1- pbinom(q=min_obs-1,size=seq_depth,prob=mutFreq*trans_eff)
     
     modelData = data.frame(cbind(mutFreq,pvals))
     colnames(modelData) = c('Mutation_frequency','p_value')
     suffix = ""
     if (min_obs != 1)
     {
       suffix = "s"
     }
     
     plot <- modelData %>%
       ggplot(aes(x=Mutation_frequency,
                  y=p_value)) +
       geom_line()+
      
       xlab('Frequency of mutation') +
       ylab(paste0('Probability of observing at least ' ,min_obs, ' mutated read',suffix))
     p <- plotly_build(plot)
     p$elementId <- NULL
     p


     #plot(mutFreq,pvals,type='l',xlab='Frequency of mutation',ylab='Probability of observing\nat least 1 read with the mutation')
     
     
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

