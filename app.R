#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(tidyverse)
library(patchwork) # for combining multiple plots into 1

#user interface
ui <- fluidPage(pageWithSidebar( 
  
  headerPanel = headerPanel("Diploid model of selection with underdominance"),
  
  sidebarPanel(
    
    HTML("<p style='font-size:14px'><B>Fitness of AA (wAA)= 1.2<BR>
         Fitness of aa (waa)= 1.0<BR>
         Fitness of Aa (wAa)= 1 - h_under</B>"),

    sliderInput(inputId = "h_under", label = "Underdominance coefficient (h_under)", value = -0.2, 
                min = -1, max = 0, step = 0.01),
    
    sliderInput(inputId = "p_0", label = "Initial p", value = 0.32, 
                min = 0, max = 1, step = 0.02),
    
    sliderInput(inputId = "gen", label = "Number of generations", value = 50, 
                min = 5, max = 1000, step = 5),
  
    HTML("<p style='font-size:8px'>
          Modified from: Copyright (c) 2017 Silas Tittes, MIT License, <A href='https://github.com/silastittes/shiny_popgen'>https://github.com/silastittes/shiny_popgen</a></p>")
    
  ), 
  
  mainPanel =  mainPanel(
    plotOutput(outputId = 'viz'),
    
    selectInput("select", label = "Plot options", 
                choices = list("Mean fitness by p" = 1, "Change in p by p" = 2,
                               "Allele frequency over time" = 3, "Genotypes over time" = 4, "Mean fitness over time" = 5), selected = 4)
    
    
  )
))

#back end code and response to user input
server <- function(input, output){
  
  output$viz <- renderPlot({
    
    p <- seq(0, 1, length.out = 1000)
    #parameters
    wAA = 1.2
    h_under = input$h_under
    waa = 1.0
    p_0 = input$p_0
    gen = (input$gen)+1
    t = seq(0, gen, length.out = 1000)
    
    # calculate the fitness of the heterozygote using the dominance coefficient,
    wAa <- 1 + h_under
      
    # make a plot of the fitness landscape
    schematic_plot <- data.frame(Genotype = c("aa", "aA", "Aa", "AA"),
                            Fitness = c(waa, rep(wAa, 2), wAA)) %>% 
                    ggplot(aes(x=Genotype, y=Fitness, fill=Genotype)) +
                      geom_bar(stat="identity") +
                      scale_fill_manual(values = c("AA" = "firebrick", "Aa" = "orange", "aA" = "orange", "aa" = "gold")) +
                      theme_classic()+theme(text = element_text(size=20),
                                            legend.position="none")
    

    
    #if(input$w_plot){
    if(input$select == 1){
      W <- p^2*wAA + 2*p*(1-p)*wAa + (1-p)^2*waa
      selected_plot <- data.frame(p=p, W=W) %>%
        ggplot(aes(x = p, y = W)) +
        geom_line(color="forestgreen", size=2) +
        xlab("Mean Fitness") +
        xlab("Allele frequency, p") +
#        xlim(0, 1) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(min(waa,wAa,wAA)-0.001, max(waa,wAa,wAA)+0.001)) + 
#        theme(panel.background =  
#                element_rect(fill =  rgb(30, 144, 255, 25, 
#                                         maxColorValue = 255)),
#              text = element_text(size=16, family= "Times"))
        theme_classic()+theme(text = element_text(size=20))
      
      #} else if(input$delta_plot){
    } else if(input$select == 2){
      W <- p^2*wAA + 2*p*(1-p)*wAa + (1-p)^2*waa
      delta_p <- (p^2*wAA + p*(1-p)*wAa) / W - p
      setlim <- (max(waa,wAa,wAA)-min(waa,wAa,wAA))/(wAA + wAa + waa)
      setlimround <- round(10*setlim)/10+0.1
      selected_plot <- data.frame(p=p, delta_p=delta_p) %>%
        ggplot(aes(x = p, y = delta_p)) +
        geom_line(color="firebrick", size=2) +
        geom_hline(yintercept = 0, lty = 2) +
        ylab(expression(paste(Delta,p))) +
        xlab("Allele frequency, p") +
#        xlim(0, 1) +
        scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(-setlimround,setlimround)) + 
#        theme(panel.background =  
#                element_rect(fill =  rgb(30, 144, 255, 25, 
#                                         maxColorValue = 255)),
#              text = element_text(size=16, family= "Times"))
        theme_classic()+theme(text = element_text(size=20))
      
      #} else if(input$time_plot){
    } else if(input$select == 3){
      
      p_t <- rep(NA, gen)
      p_t[1] <- p_0
      for(i in 2:gen){
        W <- p_t[i-1]^2*wAA + 2*p_t[i-1]*(1-p_t[i-1])*wAa + (1-p_t[i-1])^2*waa
        p_t[i] <- (p_t[i-1]^2 * wAA + p_t[i-1] * (1-p_t[i-1]) * wAa) / W
      }
      
      selected_plot <- data.frame(t = 1:gen, p_t = p_t) %>%
        ggplot(aes(x = t-1, y = p_t)) +
        geom_point(color="firebrick", size=2) +
        xlab("Generation") +
        ylab("Allele frequency, p") +
        annotate("text", x = 0.9*gen, y = 0.14, label = "Final p")+
        annotate("text", x = 0.9*gen, y = 0.1, label = round(p_t[gen],5))+
#        ylim(0, 1) + 
        scale_x_continuous(expand = c(0, 0), limits = c(0, input$gen+2)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) + 
#        theme(panel.background =  
#                element_rect(fill =  rgb(30, 144, 255, 25, 
#                                         maxColorValue = 255)),
#              text = element_text(size=16, family= "Times"))
        theme_classic()+theme(text = element_text(size=20))
    } else if(input$select == 4){
      
      p_t <- rep(NA, gen)
      p_t[1] <- p_0
      for(i in 2:gen){
        W <- p_t[i-1]^2*wAA + 2*p_t[i-1]*(1-p_t[i-1])*wAa + (1-p_t[i-1])^2*waa
        p_t[i] <- (p_t[i-1]^2 * wAA + p_t[i-1] * (1-p_t[i-1]) * wAa) / W
      }
      
      # simulate forward in time for each genotype
      temp_df <- data.frame(t = 1:gen, p_t = p_t)
      temp_df <- rbind(cbind(temp_df, y=(1-temp_df$p_t)^2, Genotype = "genotype aa"),
                       cbind(temp_df, y=2*temp_df$p_t*(1-temp_df$p_t), Genotype = "genotypes aA & Aa"),
                       cbind(temp_df, y=temp_df$p_t^2, Genotype = "genotype AA"))
      temp_df$Genotype <- factor(temp_df$Genotype,
                                 levels = c("genotype aa", "genotypes aA & Aa", "genotype AA"))
        
      selected_plot <-  temp_df %>%
        ggplot(aes(x = t-1, y = y, colour=Genotype)) +
        geom_point(size=1) +
        scale_color_manual(values = c("genotype AA" = "firebrick", "genotypes aA & Aa" = "orange", "genotype aa" = "gold")) +
        labs(x= "Generation",y="Allele frequency, p",
             color = "Legend") +
        scale_x_continuous(expand = c(0, 0), limits = c(0, input$gen+2)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(0, 1.01)) + 
        #        theme(panel.background =  
        #                element_rect(fill =  rgb(30, 144, 255, 25, 
        #                                         maxColorValue = 255)),
        #              text = element_text(size=16, family= "Times"))
        theme_classic()+theme(text = element_text(size=20))
    } else if(input$select == 5){
      p_t <- rep(NA, gen)
      W_t <- rep(NA, gen)
      p_t[1] <- p_0
      W_t[1] <- p_0^2*wAA + 2*p_0*(1-p_0)*wAa + (1-p_0)^2*waa
      for(i in 2:gen){
        W_t[i] <- p_t[i-1]^2*wAA + 2*p_t[i-1]*(1-p_t[i-1])*wAa + (1-p_t[i-1])^2*waa
        p_t[i] <- (p_t[i-1]^2 * wAA + p_t[i-1] * (1-p_t[i-1]) * wAa) / W_t[i]
      }
      
      selected_plot <- data.frame(t = 1:gen, W_t = W_t) %>%
        ggplot(aes(x = t-1, y = W_t)) +
        geom_point(color="forestgreen", size=2) +
        xlab("Generation") +
        ylab("Mean fitness") +
        #        ylim(0, 1) + 
        scale_x_continuous(expand = c(0, 0), limits = c(0, input$gen+2)) + 
        scale_y_continuous(expand = c(0, 0), limits = c(min(waa,wAa,wAA)-0.001, max(waa,wAa,wAA)+0.001)) + 
        #        theme(panel.background =  
        #                element_rect(fill =  rgb(30, 144, 255, 25, 
        #                                         maxColorValue = 255)),
        #              text = element_text(size=16, family= "Times"))
        theme_classic()+theme(text = element_text(size=20))
    }
    
    return(wrap_plots(schematic_plot, 
                      selected_plot, 
                      nrow = 2, ncol = 1, byrow = FALSE))
  })
}
# Run the application 
shinyApp(ui = ui, server = server)