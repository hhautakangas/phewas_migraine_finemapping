# Load packages ----
library(shiny)
library(shinydashboard)
library(data.table)
library(RColorBrewer)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(DT)

# load data
# for 3 color heatmap!
# option1

hm <- fread("data/temp_for_phewas_candidate_variants_FGR10_heatmap_3color.txt", data.table=F)
dt <- fread("data/phewas_FGR10_meta_betas_flipped_for_shiny_plot.txt", data.table=F)
colnames(dt)[ncol(dt)] <- "category_FG"
colnames(dt)[3:4] <- c("effect allele", "other allele")
# change columns for factors, for table filtering
for( ii in c(1:5,9:11)){
  dt[,ii] <- as.factor( dt[,ii])
}
# Source helpers ----


# User interface ----
ui <-ui <- dashboardPage(
  dashboardHeader(title ="PheWAS: Migraine candidate variants in FinnGen R10"),
  dashboardSidebar(sidebarMenu(
    menuItem("Candidate variants", tabName = "dashboard", icon = icon("dashboard")),
    menuItem("Functional variants", tabName = "functional_var", icon = icon("th")),
    menuItem("PIP > 0.1 variants", tabName = "pip_var", icon = icon("th"))
  )),
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "dashboard",
              
              fluidRow(
                box(
                  title = "PheWAS: Migraine candidate variants in FinnGen R10: Heatmap", width = 12, solidHeader = TRUE, status = "primary",
                  originalHeatmapOutput("ht", title = NULL)
                )),
              
              fluidRow(
                box(
                  title = "Sub-heatmap", width = 6, solidHeader = TRUE, status = "primary",
                  subHeatmapOutput("ht", title = NULL)
                ),
                box(
                  title = "Output", width = 6, solidHeader = TRUE, status = "primary",
                  HeatmapInfoOutput("ht", title = NULL)
                )
              ),
              
              fluidRow(title = "Table",
                       box(div(DT::dataTableOutput("mytable"),style = "font-size: 75%; width: 75%"), width = 12))
      ),
      
      # Second tab content
      tabItem(tabName = "functional_var",
              h2("Functional variants among the credible sets")
      ),
      # third
      tabItem(tabName = "pip_var",
              h2("Variants with PIP > 0.1 among the credible sets")
      )
    )
  )
)



# Server logic
server <- function(input, output, session){
  observe({       # <- this line is changed
    cluster_rows = FALSE
    cluster_columns = FALSE
    row_order = NULL
    column_order = NULL
    row_km = NULL
    column_km = NULL
    matrix = as.matrix(hm[,2:202])
    row.names(matrix) <- hm[,1]
    ht = Heatmap(matrix,  name = "migraine",
                 cluster_rows = cluster_rows, cluster_columns = cluster_columns,
                 row_order = row_order, column_order = column_order, 
                 row_km = row_km, column_km = column_km)
    
    makeInteractiveComplexHeatmap(input, output, session, ht, "ht")
  })
  
  observe({
    print(input$row1)
  })
  
  output$mytable = DT::renderDataTable({
    #Display table with select
    DT::datatable(dt, width = 7, rownames = FALSE, filter = 'top',
                  options = list(orderClasses = TRUE,
                                 lengthMenu = c( 10, 25, 50),
                                 pageLength = 10 ,
                                 
                                 drawCallback= JS(
                                   'function(settings) {
                                     Shiny.bindAll(this.api().table().node());}')
                  ),selection='none',escape=F)
    
    
  } )
}

# Run the app
shinyApp(ui = ui, server = server)





##
# Server logic
#server <- function(input, output){

#  output$heatmap <- renderPlot({
#    data <- hm
#    cols <- colorRampPalette(brewer.pal(8,"RdBu"))(3)
#    heatmap(as.matrix(data[,2:202]), Colv = NA, Rowv = NA, labRow =data$rsid, col=cols, cexRow=0.7, cexCol = 0.7, scale="none")

#  }, height = 800, width = 800)


#  output$info <- renderPrint({
#    nearPoints(hm, input$plot_click, threshold = 10, maxpoints = 1,
#               addDist = TRUE)
#  })

#}

