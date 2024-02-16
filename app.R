# Load packages ----
library(shiny)
library(shinydashboard)
library(BiocManager)
options(repos = BiocManager::repositories())
library(data.table)
library(RColorBrewer)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(DT)
library(ggplot2)
library(GetoptLong)

# load data
# for 3 color heatmap!
# option1
# only 181 variants
inc <- fread("data/top_config181_rsids_for_phewas.txt", data.table=F, header=F)
hm <- fread("data/temp_for_phewas_candidate_variants_FGR10_heatmap_3color2.txt", data.table=F)
#hm <- fread("data/temp_for_phewas_candidate_variants_FGR10_heatmap.txt", data.table=F)
matrix = as.matrix(hm[,2:ncol(hm)])
row.names(matrix) <- hm[,1]

dt <- fread("data/phewas_FGR10_meta_betas_flipped_for_shiny_plot2.txt", data.table=F)
dt$pval_FGR10 <- as.numeric(format(dt$pval_FGR10, digits = 3, scientific = 5))
colnames(dt)[ncol(dt)] <- "category_FG"
#colnames(dt)[5:6] <- c("effect allele", "other allele")
# change columns for factors, for table filtering
for( ii in c(1:8,13:15)){
  dt[,ii] <- as.factor( dt[,ii])
}

fm <- fread("data/temp_functional_for_phewas_candidate_variants_FGR10_heatmap_3color.txt", data.table=F)
matrix2 = as.matrix(fm[,2:ncol(fm)])
row.names(matrix2) <- fm[,1]

dtf <- fread("data/phewas_FGR10_functional_meta_betas_flipped_for_shiny_plot2.txt", data.table=F)
dtf$pval_FGR10 <- as.numeric(format(dtf$pval_FGR10, digits = 3, scientific = 5))
colnames(dtf)[ncol(dtf)] <- "category_FG"
#colnames(dtf)[5:6] <- c("effect allele", "other allele")
colnames(dtf)[2] <- "functional_variant"
# change columns for factors, for table filtering
for( ii in c(1:9,14:16)){
  dtf[,ii] <- as.factor( dtf[,ii])
}

fmh <- fread("data/temp_high_pip0.1_for_phewas_candidate_variants_FGR10_heatmap_3color.txt", data.table=F)

matrix3 = as.matrix(fmh[,2:ncol(fmh)])
row.names(matrix3) <- fmh[,1]

dth <- fread("data/phewas_FGR10_high_pip0.1_meta_betas_flipped_for_shiny_plot2.txt", data.table=F)
dth$pval_FGR10 <- as.numeric(format(dth$pval_FGR10, digits = 3, scientific = 5))
colnames(dth)[ncol(dth)] <- "category_FG"
#colnames(dth)[6:7] <- c("effect allele", "other allele")
colnames(dth)[2] <- "credible_set_variant"
# change columns for factors, for table filtering
for( ii in c(1,2,4:9,14:16)){
  dth[,ii] <- as.factor( dth[,ii])
}


# fd ids for interactive clicks
ids <- fread("data/fgids_for_phewas_app2.txt", data.table=F)

# User interface ----
ui <- ui  <- dashboardPage(
  dashboardHeader(title ="PheWAS: Migraine candidate variants in FinnGen R10"),
  dashboardSidebar(sidebarMenu(
    menuItem("Candidate variants", tabName = "dashboard", icon = icon("table")),
    menuItem("Functional variants", tabName = "functional_var", icon = icon("table")),
    menuItem("PIP > 0.1 variants", tabName = "pip_var", icon = icon("table"))
  )),
  dashboardBody(
    
    tabItems(
      # First tab content
      tabItem(tabName = "dashboard",
            h2("PheWAS: Migraine candidate variants in FinnGen R10"),
            fluidRow(title = "Table",
                      box(div(DT::dataTableOutput("mytable"),style = "font-size: 75%; width: 100%"), width = 12)),
            
            
           fluidRow(
               box(title = "Heatmap", footer= " ",
                InteractiveComplexHeatmapOutput("ht1", width1 = 700, height1 = 450, output_ui_float = TRUE,
                                                output_ui = htmlOutput("go_info")), width = 12)),
      
           fluidRow(
             box( selectInput("vars","Select locus",
                           choices = unique(dt[,1]), multiple=T, selected="CALCB" ),
   #          title = "Heatmap: selected locus",
             plotOutput("tile1"), width=12)),
      
       fluidRow(
            box(  selectInput("cat","Select category",
                              choices = unique(dt[,"category_FG"] ), multiple=T, selected="Neurological endpoints"),
                plotOutput("tile2"), width = 12))
      ),
  
    # Second tab content
      tabItem(tabName = "functional_var",
            h2("Functional variants among the credible sets"),
            
            fluidRow(title = "Table",
                     box(div(DT::dataTableOutput("mytablef"),style = "font-size: 75%; width: 90%"), width = 12)),

            
            fluidRow(
              box(title = "PheWAS: Functional variants in FinnGen R10", footer="Heatmap",
                  InteractiveComplexHeatmapOutput("ht1f", width1 = 700, height1 = 450, output_ui_float = TRUE,
                                                  output_ui = htmlOutput("go_infoF")), width = 12)),
            
            fluidRow( box(
              selectInput("varsf","Select locus",
                          choices = unique(dtf[,1]) ),
              plotOutput("tile1f"), width = 12)),
            
              fluidRow( box(
              selectInput("catf","Select category",
                          choices = unique(dtf[,"category_FG"] )),
              plotOutput("tile2f"), width = 12))
              
      ),
    # third
      tabItem(tabName = "pip_var",
            h2("Variants with PIP > 0.1 among the credible sets"),
            
            fluidRow(title = "Table",
                     box(div(DT::dataTableOutput("mytableh"),style = "font-size: 75%; width: 90%"), width = 12)),
            
            fluidRow(
              box(title = "PheWAS: PIP > 0.1 variants in FinnGen R10", footer="Heatmap",
                  InteractiveComplexHeatmapOutput("ht1h", width1 = 700, height1 = 450, output_ui_float = TRUE,
                                                  output_ui = htmlOutput("go_infoH")), width = 12)),
            
            fluidRow( box( 
              selectInput("varsh","Select locus",
                          choices = unique(dth[,1]) ),
              plotOutput("tile1h"), width = 12)),
           
             fluidRow( box( 
              selectInput("cath","Select category",
                          choices = unique(dth[,"category_FG"] )),
              selectInput("phenoh","Select phenotype",
                          choices = unique(dth[,"phenocode"]) ,multiple=T, selected="I9_CHD"),
              plotOutput("tile2h"),
              plotOutput("tile3h"), width = 12))
            
      )
  )
)
)

click_action = function(df, output) {
  output[["go_info"]] = renderUI({
    if(!is.null(df)) {
      rsid = rownames(matrix)[df$row_index]
      phenoc = colnames(matrix)[df$column_index]
      fgid = ids[ ids$rsid == rsid, "id"]
      
      oe = try(beta_f <- dt[ dt[,2] == rsid & dt[,"phenocode"] == phenoc, "beta_FGR10"] , silent = TRUE)
      if(inherits(oe, "try-error")) {
        beta_f = ""
      }
      oe = try(beta_m <- dt[ dt[,2] == rsid & dt[,"phenocode"] == phenoc, "beta_meta"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        beta_m = ""
      }
      oe = try(effect_allele <- dt[ dt[,2] == rsid & dt[,"phenocode"] == phenoc, "effect_allele"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        effect_allele = ""
      }
      oe = try(phenotype <- dt[ dt[,2] == rsid & dt[,"phenocode"] == phenoc, "phenotype"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        phenotype = ""
      }
      oe = try(catg <- dt[ dt[,2] == rsid & dt[,"phenocode"] == phenoc, "category_FG"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        catg = ""
      }
      oe = try(gene <- dt[ dt[,2] == rsid & dt[,"phenocode"] == phenoc, "nearest_genes"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        gene = ""
      }
      
      
      HTML(qq(
        "<div style='padding:5px 10px;border:1px solid black; width:600px; background-color:white;'>
        <h5>FinnGen R10 </h5>
        <p> <b>SNP</b>: <a href='https://r10.finngen.fi/variant/@{fgid}' target='_blank'>@{rsid}</a><br>
        <b>Phenocode</b>: <a href='https://r10.finngen.fi/pheno/@{phenoc}' target='_blank'>@{phenoc}</a> </p>
        <p> <b>Phenotype</b>: @{phenotype} <br>
        <b>Category</b>: @{catg} </p>
      <p> <b>Nearest genes</b>: @{gene} </p>
      <p> <b>Effect allele</b>: @{effect_allele} <br>
      <b>Log odds-ratio (FG)</b>: @{beta_f}<br>
      <b>Log odds-ratio (migraine)</b>: @{beta_m}</p>
    </div>"
      ))
    }
  })
}


click_actionF = function(df, output) {
  output[["go_infoF"]] = renderUI({
    if(!is.null(df)) {
      rsid = rownames(matrix2)[df$row_index]
      phenoc = colnames(matrix2)[df$column_index]
      fgid = ids[ ids$rsid == rsid, "id"]
      
      oe = try(beta_f <- dtf[ dtf[,2] == rsid & dtf[,"phenocode"] == phenoc, "beta_FGR10"] , silent = TRUE)
      if(inherits(oe, "try-error")) {
        beta_f = ""
      }
      oe = try(beta_m <- dtf[ dtf[,2] == rsid & dtf[,"phenocode"] == phenoc, "beta_meta"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        beta_m = ""
      }
      oe = try(effect_allele <- dtf[ dtf[,2] == rsid & dtf[,"phenocode"] == phenoc, "effect_allele"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        effect_allele = ""
      }
      oe = try(phenotype <- dtf[ dtf[,2] == rsid & dtf[,"phenocode"] == phenoc, "phenotype"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        phenotype = ""
      }
      oe = try(catg <- dtf[ dtf[,2] == rsid & dtf[,"phenocode"] == phenoc, "category_FG"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        catg = ""
      }
      oe = try(gene <- dtf[ dtf[,2] == rsid & dtf[,"phenocode"] == phenoc, "gene"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        gene = ""
      }
      
      
      HTML(qq(
        "<div style='padding:5px 10px;border:1px solid black; width:400px; background-color:white;'>
        <h5>FinnGen R10 </h5>
        <p> <b>SNP</b>: <a href='https://r10.finngen.fi/variant/@{fgid}' target='_blank'>@{rsid}</a><br>
        <b>Phenocode</b>: <a href='https://r10.finngen.fi/pheno/@{phenoc}' target='_blank'>@{phenoc}</a> </p>
        <p> <b>Phenotype</b>: @{phenotype} <br>
        <b>Category</b>: @{catg} </p>
      <p> <b>Nearest genes</b>: @{gene} </p>
      <p> <b>Effect allele</b>: @{effect_allele} <br>
      <b>Log odds-ratio (FG)</b>: @{beta_f}<br>
      <b>Log odds-ratio (migraine)</b>: @{beta_m}</p>
    </div>"
      ))
    }
  })
}


click_actionH = function(df, output) {
  output[["go_infoH"]] = renderUI({
    if(!is.null(df)) {
      rsid = rownames(matrix3)[df$row_index]
      phenoc = colnames(matrix3)[df$column_index]
      fgid = ids[ ids$rsid == rsid, "id"]
      
      oe = try(beta_f <- dth[ dth[,2] == rsid & dth[,"phenocode"] == phenoc, "beta_FGR10"] , silent = TRUE)
      if(inherits(oe, "try-error")) {
        beta_f = ""
      }
      oe = try(beta_m <- dth[ dth[,2] == rsid & dth[,"phenocode"] == phenoc, "beta_meta"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        beta_m = ""
      }
      oe = try(effect_allele <- dth[ dth[,2] == rsid & dth[,"phenocode"] == phenoc, "effect_allele"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        effect_allele = ""
      }
      oe = try(phenotype <- dth[ dth[,2] == rsid & dth[,"phenocode"] == phenoc, "phenotype"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        phenotype = ""
      }
      oe = try(catg <- dth[ dth[,2] == rsid & dth[,"phenocode"] == phenoc, "category_FG"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        catg = ""
      }
      oe = try(gene <- dth[ dth[,2] == rsid & dth[,"phenocode"] == phenoc, "nearest_genes"], silent = TRUE)
      if(inherits(oe, "try-error")) {
        gene = ""
      }
      
      
      HTML(qq(
        "<div style='padding:5px 10px;border:1px solid black; width:400px; background-color:white;'>
        <h5>FinnGen R10 </h5>
        <p> <b>SNP</b>: <a href='https://r10.finngen.fi/variant/@{fgid}' target='_blank'>@{rsid}</a><br>
        <b>Phenocode</b>: <a href='https://r10.finngen.fi/pheno/@{phenoc}' target='_blank'>@{phenoc}</a> </p>
        <p> <b>Phenotype</b>: @{phenotype} <br>
        <b>Category</b>: @{catg} </p>
      <p> <b>Nearest genes</b>: @{gene} </p>
      <p> <b>Effect allele</b>: @{effect_allele} <br>
      <b>Log odds-ratio (FG)</b>: @{beta_f}<br>
      <b>Log odds-ratio (migraine)</b>: @{beta_m}</p>
    </div>"
      ))
    }
  })
}


# Server logic
server <- function(input, output, session){
  observe({      
    cluster_rows = FALSE
    cluster_columns = FALSE
    row_order = NULL
    column_order = NULL
    row_km = NULL
    column_km = NULL
    matrix = as.matrix(hm[,2:ncol(hm)])
    row.names(matrix) <- hm[,1]
    ht = Heatmap(matrix,  name = "migraine",
                 cluster_rows = cluster_rows, cluster_columns = cluster_columns,
                 row_order = row_order, column_order = column_order, 
                 row_km = row_km, column_km = column_km, show_row_names = TRUE, show_column_names = TRUE)
    
    makeInteractiveComplexHeatmap(input, output, session, ht, "ht1",click_action = click_action)
  })
  
  observe({      
    cluster_rows = FALSE
    cluster_columns = FALSE
    row_order = NULL
    column_order = NULL
    row_km = NULL
    column_km = NULL
    matrix = as.matrix(fm[,2:ncol(fm)])
    row.names(matrix) <- fm[,1]
    ht = Heatmap(matrix,  name = "migraine",
                 cluster_rows = cluster_rows, cluster_columns = cluster_columns,
                 row_order = row_order, column_order = column_order, 
                 row_km = row_km, column_km = column_km, show_row_names = TRUE, show_column_names = TRUE)
    
    makeInteractiveComplexHeatmap(input, output, session, ht, "ht1f",click_action = click_actionF)
  })
  
  observe({      
    cluster_rows = FALSE
    cluster_columns = FALSE
    row_order = NULL
    column_order = NULL
    row_km = NULL
    column_km = NULL
    matrix = as.matrix(fmh[,2:ncol(fmh)])
    row.names(matrix) <- fmh[,1]
    ht = Heatmap(matrix,  name = "migraine",
                 cluster_rows = cluster_rows, cluster_columns = cluster_columns,
                 row_order = row_order, column_order = column_order, 
                 row_km = row_km, column_km = column_km, show_row_names = TRUE, show_column_names = TRUE)
    
    makeInteractiveComplexHeatmap(input, output, session, ht, "ht1h",click_action = click_actionH)
  })
  
  
  output$mytable = DT::renderDataTable({
    #Display table with select
    DT::datatable(dt, width = 7, rownames = FALSE, filter = 'top',
                  colnames= c("locus name","candidate rsid", "locus (lead rsid)","chr" ,"position","effect allele", "other allele" ,"nearest genes", "beta (FGR10)", "beta (migraine)", "pval (FGR10)", "pval (migraine)", "phenocode (FGR10)", "phenotype", "category (FG)"),
                  options = list(orderClasses = TRUE,
                                 lengthMenu = c(5, 10, 25, 50),
                                 pageLength = 5 ,
                                 autoWidth = TRUE,
                                 columnDefs = list(list(width = '300px', targets = c(13,14))),
                                 drawCallback= JS(
                                   'function(settings) {
                                     Shiny.bindAll(this.api().table().node());}')
                                 
                  ),selection='none',escape=F)
    
    
  } )
  
  output$mytablef = DT::renderDataTable({
    #Display table with select
    DT::datatable(dtf, width = 7, rownames = FALSE, filter = 'top',
                 colnames= c("locus name","functional variant rsid", "locus (lead rsid)","chr" ,"position","effect allele", "other allele","function (VEP)" ,"gene", "beta (FGR10)", "beta (migraine)", "pval (FGR10)", "pval (migraine)", "phenocode (FGR10)", "phenotype", "category (FG)"),
                                options = list(orderClasses = TRUE,
                                 lengthMenu = c( 5,10, 25, 50),
                                 pageLength = 5 ,
                                autoWidth = TRUE,
                                 columnDefs = list(list(width = '300px', targets = c(14,15))),
                                 drawCallback= JS(
                                   'function(settings) {
                                     Shiny.bindAll(this.api().table().node());}')
                  ),selection='none',escape=F)
    
    
  } )
  
  output$mytableh = DT::renderDataTable({
    #Display table with select
    DT::datatable(dth, width = 7, rownames = FALSE, filter = 'top',
                  caption = 'Table 1: Variants with PIP > 0.1 among the credible sets.',
                  colnames= c("locus name", "credible set variant rsid","PIP credible set variant" ,"locus (lead rsid)","chr" ,"position","effect allele", "other allele", "nearest genes", "beta (FGR10)", "beta (migraine)", "pval (FGR10)", "pval (migraine)", "phenocode (FGR10)", "phenotype", "category (FG)"),
                  options = list(orderClasses = TRUE,
                                 lengthMenu = c( 5,10, 25, 50),
                                 pageLength = 5 ,
                                 autoWidth = TRUE,
                                 columnDefs = list(list(width = '200px', targets = c(14,15))),
                                 drawCallback= JS(
                                   'function(settings) {
                                     Shiny.bindAll(this.api().table().node());}')
                  ),selection='none',escape=F)
    
    
  } )
  
  
#  observeEvent(input$vars, {
 #   hm2 <- hm[ hm[,1] ==  input$vars,]
  #  if(nrow(hm2) == 1){
   #   hm2[2,]<- hm2[1,]
  #  } 
  #  m2 <- as.matrix(hm2[,2:202])
  #  row.names(m2) <- hm2[,1]
  #  cluster_rows = FALSE
  #  cluster_columns = FALSE
  #  row_order = NULL
  #  column_order = NULL
   # row_km = NULL
  #  column_km = NULL
    
#    ht2 = Heatmap(m2, name = "migraine",
 #                 cluster_rows = cluster_rows, cluster_columns = cluster_columns,
  #                row_order = row_order, column_order = column_order, 
  #                row_km = row_km, column_km = column_km)
  #  makeInteractiveComplexHeatmap(input, output, session, ht2, "ht2")
#  })
  
  output$tile1 <- renderPlot({
    hm2 <- dt[ dt[,1] %in% input$vars,]
    ggplot(hm2,aes(x=phenocode, y=rsid_candidate, fill=beta_FGR10))+
      geom_tile(color = "white",
                lwd = 1.5,
                linetype = 1) +
     # scale_fill_gradientn(colors = hcl.colors(20, "RdBu"))+
      scale_fill_gradient2(low="red",high="Blue")+
      coord_fixed()+
      theme_classic()+
      xlab("")+
      ggtitle(paste(input$vars, collapse=", "))+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$tile2 <- renderPlot({
    hm2 <- dt[ dt[,"category_FG"] %in%  input$cat,]
    ggplot(hm2,aes(x=phenocode, y=rsid_candidate, fill=beta_FGR10))+
      geom_tile(color = "white",
                lwd = 1.5,
                linetype = 1) +
      scale_fill_gradient2(low="red",high="Blue")+
 #     scale_fill_gradientn(colors = hcl.colors(20, "RdBu"))+
      coord_fixed()+
      theme_classic()+
      xlab("")+
      ggtitle(paste(input$cat, collapse = ", "))+
      theme(axis.text.x = element_text(angle = 45, hjust=1))
    
    
  })
  
  output$tile1f <- renderPlot({
    hm2 <- dtf[ dtf[,1] %in%  input$varsf,]
    ggplot(hm2,aes(x=phenocode, y=functional_variant, fill=beta_FGR10))+
      geom_tile(color = "white",
                lwd = 1.5,
                linetype = 1) +
      scale_fill_gradient2(low="red",high="Blue")+
     # scale_fill_gradientn(colors = hcl.colors(20, "RdBu"))+
      coord_fixed()+
      theme_classic()+
      xlab("")+
      ggtitle(input$varsf)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$tile2f <- renderPlot({
    hm2 <- dtf[ dtf[,"category_FG"] %in%  input$catf,]
    ggplot(hm2,aes(x=phenocode, y=functional_variant, fill=beta_FGR10))+
      geom_tile(color = "white",
                lwd = 1.5,
                linetype = 1) +
      scale_fill_gradient2(low="red",high="Blue")+
     # scale_fill_gradientn(colors = hcl.colors(20, "RdBu"))+
      theme_classic()+
      xlab("")+
      ggtitle(input$catf)+
      theme(axis.text.x = element_text(angle = 45, hjust=1))
    
    
  })

  output$tile1h <- renderPlot({
    hm2 <- dth[ dth[,1] %in%  input$varsh,]
    ggplot(hm2,aes(x=phenocode, y=credible_set_variant, fill=beta_FGR10))+
      geom_tile(color = "white",
                lwd = 1.5,
                linetype = 1) +
      scale_fill_gradient2(low="red",high="Blue")+
     # scale_fill_gradientn(colors = hcl.colors(20, "RdBu"))+
      coord_fixed()+
      theme_classic()+
      xlab("")+
      ggtitle(input$varsh)+
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$tile2h <- renderPlot({
    hm2 <- dth[ dth[,"category_FG"] %in% input$cath,]
    ggplot(hm2,aes(x=phenocode, y=credible_set_variant, fill=beta_FGR10))+
      geom_tile(color = "white",
                lwd = 1.5,
                linetype = 1) +
      scale_fill_gradient2(low="red",high="Blue")+
     # scale_fill_gradientn(colors = hcl.colors(20, "RdBu"))+
      coord_fixed()+
      theme_classic()+
      xlab("")+
      ggtitle(input$cath)+
      theme(axis.text.x = element_text(angle = 45, hjust=1))
    
    
  })
  
  output$tile3h <- renderPlot({
    hm2 <- dth[ dth[,"phenocode"] %in% input$phenoh,]
    ggplot(hm2,aes(y=phenocode, x=credible_set_variant, fill=beta_FGR10))+
      geom_tile(color = "white",
                lwd = 1.5,
                linetype = 1) +
      scale_fill_gradient2(low="red",high="Blue")+
      # scale_fill_gradientn(colors = hcl.colors(20, "RdBu"))+
      coord_fixed()+
      theme_classic()+
      ylab("")+
      ggtitle(paste(input$phenoh, collapse=", "))+
      theme(axis.text.x = element_text(angle = 45, hjust=1))
    
    
  })
  
  
  output$heatmap <- renderPlot({
        data <- hm
        cols <- colorRampPalette(brewer.pal(8,"RdBu"))(3)
        heatmap(as.matrix(data[,2:202]), Colv = NA, Rowv = NA, labRow =data$rsid, col=cols, cexRow=0.7, cexCol = 0.7, scale="none")
    
      }, height = 800, width = 800)
    
  
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

