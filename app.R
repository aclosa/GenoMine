# shiny app to DE analsys and exploration data

# libraries
library(shiny)
library(bs4Dash)
library(summarytools)
library(leaflet)
library(DT)
library(dplyr)
library(shinycssloaders)
library(DESeq2)
library(limma)
library(edgeR)
library(RColorBrewer)
library(plotly)
library(shinyWidgets)
library(rlang)
library(BiocParallel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(gridExtra)
library(shinyjqui)
# options

options(shiny.maxRequestSize=1000*1024^2, spinner.type = 5, spinner.color = "#3D9970")

# variables

de.model = c("DESeq2","Limma")
norm.method = c("TMM","TMMwsp", "RLE", "upperquartile", "none")
go.opt = c("BP","MF","CC","ALL")

# functions
bs_card <- function(font_color, bg_color, header_text, title_text, body_text, switchtab, go_path) {
  
  shiny::div(
    class=glue::glue("card text-{font_color} bg-{bg_color} mb-3"), 
    style="max-width: 18rem;margin:10px", 
    
    shiny::div(h3(
      class="card-header", 
      header_text
    )
    ), 
    
    shiny::div(
      class="card-body", 
      
      shiny::h2(
        class="card-title", 
        title_text
      ), 
      
      shiny::p(
        class="card-text", 
        body_text
      ),
      
      shiny::actionButton(
        switchtab,
        go_path
      )
      
    )
    
  )
  
}

# app

shinyApp(
  
  ui = dashboardPage(scrollToTop = TRUE, fullscreen = TRUE,
                     dashboardHeader(title = dashboardBrand(title = "GENOMINE",
                                                            color = "olive",
                                                            image=""),
                                     fixed = TRUE, 
                                     border = TRUE),
                     dashboardSidebar(width = 500, status = "olive",
                                      sidebarMenu(id = "sidebar",
                                                  menuItem("Home",tabName = "home",icon=icon("home")),
                                                  menuItem("Data",tabName = "data",icon=icon("file")),
                                                  menuItem("Clustering",tabName = "cluster",icon=icon("box")),
                                                  menuItem("DE",tabName = "difexp",icon=icon("filter")),
                                                  menuItem("PathWay",tabName = "pathway",icon=icon("wrench")),
                                                  menuItem("User Guidelines",tabName = "userInfo",icon = icon("list")),
                                                  menuItem("Contact",tabName = "contact",icon = icon("user"))
                                      )
                     ),
                     dashboardBody(status = "olive",tabItems(
                       tabItem(tabName = "home",
                               fluidRow(
                                 div(class = "row row-cols-1 row-cols-md-2 g-4",style="max-width: 40rem; margin: auto",
                                     bs_card(
                                       font_color = "black", 
                                       bg_color = "light", 
                                       header_text = tags$b("STEP 1"), 
                                       title_text = "Upload files:", 
                                       body_text = "Upload your counts and clinical information tables.",
                                       switchtab = 'switchData',
                                       go_path = "Go Data"
                                     ),
                                     bs_card(
                                       font_color = "black", 
                                       bg_color = "light", 
                                       header_text = tags$b("STEP 2"), 
                                       title_text = "Clustering:", 
                                       body_text = "Explore your data with clustering visualizations.",
                                       switchtab = 'switchClustering',
                                       go_path = "Go Clustering"
                                     ),
                                     bs_card(
                                       font_color = "black", 
                                       bg_color = "light", 
                                       header_text = tags$b("STEP 3"), 
                                       title_text = "Diferential Expression:", 
                                       body_text = "Diferential expression analysis visualitzation.",
                                       switchtab = 'switchDE',
                                       go_path = "Go DE"
                                     ),
                                     bs_card(
                                       font_color = "black", 
                                       bg_color = "light", 
                                       header_text = tags$b("STEP 4"), 
                                       title_text = "Pathway Analysis:", 
                                       body_text = "Gene Onthoplogy enrichment and pathway analysis.",
                                       switchtab = 'switchPath',
                                       go_path = "Go PathWay"
                                     )
                                 ),
                               )),
                       tabItem(tabName = "data",
                               fluidRow( tags$head(tags$style(".progress-bar{background-color:#3D9970;}")),
                                 box(title = tagList(h1("Data files",style="display: inline-block; position: relative; left: 50%; transform: translateX(-50%)")), width = 12,collapsible=FALSE,
                                     
                                     column(6, align= "center",fileInput("exp.file","Upload a counts expression file (tsv)",buttonLabel = "Upload...",multiple = FALSE),
                                            style="min-width:200px; float:left; margin: auto; margin-right:6px; align:center"),
                                     column(6, align = "center",fileInput("info.file","Upload a clinical info file (tsv)",buttonLabel = "Upload...",multiple = FALSE),
                                            tags$hr(),
                                            #checkboxInput("header", "Header", TRUE),
                                            prettySwitch("header", "Header", TRUE),
                                            style="min-width:200px; float:left; margin: auto; align: center")
                                 ),
                                 tabBox(
                                   title = tagList(shiny::icon("tasks"),"Tables"), side = "right",width = 12, status = "olive",solidHeader = FALSE,type = "tabs",
                                   tabPanel("Summary Info", shinycssloaders::withSpinner(htmlOutput("profileSummary"))),
                                   tabPanel( "Table Info", shinycssloaders::withSpinner(DT::dataTableOutput("contents"))),
                                   tabPanel( "Counts", shinycssloaders::withSpinner(DT::dataTableOutput("counts")))
                                 )
                               )
                       ),
                       
                       tabItem( tabName = "cluster",
                          fluidRow(
                         box(title = tagList(shiny::icon("gear"),"Parameters"), width = 4,
                             selectInput("norm.method", label = "",""),
                             numericInput("lcpm.filt",label="Average lcpm filter",value = 0),
                             selectInput("covariable", label = "Select covariable",""),
                             actionButton("run","Run")
                           
                         ),
                         tabBox(
                           title = tagList(shiny::icon("tasks"),"Outputs"), side = "right", width = 8, status = "olive",solidHeader = FALSE,type = "tabs",
                           tabPanel("MDS",
                                    br(),
                                    DT::dataTableOutput("scores"),
                                    shinycssloaders::withSpinner(plotlyOutput("mds.plot",inline = TRUE))
                           )
                       ))),
                       
                       tabItem( tabName = "difexp",
                                fluidRow(
                                box(title = tagList(shiny::icon("gear"),"Parameters"), width = 4,
                                    selectInput("model", label = "",""),
                                    numericInput("count.filt",label="Count filter",value = 10),
                                    selectInput("covariableDE", label = "Select covariable",""),
                                    actionButton("runDE","Run")
                                    
                                ),
                                tabBox(
                                  title = tagList(shiny::icon("tasks"),"Outputs"), side = "right", width = 8, status = "olive",solidHeader = FALSE,type = "tabs",
                                  tabPanel("MA plot",
                                           br(),
                                           fluidRow(valueBoxOutput("modelBox"),
                                                    valueBoxOutput("filterBox"),
                                                    valueBoxOutput("covariableBox")
                                                    ),
                                           shinycssloaders::withSpinner(plotOutput("ma.plot",inline = TRUE)),
                                           DT::dataTableOutput("DE.table")
                                  ),
                                  tabPanel("DE genes",
                                           br(),
                                           selectInput("genes",label="Select gene",""),
                                           shinycssloaders::withSpinner(plotlyOutput("box.plot",inline = TRUE)),
                                           DT::dataTableOutput("DE.table2")
                                           )
                                )
                       )),
                       
                       tabItem(tabName = "pathway",
                               fluidRow( tags$head(tags$style(".progress-bar{background-color:#3D9970;}")),
                                         box(title = tagList(shiny::icon("gear"),"Parameters"), width = 4,collapsible=TRUE,
                                             
                                             fileInput("path.file","Upload DE results file (csv)",buttonLabel = "Upload...",multiple = FALSE),
                                             splitLayout(cellWidths = c("50%", "50%"),
                                                         numericInput("logFC",label="|LogFC|",value="0.5"),
                                                         #helpText("Write drug name tested"),
                                                         numericInput("padj",label="Padj",value = "0.05"),
                                                         #helpText("Write concentration from drug tested (uM)")
                                             ),
                                             selectInput("goOpt",label="Select Gene Onthology",choices = ""),
                                             actionButton("runGO","Run")
                                         ),
                                         tabBox(
                                           title = tagList(shiny::icon("tasks"),"Outputs"), side = "right", width = 8, status = "olive",solidHeader = FALSE,type = "tabs",
                                           tabPanel("GO enrichment",
                                                    br(),
                                                    shinycssloaders::withSpinner(plotOutput("go.plot",width = '800px',height = '1000px',inline = TRUE)),
                                                    DT::dataTableOutput("go.table")
                                           )
                                         )
                               )
                         
                       ),
                       
                       tabItem(tabName = "contact",
                               box(title = h4(tagList(shiny::icon("globe"),strong("Genome Sciences and Cancer Division"))),collapsible = FALSE,width = 12,
                                   p(
                                     h6("The John Curtin School of Medical Research"),
                                     h6("131 Garran Road"),
                                     h6("The Australian National University"),
                                     h6("Acton ACT 2601")),
                                   leafletOutput("jcsmrMap")),
                               hr(),
                               box(title = h4(tagList(shiny::icon("bell"),strong("Contact Us"))), collapsible = FALSE, width = 12,
                                   p(em(tagList(shiny::icon("user"),"Contact : adria.closamosquera@anu.edu.au"))),
                                   p(a(tagList(shiny::icon("github"),"GitHub Drug Screening"),href="https://github.com/comprna/drug_screening"))
                               )
                       )
                     )),
  ),
  
  
  
  server = function(input, output, session) {
    
    # observe Data
    observeEvent(input$switchData, {
      newtab <- switch(input$sidebar, "home" = "data","data" = "home")
      updateTabItems(session, "sidebar", newtab)
    })
    # observe Cluster
    observeEvent(input$switchClustering, {
      newtab <- switch(input$sidebar, "home" = "cluster","cluster" = "home")
      updateTabItems(session, "sidebar", newtab)
    })
    # observe DE
    observeEvent(input$switchDE, {
      newtab <- switch(input$sidebar, "home" = "difexp","difexp" = "home")
      updateTabItems(session, "sidebar", newtab)
    })
    # observe DE
    observeEvent(input$switchPath, {
      newtab <- switch(input$sidebar, "home" = "pathway","pathway" = "home")
      updateTabItems(session, "sidebar", newtab)
    })
    
    # map location
    output$jcsmrMap = renderLeaflet({
      m = leaflet() %>%
        addTiles() %>%  # Add default OpenStreetMap map tiles
        addMarkers(lng=149.1149376416087, lat=-35.28187528455249, popup="JCSMR: The John Curtin School of Medical Research")
      m
    })
    
    # laoding files and data summary
    
    info.data = reactive({
      req(input$info.file)
      inFile = input$info.file
      if (is.null(inFile)){
        return(NULL)
      } else {
        info.df = read.table(inFile$datapath, header = input$header, sep = "\t", row.names = 1)
      }
      info.df
    })
    
    info.cov = reactive({
      req(input$info.file)
      inFile = input$info.file
      if (is.null(inFile)){
        return(NULL)
      } else {
        info.df = read.table(inFile$datapath, header = input$header, sep = "\t",row.names = 1)
      }
      colnames(info.df)
    })
    
    counts.data = reactive({
      req(input$exp.file)
      inFile = input$exp.file
      if (is.null(inFile)){
        return(NULL)
      } else {
        counts.df = read.table(inFile$datapath, header = input$header, sep = "\t")
      }
      counts.df
    })
    
    path.data = reactive({
      req(input$path.file)
      inFile = input$path.file
      if (is.null(inFile)){
        return(NULL)
      } else {
        path.df = read.csv(inFile$datapath, header = T, sep = ",", row.names = 1)
      }
      path.df
    })
    
    output$profileSummary <- renderUI({
      
      inFile <- input$info.file
      
      if (is.null(inFile)) {
        return(NULL)
      } else {
        tmp <- read.table(inFile$datapath, header = input$header, sep = "\t")  
        print(dfSummary(tmp),                                    
              method = "render",                                 
              Data.frame = inFile[[1]])                          
      }
    })
    
    output$contents = DT::renderDataTable({
      
      infodf = info.data()
      #lst.syn.df = as.data.frame(do.call(rbind,all.synergy))
      infodf %>% mutate_if(is.numeric,round,digits = 4) %>%
        datatable(extensions = 'Buttons',
                  options = list(dom = 'Blfrtip',
                                 buttons = c('copy', 'csv', 'excel', 'print'),
                                 scrollX = TRUE,
                                 lengthMenu = list(c(10,25,50,-1),
                                                   c(10,25,50,"All"))))
    })
    
    output$counts = DT::renderDataTable({
      
      countsdf = counts.data()
      #lst.syn.df = as.data.frame(do.call(rbind,all.synergy))
      countsdf %>% mutate_if(is.numeric,round,digits = 4) %>%
        datatable(extensions = 'Buttons',
                  options = list(dom = 'Blfrtip',
                                 buttons = c('copy', 'csv', 'excel', 'print'),
                                 scrollX = TRUE,
                                 lengthMenu = list(c(10,25,50,-1),
                                                   c(10,25,50,"All"))))
    })
    
    # plotting data normalitzation for clustering
    
     runModel = eventReactive(input$run, {
      
      counts = counts.data()
      avgexp = aveLogCPM(counts)
      mask = avgexp > input$lcpm.filt
      m.dge <- DGEList(counts, genes=as.data.frame(rownames(counts)))
      
      dge.filt = m.dge[mask, , keep.lib.sizes=FALSE]
      
      logCPM = edgeR::cpm(dge.filt, log=TRUE, prior.count = 0.5)
      
      x1 = calcNormFactors(dge.filt, method = input$norm.method)
      x2 = x1
      x2$samples$norm.factors = 1
      x2$counts[,1] = ceiling(x2$counts[,1]*0.05)
      x2$counts[,2] = x2$counts[,2]*5
      
      lcpm = edgeR::cpm(x2, log=TRUE)
      
      x2 = calcNormFactors(x2, method = input$norm.method)  
      lcpm = edgeR::cpm(x2, log=TRUE)
      
      d <- dist(t(lcpm))
      fit <- cmdscale(d, eig=TRUE, k=2)
      scores <- data.frame(fit$points[,1],fit$points[,2],colnames(lcpm))
      
      #pol1.info = read.table("/media/adria/Seagate Exp/Nadine/info/info_pol1_batch2.txt",header = T, sep = "\t")
      infodf = info.data()
      scores = cbind(scores, infodf[,colnames(infodf) == input$covariable])
      colnames(scores) = c("Coordinate_1","Coordinate_2","Sample_ID",input$covariable)
      return(scores)
     })
      
     output$scores = DT::renderDataTable({
       
       runModel() %>% mutate_if(is.numeric,round,digits = 4) %>%
        datatable(extensions = 'Buttons',
                  options = list(dom = 'Blfrtip',
                                 buttons = c('copy', 'csv', 'excel', 'print'),
                                 scrollX = TRUE,
                                 lengthMenu = list(c(10,25,50,-1),
                                                   c(10,25,50,"All"))))
    })
    
    
     output$mds.plot = renderPlotly({
      
      plot_ly(data = runModel(), x = ~Coordinate_1, y = ~Coordinate_2, color = ~base::get(input$covariable)) %>%
        layout(data = runModel(),
               title = "",
               xaxis = list(title = "Coordinate 1"),
               yaxis = list(title = "Coordinate 2"))
    })
    
    # Run DE expression models
     
     # render valueBox
     renderModel = eventReactive(input$runDE, {
       valueBox(
         subtitle = "Model", value = HTML(paste0(h4(input$model))), icon = icon("box"),
         color = "info", gradient = TRUE
       )
     })
     renderFilter = eventReactive(input$runDE, {
       valueBox(
         subtitle = "Filter", value = HTML(paste0(h4(input$count.filt))), icon = icon("filter"),
         color = "teal", gradient = TRUE
       )
     })
     renderCov = eventReactive(input$runDE, {
       valueBox(
         subtitle = "Covariable", value = HTML(paste0(h4(input$covariableDE))), icon = icon("layer-group"),
         color = "purple", gradient = TRUE
       )
     })
       
    output$modelBox = renderInfoBox({
      renderModel()
    })
    output$filterBox = renderInfoBox({
      renderFilter()
    })
    output$covariableBox = renderInfoBox({
      renderCov()
    })
     
    runDE = eventReactive(input$runDE, {
      counts = counts.data()
      infoDE = info.data()
      ddsTxi.CX <- DESeqDataSetFromMatrix(round(as.matrix(counts)),
                                          colData = infoDE,
                                          design = as.formula(paste("~",input$covariableDE)))
      # prefiltering
      keep <- rowSums(counts(ddsTxi.CX)) >= input$count.filt
      ddsTxi.CX <- ddsTxi.CX[keep,]
      
      # deseq analysis
      register(MulticoreParam(6))
      ddsTxi.CX  <- DESeq(ddsTxi.CX,parallel = TRUE) 
      #resTC <- results(ddsTxi.CX)
      #resTC$GeneID = rownames(resTC)
      
      #head(resTC[order(resTC$padj),], 4)
      
      # normal counts
      #vsd <- vst(ddsTxi.CX , blind=FALSE)
      #rld <- rlog(ddsTxi.CX , blind=FALSE)
      #return(as.data.frame(assay(vsd)))
      #return(resTC)
      return(ddsTxi.CX)
    })
    
    # ma plot and table 
    output$ma.plot = renderPlot({
      
      resTC <- results(runDE())
      resTC$GeneID = rownames(resTC)
      
      DESeq2::plotMA(resTC, ylim=c(-5,5))
      
    }, height = 400,600)
    
    output$DE.table = DT::renderDataTable({
      
      resTC <- results(runDE())
      resTC$GeneID = rownames(resTC)
      
      topGenes <- head(order(resTC$padj),table(resTC$padj < 0.05)[2])
      dt.topGenes = resTC[topGenes,]
      #datatable(signif(as.data.frame(dt.topGenes),digits = 4), options = list(pageLength = 10))
      signif(as.data.frame(dt.topGenes[,c(1:6)]),digits = 4) %>%
        datatable(extensions = 'Buttons',
                  options = list(dom = 'Blfrtip',
                                 buttons = c('copy','csv','excel','print'),
                                 scrollX = TRUE,
                                 lengthMenu = list(c(10,25,50,-1),
                                                   c(10,25,50,"All"))))
    })
    
    # boxplot per gene and table
    output$box.plot = renderPlotly({
      
      vsd <- vst(runDE() , blind=FALSE)
      vsd.df = assay(vsd)
      
      gene.x = vsd.df[rownames(vsd.df) == input$genes,]
      cov.x = info.data()[,colnames(info.data()) == input$covariableDE]
      
      #gene.x = vsd.df[rownames(vsd.df) == "TP53",]
      #cov.x = info.data()[,colnames(info.data()) == input$covariableDE]
      
      #gene.df = data.frame(base::get(input$genes),base::get(input$covariableDE))
      gene.df = data.frame(gene.x,cov.x)
      
      gene.df %>%
        plot_ly(
          #x = ~base::get(input$covariableDE),
          #y = ~base::get(input$genes),
          #split = ~base::get(input$covariableDE),
          x = ~cov.x,
          y = ~gene.x,
          split = ~cov.x,
          type = 'violin',
          box = list(
            visible = T
          ),
          meanline = list(
            visible = T
          )
        ) %>%
        layout(
          xaxis = list(
            title = input$covariableDE
          ),
          yaxis = list(
            title = input$genes
          )
        )
    })
    
    output$DE.table2 = DT::renderDataTable({
      
      resTC <- results(runDE())
      resTC$GeneID = rownames(resTC)
      
      topGenes <- head(order(resTC$padj),table(resTC$padj < 0.05)[2])
      dt.topGenes = resTC[topGenes,]
      #datatable(signif(as.data.frame(dt.topGenes),digits = 4), options = list(pageLength = 10))
      signif(as.data.frame(dt.topGenes[,c(1:6)]),digits = 4) %>%
        datatable(extensions = 'Buttons',
                  options = list(dom = 'Blfrtip',
                                 buttons = c('copy','csv','excel','print'),
                                 scrollX = TRUE,
                                 lengthMenu = list(c(10,25,50,-1),
                                                   c(10,25,50,"All"))))
    })
    
    # go enrichment analysis 
    
    runGO = eventReactive(input$runGO, {
      
      path.df = path.data()
      
      list.genes = rownames(path.df[abs(path.df$log2FoldChange) > 0.5 & path.df$padj < 0.05,])
      
      gene.df <- bitr(list.genes, fromType = "SYMBOL",
                      toType = c("ENSEMBL", "ENTREZID"),
                      OrgDb = org.Hs.eg.db)
      gene.uni = bitr(rownames(path.df), fromType = "SYMBOL",
                      toType = c("ENSEMBL","ENTREZID"),
                      OrgDb = org.Hs.eg.db)
      
      enrich.list = path.df$log2FoldChange*(-log10(path.df$pvalue))
      names(enrich.list) = rownames(path.df)
      
      ego2 = gseGO(geneList = sort(enrich.list,decreasing = TRUE),
                   OrgDb = org.Hs.eg.db,
                   keyType = "SYMBOL",
                   ont = input$goOpt,
                   minGSSize = 10,
                   maxGSSize = 500,
                   pvalueCutoff = 1,
                   verbose = FALSE)
      return(ego2)
    })
    
    output$go.plot = renderPlot({
      
      go.model = runGO()
      
      d1 = enrichplot::dotplot(go.model,color="pvalue",showCategory = 30)
      ema1 = enrichplot::pairwise_termsim(go.model)
      ema1 = emapplot(ema1,showCategory = 30,color="pvalue")
      
      ema2 = enrichplot::pairwise_termsim(go.model)
      tr1 = enrichplot::treeplot(ema2,color = "pvalue",showCategory = 30)
      
      ga = grid.arrange(d1,ema1,tr1,ncol=1,nrow=3)
      
      grid.draw(ga)
      
    }, height = 1200,600)
    
    output$go.table = DT::renderDataTable({
      
      go.model = runGO()
      
      dt.ego2 = as.data.frame(go.model)
      dt.ego2$enrichmentScore = round(dt.ego2$enrichmentScore,digits = 4)
      dt.ego2$NES = round(dt.ego2$NES,digits = 4)
      dt.ego2$pvalue = round(dt.ego2$pvalue,digits = 4)
      dt.ego2$p.adjust = round(dt.ego2$p.adjust,digits = 4)
      dt.ego2$qvalues = round(dt.ego2$qvalues,digits = 4)
      
      dt.ego2 %>%
        datatable(extensions = 'Buttons',
                  options = list(dom = 'Blfrtip',
                                 buttons = c('copy','csv','excel','print'),
                                 scrollX = TRUE,
                                 lengthMenu = list(c(10,25,50,-1),
                                                   c(10,25,50,"All"))))
    })
  
    
    # observe options keep at the end of script
    observe({
      
      updateSelectInput(session,"goOpt",
                        choices = go.opt,
                        selected = go.opt[1])
      
      updateSelectInput(session,"model",
                        label = "Model",
                        choices= de.model,
                        selected = de.model[1])
      
      updateSelectInput(session,"norm.method",
                        label = "Norm Method",
                        choices= norm.method,
                        selected = norm.method[1])
    
      updateSelectInput(session,"covariable",
                       choices= info.cov(),
                       selected = info.cov()[1])
      
      updateSelectInput(session,"covariableDE",
                        choices= info.cov(),
                        selected = info.cov()[1])
      
      updateSelectInput(session,"genes",
                        #choices= colnames(info.data()),
                        choices= results(runDE(), tidy = TRUE)$row,
                        selected = results(runDE(), tidy = TRUE)$row[1])
                       })
    
  }
  
)
