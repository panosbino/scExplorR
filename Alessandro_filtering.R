library(shiny)
library(bslib)
library(biomaRt)
library(shinythemes)
library(shinyWidgets)
library(Seurat)
library(shinyFiles)
library(fs)

source("~/scExplorR/filtering_functions/Utils.R")



# Set the maximum request size to 300 MB
options(shiny.maxRequestSize = 300 * 1024^2)




ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css")
  ),

  theme = shinytheme("sandstone"),
  titlePanel(
    div(
      style = "text-align: center; position: relative;",
      HTML('<i class="fas fa-bomb" style="position: absolute; left: 10px;"></i>'),
      h1("Single cell EXPLODER!!!", style = "color: #2c3e50; display: inline-block;"),
      HTML('<i class="fas fa-bomb" style="position: absolute; right: 10px;"></i>'),
      p("For you who can't use R eheh", style = "font-size: 18px; color: #7f8c8d;")
    )
  ),

  mainPanel(
    tabsetPanel(
      type = "pills",

      tabPanel("Filtering",
               div(
                 style = "padding: 20px; border: 1px solid #ddd; border-radius: 10px;",
                 fluidRow(
                   column(
                     12,
                     shinyDirButton("directory", "Folder select", "Please select a folder"),
                     br(), br(),
                     selectizeInput("ensembl", "Ensembl version", choices = get_possible_ensembl_versions()),
                     selectizeInput("organism", "Organism", choices = get_possible_organisms(109)),
                     actionBttn("btn_filter_go", "Explode!", style = "jelly", color = "danger"),
                     br(), br(),
                     plotOutput("plot_output"),
                     br(), br(),
                     numericInput('gene_max', 'Filter gene max', 3, min = 0, max = 100000),
                     numericInput('gene_min', 'Filter gene min', 3, min = 0, max = 100000),
                     numericInput('count_max', 'Filter UMI max', 3, min = 0, max = 100000),
                     numericInput('count_min', 'Filter UMI min', 3, min = 0, max = 100000),
                     numericInput('mt_max', 'Filter mito max', 3, min = 0, max = 100000),
                     br(), br(),
                     actionBttn("btn_filter_go_again", "Explode...again!", style = "jelly", color = "danger"),
                     plotOutput("plot_output_again"),
                     downloadBttn("button_download_plot", "Download Plot", style = "jelly", color = "success"),
                     br(), br(),
                     downloadBttn("button_download_data", "Download filtered data", style = "jelly", color = "success")
                   )
                 )
               )
      ),tabPanel("Normalization",
                 div(
                   style = "padding: 20px; border: 1px solid #ddd; border-radius: 10px;",
                   fluidRow(
                     column(
                       12,
                       fileInput("file_input_filtered", "File Input:", buttonLabel = "Browse...", placeholder = "No file selected"),
                       selectInput("method", "Method", choices = c("cpm", "log", "scran", "asinh", "log_geom")),
                       actionBttn("btn_norm_go", "Explode!", style = "jelly", color = "danger"),
                       br(), br(),
                       downloadBttn("button_download_norm", "Download", style = "jelly", color = "success")
                     )
                   ),
                   fluidRow(
                     column(
                       12,
                       verbatimTextOutput("text_output_three"),
                       plotlyOutput("plot_output_PCA"),
                       plotOutput("plot_output_UMAP")
                     )
                   )
                 )
      ),
      tabPanel("Clustering",
               div(
                 style = "padding: 20px; border: 1px solid #ddd; border-radius: 10px;",
                 fluidRow(
                   column(
                     12,
                     sliderInput("clustering_resolution", "Clustering resolution range", min = 0, max = 2, value = c(0, 2), step = 0.1),
                     actionBttn("btn_cluster_go", "Explode!", style = "jelly", color = "danger")
                   )
                 ),
                 fluidRow(
                   column(
                     12,
                     plotOutput("clustered_UMAP")
                   )
                 )
               )
      )
    )
  )
)






############################

server <- function(input, output, session) {




#create reactuve value object
  reactive_sobj <- reactiveValues(sobj = NULL)


  volumes <- c(Home = fs::path_home(),
               "R Installation" = R.home(),
               getVolumes()())

  #set directory selector
  shinyDirChoose(input,
                 "directory",
                 roots = volumes,
                 session = session,
                 restrictions = system.file(package = "base"),
                 allowDirCreate = FALSE)

  #extract directory path of the submitted folder
  getdata <- reactive({
    req(input$directory)
    parseDirPath(volumes, input$directory)
  })




#observe event for plotting data
  observeEvent(input$btn_filter_go, {
    req(input$directory, input$ensembl, input$organism)

    data_rep <- getdata()
    ensembl_version <- input$ensembl
    organism <- input$organism

      sobj <- load_and_annotate_recipe_1(data_rep, organism, ensembl_version)

      reactive_sobj$sobj <- sobj

      output$plot_output <- renderPlot({
        show_QC_plots(reactive_sobj$sobj)


        })

    })



  observeEvent(input$btn_filter_go_again, {
    req(reactive_sobj$sobj, input$gene_min, input$gene_max, input$mt_max, input$count_max, input$count_min)




    output$plot_output_again <- renderPlot({
      preview_filtered_sobj(reactive_sobj$sobj, min_genes = input$gene_min, max_genes = input$gene_max, max_pct_mito = input$mt_max ,
                            max_UMIs = input$count_max, min_UMIs = input$count_min)


    })
  })



  output$button_download_plot <- downloadHandler(
    filename = function() {
      paste("QC_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = last_plot())
    }
  )

  output$button_download_data <- downloadHandler(
    filename = function() {
      "processed_data.rds"
    },
    content = function(file) {
      processed_data <- commit_filtered_sobj(
        reactive_sobj$sobj,
        min_genes = input$gene_min,
        max_genes = input$gene_max,
        max_pct_mito = input$mt_max,
        max_UMIs = input$count_max,
        min_UMIs = input$count_min
      )
      saveRDS(processed_data, file)
    }
  )
  observeEvent(input$btn_norm_go, {
    req(input$file_input_filtered)
    req(input$method)

    print(input$file_input_filtered$name)
    print(input$method)

    # Read the uploaded file
    sc_filtered_object <- readRDS(input$file_input_filtered$datapath)

    print("File loaded successfully")

    # Normalize the data
    normalized_seurat_object <- normalize_all(sc_input_object = sc_filtered_object, method = input$method)
    print("Normalization complete")

    # Perform dimensionality reduction
    rv$dim_reduced_seurat_obj <- dim_reduction(normalized_seurat_object)
    print("Dimensionality reduction complete")

    # Generate PCA and UMAP plots
    PCA_plot <- plot_PCA(dim_reduced_seurat_obj = rv$dim_reduced_seurat_obj)
    UMAP_plot <- plot_UMAP(dim_reduced_seurat_obj = rv$dim_reduced_seurat_obj)

    # Render the plots
    output$plot_output_PCA <- renderPlotly({
      PCA_plot
    })

    output$plot_output_UMAP <- renderPlot({
      UMAP_plot
    })
  })

  observeEvent(input$btn_cluster_go, {
    req(rv$dim_reduced_seurat_obj)

    # Perform clustering
    clustering_res <- cluster(rv$dim_reduced_seurat_obj, clustering_resolution = input$clustering_resolution)

    # Assuming clustered_UMAP is the UMAP plot
    output$clustered_UMAP <- renderPlot({
      clustering_res$UMAP
    })
  })
}





shinyApp(ui = ui, server = server)





