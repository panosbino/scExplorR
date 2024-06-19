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
                     shinyDirButton("directory_filt", "Folder select", "Please select a folder"),
                     verbatimTextOutput("filepaths_filt"),
                     br(), br(),
                     selectizeInput("ensembl", "Ensembl version", choices = get_possible_ensembl_versions()),
                     selectizeInput("organism", "Organism", choices = get_possible_organisms(109)),
                     sliderInput("range", "Genes filtering range",min = 0, max = 100, value = c(0,100)),
                     sliderInput("range", "Counts filtering range",min = 0, max = 100, value = c(0,100)),
                     sliderInput("range", "Mitocondria filtering range",min = 0, max = 100, value = c(0,100)),
                     actionBttn("btn_filter_go", "Explode!", style = "jelly", color = "danger"),
                     #progressBar(id = "progress"),
                     br(), br(),
                     plotOutput("plot_output"),
                     downloadBttn("button_download_plot", "Download Plot", style = "jelly", color = "success"),
                     br(), br(),
                     downloadBttn("button_download_data", "Download filtered data", style = "jelly", color = "success")
                   )
                 ),
                 fluidRow(
                   column(
                     12,
                     verbatimTextOutput("file_output")
                   )
                 )
               )
      ),

      tabPanel("Normalisation",
               div(
                 style = "padding: 20px; border: 1px solid #ddd; border-radius: 10px;",
                 fluidRow(
                   column(
                     12,
                     shinyDirButton("directory_norm", "Folder select", "Please select a folder"),
                     verbatimTextOutput("filepaths_norm"),
                     br(), br(),
                     selectInput("method", "Method", choices = c("Cpm", "Log", "Scran", "Asinh", "Log_geom")),
                     actionBttn("btn_norm_go", "Explode!", style = "jelly", color = "danger"),
                     br(), br(),
                     downloadBttn("button_download_norm", "Download", style = "jelly", color = "success")
                   )
                 ),
                 fluidRow(
                   column(
                     12,
                     verbatimTextOutput("text_output_norm"),
                     plotOutput("plot_output_norm")
                   )
                 )
               )
      ),

      tabPanel("Visualization",
               div(
                 style = "padding: 20px; border: 1px solid #ddd; border-radius: 10px;",
                 fluidRow(
                   column(
                     12,
                     fileInput("file_input_three", "File Input:", buttonLabel = "Browse...", placeholder = "No file selected"),
                     textInput("target_three", "Target:", placeholder = "Enter target..."),
                     actionBttn("btn_three_go", "Submit", style = "jelly", color = "danger")
                   )
                 ),
                 fluidRow(
                   column(
                     12,
                     verbatimTextOutput("text_output_three"),
                     plotOutput("plot_output_three")
                   )
                 )
               )
      )
    )
  )
)




############################

server <- function(input, output, session) {

  volumes <- c(Home = fs::path_home(),
               "R Installation" = R.home(),
               getVolumes()())

  shinyDirChoose(input,
                 "directory_filt",
                 roots = volumes,
                 session = session,
                 restrictions = system.file(package = "base"),
                 allowDirCreate = FALSE)

  getdata <- reactive({
    req(input$input$directory_filt)
    parseDirPath(volumes, input$directory_filt)
  })

  output$filepaths_filt <- renderPrint({
    if (is.null(input$directory_filt)) {
      "No directory selected"
    } else {
      getdata()
    }
  })



  observeEvent(input$btn_filter_go, {
    req(input$directory_filt, input$ensembl, input$organism)

    data_rep <- getdata()
    ensembl_version <- input$ensembl
    organism <- input$organism

      # Perform filtering or other operations here
      sobj <- load_and_annotate_recipe_1(data_rep, organism, ensembl_version)

      output$plot_output <- renderPlot({
        show_QC_plots(sobj)
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
      paste("filtered_data", Sys.Date(), ".rds", sep = "")
    },
    content = function(file) {
      saveRDS(getdata(), file)
    }
  )




  #Normalization
  volumes <- c(Home = fs::path_home(),
               "R Installation" = R.home(),
               getVolumes()())
  shinyDirChoose(input, "directory",
                 roots = volumes,
                 session = session,
                 restrictions = system.file(package = "base"),
                 allowDirCreate = FALSE)


  output$filepaths_norm <- renderPrint({
    if (is.null(input$directory)) {
      "No directory selected"
    } else {
      getdata()
    }
  })


  getdata <- reactive({
    req(input$directory_norm)
    parseDirPath(volumes, input$directory_norm)
  })

  observeEvent(input$btn_filter_go, {
    req(input$directory_norm, input$method)

    path_to_filtered_seurat_obj <- getdata()
    normalization_method <- input$method

    # Perform Normalization or other operations here
    normalize_and_plot_main(path_to_filtered_seurat_obj, normalization_method = "cpm")

    output$plot_output <- renderPlot({
      plot_umap
    })

  })


}

shinyApp(ui = ui, server = server)





