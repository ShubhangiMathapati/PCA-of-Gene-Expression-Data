library(shiny)
library(DT)  
library(shinyWidgets)  
library(plotly)  


custom_css <- "
  .title-box {
    background-color: #f8f9fa;  
    padding: 15px;
    border-radius: 8px;
    text-align: center;
    font-weight: bold;
    font-size: 22px;
    box-shadow: 2px 2px 5px rgba(0,0,0,0.1);
  }
  
  .sidebar-custom {
    background-color: #f8f9fa;  
    padding: 15px;
    border-radius: 8px;
    box-shadow: 2px 2px 5px rgba(0,0,0,0.1);
  }
  
  .tab-panel-custom .nav-tabs {
    background-color: #f8f9fa;  
    border-radius: 8px;
    padding: 5px;
  }

  .btn-primary {
    background-color: #007bff;  
    border: none;
  }
  
  .btn-primary:hover {
    background-color: #0056b3;  
  }
"

ui <- fluidPage(
  
  # Applying custom CSS
  tags$head(tags$style(HTML(custom_css))),
  
  # Center-aligning Title inside a styled box
  div(class = "title-box", "PCA Analysis of Gene Expression Data"),
  
  sidebarLayout(
    sidebarPanel(
      class = "sidebar-custom",  
      fileInput("file", "Upload GCT File", accept = c(".gct")),
      selectInput("normalization", "Select Normalization Method:",
                  choices = c("Log2 Transformation", "Z-score Normalization", "Quantile Normalization")),
      pickerInput("filter_group", "Filter by Tissue Type:",  
                  choices = c("All", "Pancreatic Islets", "Lung", "Pancreatic LN"),  
                  selected = "All", multiple = FALSE),
      actionButton("run_pca", "Run PCA", class = "btn-primary"),
      br(), br(),
      downloadButton("downloadData", "Download Processed Data", class = "btn-success"),
      downloadButton("download_pca_plot", "Download PCA Plot", class = "btn-success")
    ),
    
    mainPanel(
      div(
        class = "tab-panel-custom",  
        tabsetPanel(
          type = "tabs",
          tabPanel("PCA Plot", 
                   textOutput("variance_explained"),  
                   plotlyOutput("pca_plot", width = "100%", height = "500px")),  
          tabPanel("Scree Plot", 
                   plotOutput("scree_plot"),
                   downloadButton("download_scree_plot", "Download Scree Plot", class = "btn-success")),
          tabPanel("Data Preview", tableOutput("data_preview")),
          tabPanel("Metadata", dataTableOutput("metadata_table")),  
          
          # "About" Section
          tabPanel("About",
                   h3("Normalization Methods"),
                   p(strong("1. Log2 Transformation:"), " Converts expression values to a logarithmic scale, reducing data skewness and making variance more stable."),
                   p(strong("2. Z-score Normalization:"), " Standardizes data by subtracting the mean and dividing by the standard deviation, ensuring each gene has zero mean and unit variance."),
                   p(strong("3. Quantile Normalization:"), " Adjusts distributions of gene expression values to be similar across samples, minimizing systematic differences."),
                   
                   h3("Data preprocessing mainstep"),
                   p("First 60 lines of metadata are skipped for the PCA analysis"),
                   p("All the rows/columns with 0 and constant values are removed, "),
                   
                   h3("Principal Component Analysis (PCA)"),
                   p("PCA is a dimensionality reduction technique that transforms high-dimensional data into a set of principal components (PCs) that capture the most variance."),
                   p("In the PCA plot:"),
                   tags$ul(
                     tags$li("Samples closer together have similar gene expression profiles."),
                     tags$li("PC1 captures the most variance, followed by PC2."),
                     tags$li("Different tissue types may cluster separately, indicating distinct expression patterns.")
                   ),
                   p("Scree plots help determine how many principal components explain most of the variance in the data.")
          )
        )
      )
    )
  )
)
