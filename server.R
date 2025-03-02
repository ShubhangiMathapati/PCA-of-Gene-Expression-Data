library(shiny)
library(ggplot2)
library(plotly)
library(DT)  
library(preprocessCore)  
library(shinyWidgets)

options(shiny.maxRequestSize = 50 * 1024^2)  

server <- function(input, output, session) {
  
  dataset <- reactive({
    req(input$file)
    
    # If a file other than .gct file is uploaded it will throw error saying invalid file
    validate(need(grepl("\\.gct$", input$file$name), "Please upload a valid GCT file."))
    
    metadata_raw <- readLines(input$file$datapath, n = 60)
    
    metadata_list <- lapply(metadata_raw, function(line) {
      parts <- unlist(strsplit(line, "\t"))
      if (length(parts) > 1) {
        return(data.frame(Key = parts[1], Value = paste(parts[-1], collapse = " ")))
      } else {
        return(NULL)
      }
    })
    
    metadata_df <- do.call(rbind, metadata_list)
    
    df <- read.table(input$file$datapath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, skip = 60)
    sample_names <- colnames(df)[-1]
    
    rownames(df) <- df[, 1]
    df <- df[, -1]
    df <- as.data.frame(lapply(df, as.numeric))
    
    return(list(data = df, samples = sample_names, metadata = metadata_df))
  })
  
  normalized_data <- reactive({
    req(dataset())
    df <- dataset()$data
    
    df <- df[, apply(df, 2, var, na.rm = TRUE) > 0]   # Remove genes with zero variance
    df[df <= 0] <- NA                               # Replace values less than or equal to 0 as NA
    df <- df[, colSums(!is.na(df)) > 0]  
    
    original_colnames <- colnames(df)
    original_rownames <- rownames(df)
    
    if (input$normalization == "Log2 Transformation") {
      df <- log2(df + 1)    # Converting expression values to log scale
    } else if (input$normalization == "Z-score Normalization") {
      df <- scale(df)   # Centers each gene
    } else if (input$normalization == "Quantile Normalization") {
      df <- normalize.quantiles(as.matrix(df)) # Makes gene expression distributions uniform across samples
      df <- as.data.frame(df)
      colnames(df) <- original_colnames
      rownames(df) <- original_rownames
    }
    
    # If log2(0) or scale() generates NA or Inf, we replace them with the mean of that gene.
    df <- as.data.frame(apply(df, 2, function(x) {
      x[is.infinite(x) | is.na(x)] <- mean(x, na.rm = TRUE)
      return(x)
    }))
    
    df <- df[, apply(df, 2, var, na.rm = TRUE) > 0]
    
    return(df)
  })
  
  pca_result <- reactive({
    req(normalized_data())
    df <- normalized_data()
    
    df_t <- t(df)
    colnames(df_t) <- rownames(df)
    rownames(df_t) <- dataset()$samples
    
    df_t <- scale(df_t)
    
    df_t <- df_t[, apply(df_t, 2, var, na.rm = TRUE) > 0]
    
    prcomp(df_t, center = TRUE, scale. = FALSE)
  })
  
  output$variance_explained <- renderText({
    req(pca_result())
    pca <- pca_result() 
    
    # Calculating the variance here 
    variance_explained <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 2)
    paste("Variance Explained: PC1 =", variance_explained[1], "% | PC2 =", variance_explained[2], "%")
  })
  
  output$pca_plot <- renderPlotly({
    req(pca_result())
    pca <- pca_result()
    
    pca_df <- data.frame(     # Creating PCA dataframe here
      Sample = dataset()$samples,
      PC1 = pca$x[,1],
      PC2 = pca$x[,2],
      Group = rep(c("Pancreatic Islets", "Lung", "Pancreatic LN"), length.out = nrow(pca$x))
    )
    
    # Defining consistent colors here so that when every time one tissue type is selected
    # it will show different colors assigned to respective tissue type
    color_map <- c("Pancreatic Islets" = "green", "Lung" = "red", "Pancreatic LN" = "blue")
    
    # Placing a tissue type filter so that individual tissue types can be visualized if needed
    if (input$filter_group != "All") {
      pca_df <- pca_df[pca_df$Group == input$filter_group, ]
    }
    
    # Showing PC1 and PC2 scores while hovering over the dots
    pca_df$HoverText <- paste(
      "PC1:", round(pca_df$PC1, 3),
      "<br>PC2:", round(pca_df$PC2, 3)
    )
    
    # Creating the plot
    p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, text = HoverText)) +
      geom_point(size = 3) +
      scale_color_manual(values = color_map) +
      labs(title = "PCA of Gene Expression Data (Normalized)",
           x = "Principal Component 1",
           y = "Principal Component 2",
           color = "Tissue Type") +
      theme_minimal()
    
    # Converting the plot to plotly with hover text to be able to add the hovering functionality
    ggplotly(p, tooltip = "text")
  })
  
  output$scree_plot <- renderPlot({      # function to create a Scree plot
    req(pca_result())
    pca <- pca_result()
    variance_explained <- pca$sdev^2 / sum(pca$sdev^2)
    
    ggplot(data = data.frame(PC = 1:length(variance_explained), Variance = variance_explained),
           aes(x = PC, y = Variance)) +
      geom_bar(stat = "identity", fill = "steelblue") +
      geom_line(color = "red") +
      geom_point(size = 3) +
      labs(title = "Scree Plot", x = "Principal Components", y = "Variance Explained") +
      theme_minimal()
  })
  
  output$data_preview <- renderTable({     # Displays the processed data for PCA plotting
    head(normalized_data(), 10)
  })
  
  output$metadata_table <- renderDataTable({    # Extracts metadata to display in metadata tab
    req(dataset())
    datatable(dataset()$metadata, options = list(pageLength = 10, autoWidth = TRUE))
  })
  
  output$download_pca_plot <- downloadHandler(    # allows to download the PCA plot as a .png file
    filename = function() { "pca_plot.png" },
    content = function(file) {
      req(pca_result())
      
      pca <- pca_result()
      pca_df <- data.frame(Sample = dataset()$samples, 
                           PC1 = pca$x[,1], 
                           PC2 = pca$x[,2],
                           Group = rep(c("Pancreatic Islets", "Lung", "Pancreatic LN"), 
                                       length.out = nrow(pca$x)))
      
      color_map <- c("Pancreatic Islets" = "green", "Lung" = "red", "Pancreatic LN" = "blue")
      
      if (input$filter_group != "All") {
        pca_df <- pca_df[pca_df$Group == input$filter_group, ]
      }
      
      p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
        geom_point(size = 3) +
        scale_color_manual(values = color_map) +
        labs(title = "PCA of Gene Expression Data",
             x = "PC 1",
             y = "PC 2",
             color = "Tissue Type") +
        theme_minimal()
      
      ggsave(file, plot = p, width = 7, height = 5, dpi = 300)
    }
  )
  
  output$downloadData <- downloadHandler(       # Allows to download processed csv file.
    filename = function() { "processed_data.csv" },
    content = function(file) {
      write.csv(normalized_data(), file, row.names = TRUE)
    }
  )
  
  output$download_scree_plot <- downloadHandler(    # allows to download the scree plot as a .png file
    filename = function() { "scree_plot.png" },
    content = function(file) {
      req(pca_result())
      pca <- pca_result()
      variance_explained <- pca$sdev^2 / sum(pca$sdev^2)
      
      # Creating the Scree plot
      p <- ggplot(data = data.frame(PC = 1:length(variance_explained), Variance = variance_explained),
                  aes(x = PC, y = Variance)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        geom_line(color = "red") +
        geom_point(size = 3) +
        labs(title = "Scree Plot", x = "Principal Components", y = "Variance Explained") +
        theme_minimal()
      
      # Saving the plot as PNG
      ggsave(file, plot = p, width = 7, height = 5, dpi = 300)
    }
  )
}
