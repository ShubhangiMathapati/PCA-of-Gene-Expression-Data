# PCA-of-Gene-Expression-Data
This repository contains code for creating an RShiny dashboard which shows the PCA of a Gene Expression Dataset.
Packages to be installed are - shiny, data.table, ggplot2, plotly, preprocessCore, DT, shinyWidgets.
Place the raw gct file, app.R, ui.R, server.R files in one same directory. Create a folder named "www" in this directory and inside the folder place style.css file which helps in styling the dashboard.
Run the app.R script which will automatically run the other two scripts ui.R and server.R
Once the dashboard is available,
1. Upload your raw .gct file
2. Select the normalization method and click "Run PCA button".
3. Hover over the points to know the exact point value.
4. You can view individual tissue type by selecting the "Filter for tissue type" button.
5. Click on "Scree Plot" tab to view the Scree plot. It changes according to the normalization method selected.
6. Click on the "Metadata" tab for the Metadata information.
7. Click on the "Data Preview" tab to view the processed data for PCA plot.
8. Click on "Download Processed Data" button to download the csv file of the processed data for PCA plot.
9. Click on "Download PCA Plot" button to download the plot in .png format.

APPROACH- 
As the data had many zero values, I removed the rows where zero values were present. Then, zero variance still present are replaced with NA for log2 transformation and then all the 0 values or NA values still present are removed once again to ensure smooth PCA. 
I have also used imputation here as there were mane 0 values but also there were rows with non zero and zero values both hence, they had to be imputated.
These were the key decisions made on the data received. 
Metadata mapping was a bit tedious here hence, did not add the metadata information to view on hovering over the points.
For better viewing, I have added tissue type filter to view individual tissue type filter.
