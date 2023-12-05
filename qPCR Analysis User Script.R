################################################################################
## Title: qPCR Analysis User Script
## Creator: Eric Canton Dominguez
################################################################################

## This scrip uses 3 main functions: read_pdf, ddct_analysis and ddct_plot. 
## There is also the complete_ddct_analysis function, that combines the 3 functions above.

# 1. Read_pdf: Use the read_pdf function to extract the qPCR data from a pdf 
# file and export it to a .xslx file.
# 2. ddct_analysis: Use the ddct_analysis function to quicklyperformthe ddct
# method from an excel file or the result from the pdf read by the read_pdf function
# 3. ddct_plot: Use the ddct_plot to generate a plot comparing the values of 
# the different conditions

# To use this function it is required to have Java 6 and JDK previously 
# installed. 

################################################################################
################################ INSTRUCTIONS ##################################
################################################################################

### This is the only section the user must fill ###

# Name of the experiment (This will be used to name the files).
# exp_name <- "Name_of_the_experiment"
exp_name <- "Exp010102"

# Path where the results will be stored
# result_path <- "path_to_store_results"
result_path = "~/Downloads/qPCR_data/Results"

# Housekeeping genes
# housekeeping_genes = c("ACTB", "...",...)
housekeeping_genes = c("ACTB")

# Number of wells used for each gene
# wells = number of wells (15)
wells = 15

# Names of all the conditions ("Groups") of the analysis
# analized_groups = c("UNT", "CAF", "...",...)
analized_groups = c("UNT", "CAF")

# (Only needed if you *won't*  use the read_pdf function)
# qPCR excel file
# Result = "path_of_the_excel_file"
Excel_path = "~/Downloads/qPCR_data/Exp010102_qPCR_data.xlsx"

# (Only needed if you will use the read_pdf function)
# Path of the LightCycler 480 PDF Report 
# pdf_path <- "path_of_the_file"
pdf_path <- "~/Downloads/qPCR_data/qPCR_RawData_LauraS.pdf"


# (Only needed if you will use the ddct_analysis function)
# Variable used as the control variable for the ddct analysis 
# control_variable = "UNT"
control_variable = "UNT"


# (Only needed if you will use the ddct_plot function)
# Genes wanted to be plotted 
# Genes_of_interest = c("NNMT1", "SERPINE1", "SNAI2", "THBS1")
Genes_of_interest = c("SNAI2", "SERPINE1", "NNMT1", "THBS1")

# (Only needed if you will use the ddct_plot function)
# Title of the plot generated
# title = "Title"
title = "ddct results"


################################################################################
########################## 1. Install all Packages #############################
################################################################################

# All the packages must be installed prior to the execution of the function.

if(!require(tabulizer)){
  devtools::install_github("ropensci/tabulizer")
  library(tabulizer) # Extract the tables from the pdf file 
}

if(!require(pdftools)){
  install.packages("pdftools")
  library(pdftools) # Extract the text from the pdf file 
}

if(!require(stringr)){
  install.packages("stringr")
  library(stringr) # Recognize patterns
}

if(!require(openxlsx)){
  install.packages("openxlsx")
  library(openxlsx) # Save the results in an excel file
}

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr) # Utility
}


if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2) # Plots
}

if(!require(ggpubr)){
  install.packages("ggpubr")
  library(ggpubr) # Plot statistics
}

if(!require(readxl)){
  install.packages("readxl")
  library(readxl) # Read excel files
}

source("qPCR Analysis Functions.R") # All of the functions we will need 

################################################################################
############################# 2. Execute functions #############################
################################################################################

# Read the pdf file of the qPCR
Results <- read_pdf(pdf_path = pdf_path, 
                    result_path = result_path, 
                    exp_name = exp_name, 
                    wells = wells, 
                    analized_groups = analized_groups, 
                    housekeeping_genes = housekeeping_genes)

# Read the excel fileof the qPCR
Results <- read_excel(Excel_path = Excel_path)

# Perform a ddct analysis of the file
ddct_results <- ddct_analysis(data = Results,
                              result_path = result_path, 
                              exp_name = exp_name, 
                              housekeeping_genes = housekeeping_genes,
                              control_variable = control_variable)

# Plot the results of the ddct analysis
plot_result <- ddct_plot(data = Results,
                      ddct_values = ddct_results,
                      Genes_of_interest = Genes_of_interest,
                      title = title,
                      result_path = result_path)

# Perform all functions at one from the pdf file 
complete_ddct_analysis <- complete_ddct_analysis()

