################################################################################
## Title: qPCR Analysis User Script
## Creator: Eric Canton Dominguez
################################################################################

## This script contains these functions: 

# 1. Read_pdf: Use the read_pdf function to extract the qPCR data from a pdf 
# file and export it to a .xlsx file.

# To use this function it is required to have Java 6 and JDK previously ########
# installed. ###################################################################

# 2. Read_excel: Use the read_excel function to import the data directly from  
# a .xlsx file.

# 3. ddct_analysis: Use the ddct_analysis function to quickly perform the ddct
# method from an excel file or the result from the pdf read by the read_pdf function

# 4. ddct_plot: Use the ddct_plot to generate a plot comparing the values of 
# the different conditions

## 5. complete_ddct_analysis: Combination of the read_pdf, ddct_analysis and ddct_plot functions.

# 6. qPCR_experiments: Combine 2 or more experiments and make a plot with their ddct values

################################################################################
################################ INSTRUCTIONS ##################################
################################################################################

### This is the only section the user must fill ###

# Name of the experiment (This will be used to name the files).
# exp_name <- "Name_of_the_experiment"
exp_name <- "2023-10-10"

# Path where the results will be stored
# result_path <- "path_to_store_results"
result_path = "~/Downloads/qPCR_data/Results"

# Housekeeping genes
# housekeeping_genes = c("ACTB", "...",...)
housekeeping_genes = c("HPRT1")

# Number of wells used for each gene
# wells = number of wells (15)
wells = 18

# Names of all the conditions ("Groups") of the analysis
# If one group contains another ("SPP1" and "SPP12") put the longest one first.
# If there are special characters in a group ("SPP1+EPZ) put "\\" before the 
# special character ("SPP1\\+EPZ")
# analized_groups = c("UNT", "CAF", "...",...)
analized_groups = c("UNT", "CAF") #Hacer pruebas con "\\"

# (Only needed if you *won't*  use the read_pdf function)
# qPCR excel file
# Result = "path_of_the_excel_file"
Excel_path = "~/Downloads/qPCR_data/Results/2023-01-20_qPCR_data.xlsx"

# (Only needed if you will use the read_pdf function)
# Path of the LightCycler 480 PDF Report 
# pdf_path <- "path_of_the_file"
pdf_path <- "~/Downloads/qPCR_data/2023-10-10 Laura S_EXAMPLE2.PDF"


# (Only needed if you will use the ddct_analysis function)
# Variable used as the control variable for the ddct analysis 
# control_variable = "UNT"
control_variable = "UNT"


# (Only needed if you will use the ddct_plot function)
# Genes wanted to be plotted 
# Genes_of_interest = c("NNMT1", "SERPINE1", "SNAI2", "THBS1")
Genes_of_interest = c("MRAS", "TUBB6")

# Names of all the conditions ("Groups") to plot
# Groups_of_interest = c("UNT", "CAF", "...",...)
Groups_of_interest = c("UNT", "CAF")


# (Only needed if you will use the ddct_plot function)
# Title of the plot generated
# title = "Title"
title = "ddct results"


# (Only needed if you will use the qPCR_experiments function)
# List of all the experiments to mix
# exp_paths = list(~/Exp/Path/Exp1.xlsx", 
#                 "~/Exp/Path/Exp2.xlsx")
exp_paths = list("~/Downloads/qPCR_data/Exp10102023_qPCR_ddct.xlsx",
                 "~/Downloads/qPCR_data/Exp10102023_qPCR_ddct.xlsx", 
                 "~/Downloads/qPCR_data/Exp10102023_qPCR_ddct.xlsx")

# (Only needed if you will use the qPCR_experiments function)
# List of all the experiments to mix
# exp_paths = list("2023-11-10", "2023-12-10", "2023-12-20")
exp_names <- c("Experimento1", "Experimento2", "Experimento3")


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

if(!require(ggsignif)){
  install.packages("ggsignif")
  library(ggsignif) # Plot statistics
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


# From Excel to a HT-qPCR object


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
                      result_path = result_path,
                      Groups_of_interest = Groups_of_interest,
                      exp_name = exp_name)

# Perform all functions at one from the pdf file 
complete_ddct_analysis <- complete_ddct_analysis()

# Join multiple experiments and plot them
combine_qPCRs <- qPCR_experiments(exp_paths = exp_paths,
                              exp_names = exp_names,
                              Genes_of_interest = Genes_of_interest,
                              title = title,
                              result_path = result_path,
                              Groups_of_interest = Groups_of_interest,
                              exp_name = exp_name)

