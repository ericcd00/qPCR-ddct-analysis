################################################################################
## Title: qPCR Analysis script
## Creator: Eric Canton Dominguez
################################################################################

## This scrip uses 3 main functions: read_pdf, ddct_analysis and ddct_plot. 
## There is also the complete_ddct_analysis function, that combines the 3 functions above

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
Genes_of_interest = c("NNMT1", "SERPINE1", "SNAI2", "THBS1")

# (Only needed if you will use the ddct_plot function)
# Title of the plot generated
# title = "Title"
title = "Barplot comparativo ddct"


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
  library(dplyr)
}

if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)
}

################################################################################
################################## Function ####################################
################################################################################

read_pdf <- function(pdf_path,
                     result_path,
                     exp_name,
                     wells,
                     analized_groups,
                     housekeeping_genes) {
  
  ################################################################################
  ################################ Sanity Checks #################################
  ################################################################################
  
  # Sanity check 1: Paths  
  if (!file.exists(pdf_path)) {
    stop("This PDF file doesn't exist")
  }
  
  if (!file.exists(result_path)) {
    stop("This path doesn't exist")
  }
  
  # Sanity check 2: Mandatory parameters
  if (missing(pdf_path) || missing(result_path) || missing(exp_name) ||
      missing(wells) || missing(analized_groups) || missing(housekeeping_genes)) {
    stop("There are missing parameters. Check if all the parameters are filled.")
  }
  
  # Sanity check 3: Parameter values
  if (!wells == round(wells) || wells <= 0) {
    stop("The parameter *wells* must be a positive integer number.")
  }
  
  if (!is.character(analized_groups) || length(analized_groups) == 0) {
    stop("The parameter *analized_groups* must be a string vector.")
  }
  
  if (!is.character(housekeeping_genes) || length(housekeeping_genes) == 0) {
    stop("The parameter *housekeeping_genes* must be a string vector.")
  }
  
  
  ################################################################################
  #################################### Code ######################################
  ################################################################################
  
  
  # Extract the tables into a list of dataframes
  Dataframes <- extract_tables(pdf_path, method = "stream", 
                               output= "data.frame")
  
  # Extract the text from the pdf 
  pdf_text <- pdf_text(pdf_path)
  pdf_text_combined <- paste(pdf_text, collapse = "\n")
  
  # Use the pattern that always repeats before and after the gene name to 
  # extract that name
  
  pattern <- "Abs Quant/2nd Derivative Max for (.*?) \\(Abs Quant/2nd Derivative Max\\)"
  
  pattern_string <- str_extract_all(pdf_text_combined, pattern)
  
  # Keep the gene name without the pattern
  gene_names <- strsplit(pattern_string[[1]], 
                         split = " ", 
                         fixed= TRUE)
  
  gene_names <- unlist(lapply(gene_names, function(df){
    df[6]
  }))
  
  
  ## Keep the tables that contain CT values
  #Delete the unnecessary tables
  delete_dataframes <- function(list) {
    
    # Function that filters the list of dataframes
    filtered_list <- lapply(list, function(df) {
      if ("CP" %in% colnames(df)) {
        return(df)
      }
    })
    
    filtered_list <- filtered_list[!sapply(filtered_list, is.null)]
    
    return(filtered_list)
  }
  
  CT_List <- delete_dataframes(Dataframes)
  
  
  ## Join the dataframes that are split in two
  # Due to the placement of the tables in the file, we need to use the number of
  # wells to join the tables that were split
  
  combine_and_delete <- function(list, ind) {
    list[[ind]] <- rbind(list[[ind]], list[[ind + 1]])
    list[[ind + 1]] <- NULL
    list <- list[!sapply(list, is.null)]
    return(list)
  }
  
  
  # Verify that the lenght of the list is the same as the vector of gene names
  i <- 1
  while (i <= length(CT_List) - 2) {
    current_df <- CT_List[[i]]
    next_df <- CT_List[[i + 1]]
    
    # Verify if the dataframe has the same number of rows than wells
    if (nrow(current_df) != wells) {
      CT_List <- combine_and_delete(CT_List, i)
    } else {
      i <- i + 1
    }
  }
  
  # Delete possible 0 row dataframes
  CT_List <- Filter(function(df) {
    nrow(df) > 0
  }, CT_List)
  
  if (length(CT_List) == length(gene_names)) {
    
    # Asign gene names to the names of the dataframes
    names(CT_List) <- gene_names
  } else {
    # error handling (stop)
    stop("The number of dataframes is not the same as the amount of genes.")
  }
  
  ## List that contains the dataframes with the necessary values obtained from the pdf
  Result_List <- lapply(CT_List, function(df) {
    
    # Column with the position of the well of each sample
    pos_well <- c(1:wells)
    df <- cbind(df, Well = pos_well)
    
    # Column with the groups submitted by the user ( mirar which() )
    df$Group <- apply(df, 1, function(row) {
      matching_group <- analized_groups[str_detect(row["Name"], analized_groups)]
      if (length(matching_group) > 0) {
        return(matching_group)
      } else {
        return("No match")
      }
    })
    
    
    # Column with the replicas of each group
    df$Rep <- ave(seq_along(df$Name), df$Name, df$Group, FUN = seq_along)
    
    return(df)
  })
  
  # Check if any of the groups introduced by the user is not in the pdf file
  missing <- analized_groups[!analized_groups %in% 
                               unlist(lapply(Result_List, function(df) unique(df$Group)))]
  if (length(missing) > 0) {
    warning(paste("Some groups were not found in the pdf file:", paste(missing, collapse = ", ")))
  }
  
  # Add a new column "HK". It contais information about the gene being a housekeeping.
  i = 1
  for (i in 1:(length(Result_List))) {
    
    is_hk <- housekeeping_genes[str_detect(gene_names[i], housekeeping_genes)]
    if (length(is_hk) > 0) {
      Result_List[[i]]$HK <- "Yes"
    }  else {
      Result_List[[i]]$HK <- "No"
    }
    
  }
  
  # Reorganize dataframes and delete H2O samples
  Result_List <- lapply(Result_List, function(df) {
    order <- c("Pos", "Group", "Name", "Rep", "CP", "Status", "HK")
    df <- df[, order]
    df <- subset(df, Name != "H2O")
    
    return(df)
  })
  
  hs <- createStyle(textDecoration = "BOLD", fontColour = "white", fontSize = 12, 
                    fontName = "Arial Narrow", fgFill = "slateblue")
  
  write.xlsx(x = Result_List, 
             file = paste0(result_path, "/", exp_name, "_qPCR_data.xlsx"),
             headerStyle = hs,
             rowNames = FALSE)
  
  save(Result_List, 
       result_path,
       exp_name,
       wells,
       analized_groups,
       housekeeping_genes, 
       file = paste0(result_path, "/", exp_name, ".RData"))
  
  return(Result_List)
  
}


ddct_analysis <- function(data,
                          result_path,
                          exp_name,
                          housekeeping_genes,
                          control_variable){
  
  ################################################################################
  ################################ Sanity Checks #################################
  ################################################################################
  
  # Sanity check 1: Paths
  if (!file.exists(result_path)) {
    stop("This path doesn't exist")
  }
  
  # Sanity check 2: Mandatory parameters
  if (missing(data) || missing(result_path) || missing(exp_name) ||
      missing(housekeeping_genes) || missing(control_variable)) {
    stop("There are missing parameters. Check if all the parameters are filled.")
  }
  
  # Sanity check 3: Parameter values
  if (!wells == round(wells) || wells <= 0) {
    stop("The parameter *wells* must be a positive integer number.")
  }
  
  if (!is.character(housekeeping_genes) || length(housekeeping_genes) == 0) {
    stop("The parameter *housekeeping_genes* must be a string vector.")
  }
  
  if (!is.character(control_variable) || length(control_variable) == 0) {
    stop("The parameter *control_variable* must be a string vector.")
  }
  
  ################################################################################
  #################################### Code ######################################
  ################################################################################
  
  
  sample_means <- lapply(data, function(df) {
    means <- aggregate(CP ~ Name, df, mean)
    return(means)
  })
  
  results <- lapply(data, function(df) {
    
    merged_df <- merge(df, sample_means[[housekeeping_genes]], by = "Name", all.x = TRUE)
    names(merged_df)[names(merged_df) == 'CP.x'] <- 'CP'
    names(merged_df)[names(merged_df) == 'CP.y'] <- 'CP.HK'
    
    merged_df$delta1 <- merged_df$CP - merged_df$CP.HK
    
    meanControl <- aggregate(delta1 ~ Group, merged_df, mean)
    meanControl <- subset(meanControl, Group == control_variable)
    mean <- meanControl$delta1
    
    merged_df$delta2 <- merged_df$delta1 - mean
    
    merged_df$ddct <- 2^(-merged_df$delta2)
    
    return(merged_df)
  })
  
  final <- lapply(results, function(df){
    merged_df <- aggregate(ddct ~ Group, df, mean)
    sd <- aggregate(ddct ~ Group, df, sd)
    merged_df$sd <- sd$ddct
    merged_df$ddct <- round(merged_df$ddct, 2)
    merged_df$sd <- round(merged_df$sd, 2)
    return(merged_df)
    
  })
  
  hs <- createStyle(textDecoration = "BOLD", fontColour = "white", fontSize = 14, 
                    fontName = "Arial Narrow", fgFill = "turquoise4")
  
  write.xlsx(x = final, 
             file = paste0(result_path, "/", exp_name, "_qPCR_ddct.xlsx"),
             headerStyle = hs,
             rowNames = FALSE)
  
  return(final)
}


ddct_plot <- function(data,
                      ddct_values, 
                      Genes_of_interest, 
                      title,
                      result_path){
  
  datos <-  lapply(data, function(df) {
    df[, c("Group", "CP"), drop = FALSE]
  })
  datos <- bind_rows(datos, .id = "ID")
  datos <- datos %>% filter(ID %in% Genes_of_interest)
  
  df_combined <- bind_rows(ddct_values, .id = "ID")
  
  df_filtrado <- df_combined %>% filter(ID %in% Genes_of_interest)
  
  plot <- ggplot(df_filtrado, aes(x = ID, y=ddct, fill=Group)) + 
    geom_bar(stat="identity", position = position_dodge(.6), 
             width = .4) +
    
    geom_errorbar(aes(ymin = ddct - sd, ymax = ddct + sd), 
                  position = position_dodge(.6), width = .2) +
    
    ## geom_point(data = datos, aes(ID, CP), position = position_dodge(.6)) +
    
    ggtitle(title) +
    xlab("") +
    ylab("Relative mRNA levels") +
    scale_fill_grey(start = .1, end = .5) +
    
    scale_y_continuous(expand = c(0, 0)) +
    
    geom_hline(yintercept=0) +
    theme_bw() + 
    
    theme(plot.title = element_text(hjust = 0.5),
          panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.title = element_blank(),
          axis.title.y = element_text(face="bold"),
          axis.line = element_line(colour = "black", linewidth = 1),
          axis.text.x = element_text(face="bold", size=12, angle=45, 
                                     vjust = 0.6, colour="Black"),
          axis.text.y = element_text(face="bold", size="12"))
  
  print(plot)
  
  ggsave(file = paste0(result_path, "/", exp_name, "_ddct_plot.pdf"),
         plot= last_plot(),
         device = pdf,
         path=result_path)          
  
}


complete_ddct_analysis <- function() {
  
  Results <- read_pdf(pdf_path = pdf_path, 
                      result_path = result_path, 
                      exp_name = exp_name, 
                      wells = wells, 
                      analized_groups = analized_groups, 
                      housekeeping_genes = housekeeping_genes)
  
  ddct_results <- ddct_analysis(data = Results,
                                result_path = result_path, 
                                exp_name = exp_name, 
                                housekeeping_genes = housekeeping_genes,
                                control_variable = control_variable)
  
  plotting <- ddct_plot(data = Results,
                        ddct_values = ddct_results,
                        Genes_of_interest = Genes_of_interest,
                        title = title,
                        result_path = result_path)
  
  Complete_list <- list(Results,
                        ddct_results)
  
  names(Complete_list) <- c("qPCR Results",
                            "DDCT results")
  return(Complete_list)
  
}


################################################################################
############################## Execute function ################################
################################################################################

Results <- read_pdf(pdf_path = pdf_path, 
                    result_path = result_path, 
                    exp_name = exp_name, 
                    wells = wells, 
                    analized_groups = analized_groups, 
                    housekeeping_genes = housekeeping_genes)


ddct_results <- ddct_analysis(data = Results,
                              result_path = result_path, 
                              exp_name = exp_name, 
                              housekeeping_genes = housekeeping_genes,
                              control_variable = control_variable)


plotting <- ddct_plot(data = Results,
                      ddct_values = ddct_results,
                      Genes_of_interest = Genes_of_interest,
                      title = title,
                      result_path = result_path)


complete_ddct_analysis <- complete_ddct_analysis()
