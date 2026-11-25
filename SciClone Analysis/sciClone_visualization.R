
install.packages("ggplot2")
install.packages("magick")
install.packages("readr")
install.packages("dplyr")
library(ggplot2)
library(magick)
library(readr)
library(dplyr)


# Creating animations of SciClone clustering data with imagemagick
# https://github.com/genome/sciclone/blob/master/README.md
# Read in files
input_raw <- read_tsv("./SciClone Analysis/sciclone_results/tables/sciclone_total.tsv")
input_clean <- read_tsv("./SciClone Analysis/sciclone_results/tables/sciclone_clustered.tsv")
print(input_clean)

# Function definition for plotting patient VAFs by frame
# Parameters:
# patients -> a list of only unique patient IDs from raw input files
create_plot <- function(patient_id){
  patient_data <- input_raw %>%
    filter(Patient_ID == patient_id)
  
  vaf_ratio <- sum(patient_data$Metastasis_VAF, na.rm = TRUE) / 
               sum(patient_data$Primary_VAF, na.rm = TRUE)
  
  p <- ggplot(
    data = patient_data,
    mapping = aes(
      x = Primary_VAF, 
      y = Metastasis_VAF, 
      color = factor(cluster)
    ) 
  ) + 
      geom_point(size = 6) + 
      geom_abline(
        slope = vaf_ratio, intercept = 0, linetype = "dashed"
      ) +
      labs(
        title= paste("Sciclone Clustering: ", patient_id), 
        subtitle = paste0("VAF Ratio = ", round(vaf_ratio,2)),
        x = "Primary VAF (%)", 
        y = "Metastasis VAF (%)", 
        color = "Cluster"
      ) + 
      coord_cartesian(
        xlim = c(0,100),
        ylim = c(0,100)
      ) + 
      theme_minimal()+
      theme(
        plot.title = element_text(size = 18, face = "bold"
        )
  )
  
  return(p)
}

# Call the function to animate
patients <- unique(input_raw$Patient_ID)
frames <- lapply(patients, create_plot)


test_plot <- create_plot("EV-070")
test_plot

frame_files <- sapply(
  seq_along(frames),
  function(i) {
    
    file <- tempfile(
      fileext = ".png"
    )
    
    ggsave(
      filename = file,
      plot = frames[[i]],
      width = 7,
      height = 6,
      dpi = 150
    )
    
    file
  }
)
gif <- image_read(frame_files)
gif <- image_animate(gif, fps = 1)
 
image_write(gif,"./SciClone Analysis/sciclone_results/plots/sciclone_animation.gif")


annotated_mutations <- read_tsv("./SciClone Analysis/annotated_mutations.tsv")
hist(as.factor(annotated_mutations$Hugo_Symbol))
ggplot(annotated_mutations, aes(x = Hugo_Symbol)) + geom_bar(fill = "green") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggplot(head(count(annotated_mutations, Hugo_Symbol) %>% arrange(desc(n)), 10), aes(reorder(Hugo_Symbol, n), n)) + geom_col(fill="blue") + labs(title = "Top 10 Hugo Symbols") + theme(axis.text.x = element_text(size=12))
