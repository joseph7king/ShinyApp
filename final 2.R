
library(curatedMetagenomicData); packageVersion("curatedMetagenomicData")
library(tidyverse); packageVersion("tidyverse")                                                 
library(phyloseq); packageVersion("phyloseq")  
library(LinDA); packageVersion("LinDA")  
library(topicmodels); packageVersion("topicmodels")     
library(tidytext); packageVersion("tidytext")            
library(ldatuning); packageVersion("ldatuning")   
library(cowplot); packageVersion("cowplot")  
library(dplyr); packageVersion("dplyr")  
library(viridis); packageVersion("viridis")  
library(glue); packageVersion("glue")  
 

# List all available datasets in curatedMetagenomicData
available_datasets <- curatedMetagenomicData::sampleMetadata

# Print the number of available datasets
cat("Number of datasets available in curatedMetagenomicData:", length(available_datasets), "\n")

# Print a brief summary of the datasets
summary_df <- data.frame(Dataset = names(available_datasets))
summary(summary_df)
 

meta_df <- curatedMetagenomicData::sampleMetadata

mydata <- meta_df %>%
  dplyr::filter(study_name == "FengQ_2015") %>%
  dplyr::select(sample_id, subject_id, body_site, disease, age, gender, BMI) %>%
  dplyr::mutate(case_status = ifelse(disease == "CRC", "CRC", 
                                     ifelse(disease == "healthy", "Control", "Other"))) %>%
  dplyr::filter(case_status != "Other") %>%
  na.omit(.)

keep_id <- mydata$sample_id

crc_se_sp <- sampleMetadata %>%
  dplyr::filter(sample_id %in% keep_id) %>%
  curatedMetagenomicData::returnSamples("relative_abundance", counts = TRUE)

crc_sp_df <- as.data.frame(as.matrix(assay(crc_se_sp)))

crc_sp_df <- crc_sp_df %>%
  rownames_to_column(var = "Species")

otu_tab <- crc_sp_df %>%
  tidyr::separate(Species, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|") %>%
  dplyr::select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus) %>%
  dplyr::mutate(Species = gsub("s__", "", Species)) %>%
  tibble::column_to_rownames(var = "Species")

tax_tab <- crc_sp_df %>%
  dplyr::select(Species) %>%
  tidyr::separate(Species, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|") %>%
  dplyr::mutate(Kingdom = gsub("k__", "", Kingdom),
                Phylum = gsub("p__", "", Phylum),
                Class = gsub("c__", "", Class),
                Order = gsub("o__", "", Order),
                Family = gsub("f__", "", Family),
                Genus = gsub("g__", "", Genus),
                Species = gsub("s__", "", Species)) %>%
  dplyr::mutate(spec_row = Species) %>%
  tibble::column_to_rownames(var = "spec_row")

rownames(mydata) <- NULL
mydata <- mydata %>%
  tibble::column_to_rownames(var = "sample_id")

(ps <- phyloseq(sample_data(mydata),
                otu_table(otu_tab, taxa_are_rows = TRUE),
                tax_table(as.matrix(tax_tab))))

ps
 

sample_data(ps)$read_count <- sample_sums(ps)
#ps <- subset_samples(ps, read_count > 10^7)

minTotRelAbun <- 1e-5           
x <- taxa_sums(ps)
keepTaxa <- (x / sum(x)) > minTotRelAbun
(ps <- prune_taxa(keepTaxa, ps))
 

head(sample_data(ps))
 

table(sample_data(ps)$case_status)
 

rm(list= ls()[!(ls() %in% c("ps"))])
 

(ps_g <- tax_glom(ps, taxrank = "Genus"))
 

taxa_names(ps_g) <- tax_table(ps_g)[, 6]
head(taxa_names(ps_g))
 

set.seed(42)

count_matrix <- data.frame(t(data.frame(otu_table(ps_g))))

# Remove terms and documents with zero counts
count_matrix <- count_matrix[rowSums(count_matrix) > 0, ]
count_matrix <- count_matrix[, colSums(count_matrix) > 0]

result <- FindTopicsNumber(
  count_matrix,
  topics = seq(from = 3, to = 50, by = 3),
  metrics = c("CaoJuan2009", "Arun2010"),
  method = "VEM",
  control = list(seed = 42),
  mc.cores = 4,
  verbose = TRUE
)
 

my_plot <- function (values)  # updating to allow plotting via cowplot::plot_grid
{
  if ("LDA_model" %in% names(values)) {
    values <- values[!names(values) %in% c("LDA_model")]
  }
  columns <- base::subset(values, select = 2:ncol(values))
  values <- base::data.frame(values["topics"], base::apply(columns, 
                                                           2, function(column) {
                                                             scales::rescale(column, to = c(0, 1), from = range(column))
                                                           }))
  values <- reshape2::melt(values, id.vars = "topics", na.rm = TRUE)
  values$group <- values$variable %in% c("Griffiths2004", "Deveaud2014")
  values$group <- base::factor(values$group, levels = c(FALSE, TRUE), labels = c("Minimize", "Maximize"))
  p <- ggplot(values, aes_string(x = "topics", y = "value", group = "variable"))
  p <- p + geom_line()
  p <- p + geom_point(aes_string(shape = "variable"), size = 3)
  p <- p + guides(size = FALSE, shape = guide_legend(title = "Metrics:"))
  p <- p + scale_x_continuous(breaks = values$topics)
  p <- p + labs(
    title = "Evaluation of Topic Models: Number of Topics vs Metrics",
    x = "\nNumber of Topics",
    y = NULL
  )
  p <- p + facet_grid(group ~ .)
  p <- p + theme_bw() %+replace% theme(panel.grid.major.y = element_blank(), 
                                       panel.grid.minor.y = element_blank(), 
                                       panel.grid.major.x = element_line(colour = "grey70"), 
                                       panel.grid.minor.x = element_blank(), 
                                       legend.key = element_blank(), 
                                       strip.text.y = element_text(angle = 90),
                                       legend.position = "bottom")
  return(p)
}

column_names <- c("topics","CaoJuan2009" ,"Arun2010")
(p1 <- my_plot(result %>% select(all_of(column_names))))
 
# Min-max normalization for CaoJuan2009
result$CaoJuan2009_norm <- (result$CaoJuan2009 - min(result$CaoJuan2009)) / (max(result$CaoJuan2009) - min(result$CaoJuan2009))

# Min-max normalization for Arun2010
result$Arun2010_norm <- (result$Arun2010 - min(result$Arun2010)) / (max(result$Arun2010) - min(result$Arun2010))

# Calculate the mean of the normalized metrics for each row
result$mean_metric <- rowMeans(result[, c("CaoJuan2009_norm", "Arun2010_norm")])

# Identify the row with the lowest mean value
optimal_topic <- result[which.min(result$mean_metric), ]

# Print the optimal number of topics and corresponding metrics
optimal_topic

#   topics CaoJuan2009 Arun2010 CaoJuan2009_norm Arun2010_norm mean_metric
# 18	0.0390238	100.0284	0	0.1419839	0.07099197
 

set.seed(224)

lda_k <- LDA(count_matrix, k = optimal_topic$topics, method = "VEM", control = list(seed = 42))

b_df <- data.frame(tidy(lda_k, matrix = "beta"))

g_df <- data.frame(tidy(lda_k, matrix = "gamma")) %>%
  arrange(document, topic)

head(b_df)
 

lib_size_df <- data.frame(sample_sums(ps_g)) %>%
  dplyr::rename("read_count" = "sample_sums.ps_g.") %>%
  rownames_to_column(var = "document")

tm_df <- left_join(lib_size_df, g_df) %>%
  mutate(topic_count = read_count * gamma,
         topic_count = round(topic_count, 0)) %>%
  dplyr::select(-read_count, -gamma) %>%
  pivot_wider(names_from = topic, values_from = topic_count) %>%
  dplyr::rename_with(~ paste0("Topic_", .), -document) %>%
  column_to_rownames(var = "document") %>%
  t(.) %>%
  data.frame(.)

(ps_topic_g <- phyloseq(
  sample_data(ps_g),
  otu_table(tm_df, taxa_are_rows = TRUE)))

 



ps_otu_df <- data.frame(otu_table(ps_topic_g))
ps_otu_df <- ps_otu_df[-nrow(ps_otu_df),]
ps_otu_df[is.na(ps_otu_df)] <- 0
linda <- linda(otu.tab = ps_otu_df, meta = data.frame(sample_data(ps_topic_g)), 
               formula = '~ case_status + age + gender + BMI', alpha = 0.05, winsor.quan = 0.97,
               prev.cut = 0.3, lib.cut = 0, n.cores = 1)
linda
 

# Process and arrange data by absolute value of log2FoldChange
linda_res_df <- data.frame(linda$output$case_statusCRC) %>%
  dplyr::arrange(padj) %>%
  rownames_to_column(var = "Topic")

fdr_linda <- linda_res_df %>%
  mutate(reject = ifelse(padj < 0.05, "Yes", "No")) %>%
  separate(Topic, into = c("t", "t_n"), remove = FALSE, convert = TRUE) %>%
  mutate(Topic = gsub("_", " ", Topic)) %>%
  arrange(desc(abs(log2FoldChange))) %>%  # Order by absolute value
  mutate(Topic = factor(Topic, levels = rev(Topic)))  # Reverse levels for top-down order

# Create the plot
p2 <- ggplot(data = fdr_linda, aes(x = Topic, y = log2FoldChange, fill = reject)) +
  geom_col(alpha = .8) +
  labs(
    title = "Differential Abundance of Bacterial Genera in Colorectal Cancer",
    y = "\nLog2 Fold-Change",
    x = ""
  ) +
  theme(
    axis.text.x = element_text(color = "black", size = 8),
    axis.text.y = element_text(color = "black", size = 8),
    axis.title.y = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.position = "none"
  ) +
  coord_flip() +
  scale_fill_manual(values = c("grey", "dodgerblue")) +
  geom_hline(yintercept = 0, linetype = "dotted")

p2
 
  
top1_topic = head(fdr_linda,1)$t_n
top1_topic
 

  
# Improved Data Visualization for b_plot
b_plot <- b_df %>%
  dplyr::filter(topic == top1_topic) %>%
  arrange(desc(beta)) %>%
  slice_head(n = 20) %>%
  arrange(desc(beta)) %>%
  mutate(term = factor(term, levels = rev(term)))  # Reverse the order for better display

# Enhanced Plot p3
p3 <- ggplot(data = b_plot, aes(x = beta, y = term, color = beta)) +
  geom_point(aes(size = beta), alpha = 0.7) +
  scale_color_viridis_c(direction = -1) +  # Reverse color scale direction
  labs(
    title = glue("Top 20 Terms in Topic {top1_topic}"),
    subtitle = "Terms are ordered by their contribution to the topic",
    x = "\nBeta (Topic-Term Probability)",
    y = NULL,
    size = "Beta"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",  # Remove legend
    axis.text.y = element_text(size = 10, color = "black"),  # Enhance y-axis labels
    axis.text.x = element_text(size = 10, color = "black"),  # Enhance x-axis labels
    axis.title.x = element_text(size = 12, face = "bold"),  # Enhance x-axis title
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Center and bold title
    plot.subtitle = element_text(size = 12, hjust = 0.5)  # Center subtitle
  )

# Display the plot
p3
 

  
linda_g_all <- linda(otu.tab = data.frame(otu_table(ps_g)), meta = data.frame(sample_data(ps_g)), 
                     formula = '~ case_status + age + gender + BMI', alpha = 0.05, winsor.quan = 0.97,
                     prev.cut = 0, lib.cut = 0, n.cores = 1)
 
  
data.frame(linda_g_all$output$case_statusCRC) %>%
  dplyr::arrange(padj) %>%
  rownames_to_column(var = "Genus")
 

  
linda_res_g_df <- data.frame(linda_g_all$output$case_statusCRC) %>%
  dplyr::arrange(padj) %>%
  rownames_to_column(var = "Genus")

keep_g <- b_plot$term

# Process and arrange data for p4 plot
log2_df <- linda_res_g_df %>%
  filter(Genus %in% keep_g) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  mutate(Genus = factor(Genus, levels = rev(Genus)))  # Reverse the order for correct plotting

# Enhanced Plot without the Legend
p4 <- ggplot(data = log2_df, aes(x = Genus, y = log2FoldChange, fill = abs(log2FoldChange))) +
  geom_col(alpha = 0.8) +
  scale_fill_viridis_c(direction = -1) +  # Use viridis color scale and reverse for darker colors for larger values
  labs(
    title = "Log2 Fold-Change in Genus Abundance",
    subtitle = "Comparing case_statusCRC",
    y = "Log2 Fold-Change",
    x = "",
    fill = "Log2 Fold-Change"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(color = "black", size = 10),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.title.x = element_text(size = 12, face = "bold"),
    legend.position = "none",  # Remove the legend
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5)
  ) +
  coord_flip() +  # Flip coordinates for better readability
  geom_hline(yintercept = 0, linetype = "dotted", color = "red")  # Highlight the zero line

# Display the plot
p4
 

  
# Combine plots
combined_plot <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 2, nrow = 2, scale = 0.9, labels = c("A", "B", "C", "D")) +
  theme(plot.background = element_rect(fill = "white"),
        text = element_text(family = "sans"))  # Use 'sans' font family

combined_plot

ggsave("combined_plot4.png", plot = combined_plot, width = 14, height = 10, units = "in", dpi = 600)
 
ui <- fluidPage(
  titlePanel("Combined Topic Modeling and Differential Abundance Plot"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("This app displays the combined plot showing topic modeling and differential abundance analysis."),
      hr(),
      h4("Instructions:"),
      p("The plot below visualizes different analyses:"),
      p("A - Evaluation of Topic Models: Number of Topics vs Metrics"),
      p("B - Differential Abundance of Bacterial Genera in Colorectal Cancer"),
      p("C - Top 20 Terms in Topic X"),
      p("D - Log2 Fold-Change in Genus Abundance"),
      hr()
    ),
    
    # Main panel to display the plot
    mainPanel(
      plotOutput("combinedPlot", height = "800px")
    )
  )
)

# Define server logic
server <- function(input, output) {
  output$combinedPlot <- renderPlot({
    # Combined plot generated from your code
    combined_plot
  })
}

# Run the application
shinyApp(ui = ui, server = server)




