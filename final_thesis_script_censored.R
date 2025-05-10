## SCRIPT DIFFERENTIAL MIRNA EXPRESSION ANALYSIS AND PATHWAY ENRICHMENT
    # including visualization / statistical analysis

setwd("/.../Thesis Data")

#Load Clinical Data Script
source("/.../Clinical_Metadata")


##### 1. Differential Expression Analysis #####

  # 1.1 Design Matrix construction + covariates

    seq_files_df_with_clinical$chemo_radiotherapy[seq_files_df_with_clinical$chemo_radiotherapy %in% "Other"] <- "Chemoradiation <=2 months"
    seq_files_df_with_clinical$RT_group <- seq_files_df_with_clinical %>% {ifelse(.$tot_straledose < 47, "Low", "High")}
    clinical_effect <- ({paste(seq_files_df_with_clinical$chemo_radiotherapy, 
                               seq_files_df_with_clinical$RT_group, seq_files_df_with_clinical$before_after_progression, sep="_")} 
                        %>% {sub("(Baseline|Chemo effect)_(High|Low)", "\\1", .)})
    clinical_effect[clinical_effect %in% c("Chemoradiation <=2 months_High_1", "Chemoradiation <=2 months_Low_1")] <- "Chemoradiation <=2 months_1"
    clinical_effect <- factor(clinical_effect)
    seq_files_df_with_clinical$age_at_t0 <- as.numeric(difftime(seq_files_df_with_clinical$ink_1, seq_files_df_with_clinical$f_dato, units = "days")) / 365.25
    median_age <- median(seq_files_df_with_clinical$age_at_t0)
    seq_files_df_with_clinical$age_factor <- (seq_files_df_with_clinical$age_at_t0 - median_age) / 10

        # Define levels of matrix
        levels(clinical_effect)
        levels(clinical_effect) <- c("Baseline", "ChE", "ChR_LE2_1", "ChR_LE2_H_0", "ChR_LE2_L_0", "ChR_GR2_H_0", "ChR_GR2_H_1", "ChR_GR2_L_0", "ChR_GR2_L_1")
        design_clinical <- model.matrix(~ 0 + clinical_effect 
                                        + factor(seq_files_df_with_clinical$library_prep_plate) 
                                        + seq_files_df_with_clinical$age_factor)
        rownames(design_clinical) <- seq_files_df_with_clinical$ID_sample
        colnames(design_clinical) <- c("Baseline", "ChE", "ChR_LE2_1", "ChR_LE2_H_0", "ChR_LE2_L_0", "ChR_GR2_H_0", "ChR_GR2_H_1", "ChR_GR2_L_0", "ChR_GR2_L_1", "prep_plate2", "prep_plate3", "prep_plate4", "age_factor")
        dim(design_clinical)
        head(design_clinical)
        
    # 1.2 Expression Matrix and Model Fitting
        
      library(limma)
      library(edgeR)
      library(ggplot2)
        
        mature_expr_matrix_controls <- read.table("/.../Expression.Matrix.csv", row.names = 1, sep = "\t", header = TRUE)
        
        # filter out controls and unmatched samples
        control_samples <- grep("Control", colnames(mature_expr_matrix_controls), value = TRUE)
        non_control_samples <- colnames(mature_expr_matrix_controls)[!colnames(mature_expr_matrix_controls) %in% control_samples]
        mature_expr_matrix <- mature_expr_matrix_controls[, non_control_samples]
        colnames(mature_expr_matrix) <- gsub("^X", "", colnames(mature_expr_matrix))
        setdiff(colnames(mature_expr_matrix), rownames(design_clinical))
        mature_expr_matrix_filtered <- mature_expr_matrix[, !colnames(mature_expr_matrix) %in% "unmatched_sample_number"]
        
        # Voom transformation and DGE
        dge <- DGEList(counts = mature_expr_matrix_filtered)
        dge <- calcNormFactors(dge) 
        v_clinical <- voom(dge, design_clinical, plot = TRUE)
        
        # Filter out low expressed miRNAs under 4 reads per million
        v_clinical_rpm1 <- v_clinical %>% {mean.rpm <- rowSums(2**.$E) / ncol(.); .[mean.rpm >= 4, ]}
        dge_rpm1 <- dge [rownames(v_clinical_rpm1), ]
        v_rpm1 <- voom(dge_rpm1, design_clinical, plot = TRUE)
        
        # Duplicate Correlation
        dupcor_clinical <- duplicateCorrelation(v_rpm1, design_clinical, block = seq_files_df_with_clinical$studieid)
        fit_clinical <- lmFit(v_rpm1, design_clinical , block=seq_files_df_with_clinical$studieid, correlation=dupcor_clinical$consensus) 
        
    # 1.3 Create Contrasts
        
        contrast_matrix_clinical <- makeContrasts(
          
                # 1. Baseline vs Treatment/Timepoints
                BaseVsChE = ChE - Baseline,
                BaseVsLE2_H = ChR_LE2_H_0 - Baseline,
                BaseVsLE2_L = ChR_LE2_L_0 - Baseline,
                BaseVsGR2_H = 0.5 * (ChR_GR2_H_0 + ChR_GR2_H_1) - Baseline,
                BaseVsGR2_L = 0.5 * (ChR_GR2_L_0 + ChR_GR2_L_1) - Baseline,
                BaseVsH = (1/3) * (ChR_LE2_H_0 + ChR_GR2_H_0 + ChR_GR2_H_1) - Baseline,
                BaseVsL = (1/3) * (ChR_LE2_L_0 + ChR_GR2_L_0 + ChR_GR2_L_1) - Baseline,
                
                # 2. Chemo vs Dose effect
                ChEVsH = (1/3) * (ChR_LE2_H_0 + ChR_GR2_H_0 + ChR_GR2_H_1) - ChE,
                ChEVsL = (1/3) * (ChR_LE2_L_0 + ChR_GR2_L_0 + ChR_GR2_L_1) - ChE,
                
                # 3. Progression vs No Progression
                ProgVsNoprog_tot = (1/3) * (ChR_LE2_1 + ChR_GR2_L_1 + ChR_GR2_H_1) -
                  0.25 * (ChR_LE2_L_0 + ChR_LE2_H_0 + ChR_GR2_L_0 + ChR_GR2_H_0),
                
                ProgVsNoprog_ChE = (1/3) * (ChR_LE2_1 + ChR_GR2_L_1 + ChR_GR2_H_1) -
                  (1/5) * (ChR_LE2_L_0 + ChR_LE2_H_0 + ChR_GR2_L_0 + ChR_GR2_H_0 + ChE),
                
                ProgVsNoprog_GR2 = 0.5 * (ChR_GR2_H_1 + ChR_GR2_L_1) - 0.5 * (ChR_GR2_H_0 + ChR_GR2_L_0),
                ProgVsNoprog_LE2 = ChR_LE2_1 - 0.5 * (ChR_LE2_H_0 + ChR_LE2_L_0),
                
                # 4. Baseline vs Progression / Noprog
                BaseVsProg_tot = (1/3) * (ChR_LE2_1 + ChR_GR2_L_1 + ChR_GR2_H_1) - Baseline,
                BaseVsNoprog_tot = 0.25 * (ChR_LE2_L_0 + ChR_LE2_H_0 + ChR_GR2_L_0 + ChR_GR2_H_0) - Baseline,
                Base_Vs_Prog_H = ChR_GR2_H_1 - Baseline,
                BaseVsNoprog_LE2 = 0.5 * (ChR_LE2_L_0 + ChR_LE2_H_0) - Baseline,
                BaseVsNoprog_GR2 = 0.5 * (ChR_GR2_L_0 + ChR_GR2_H_0) - Baseline,
                
                # 5. Chemo vs Progression / Noprog
                ChEVsProg_tot = (1/3) * (ChR_LE2_1 + ChR_GR2_L_1 + ChR_GR2_H_1) - ChE,
                ChEVsNoprog_tot = 0.25 * (ChR_LE2_L_0 + ChR_LE2_H_0 + ChR_GR2_L_0 + ChR_GR2_H_0) - ChE,
                
                # 6. Interaction contrasts
                ProgVsNoprog_HvsL = (ChR_GR2_H_1 - ChR_GR2_H_0) - (ChR_GR2_L_1 - ChR_GR2_L_0),
                Prog_LE2vsGR2 = ChR_LE2_1 - 0.5 * (ChR_GR2_H_1 + ChR_GR2_L_1),
                ProgVsNoprog_VsBase = ((1/3) * (ChR_LE2_1 + ChR_GR2_L_1 + ChR_GR2_H_1) - Baseline) -
                  (0.25 * (ChR_LE2_L_0 + ChR_LE2_H_0 + ChR_GR2_L_0 + ChR_GR2_H_0) - Baseline),
                
                levels = design_clinical
              )
        
        # Fit Contrasts
        fit2_clinical <- contrasts.fit(fit_clinical, contrast_matrix_clinical)
        fit2_clinical <- eBayes(fit2_clinical)
        
    # 1.4 Extract Differential Expression Results
        
        results_list_clinical <- list()
        
        for (contrast in colnames(contrast_matrix_clinical)) {
          top_table <- topTable(fit2_clinical, coef = contrast, number = Inf)
          results_list_clinical[[contrast]] <- top_table
        }
        
        names(results_list_clinical)
        
    # 1.5 Visualization of DEA Results
        
        # 1.5.1 Volcano Plots
        
        output_dir <- "/.../Volkano Plots/"
        
        volcano_objects <- list()
        
        volcano_plot_clinical <- function(df, contrast_name) {
          ggplot(df, aes(x = logFC, y = -log10(adj.P.Val))) +
            geom_point(aes(color = adj.P.Val < 0.05), alpha = 0.6, size = 0.7) +
            scale_color_manual(values = c("gray", "red")) +
            theme_minimal() +
            labs(
              title = NULL,
              x = "Log Fold Change (logFC)",
              y = "-log10(Adjusted P-value)"
            ) +
            theme_minimal(base_size = 10) +
            theme(legend.position = "right")
        }
        
        for (contrast in names(results_list_clinical)) {
          df <- results_list_clinical[[contrast]]
          plot <- volcano_plot_clinical(df, contrast)
          volcano_objects[[contrast]] <- plot
          ggsave(
            filename = paste0(output_dir, contrast, "_volcano.png"),
            plot = plot,
            device = "png",
            width = 2.5,
            height = 3,
          )
        }
        
        
        # 1.5.2 Heatmaps
        
        library(ComplexHeatmap)
        library(grid)
        library(circlize)
        library(RColorBrewer)
        
        output_dir <- "/.../Heatmaps/"
        
        column_split_vector <- rep(NA, 291)
        
        column_split_vector[clinical_effect == "Baseline"] <- "Baseline"  
        column_split_vector[clinical_effect == "ChE"] <- "Chemo Effect"  
        column_split_vector[clinical_effect %in% c("ChR_LE2_1", "ChR_LE2_H_0", "ChR_LE2_L_0")] <- "<= 2 months"  
        column_split_vector[clinical_effect %in% c("ChR_GR2_H_0", "ChR_GR2_H_1", "ChR_GR2_L_0", "ChR_GR2_L_1")] <- "> 2 months"
        
        column_split_vector <- factor(column_split_vector, 
                                      levels = c("Baseline", 
                                                 "Chemo Effect", 
                                                 "<= 2 months", 
                                                 "> 2 months"))
        
        
        all_sig_miRNAs <- list()
        all_up_miRNAs <- list()
        all_down_miRNAs <- list()
        heatmap_objects <- list()
        
        for (contrast in colnames(contrast_matrix_clinical)) {
          sig_miRNA <- results_list_clinical[[contrast]] %>%
            filter(adj.P.Val <= 0.05)
          
          up_miRNA <- sig_miRNA %>% filter(logFC > 0) %>% arrange(desc(logFC))
          down_miRNA <- sig_miRNA %>% filter(logFC < 0) %>% arrange(logFC)
          
          sig_miRNA_ordered <- bind_rows(up_miRNA, down_miRNA)
          
          
          if (nrow(sig_miRNA) > 0) {
            all_sig_miRNAs[[contrast]] <- rownames(sig_miRNA)
            
            heatmap_clinical_data <- v_clinical$E[rownames(sig_miRNA_ordered), ]
            heatmap_clinical_data[is.na(heatmap_clinical_data)] <- 0
            heatmap_clinical_matrix <- as.matrix(heatmap_clinical_data)
            
            colour_mapping <- 
              c(
                "Baseline" = "red",
                "ChE" = "darkred",
                "ChR_LE2_1" = "darkorange2",
                "ChR_LE2_H_0" = "goldenrod",
                "ChR_LE2_L_0" = "darkkhaki",
                "ChR_GR2_H_0" = "darkgrey",
                "ChR_GR2_H_1" = "steelblue",
                "ChR_GR2_L_0" = "dodgerblue4",
                "ChR_GR2_L_1" = "blueviolet"
              )
            
            ha = HeatmapAnnotation(
              bar = clinical_effect,
              col = list(
                bar = colour_mapping,
                split = c(
                  "Baseline" = "seagreen",
                  "Chemo Effect" = "lightgreen",
                  "<= 2 months" = "cadetblue",
                  "> 2 months" = "cadetblue3"
                )),
              annotation_legend_param = list(
                title_gp = gpar(fontsize = 8, fontface = "bold"),
                labels_gp = gpar(fontsize = 8)
              ),
              split = column_split_vector
            )
            
            regulation_status <- ifelse(sig_miRNA_ordered$logFC > 0, "Up", "Down")
            
            row_ha <- rowAnnotation(
              Regulation = regulation_status,
              col = list(Regulation = c("Up" = "firebrick3", "Down" = "steelblue")),
              show_annotation_name = FALSE,
              annotation_legend_param = list(title = "Regulation", labels_gp = gpar(fontsize = 8))
            )
            
            
            heatmap_clinical <- Heatmap(
              t(scale(t(heatmap_clinical_matrix), scale = FALSE)),
              name = " ",
              heatmap_legend_param = list(direction = "horizontal",
                                          legend_position = "bottom",
                                          legend_margin = unit(c(0, 0, 0, 1), "cm"),
                                          title = "Expression",
                                          title_position = "topcenter",
                                          labels_gp = gpar(fontsize = 8)),
              top_annotation = ha,
              left_annotation = row_ha,
              row_title = " ",
              column_title = " ",
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              column_order = colnames(heatmap_clinical_matrix)[order(clinical_effect)],
              column_split = column_split_vector,
              row_dend_reorder = TRUE,
              show_row_names = TRUE,
              show_column_names = FALSE,
              row_names_gp = gpar(fontsize = 8), 
              col = colorRampPalette(c("blue", "gray", "red"))(100),
            )
            
            heatmap_objects[[contrast]] <- heatmap_clinical
          }
        }
        
    # 1.6 Significant miRNAs - grouping for later visualization
        
        output_dir <- "/.../Mirna Plots/"
        
            for (contrast in colnames(contrast_matrix_clinical)) {sig_miRNA <- results_list_clinical[[contrast]] %>% {.[.$adj.P.Val <= 0.05, ]}
            if (nrow(sig_miRNA) > 0) {
              all_sig_miRNAs[[contrast]] <- rownames(sig_miRNA)
              up_miRNA <- rownames(sig_miRNA[sig_miRNA$logFC > 0,])
              down_miRNA <- rownames(sig_miRNA[sig_miRNA$logFC < 0,])
              all_up_miRNAs[[contrast]] <- up_miRNA
              all_down_miRNAs[[contrast]] <- down_miRNA
            }
            }
        
        # Create vector containing all significant miRNAs
        
        all_sig_miRNAs_vector <- unique(unlist(all_sig_miRNAs))
        sig_miRNA_data <- v_rpm1$E[all_sig_miRNAs_vector, ]
        sig_miRNA_data_df <- as.data.frame(sig_miRNA_data)
        sig_miRNA_data_df$miRNA_ID <- rownames(sig_miRNA_data_df)
        
        # add columns to clinical data for x axis of future plots 
        
        seq_files_df_with_clinical$timepoint <- seq_files_df_with_clinical %>% {paste(.$chemo_radiotherapy, .$before_after_progression, sep="_")} %>% {sub("(Baseline|Chemo effect)_(High|Low)", "\\1", .)} %>% {sub(".*_1$", "Progression", .)}
        seq_files_df_with_clinical$days_since_t0 <- as.numeric(seq_files_df_with_clinical$provedato_corr - seq_files_df_with_clinical$ink_1)
        
        
        library(reshape2)
        library(dplyr)
        
        # create a new dataframe for the significant miRNAs and merge with clinical metadata
        
        sig_miRNA_long <- melt(sig_miRNA_data_df, id.vars = "miRNA_ID", variable.name = "ID_sample", value.name = "Expression")
        
        clinical_metadata_subset <- seq_files_df_with_clinical[, c(
          "ID_sample",
          "id_person_corr",
          "timepoint",
          "days_since_t0",
          "before_after_progression",
          "RT_group"
        )]
        
        sig_miRNA_long <- merge(sig_miRNA_long, clinical_metadata_subset, by = "ID_sample")
        miRNAs <- unique(sig_miRNA_long$miRNA_ID)
        
        # 1.6.1 UpSet
        
            comb_matrix <- make_comb_mat(all_sig_miRNAs)
            UpSet(comb_matrix, 
                  comb_order = order(comb_size(comb_matrix), decreasing = TRUE),
                  top_annotation = HeatmapAnnotation(
                    "Intersection Size" = anno_barplot(
                      comb_size(comb_matrix),
                      gp = gpar(fill = "skyblue"),
                      axis_param = list(gp = gpar(fontsize = 10)) 
                    ),
                    annotation_name_side = "left", 
                    annotation_name_gp = gpar(fontsize = 10),
                    annotation_name_rot = 0
                  ),
                  right_annotation = rowAnnotation(
                    "Set Size" = anno_barplot(
                      set_size(comb_matrix),
                      gp = gpar(fill = "lightgrey"),
                      axis_param = list(gp = gpar(fontsize = 10)) 
                    ),
                    annotation_name_gp = gpar(fontsize = 10)
                  ),
                  column_title_gp = gpar(fontsize = 10),
                  row_names_gp = gpar(fontsize = 10)    
            )
        
        # 1.6.2 Grouping of miRNAs
        
          # GROUP 1 - Regulation at progression
          prog_results <- topTable(fit2_clinical, coef = "ProgVsNoprog_GR2", adjust.method = "BH", number = Inf)
          sig_miRNAs_prog <- prog_results[prog_results$adj.P.Val < 0.05, ]
          sig_miRNAs_prog$miRNA_ID <- rownames(sig_miRNAs_prog)
          sig_miRNAs_prog <- sig_miRNAs_prog[order(sig_miRNAs_prog$P.Value), ]
          up_prog_miRNAs <- rownames(sig_miRNAs_prog[sig_miRNAs_prog$logFC > 0, ])
          down_prog_miRNAs <- rownames(sig_miRNAs_prog[sig_miRNAs_prog$logFC < 0, ])
          
          prog_miRNAs_list <- list(upregulated = up_prog_miRNAs, downregulated = down_prog_miRNAs)
          print(prog_miRNAs_list)
          
          
          # GROUP 2 - Interaction contrast
          interaction_results <- topTable(fit2_clinical, coef = "ProgVsNoprog_HvsL", adjust.method = "BH", number = Inf)
          sig_miRNAs_interaction <- interaction_results[interaction_results$adj.P.Val < 0.05, ]
          sig_miRNAs_interaction$miRNA_ID <- rownames(sig_miRNAs_interaction)
          sig_miRNAs_interaction <- sig_miRNAs_interaction[order(sig_miRNAs_interaction$P.Value), ]
          up_interact_miRNAs <- rownames(sig_miRNAs_interaction[sig_miRNAs_interaction$logFC > 0, ])
          down_interact_miRNAs <- rownames(sig_miRNAs_interaction[sig_miRNAs_interaction$logFC < 0, ])
          
          # GROUP 3 - Base versus no progression
          noprog_results <- topTable(fit2_clinical, coef = "BaseVsNoprog_tot", adjust.method = "BH", number = Inf)
          sig_miRNAs_noprog <- noprog_results[noprog_results$adj.P.Val < 0.05, ]
          sig_miRNAs_noprog$miRNA_ID <- rownames(sig_miRNAs_noprog)
          sig_miRNAs_noprog <- sig_miRNAs_noprog[order(sig_miRNAs_noprog$P.Value), ]
          up_noprog_miRNAs <- rownames(sig_miRNAs_noprog[sig_miRNAs_noprog$logFC > 0, ])
          down_noprog_miRNAs <- rownames(sig_miRNAs_noprog[sig_miRNAs_noprog$logFC < 0, ])
          
          # GROUP 4 - remainder
          groups_significant_miRNAs <- unique(c(rownames(sig_miRNAs_prog), rownames(sig_miRNAs_interaction), rownames(sig_miRNAs_noprog)))
          all_remaining_miRNAs <- setdiff(miRNAs, groups_significant_miRNAs)
        
        
    # 1.7 Expression Boxplots per Group
        
        # all significant miRNAs
        
        ggplot(sig_miRNA_long,
               aes(x = timepoint, y = Expression, fill = timepoint)) +
          geom_boxplot() +
          facet_wrap(~ miRNA_ID, scales = "free_y", ncol = 6) +
          labs(title = " ",
               x = " ",
               y = "Expression Level",
               fill = "Timepoint") +
          theme_minimal() +
          theme(
            legend.position = "bottom",
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            strip.text = element_text(size = 8)
          ) +
          scale_fill_manual(
            values = c(
              "Baseline_0" = "lightsalmon",
              "Chemo effect_0" = "lightgoldenrod",
              "Chemoradiation <=2 months_0" = "lightgreen",
              "Chemoradiation >2 months_0" = "lightblue",
              "Progression" = "lightpink3"
            ),
            labels = c(
              "Baseline_0" = "Baseline",
              "Chemo effect_0" = "Chemo effect",
              "Chemoradiation <=2 months_0" = "Chemoradiation \u2264 2 months",
              "Chemoradiation >2 months_0" = "Chemoradiation > 2 months",
              "Progression" = "Progression"
            )
          )
        
        
        # GROUP 1 - Box plot of miRNAs up- or downregulated at progression (greater than 2 months)
        
        up_miRNAs_long <- sig_miRNA_long %>% filter(miRNA_ID %in% up_prog_miRNAs)
        up_miRNAs_long <- up_miRNAs_long %>% left_join(sig_miRNAs_prog[, c("adj.P.Val", "miRNA_ID")], by = "miRNA_ID")
        up_miRNAs_long <- up_miRNAs_long %>% mutate(miRNA_ID = factor(miRNA_ID, levels = sig_miRNAs_prog$miRNA_ID))
        
        up_miRNAs_long <- up_miRNAs_long[order(up_miRNAs_long$adj.P.Val), ]
        
        plot_up_miRNAs_progression <- ggplot(up_miRNAs_long,
                                             aes(x = timepoint, y = Expression, fill = timepoint)) +
          geom_boxplot() +
          facet_wrap( ~ miRNA_ID, scales = "free_y", ncol = 7) +
          labs(title = " ",
               x = " ",
               y = " ") +
          theme_minimal() +
          theme(legend.position = "right", axis.text.x = element_blank()) +
          scale_fill_manual(
            values = c(
              "Baseline_0" = "lightsalmon",
              "Chemo effect_0" = "lightgoldenrod",
              "Chemoradiation <=2 months_0" = "lightgreen",
              "Chemoradiation >2 months_0" = "lightblue",
              "Progression" = "lightpink3"
            ),
            labels = c(
              "Baseline_0" = "Baseline",
              "Chemo effect_0" = "Chemo effect",
              "Chemoradiation <=2 months_0" = "Chemoradiation \u2264 2 months",
              "Chemoradiation >2 months_0" = "Chemoradiation > 2 months",
              "Progression" = "Progression"
            )
          )
        
        
        down_miRNAs_long <- sig_miRNA_long %>% filter(miRNA_ID %in% down_prog_miRNAs)
        down_miRNAs_long <- down_miRNAs_long %>% left_join(sig_miRNAs_prog[, c("adj.P.Val", "miRNA_ID")], by = "miRNA_ID")
        down_miRNAs_long <- down_miRNAs_long %>% mutate(miRNA_ID = factor(miRNA_ID, levels = sig_miRNAs_prog$miRNA_ID))
        
        down_miRNAs_long <- down_miRNAs_long[order(down_miRNAs_long$adj.P.Val), ]
        
        plot_down_miRNAs_progression <- ggplot(down_miRNAs_long,
                                               aes(x = timepoint, y = Expression, fill = timepoint)) +
          geom_boxplot() +
          facet_wrap( ~ miRNA_ID, scales = "free_y", ncol = 7) +
          labs(title = " ",
               x = "Time Point",
               y = "Expression Level") +
          theme_minimal() +
          theme(legend.position = "right", axis.text.x = element_blank()) +
          scale_fill_manual(
            values = c(
              "Baseline_0" = "lightsalmon",
              "Chemo effect_0" = "lightgoldenrod",
              "Chemoradiation <=2 months_0" = "lightgreen",
              "Chemoradiation >2 months_0" = "lightblue",
              "Progression" = "lightpink3"
            ),
            labels = c(
              "Baseline_0" = "Baseline",
              "Chemo effect_0" = "Chemo effect",
              "Chemoradiation <=2 months_0" = "Chemoradiation \u2264 2 months",
              "Chemoradiation >2 months_0" = "Chemoradiation > 2 months",
              "Progression" = "Progression"
            )
          )
        
        
        # GROUP 2 - Box plot of miRNAs up- or downregulated in the interaction contrast
        
        up_miRNAs_interact_long <- sig_miRNA_long %>% filter(miRNA_ID %in% up_interact_miRNAs)
        up_miRNAs_interact_long <- up_miRNAs_interact_long %>% left_join(sig_miRNAs_interaction[, c("adj.P.Val", "miRNA_ID")], by = "miRNA_ID")
        up_miRNAs_interact_long <- up_miRNAs_interact_long %>% mutate(miRNA_ID = factor(miRNA_ID, levels = sig_miRNAs_interaction$miRNA_ID))
        
        up_miRNAs_interact_long <- up_miRNAs_interact_long[order(up_miRNAs_interact_long$adj.P.Val), ]
        
        plot_up_miRNAs_interaction <- ggplot(up_miRNAs_interact_long,
                                             aes(x = timepoint, y = Expression, fill = RT_group)) +
          geom_boxplot() +
          facet_wrap( ~ miRNA_ID, scales = "free_y", ncol = 3) +
          scale_x_discrete(
            labels = c(
              "Baseline_0" = "Baseline",
              "Chemo effect_0" = "Chemo effect",
              "Chemoradiation <=2 months_0" = "Chemoradiation \u2264 2 months",
              "Chemoradiation >2 months_0" = "Chemoradiation > 2 months",
              "Progression" = "Progression")) +
          labs(title = " ",
               x = "Time Point",
               y = "Expression Level") +
          theme_minimal() +
          theme(
            legend.position = "right",
            axis.text.x = element_text(
              size = 10,
              angle = 90,
              vjust = 0.5,
              hjust = 1
            )
          )
        
        
        down_miRNAs_interact_long <- sig_miRNA_long %>% filter(miRNA_ID %in% down_interact_miRNAs)
        down_miRNAs_interact_long <- down_miRNAs_interact_long %>% left_join(sig_miRNAs_interaction[, c("adj.P.Val", "miRNA_ID")], by = "miRNA_ID")
        down_miRNAs_interact_long <- down_miRNAs_interact_long %>% mutate(miRNA_ID = factor(miRNA_ID, levels = sig_miRNAs_interaction$miRNA_ID))
        
        down_miRNAs_interact_long <- down_miRNAs_interact_long[order(down_miRNAs_interact_long$adj.P.Val), ]
        
        plot_down_miRNAs_interaction <- ggplot(down_miRNAs_interact_long,
                                               aes(x = timepoint, y = Expression, fill = RT_group)) +
          geom_boxplot() +
          facet_wrap( ~ miRNA_ID, scales = "free_y") +
          scale_x_discrete(
            labels = c(
              "Baseline_0" = "Baseline",
              "Chemo effect_0" = "Chemo effect",
              "Chemoradiation <=2 months_0" = "Chemoradiation \u2264 2 months",
              "Chemoradiation >2 months_0" = "Chemoradiation > 2 months",
              "Progression" = "Progression")) +
          labs(title = " ",
               x = "Time Point",
               y = "Expression Level") +
          theme_minimal() +
          theme(
            legend.position = "right",
            axis.text.x = element_text(
              size = 8,
              angle = 90,
              vjust = 0.5,
              hjust = 1
            )
          )
        
        # GROUP 3 - Box plot of miRNAs up- or downregulated at no progression compared to baseline
        
        up_miRNAs_noprog_long <- sig_miRNA_long %>% filter(miRNA_ID %in% up_noprog_miRNAs)
        up_miRNAs_noprog_long <- up_miRNAs_noprog_long %>% left_join(sig_miRNAs_noprog[, c("adj.P.Val", "miRNA_ID")], by = "miRNA_ID")
        up_miRNAs_noprog_long <- up_miRNAs_noprog_long %>% mutate(miRNA_ID = factor(miRNA_ID, levels = sig_miRNAs_noprog$miRNA_ID))
        
        up_miRNAs_noprog_long <- up_miRNAs_noprog_long[order(up_miRNAs_noprog_long$adj.P.Val), ]
        
        plot_up_miRNAs_noprogression <- ggplot(up_miRNAs_noprog_long,
                                               aes(x = timepoint, y = Expression, fill = timepoint)) +
          geom_boxplot() +
          facet_wrap( ~ miRNA_ID, scales = "free_y") +
          labs(title = " ",
               x = " ",
               y = "Expression Level") +
          theme_minimal() +
          theme(legend.position = "right", axis.text.x = element_blank()) +
          scale_fill_manual(
            values = c(
              "Baseline_0" = "lightsalmon",
              "Chemo effect_0" = "lightgoldenrod",
              "Chemoradiation <=2 months_0" = "lightgreen",
              "Chemoradiation >2 months_0" = "lightblue",
              "Progression" = "lightpink3"
            ),
            labels = c(
              "Baseline_0" = "Baseline",
              "Chemo effect_0" = "Chemo effect",
              "Chemoradiation <=2 months_0" = "Chemoradiation \u2264 2 months",
              "Chemoradiation >2 months_0" = "Chemoradiation > 2 months",
              "Progression" = "Progression"
            ))
        
        down_miRNAs_noprog_long <- sig_miRNA_long %>% filter(miRNA_ID %in% down_noprog_miRNAs)
        down_miRNAs_noprog_long <- down_miRNAs_noprog_long %>% left_join(sig_miRNAs_noprog[, c("adj.P.Val", "miRNA_ID")], by = "miRNA_ID")
        down_miRNAs_noprog_long <- down_miRNAs_noprog_long %>% mutate(miRNA_ID = factor(miRNA_ID, levels = sig_miRNAs_noprog$miRNA_ID))
        
        down_miRNAs_noprog_long <- down_miRNAs_noprog_long[order(down_miRNAs_noprog_long$adj.P.Val), ]
        
        plot_down_miRNAs_noprogression <- ggplot(down_miRNAs_noprog_long,
                                                 aes(x = timepoint, y = Expression, fill = timepoint)) +
          geom_boxplot() +
          facet_wrap( ~ miRNA_ID, scales = "free_y", ncol = 6) +
          labs(title = " ",
               x = " ",
               y = "Expression Level") +
          theme_minimal() +
          theme(legend.position = "right", axis.text.x = element_blank()) +
          scale_fill_manual(
            values = c(
              "Baseline_0" = "lightsalmon",
              "Chemo effect_0" = "lightgoldenrod",
              "Chemoradiation <=2 months_0" = "lightgreen",
              "Chemoradiation >2 months_0" = "lightblue",
              "Progression" = "lightpink3"
            ),
            labels = c(
              "Baseline_0" = "Baseline",
              "Chemo effect_0" = "Chemo effect",
              "Chemoradiation <=2 months_0" = "Chemoradiation \u2264 2 months",
              "Chemoradiation >2 months_0" = "Chemoradiation > 2 months",
              "Progression" = "Progression"
            ))
        
        
        # GROUP 4 - Box plot remaining miRNAs (including best pval and contrast)
        
        results_list_clinical <- lapply(results_list_clinical, function(x) {x$miRNA_ID <- rownames(x); x}) 
        results_list_long <- melt(results_list_clinical, id.vars=colnames(results_list_clinical[[1]]))
        best_pvals <- results_list_long %>%
          group_by(miRNA_ID) %>%
          summarize(
            best_pvalue = min(adj.P.Val, na.rm = TRUE), 
            best_contrast = L1[which.min(adj.P.Val)] 
          )
        
        
        remaining_miRNAs_long <- sig_miRNA_long %>% filter(miRNA_ID %in% all_remaining_miRNAs)
        remaining_miRNAs_long <- remaining_miRNAs_long %>% left_join(best_pvals, by = "miRNA_ID")
        
        
        plot_remaining_miRNAs_pval <- ggplot(remaining_miRNAs_long,
                                             aes(x = timepoint, y = Expression, fill = timepoint)) +
          geom_boxplot() +
          geom_text(
            aes(
              x = 3,  
              y = max(Expression) - 0.4,  
              label = paste("p =", signif(best_pvalue, digits = 1), "\n", " ", best_contrast)
            ),
            inherit.aes = TRUE,  
            size = 2.5, 
            color = "darkgray",
            hjust = 0.7, 
            vjust = 1,
            fontface = "plain"
          ) +
          facet_wrap( ~ miRNA_ID, scales = "free_y", ncol = 4) +
          scale_x_discrete(
            labels = c(
              "Baseline_0" = "Baseline",
              "Chemo effect_0" = "Chemo effect",
              "Chemoradiation <=2 months_0" = "\u2264 2 months",
              "Chemoradiation >2 months_0" = "> 2 months",
              "Progression" = "Progression")) +
          scale_fill_manual(
            values = c(
              "Baseline_0" = "lightsalmon",
              "Chemo effect_0" = "lightgoldenrod",
              "Chemoradiation <=2 months_0" = "lightgreen",
              "Chemoradiation >2 months_0" = "lightblue",
              "Progression" = "lightpink3"
            )) +
          labs(title = " ",
               x = "Time Point",
               y = "Expression Level") +
          theme_minimal() +
          theme(legend.position = "right", axis.text.x = element_text(angle = 45, hjust = 1))
        
    # 1.8 miRNA Expression over time - single plots
        
        sig_miRNA_long$before_after_progression <- factor(sig_miRNA_long$before_after_progression, levels = c(0, 1))
        sig_miRNA_long$RT_group <- factor(sig_miRNA_long$RT_group, levels = c("High", "Low"))
        
        mirna_plots <- list()
        
        for (current_miRNA in miRNAs) {
          
          mirna_plots[[current_miRNA]] <- ggplot(
            sig_miRNA_long %>% filter(miRNA_ID == current_miRNA),
            aes(x = days_since_t0,
                y = Expression,
                color = factor(before_after_progression))) +
            geom_point(size = 0.3,
                       alpha = 0.7) +
            geom_line(aes(group = id_person_corr), linewidth = 0.2, alpha = 0.5) +
            facet_wrap( ~ RT_group, labeller = 
                          labeller(RT_group = c("High" = "60 Gy",
                                                "Low" = "45 Gy")),
                        ncol = 1) +
            labs(title = paste(" ", current_miRNA, " "),
                 x = "Days since t0",
                 y = "Expression Level",
                 color = "Progression Status") +
            scale_color_manual(
              values = c("0" = "gray", "1" = "red"),
              labels = c("Not Progressed", "Progressed")) +
            theme_minimal(base_size = 10) +
            theme(
              axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1),
            ) +
            theme(legend.position = "bottom")
        }
        
        for (current_miRNA in miRNAs) {
          
          ggsave(
            paste0(output_dir, "expression_", current_miRNA, "by_dosage.png"), 
            plot = mirna_plots[[current_miRNA]], 
            width = 3, 
            height = 4
          )
        }
        
        # 1.8.1  grid of 8 miRNAS over time (for Results contrast B)
        
        library(patchwork)
        
        miRNA_grid_plots <- mirna_plots[c("hsa-miR-xxx", "hsa-miR-xxx", "hsa-miR-xxx", "hsa-miR-xxx", 
                                          "hsa-miR-xxx", "hsa-miR-xxx", "hsa-miR-xxx", "hsa-miR-xxx")]
        
        miRNA_grid_plots <- lapply(miRNA_grid_plots, function(p) {
          p + theme(
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = "none"
          )
        })
        
        miRNA_grid_plots[[4]] <- miRNA_grid_plots[[4]] +
          labs(y = "Expression Level") +
          theme(axis.title.y = element_text(size = 10, angle = 90, margin = margin(r = 10)))
        
        miRNA_grid_plots[[8]] <- miRNA_grid_plots[[8]] +
          labs(x = "Days since t0") +
          theme(axis.title.x = element_text(size = 10, margin = margin(b = 20)))
        
        
        miRNA_grid <- (miRNA_grid_plots[[1]] + miRNA_grid_plots[[2]] + miRNA_grid_plots[[3]]) /
          (miRNA_grid_plots[[4]] + miRNA_grid_plots[[5]] + miRNA_grid_plots[[6]]) /
          (miRNA_grid_plots[[7]] + miRNA_grid_plots[[8]] + plot_spacer()) +
          plot_layout(widths = c(1, 1, 1), heights = c(1.5, 1.5, 1.5), guides = "collect") &
          theme(
            legend.position = "bottom",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            axis.text.x = element_text(margin = margin(b = 10), size = 8, angle = 45, vjust = 0.5, hjust = 0.5),
            axis.text.y = element_text(size = 8),
            plot.title = element_text(hjust = 0.5, size = 10, face = "plain")
          )
        
        miRNA_grid <- miRNA_grid + plot_annotation(tag_levels = 'A')
        
    # 1.9 Extracting Results (tables, additional information)
        
        sig_miRNAs_table <- data.frame()
        
        for (contrast in names(all_sig_miRNAs)) {
          sig_miRNAs <- all_sig_miRNAs[[contrast]]
          sig_data <- results_list_clinical[[contrast]][rownames(results_list_clinical[[contrast]]) %in% sig_miRNAs,]
          
          sig_miRNAs_table <- rbind(
            sig_miRNAs_table,
            data.frame(
              Contrast = contrast,
              miRNA = rownames(sig_data),
              logFC = sig_data$logFC,
              adj.P.Val = sig_data$adj.P.Val
            )
          )
        }
        
        sig_miRNAs_table <- sig_miRNAs_table[order(sig_miRNAs_table$`adj.P.Val`), ]
        
        rownames(sig_miRNAs_table) <- NULL
        sig_miRNAs_table <- sig_miRNAs_table[, c("miRNA", "logFC", "adj.P.Val", "Contrast")]
        
        
            write.csv(sig_miRNAs_table, "/.../top_significant_miRNAs_table.csv", row.names = FALSE)       
            
            sig_miRNAs_table_progvsnoprog_GR2 <- sig_miRNAs_table %>% filter(Contrast == "ProgVsNoprog_GR2")   
            sig_miRNAs_table_progvsnoprog_GR2 <- sig_miRNAs_table_progvsnoprog_GR2[, c("miRNA", "logFC", "adj.P.Val")]
            
            sig_miRNAs_table_ProgVsNoprog_HvsL <- sig_miRNAs_table %>% filter(Contrast == "ProgVsNoprog_HvsL")   
            sig_miRNAs_table_ProgVsNoprog_HvsL <- sig_miRNAs_table_ProgVsNoprog_HvsL[, c("miRNA", "logFC", "adj.P.Val")]
            
            sig_miRNAs_table_BaseVsNoprog_tot <- sig_miRNAs_table %>% filter(Contrast == "BaseVsNoprog_tot")   
            sig_miRNAs_table_BaseVsNoprog_tot <- sig_miRNAs_table_BaseVsNoprog_tot[, c("miRNA", "logFC", "adj.P.Val")]
            
            sig_miRNAs_table_remaining <- sig_miRNAs_table %>% filter(!Contrast %in% c("ProgVsNoprog_GR2", "ProgVsNoprog_HvsL", "BaseVsNoprog_tot"))  
            sig_miRNAs_table_remaining <- sig_miRNAs_table_remaining[, c("miRNA", "logFC", "adj.P.Val", "Contrast")]
            
            
            write.csv(sig_miRNAs_table_progvsnoprog_GR2, "/.../top_prog_miRNAs_table.csv", row.names = FALSE) 
            write.csv(sig_miRNAs_table_ProgVsNoprog_HvsL, "/.../top_interaction_miRNAs_table.csv", row.names = FALSE) 
            write.csv(sig_miRNAs_table_BaseVsNoprog_tot, "/.../top_noprog_miRNAs_table.csv", row.names = FALSE) 
            write.csv(sig_miRNAs_table_remaining, "/.../top_remaining_miRNAs_table.csv", row.names = FALSE)
        
        # 1.9.1   extracting sample data, patient data, normalized age, etc.
            
            seq_files_df_with_clinical$age_at_t0 <- as.numeric(difftime(seq_files_df_with_clinical$ink_1, seq_files_df_with_clinical$f_dato, units = "days")) / 365.25    
            unique_patients <- seq_files_df_with_clinical %>% distinct(id_person_corr, .keep_all = TRUE)     
            median(unique_patients$age_at_t0)   
            min(unique_patients$age_at_t0) 
            max(unique_patients$age_at_t0) 
            
            age_stats <- unique_patients %>% group_by(RT_group == "Low") %>%
              summarize(
                median_age = median(age_at_t0, na.rm = TRUE),
                max_age = max(age_at_t0, na.rm = TRUE),
                min_age = min(age_at_t0, na.rm = TRUE),
                count = n()
              )
            
            unique_patients[unique_patients$age_at_t0 >= 70, ] %>% filter(RT_group == "Low")
            
            unique_patients %>% group_by(sex, RT_group) %>% summarize(count = n())
            
            unique_patients %>% group_by(stage_tnm7_taha, RT_group) %>% summarize(count = n())
            
            #normalized age =^ median age, 1 = +10yrs, -1 = =10yrs  value/10
            median_age <- median(seq_files_df_with_clinical$age_at_t0)
            seq_files_df_with_clinical$age_factor <- (seq_files_df_with_clinical$age_at_t0 - median_age) / 10
            
            #optional - SD
            mean_age <- mean(seq_files_df_with_clinical$age_at_t0, na.rm = TRUE)
            sd_age <- sd(seq_files_df_with_clinical$age_at_t0, na.rm = TRUE)
            seq_files_df_with_clinical$age_normalized <- (seq_files_df_with_clinical$age_at_t0 - mean_age) / sd_age
            seq_files_df_with_clinical$age_deviation <- seq_files_df_with_clinical$age_normalized * sd_age
            
            #### THIS SECTION WAS MAINLY UNUSED BUT LEFT FOR FUTURE REFERENCE ####
            
            seq_files_df_with_clinical$kur1_sample <-
              seq_files_df_with_clinical$provedato_corr > seq_files_df_with_clinical$kur1 &
              seq_files_df_with_clinical$provedato_corr <= seq_files_df_with_clinical$kur2
            
            seq_files_df_with_clinical$kur2_sample <-
              seq_files_df_with_clinical$provedato_corr > seq_files_df_with_clinical$kur2 &
              seq_files_df_with_clinical$provedato_corr <= seq_files_df_with_clinical$kur3
            
            seq_files_df_with_clinical$kur3_sample <-
              seq_files_df_with_clinical$provedato_corr > seq_files_df_with_clinical$kur3 &
              seq_files_df_with_clinical$provedato_corr <= seq_files_df_with_clinical$kur4
            
            seq_files_df_with_clinical$kur4_sample <-
              seq_files_df_with_clinical$provedato_corr > seq_files_df_with_clinical$kur4
            
            
            
            seq_files_df_with_clinical$days_since_kur1 <- ifelse(seq_files_df_with_clinical$kur1_sample,
                                                                 as.numeric(difftime(seq_files_df_with_clinical$provedato_corr, seq_files_df_with_clinical$kur1, units = "days")),
                                                                 NA)
            
            seq_files_df_with_clinical$days_since_kur2 <- ifelse(seq_files_df_with_clinical$kur2_sample,
                                                                 as.numeric(difftime(seq_files_df_with_clinical$provedato_corr, seq_files_df_with_clinical$kur2, units = "days")),
                                                                 NA)
            
            seq_files_df_with_clinical$days_since_kur3 <- ifelse(seq_files_df_with_clinical$kur3_sample,
                                                                 as.numeric(difftime(seq_files_df_with_clinical$provedato_corr, seq_files_df_with_clinical$kur3, units = "days")),
                                                                 NA)
            
            seq_files_df_with_clinical$days_since_kur4 <- ifelse(seq_files_df_with_clinical$kur4_sample,
                                                                 as.numeric(difftime(seq_files_df_with_clinical$provedato_corr, seq_files_df_with_clinical$kur4, units = "days")),
                                                                 NA)
            
            
            mean(seq_files_df_with_clinical$days_since_kur1, na.rm = TRUE)
            mean(seq_files_df_with_clinical$days_since_kur2, na.rm = TRUE)
            mean(seq_files_df_with_clinical$days_since_kur3, na.rm = TRUE)
            mean(seq_files_df_with_clinical$days_since_kur4, na.rm = TRUE)
            
            kur_cols <- c("kur1_sample", "kur2_sample", "kur3_sample", "kur4_sample")
            seq_files_df_with_clinical$chemo_nr <- apply(seq_files_df_with_clinical[, kur_cols], 1, function(row) {
              chemo_num <- which(row == TRUE)
              if (length(chemo_num) == 1) {
                paste("chemo", chemo_num)
              } else {
                NA  
              }
            })
            
            # summary table 
            sig_miRNAs_table_progvsnoprog_GR2 <- sig_miRNAs_table_progvsnoprog_GR2 %>% mutate( Contrast = "A")
            sig_miRNAs_table_ProgVsNoprog_HvsL <- sig_miRNAs_table_ProgVsNoprog_HvsL %>% mutate( Contrast = "B")
            sig_miRNAs_table_BaseVsNoprog_tot <- sig_miRNAs_table_BaseVsNoprog_tot %>% mutate( Contrast = "C")
            
            summary_table_pre <- bind_rows(sig_miRNAs_table_progvsnoprog_GR2, sig_miRNAs_table_ProgVsNoprog_HvsL, sig_miRNAs_table_BaseVsNoprog_tot)
            summary_table_select <- summary_table_pre %>% select(miRNA, Contrast, logFC, adj.P.Val) %>% arrange(miRNA, Contrast, adj.P.Val)
            
            print(summary_table_select)
            write.csv(summary_table_select, "/.../miRNA_summary_table.csv", row.names = FALSE)
            
            
        # 1.9.2 Sample Overview Plot & data (+replicates)
            
            seq_files_df_with_clinical$replicate_number <- ave(
              seq_files_df_with_clinical$ID_sample,
              interaction(
                seq_files_df_with_clinical$id_person_corr,
                seq_files_df_with_clinical$timepoint
              ),
              FUN = seq_along
            )
            
            timepoint_levels <- c(
              "Baseline_0" = 1,
              "Chemo effect_0" = 2,
              "Chemoradiation <=2 months_0" = 3,
              "Chemoradiation >2 months_0" = 4,
              "Progression" = 5
            )
            
            seq_files_df_with_clinical$timepoint_numeric <-
              timepoint_levels[seq_files_df_with_clinical$timepoint]
            
            n_samples <- seq_files_df_with_clinical %>%
              group_by(timepoint, RT_group) %>%
              summarise(
                count_samples = n(),
                days_since_t0 = mean(days_since_t0, na.rm = TRUE),
                .groups = "drop"
              )
            
            
            ggplot(seq_files_df_with_clinical,
                   aes(x = days_since_t0, y = timepoint)) +
              geom_boxplot(alpha = 0.6, fill = "lightblue", outlier.shape = NA) +
              geom_point(aes(y = timepoint_numeric + 0.1 * as.numeric(replicate_number)),
                         size = 2,
                         alpha = 0.3) +
              geom_line(aes(y = timepoint_numeric + 0.1 * as.numeric(replicate_number), 
                            group = id_person_corr), alpha = 0.2) +
              geom_text(data = n_samples,
                        x = 1050,
                        aes(y = timepoint, label = paste0("n = ", count_samples))) +
              coord_cartesian(xlim = c(0, 1100)) +
              facet_wrap( ~ RT_group, scales = "free_y", ncol = 1) +
              labs(
                title = "",
                x = "",
                y = "",
                color = " "
              ) +
              scale_y_discrete(
                labels = c(
                  "Baseline_0" = "Baseline",
                  "Chemo effect_0" = "Chemo effect",
                  "Chemoradiation <=2 months_0" = "\u2264 2 months",
                  "Chemoradiation >2 months_0" = "> 2 months",
                  "Progression" = "Progression"
                )
              ) +
              theme_minimal(base_size = 10) +
              theme(
                axis.text.x = element_text(
                  angle = 45,
                  hjust = 1,
                  size = 10
                ),
                axis.text.y = element_text(size = 10),
                legend.position = "bottom",
                legend.title = element_text(size = 10),
                legend.text = element_text(size = 9),
                strip.text = element_text(size = 12, face = "bold")
              )
            
            # Replicate table
                replicates_table <- seq_files_df_with_clinical %>%
                  group_by(RT_group, timepoint_numeric) %>%
                  summarise(
                    range_replicates = paste0(min(replicate_number), "-", max(replicate_number)),
                    range_samples_per_patient = paste0(min(table(id_person_corr)), "-", max(table(id_person_corr))),
                    total_over1 = sum(table(id_person_corr) > 1),
                    total_samples = n()
                  )
                
                write.csv(replicates_table, "/.../replicates_table.csv")
                
            # More Information mainly for extraction and in-text mentions (added later on)
                
                seq_files_df_with_clinical %>%
                  group_by(timepoint_numeric) %>%
                  summarize(
                    count = n(),  
                    unique_patients = n_distinct(id_person_corr)  
                  )
                
                seq_files_df_with_clinical %>%
                  group_by(visit_corr, chemo_radiotherapy) %>%
                  summarize(
                    unique_patients = n_distinct(id_person_corr)
                  )
                
                
                seq_files_df_with_clinical %>%
                  filter(timepoint_numeric == unique(timepoint_numeric)[5]) %>%  
                  select(timepoint_numeric, id_person_corr, days_since_t0, days_since_kur4, days_since_kur3, days_since_kur2, chemo_radiotherapy) %>% print(n = 62)
                
                
                unique_miRNAs_numbers <- setdiff(rownames(sig_miRNAs_prog), rownames(sig_miRNAs_noprog))
                unique_miRNAs_numbers_2 <- setdiff(unique_miRNAs_numbers, rownames(sig_miRNAs_prog))
                
                overlapping_miRNAs <- Reduce(intersect, list(
                  rownames(sig_miRNAs_prog), 
                  rownames(sig_miRNAs_noprog), 
                  rownames(sig_miRNAs_interaction)
                ))
                
                
                seq_files_df_with_clinical %>%
                  group_by(chemo_radiotherapy, before_after_progression) %>%
                  summarise(count = n(), .groups = "drop")                
          
##### 2. Overrepresentation/Enrichment Analysis #####
                
    library(gprofiler2)
    library(dplyr)
    library(readr)

    output_dir <- "/.../Enrichment Analysis/"


    # 2.1 read in data, define miRNAs

    targetscan_scores_file <- read_tsv("/.../Predicted_Targets_Context_Scores.default_predictions.txt")
    colnames(targetscan_scores_file)
    miRNAs_targetscan <- unique(c(up_prog_miRNAs, down_prog_miRNAs, up_noprog_miRNAs, down_noprog_miRNAs, up_interact_miRNAs, down_interact_miRNAs))  
    
    
    # 2.2 Loop Targets per miRNA and run Enrichment Analysis 
    
      targets_by_miRNA <- list()
          
          for (miR in miRNAs_targetscan) {
            
            targets <- targetscan_scores_file %>%
              filter(miRNA == miR) %>%
              select(`Gene Symbol`, `context++ score`, `weighted context++ score`)
            
            gene_symbols <- unique(targets$`Gene Symbol`)
            
            if (length(gene_symbols) > 0) {
              targets_by_miRNA[[miR]] <- gene_symbols
            } else {
              warning(paste("No targets found for", miR))
            }
          }
          
          print(head(targets_by_miRNA))
          no_target_miRNAs <- setdiff(miRNAs_targetscan, names(targets_by_miRNA))
          length(no_target_miRNAs)

    
    miR_gost <- list()
    
          for (miR in names(targets_by_miRNA)) {
            target_genes <- targets_by_miRNA[[miR]]
            
            print(paste("miRNA:", miR, " | Number of Target Genes:", length(target_genes)))
            print(head(target_genes)) 
            
            if (length(target_genes) > 0) {
              
              gost_result <- gost(
                query = target_genes,
                organism = "hsapiens",
                ordered_query = FALSE,
                multi_query = FALSE,
                significant = TRUE,
                exclude_iea = FALSE,
                measure_underrepresentation = FALSE,
                evcodes = FALSE,
                user_threshold = 0.05,
                correction_method = "g_SCS",
                domain_scope = "custom",
                custom_bg = targetscan_scores_file$`Gene Symbol`,
                numeric_ns = "",
                sources = c("GO:BP", "GO:MF", "GO:CC", "REAC"),
                as_short_link = FALSE,
                highlight = FALSE
              )
              
              if (!is.null(gost_result$result)) {
                miR_gost[[miR]] <- gost_result$result
              } else {
                print(paste("No enrichment found for miRNA:", miR))
              }
            }
          }
          
          print(miR_gost[["hsa-miR-xxx"]])
          length(miR_gost)
          
      # 2.3 Heatmap of Enrichment Results (Top 10)
          
          # Filter top terms for each source
          gost_df_heatmap <- bind_rows(miR_gost, .id = "miRNA")
          
          top_results_gost <- bind_rows(miR_gost, .id = "miRNA") %>%
            group_by(source, term_name) %>%
            summarise(min_p = min(p_value), .groups = "drop") %>% 
            group_by(source) %>%
            slice_min(min_p, n = 10)
          
          # Merge back with miRNA results
          gost_df_miR <- bind_rows(miR_gost, .id = "miRNA") %>%
            filter(term_name %in% top_results_gost$term_name) %>%
            select(miRNA, term_name, p_value) %>%
            spread(term_name, p_value)
          
          # Convert to matrix
          gost_matrix_miR <- as.matrix(gost_df_miR[, -1])  
          rownames(gost_matrix_miR) <- gost_df_miR$miRNA
          
          # -log10 (+prevent inf values)
          gost_matrix_miR <- -log10(gost_matrix_miR)
          gost_matrix_miR[is.infinite(gost_matrix_miR)] <- max(gost_matrix_miR[is.finite(gost_matrix_miR)], na.rm = TRUE)
          
          # heatmap annotation
          heatmap_terms <- colnames(gost_matrix_miR)
          
          # Get source info for those terms
          term_sources <- gost_df_heatmap %>%
            filter(term_name %in% heatmap_terms) %>%
            distinct(term_name, source) %>%
            slice(match(heatmap_terms, term_name))
          
          
          source_colours <- c("GO:BP" = "seagreen", "GO:MF" = "lightgreen", "GO:CC" = "cadetblue", "REAC" = "cadetblue3")
          source_ha <- HeatmapAnnotation(
            Source = term_sources$source,
            col = list(Source = source_colours),
            annotation_legend_param = list(title = "Source", direction = "vertical"),
            show_legend = TRUE
          )
          
          gost_matrix_miR2 <- gost_matrix_miR %>% {.[is.na(.)] <- -10; .}
          
          cor_miR_cols <- cor(gost_matrix_miR2, method = "pearson", use = "pairwise.complete.obs")  
          dist_miR_cols <- as.dist(1 - cor_miR_cols)
          dist_miR_cols[is.na(dist_miR_cols) | is.infinite(dist_miR_cols)] <- max(dist_miR_cols, na.rm = TRUE)
          hclust_cols <- hclust(dist_miR_cols)
          
          # correlation-based distance 1-corr
          cor_miR_matrix <- cor(t(gost_matrix_miR2), method = "pearson", use = "pairwise.complete.obs")  
          dist_miR_matrix <- as.dist(1 - cor_miR_matrix)
          dist_miR_matrix[is.na(dist_miR_matrix) | is.infinite(dist_miR_matrix)] <- max(dist_miR_matrix, na.rm = TRUE)
          hclust_miR <- hclust(dist_miR_matrix)
          
          # heatmap
          
          gost_miR_heatmap <- Heatmap(
            gost_matrix_miR,
            name = "-log10(p-value)",  
            cluster_rows = hclust_miR,  
            cluster_columns = hclust_cols,
            show_column_names = TRUE,
            show_row_names = TRUE,
            column_title = "Enriched Terms",
            row_title = "miRNAs",
            heatmap_legend_param = list(title = "-log10(p-value)", direction = "vertical"),
            top_annotation = source_ha,
            column_names_rot = 90,
            column_names_gp = gpar(fontsize = 8),
            column_names_max_height = unit(10, "cm"),
            row_names_gp = gpar(fontsize = 9),
            row_names_max_width = unit(6, "cm")
          )
          
          png(file.path(output_dir, "top_10_heatmap.png"), width=2500, height=2500, res=300)
          
          draw(
            gost_miR_heatmap,
            heatmap_legend_side = "right",  
            annotation_legend_side = "right",
            merge_legend = TRUE
          )
          
          dev.off()
          
          
    # 2.4 Dotplots of Enrichment Results (top 10)
          
              gost_df_dotplot <- bind_rows(miR_gost, .id = "query") %>%
                select(query, source, term_name, p_value, intersection_size, query_size, term_size, term_id)
              
              gost_df_dotplot <- gost_df_dotplot %>%
                mutate(fold_enrichment = (intersection_size / query_size) / (term_size / total_background)) %>%
                mutate(odds_ratio = (intersection_size / (query_size - intersection_size)) /
                         ((term_size - intersection_size) / (total_background - query_size - term_size + intersection_size)))
              
              gost_df_dotplot$odds_ratio <- ifelse(gost_df_dotplot$odds_ratio > 10000, 10000, gost_df_dotplot$odds_ratio)
              
              top_terms_miR <- gost_df_dotplot %>%
                group_by(source, term_name) %>%
                summarise(min_p = min(p_value), .groups = "drop") %>%
                arrange(min_p) %>%
                distinct(term_name, .keep_all = TRUE) %>%
                group_by(source) %>%
                slice_head(n = 10) %>%
                ungroup()
              
              # only include filtered top terms
              miR_df_input <- gost_df_dotplot %>% filter(term_name %in% top_terms_miR$term_name)
              
              ## GLOBAL LEGEND ##
              # global min/max for -log10(p_value) and log10(odds_ratio)
              global_color_min <- min(-log10(miR_df_input$p_value), na.rm = TRUE)
              global_color_max <- max(-log10(miR_df_input$p_value), na.rm = TRUE)
              
              global_size_min <- min(log10(miR_df_input$odds_ratio), na.rm = TRUE)
              global_size_max <- max(log10(miR_df_input$odds_ratio), na.rm = TRUE)
              
              miR_dotplots <- list()
              miR_sources <- unique(miR_df_input$source)
          
        # 2.4.1 Visualize Dotplots
              
          for (src in miR_sources) {
            
            miR_df_source <- miR_df_input %>% filter(source == src)
            
            p_miR <- ggplot(miR_df_source, aes(
              x = query, 
              y = reorder(term_name, fold_enrichment), 
              color = -log10(p_value),  
              size = log10(odds_ratio) 
            )) +
              geom_point(alpha = 0.9) +  
              scale_color_gradient(
                low = "blue", high = "red",
                limits = c(global_color_min, global_color_max)
              ) +  
              scale_size(
                range = c(2, 5),
                limits = c(global_size_min, global_size_max)
              ) +
              theme_minimal(base_size = 10) +
              labs(
                title = paste("", src),  
                x = "miRNA", 
                y = "Enriched Terms", 
                size = "Odds Ratio", 
                color = "-log10(p-value)"
              ) +
              scale_x_discrete(labels = function(x) gsub("^hsa-", "", x)) +
              theme(
                axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
              )
            
            miR_dotplots[[src]] <- p_miR
            
            ggsave(filename = paste0(output_dir, "dotplot_", src, ".png"), plot = p_miR, height = 3.2
            )
          }
          
          # Check the generated dotplots
          miR_dotplots
          
          # 2.4.2 complete dotplots for Appendix
          miR_dotplots_all <- list()
          miR_df_input_all <- gost_df_dotplot
          
          global_color_min_all <- min(-log10(miR_df_input_all$p_value), na.rm = TRUE)
          global_color_max_all <- max(-log10(miR_df_input_all$p_value), na.rm = TRUE)
          global_size_min_all <- min(log10(miR_df_input_all$odds_ratio), na.rm = TRUE)
          global_size_max_all <- max(log10(miR_df_input_all$odds_ratio), na.rm = TRUE)
          
          for (src in unique(miR_df_input_all$source)) {
            miR_df_source <- miR_df_input_all %>% filter(source == src)
            
            p_all <- ggplot(miR_df_source, aes(
              x = query,
              y = reorder(term_name, fold_enrichment),
              color = -log10(p_value),
              size = log10(odds_ratio)
            )) +
              geom_point(alpha = 0.9) +
              scale_color_gradient(low = "blue", high = "red", limits = c(global_color_min_all, global_color_max_all)) +
              scale_size(range = c(2, 5), limits = c(global_size_min_all, global_size_max_all)) +
              theme_minimal(base_size = 10) +
              labs(
                title = paste("", src, " "),
                x = "miRNA",
                y = "Enriched Terms",
                size = "Odds Ratio",
                color = "-log10(p-value)"
              ) +
              scale_x_discrete(labels = function(x) gsub("^hsa-", "", x)) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
            
            miR_dotplots_all[[src]] <- p_all
          }
          
##### 3. Additional Visualization (Grids) for finalized Thesis #####
          
output_dir <- "/.../Final Visualization Thesis/"
          
    # 3.1 Top 5 miRs per Contrast (Up and Down)
          
          top_5_miRNAs_up_p <- head(unique(up_miRNAs_long$miRNA_ID), 5)
          top5_up_p <- up_miRNAs_long %>% filter(miRNA_ID %in% top_5_miRNAs_up_p)
          
          top_5_miRNAs_down_p <- head(unique(down_miRNAs_long$miRNA_ID), 5)
          top5_down_p <- down_miRNAs_long %>% filter(miRNA_ID %in% top_5_miRNAs_down_p)
          
          top_5_miRNAs_up_i <- head(unique(up_miRNAs_interact_long$miRNA_ID), 5)
          top5_up_i <- up_miRNAs_interact_long %>% filter(miRNA_ID %in% top_5_miRNAs_up_i)
          
          top_5_miRNAs_down_i <- head(unique(down_miRNAs_interact_long$miRNA_ID), 5)
          top5_down_i<- down_miRNAs_interact_long %>% filter(miRNA_ID %in% top_5_miRNAs_down_i)   
          
          top_5_miRNAs_up_np <- head(unique(up_miRNAs_noprog_long$miRNA_ID), 5)
          top5_up_np <- up_miRNAs_noprog_long %>% filter(miRNA_ID %in% top_5_miRNAs_up_np)
          
          top_5_miRNAs_down_np <- head(unique(down_miRNAs_noprog_long$miRNA_ID), 5)
          top5_down_np <- down_miRNAs_noprog_long %>% filter(miRNA_ID %in% top_5_miRNAs_down_np)
          
          # 3.1.1 Contrast A and C
          
              top5_miRNA_list <- list(
                "Contrast A, upregulated" = top5_up_p,
                "Contrast A, downregulated" = top5_down_p,
                "Contrast C, upregulated" = top5_up_np,
                "Contrast C, downregulated" = top5_down_np
              )
              
              top5_plots <- lapply(names(top5_miRNA_list), function(name) {
                top5_df <- top5_miRNA_list[[name]]
                
                ggplot(top5_df, aes(x = timepoint, y = Expression, fill = timepoint)) +
                  geom_boxplot() +
                  facet_wrap(~ miRNA_ID, scales = "free_y", ncol = 5) +
                  labs(title = name , x = "", y = "") +
                  theme_minimal(base_size = 10) +
                  theme(legend.position = "bottom", axis.text.x = element_blank(),
                        legend.box = "vertical",
                  ) +
                  scale_fill_manual(
                    name = "Timepoint",
                    values = c(
                      "Baseline_0" = "lightsalmon",
                      "Chemo effect_0" = "lightgoldenrod",
                      "Chemoradiation <=2 months_0" = "lightgreen",
                      "Chemoradiation >2 months_0" = "lightblue",
                      "Progression" = "lightpink3"
                    ),
                    labels = c(
                      "Baseline_0" = "(1) Baseline",
                      "Chemo effect_0" = "(2) Chemo effect",
                      "Chemoradiation <=2 months_0" = "(3) Chemoradiation \u2264 2 months",
                      "Chemoradiation >2 months_0" = "(4) Chemoradiation > 2 months",
                      "Progression" = "(5) Progression"
                    )
                  )
              })
              
              top5_plots[[4]]
          
          # 3.1.2 Interaction Contrast (B)
          
              top5_interact_list <- list(
                "Contrast B, upregulated" = top5_up_i,
                "Contrast B, downregulated" = top5_down_i)
              
              top5_interact_plots <- lapply(names(top5_interact_list), function(name) {
                top5_interact_df <- top5_interact_list[[name]]
                
                ggplot(top5_interact_df,
                       aes(x = timepoint, y = Expression, fill = RT_group)) +
                  geom_boxplot() +
                  facet_wrap( ~ miRNA_ID, scales = "free_y") +
                  scale_x_discrete(
                    labels = c(
                      "Baseline_0" = "Baseline",
                      "Chemo effect_0" = "Chemo effect",
                      "Chemoradiation <=2 months_0" = "Chemoradiation \u2264 2 months",
                      "Chemoradiation >2 months_0" = "Chemoradiation > 2 months",
                      "Progression" = "Progression")
                  ) +
                  labs(title = name,
                       x = "Time Point",
                       y = "Expression Level") +
                  theme_minimal(base_size = 10) +
                  theme(
                    legend.position = "bottom",
                    axis.text.x = element_text(
                      size = 8,
                      angle = 45,
                      vjust = 1,
                      hjust = 1
                    ),
                    axis.title.x = element_blank()
                  )  
              })
              
    # 3.2 Venn Diagram
              
        library(venn) 
              
          prog_venn <- sig_miRNAs_prog$miRNA_ID
          interact_venn <- sig_miRNAs_interaction$miRNA_ID
          noprog_venn <- sig_miRNAs_noprog$miRNA_ID
              
          list_all_three <- list("Contrast A" = prog_venn,
                                     "Contrast B" = interact_venn,
                                     "Contrast C" = noprog_venn) 
              
              venn(list_all_three,
                   snames = c("Contrast A", "Contrast B", "Contrast C"),
                   ilabels = "counts",                      
                   opacity = 0.3,                      
                   plotsize = 15,                    
                   zcolor = "style",                  
                   borders = TRUE,
                   ellipse = TRUE,
                   box = FALSE,  
                   par = FALSE,
                   ilcs = 1,                
                   sncs = 1,                       
                   cex = 0.8, 
                   main = "" 
              )
              
    # 3.3 Heatmap/Volcano Plot Grids per Contrast
              
       names(volcano_objects)
       names(heatmap_objects)
              
              heatmap_grobs <- lapply(heatmap_objects, function(hm) {
                grid.grabExpr({
                  draw(hm,
                       heatmap_legend_side = "top", annotation_legend_side = "right",
                       show_heatmap_legend = TRUE, show_annotation_legend = TRUE
                  )
                })
              })
              names(heatmap_grobs)
              
      grid_p <- ((volcano_objects[[ "ProgVsNoprog_GR2" ]] + plot_spacer()) / heatmap_grobs[["ProgVsNoprog_GR2"]]) +
      plot_layout(ncol = 1, nrow = 2, widths = c(0.3, 0.3, 0.5), heights = c(0.3, 1, 1)) + plot_annotation(tag_levels = 'A')
   
      grid_np <- ((volcano_objects[[ "BaseVsNoprog_tot" ]] + plot_spacer()) / heatmap_grobs[["BaseVsNoprog_tot"]]) +
      plot_layout(ncol = 1, nrow = 2, widths = c(1, 1), heights = c(0.3, 1, 1)) + plot_annotation(tag_levels = 'A')
    
      grid_i <- ((volcano_objects[[ "ProgVsNoprog_HvsL" ]] + plot_spacer()) / heatmap_grobs[["ProgVsNoprog_HvsL"]]) +
      plot_layout(ncol = 1, nrow = 2, widths = c(1, 1), heights = c(0.3, 1, 1)) + plot_annotation(tag_levels = 'A')         
              
      
    # 3.4 Grids of Top 5 per contrast (UP/DOWN)   
      
      top5_grid_four <- top5_plots[[1]] + top5_plots[[2]] + top5_plots[[3]] + top5_plots[[4]] +
        plot_layout(heights = c(1, 1, 1, 1), widths = c(1, 1, 1, 1), guides = "collect", ncol = 1, nrow = 4) +
        plot_annotation(tag_levels = "A") &
        theme(
          legend.position = "right",
          legend.box = "vertical",
        )
      
      top5_grid_interact <- top5_interact_plots[[1]] / top5_interact_plots[[2]] +
        plot_layout(heights = c(2, 4), widths = c(1, 1, 1), guides = "collect", ncol = 1) +
        plot_annotation(tag_levels = "A") &
        theme(legend.position = "right")
      
    # 3.5 Enrichment Dotplot Grid
      
      enrichment_dotplot_grid <- miR_dotplots[[1]] / miR_dotplots[[2]] / miR_dotplots[[3]] / miR_dotplots[[4]] + plot_layout(ncol = 1, guides = "collect") +
        plot_annotation(tag_levels = "A") &
        theme(legend.position = "right",
              axis.title.y = element_blank(),
              axis.text.x = element_text(angle = 60, hjust = 1),
              axis.text.y = element_text(size = 7.5))
      
    # 3.6 Appendix Tables for Terms per miRNA / Term Counts
      
      target_counts <- sapply(targets_by_miRNA, length)
      
      target_counts_df <- data.frame(
        miRNA_ID = names(target_counts),
        num_targets = as.integer(target_counts)
      )
      
      write.csv(target_counts_df, "/.../targets_per_miRNA.csv", row.names = FALSE)
      
      # Terms Per Source Per miRNA
      
          source_term_data <- do.call(rbind, lapply(names(miR_gost), function(miRNA) {
            enrich_df <- miR_gost[[miRNA]]
            enrich_df$miRNA_ID <- miRNA
            return(enrich_df)
          }))
          
          # Count unique Terms per Source
          terms_by_source <- source_term_data %>%
            group_by(source) %>%
            summarise(unique_terms_count = n_distinct(term_name), .groups = 'drop')
          
          print(terms_by_source)
          
          all_enriched_terms <- unlist(lapply(miR_gost, function(x) x$term_id))
          terms_count <- length(unique(all_enriched_terms))
          print(terms_count)
      
          enr_source_miR <- do.call(rbind, lapply(names(miR_gost), function(miRNA) {
            enrichment_df <- miR_gost[[miRNA]]
            enrichment_df$miRNA <- miRNA
            return(enrichment_df)
          }))
      
      source_counts <- enr_source_miR %>%
        group_by(miRNA, source) %>%
        summarise(unique_terms_count = n_distinct(term_id), .groups = 'drop')
      write.csv(source_counts, "/.../term_persource_permiRNA.csv", row.names = FALSE)
      
##### END OF ANALYSIS #####
save.image(file = "/.../env_11.05.10_RData")
            