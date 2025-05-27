##############################################################################
#  app.R – R-ODAF DESeq2 pipeline                                             #
#          * All filter thresholds (CPM, FDR, etc.) are now fixed values      #
#          * Users cannot modify thresholds in the UI                        #
#          * Top / Avg gene panels use a RANGE slider only                   #
#          * Hovering shows sample names only in the Top‐gene plot           #
#          * Uses normalized counts for expression plots                     #
#          * Wrapped long gene IDs in Top gene facet labels                  #
#          * Rotated x‐axis labels 90° in Avg expression plot  
# ADD option to validate a public dataset first and understand the tool
# ADD option to ensure if the user wants to analyse TempoSeq or RNAseq since 
# current pipeline builds on RNAseq. For tempoSeq, CPM changes to 1 million instead of 
# 5 million. 5 million wipes everything away.
# Add Option to run analysis for multiple conditions,,,,
# (Not only case vs controls but also control vs low, control vs high) and....
# have a panel/drop down option to choose analysis
# Add option to add reordering data, controls first, low dose followed by subsequent dose
##############################################################################

library(shiny)
library(DT)
library(data.table)
library(dplyr)
library(tidyr)
library(DESeq2)
library(edgeR)
library(pheatmap)
library(ggplot2)
library(cowplot)
library(corrplot)
library(PCAtools)
library(plotly)
library(stringr)    # for str_wrap()

options(shiny.maxRequestSize = 100 * 1024^2)

# ─────────────────────────────────────────────────────────────────────────────
# FIXED FILTER PARAMETERS (users cannot change these)
MIN_COVERAGE   <- 5e6    # Min total reads / sample
MIN_CPM        <- 1      # Min CPM threshold
FDR_CUTOFF     <- 0.01   # FDR cutoff (alpha)
VAR_THRESHOLD  <- 20     # Outlier variance cutoff (%)
# ─────────────────────────────────────────────────────────────────────────────

############################  helper functions  ##############################
check_outliers_vst <- function(dds, meta_df, design_col, conds, var_thresh) {
  vsd_mat <- vst(dds, blind = TRUE) %>% assay()
  pca_obj <- pca(vsd_mat, metadata = meta_df, removeVar = 0.1)
  out     <- character()
  pcs     <- which(pca_obj$variance >= var_thresh)
  for (pc in pcs) {
    max_pc <- abs(max(pca_obj$rotated[, pc]) - min(pca_obj$rotated[, pc]))
    for (c in conds) {
      samp_ids <- rownames(meta_df)[meta_df[[design_col]] == c]
      vals     <- pca_obj$rotated[samp_ids, pc, drop = FALSE] %>%
        as.data.frame() %>%
        arrange(desc(.[,1]))
      dmat      <- dist(vals[,1]) %>% as.matrix()
      diag_vals <- dmat[row(dmat) == col(dmat) + 1]
      flags     <- ((pca_obj$variance[pc] * diag_vals) / max_pc) > var_thresh
      if (any(flags)) out <- c(out, rownames(vals)[which(flags)])
    }
  }
  unique(out)
}

##############################  UI  ##########################################
ui <- fluidPage(
  titlePanel("T-ROAM/R-ODAF DESeq2 Pipeline"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("counts_file", "Counts (genes × samples)", accept = c(".csv", ".txt")),
      fileInput("meta_file",   "Metadata (samples × attributes)", accept = c(".csv", ".txt")),
      textInput("design_col",  "Design column",  value = "treatment"),
      textInput("condition1",  "Control label",  value = "control"),
      textInput("condition2",  "Case label",     value = "case"),
      tags$hr(),
      tags$h5("Fixed filter settings:"),
      tags$p(strong("Min total reads / sample:"), format(MIN_COVERAGE, scientific = FALSE)),
      tags$p(strong("Min CPM threshold:"),         MIN_CPM),
      tags$p(strong("FDR cutoff (alpha):"),        FDR_CUTOFF),
      tags$p(strong("Outlier variance cutoff (%):"), VAR_THRESHOLD),
      tags$hr(),
      actionButton("runAnalysis", "Run analysis", class = "btn-primary"),
      width = 3
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Data preview",
                 h4("Counts preview"),   DTOutput("counts_preview"),
                 h4("Metadata preview"), DTOutput("meta_preview")
        ),
        
        tabPanel("PCA & outliers",
                 h4("PCA before outlier removal"), plotlyOutput("pca_before"),
                 verbatimTextOutput("outlier_log"),
                 h4("PCA after outlier removal"),  plotlyOutput("pca_after")
        ),
        
        tabPanel("DE results",
                 verbatimTextOutput("deg_summary"),
                 DTOutput("deg_table"),
                 downloadButton("download_degs", "Download DEGs"),
                 downloadButton("download_norm", "Download Norm Counts"),
                 h4("Heat‐map: all samples"),  plotOutput("heatmap_all"),
                 h4("Heat‐map: final DEGs"),   plotOutput("heatmap_degs")
        ),
        
        tabPanel("Top gene expression",
                 sliderInput("nTopGenes", "Total DEGs to consider",
                             min = 5, max = 100, step = 5, value = 20),
                 uiOutput("topRangeUI"),
                 plotlyOutput("topGenePlots")
        ),
        
        tabPanel("Avg expression",
                 sliderInput("nAvgGenes", "Total DEGs to consider",
                             min = 5, max = 100, step = 5, value = 20),
                 uiOutput("avgRangeUI"),
                 plotlyOutput("avgExprPlot")
        ),
        
        tabPanel("Summary",
                 h4("Analysis summary"),
                 verbatimTextOutput("analysis_summary"),
                 downloadButton("download_report", "Download summary (.txt)")
        )
      )
    )
  )
)

##############################  server  ######################################
server <- function(input, output, session) {
  WINDOW_SIZE <- 5
  
  ## dynamic range sliders
  output$topRangeUI <- renderUI({
    sliderInput("topRange", "Gene range",
                min = 1, max = input$nTopGenes,
                value = c(1, min(WINDOW_SIZE, input$nTopGenes)), step = 1)
  })
  output$avgRangeUI <- renderUI({
    sliderInput("avgRange", "Gene range",
                min = 1, max = input$nAvgGenes,
                value = c(1, min(WINDOW_SIZE, input$nAvgGenes)), step = 1)
  })
  observeEvent(input$nTopGenes, {
    rng <- input$topRange; if (is.null(rng)) rng <- c(1, WINDOW_SIZE)
    rng <- pmin(pmax(rng, 1), input$nTopGenes)
    updateSliderInput(session, "topRange", max = input$nTopGenes, value = rng)
  })
  observeEvent(input$nAvgGenes, {
    rng <- input$avgRange; if (is.null(rng)) rng <- c(1, WINDOW_SIZE)
    rng <- pmin(pmax(rng, 1), input$nAvgGenes)
    updateSliderInput(session, "avgRange", max = input$nAvgGenes, value = rng)
  })
  
  ## read & preview
  counts_data <- reactive({
    req(input$counts_file)
    ext <- tolower(tools::file_ext(input$counts_file$name))
    df  <- if (ext == "csv")
      read.csv(input$counts_file$datapath, row.names = 1, check.names = FALSE)
    else
      read.delim(input$counts_file$datapath, row.names = 1, check.names = FALSE)
    as.matrix(df)
  })
  meta_data <- reactive({
    req(input$meta_file)
    ext <- tolower(tools::file_ext(input$meta_file$name))
    if (ext == "csv")
      read.csv(input$meta_file$datapath, row.names = 1, check.names = FALSE)
    else
      read.delim(input$meta_file$datapath, row.names = 1, check.names = FALSE)
  })
  output$counts_preview <- renderDT({
    datatable(counts_data(), options = list(scrollX = TRUE, pageLength = 5))
  })
  output$meta_preview <- renderDT({
    datatable(meta_data(), options = list(scrollX = TRUE, pageLength = 5))
  })
  
  ## main pipeline
  pipeline <- eventReactive(input$runAnalysis, {
    withProgress(message = "Running R-ODAF analysis", value = 0, {
      # (a) align + NA→0
      cnts0 <- counts_data(); meta0 <- meta_data()
      common <- intersect(colnames(cnts0), rownames(meta0))
      validate(need(length(common) > 1, "Counts / metadata IDs don’t match."))
      cnts <- cnts0[, common, drop = FALSE]
      meta <- meta0[common, , drop = FALSE]
      cnts[is.na(cnts)] <- 0
      
      # (b) low-coverage filter (fixed)
      incProgress(0.10, detail = "Filtering low-coverage samples...")
      keep_samp <- colSums(cnts) > MIN_COVERAGE
      cnts <- cnts[, keep_samp, drop = FALSE]
      meta <- meta[keep_samp, , drop = FALSE]
      
      # (c) initial DESeq2
      incProgress(0.15, detail = "Fitting initial DESeq2...")
      design_f <- as.formula(paste0("~", input$design_col))
      dds_init <- DESeqDataSetFromMatrix(round(cnts), meta, design_f)
      dds_init <- DESeq(dds_init, quiet = TRUE)
      
      # (d) PCA outlier removal
      incProgress(0.20, detail = "PCA outlier check...")
      dds_cur <- dds_init; cnts_cur <- cnts; meta_cur <- meta
      removed <- character(); conds <- c(input$condition1, input$condition2)
      repeat {
        outs <- check_outliers_vst(dds_cur, as.data.frame(colData(dds_cur)),
                                   input$design_col, conds, VAR_THRESHOLD)
        if (!length(outs)) break
        prefixes <- unique(sub("_.*", "", outs))
        mask <- grepl(paste0("^(", paste(prefixes, collapse = "|"), ")_"),
                      colnames(cnts_cur))
        removed <- c(removed, colnames(cnts_cur)[mask])
        cnts_cur <- cnts_cur[, !mask, drop = FALSE]
        meta_cur <- meta_cur[!mask, , drop = FALSE]
        dds_cur  <- DESeqDataSetFromMatrix(round(cnts_cur), meta_cur, design_f)
        dds_cur  <- DESeq(dds_cur, quiet = TRUE)
      }
      
      # (e) 75% expression filter (fixed)
      incProgress(0.25, detail = "75% expression filter...")
      norm_counts <- counts(dds_cur, normalized = TRUE)
      CPMnorm     <- cpm(norm_counts)
      samp_tab    <- table(meta_cur[[input$design_col]])
      keep_genes  <- rownames(norm_counts)[
        sapply(rownames(norm_counts), function(g) {
          any(sapply(names(samp_tab), function(gr) {
            cols <- which(meta_cur[[input$design_col]] == gr)
            sum(CPMnorm[g, cols] >= MIN_CPM) >= 0.75 * length(cols)
          }))
        })
      ]
      
      # (f) final DESeq2 + fixed FDR
      incProgress(0.30, detail = "Final DESeq2 & R-ODAF filters...")
      res_all <- results(dds_cur[keep_genes, ],
                         contrast = c(input$design_col, input$condition2, input$condition1),
                         alpha = FDR_CUTOFF,
                         independentFiltering = FALSE,
                         cooksCutoff = FALSE,
                         pAdjustMethod = "fdr")
      res_df <- as.data.frame(res_all)
      sig_df <- subset(res_df, padj < FDR_CUTOFF)
      
      # R-ODAF extra filters
      if (nrow(sig_df) > 0) {
        Filt <- data.frame(
          low   = rep(1, nrow(sig_df)),
          quant = rep(0, nrow(sig_df)),
          spike = rep(1, nrow(sig_df)),
          row.names = rownames(sig_df), stringsAsFactors = FALSE
        )
        for (g in rownames(sig_df)) {
          gr1 <- which(meta_cur[[input$design_col]] == input$condition1)
          gr2 <- which(meta_cur[[input$design_col]] == input$condition2)
          m1  <- median(norm_counts[g, gr1]); q2 <- quantile(norm_counts[g, gr2], .75)
          m2  <- median(norm_counts[g, gr2]); q1 <- quantile(norm_counts[g, gr1], .75)
          Filt[g, "quant"] <- as.integer((m1 > q2) | (m2 > q1))
          bad_spike <- sapply(list(gr1, gr2), function(cols) {
            v <- norm_counts[g, cols]
            if (max(v) == 0) FALSE else (max(v)/sum(v)) >= 1.4 * (length(cols))^(-0.66)
          })
          Filt[g, "spike"] <- as.integer(sum(bad_spike) <= 1)
        }
        keep_final <- rownames(Filt)[rowSums(Filt) == 3]
        degs_final <- sig_df[keep_final, , drop = FALSE]
      } else {
        degs_final <- sig_df
      }
      
      incProgress(1, detail = "Done.")
      list(
        dds_init  = dds_init,
        dds_final = dds_cur,
        norm      = norm_counts,
        all_res   = res_df,
        degs      = degs_final,
        removed   = unique(removed)
      )
    })
  })
  
  ## ranked genes & window helpers
  ranked_genes <- reactive({
    req(pipeline(), nrow(pipeline()$degs) > 0)
    rownames(pipeline()$degs)[order(pipeline()$degs$padj)]
  })
  
  topGenesWin <- reactive({
    req(ranked_genes(), input$topRange)
    g   <- head(ranked_genes(), input$nTopGenes)
    idx <- seq(input$topRange[1], min(input$topRange[2], length(g)))
    g[idx]
  })
  
  topAvgGenesWin <- reactive({
    req(ranked_genes(), input$avgRange)
    g   <- head(ranked_genes(), input$nAvgGenes)
    idx <- seq(input$avgRange[1], min(input$avgRange[2], length(g)))
    g[idx]
  })
  
  ## Top-gene expression (boxplots) with sample hover
  topExprDF <- reactive({
    req(topGenesWin())
    mat <- pipeline()$norm[topGenesWin(), , drop = FALSE]
    df  <- as.data.frame(mat); df$gene <- rownames(df)
    long <- pivot_longer(df, -gene, names_to = "sample", values_to = "expression")
    meta <- as.data.frame(colData(pipeline()$dds_final)); meta$sample <- rownames(meta)
    left_join(long, meta, by = "sample")
  })
  output$topGenePlots <- renderPlotly({
    req(topExprDF())
    p <- ggplot(topExprDF(),
                aes_string(x = input$design_col, y = "expression",
                           color = input$design_col, text = "sample")) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      facet_wrap(~ gene, ncol = WINDOW_SIZE, scales = "free_y",
                 labeller = labeller(gene = function(x) str_wrap(x, width = 10))) +
      labs(x = input$design_col, y = "Normalized count") +
      theme_minimal() + theme(strip.text = element_text(size = 8))
    ggplotly(p, tooltip = c("text", "gene", "expression"))
  })
  
  ## Avg expression (bars)
  avgExprDF <- reactive({
    req(topAvgGenesWin())
    mat  <- pipeline()$norm[topAvgGenesWin(), , drop = FALSE]
    meta <- as.data.frame(colData(pipeline()$dds_final))
    do.call(rbind, lapply(rownames(mat), function(g) {
      mu <- tapply(mat[g, ], meta[[input$design_col]], mean)
      data.frame(gene = g, condition = names(mu),
                 mean_expression = as.numeric(mu), stringsAsFactors = FALSE)
    }))
  })
  output$avgExprPlot <- renderPlotly({
    req(avgExprDF())
    p <- ggplot(avgExprDF(),
                aes_string(x = "gene", y = "mean_expression", fill = "condition")) +
      geom_col(position = position_dodge()) +
      labs(x = "Gene", y = "Mean normalized count") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            axis.title.x = element_text(margin = margin(t = 10)))
    ggplotly(p, tooltip = c("gene", "condition", "mean_expression"))
  })
  
  ## PCA & outliers
  output$pca_before <- renderPlotly({
    req(pipeline())
    vsd <- vst(pipeline()$dds_init, blind = TRUE) |> assay()
    md  <- as.data.frame(colData(pipeline()$dds_init))
    plotly_build(
      biplot(pca(vsd, metadata = md, removeVar = 0.1),
             x = "PC1", y = "PC2", colby = input$design_col) +
        coord_fixed()
    )
  })
  output$outlier_log <- renderPrint({
    req(pipeline())
    if (length(pipeline()$removed) == 0) cat("No outliers removed.")
    else cat("Removed samples:", paste(pipeline()$removed, collapse = ", "))
  })
  output$pca_after <- renderPlotly({
    req(pipeline())
    vsd2 <- vst(pipeline()$dds_final, blind = TRUE) |> assay()
    md2  <- as.data.frame(colData(pipeline()$dds_final))
    plotly_build(
      biplot(pca(vsd2, metadata = md2, removeVar = 0.1),
             x = "PC1", y = "PC2", colby = input$design_col) +
        coord_fixed()
    )
  })
  
  ## DE results & heatmaps
  output$deg_summary <- renderPrint({
    req(pipeline())
    cat("Total genes tested:", nrow(pipeline()$all_res), "\n")
    cat("Final DEGs (FDR <", FDR_CUTOFF, "):", nrow(pipeline()$degs), "\n")
  })
  output$deg_table <- renderDT({
    req(pipeline())
    datatable(pipeline()$degs, options = list(scrollX = TRUE, pageLength = 5))
  })
  output$heatmap_all <- renderPlot({
    req(pipeline())
    corrplot(cor(pipeline()$norm), method = "color", tl.cex = 0.7)
  })
  output$heatmap_degs <- renderPlot({
    req(pipeline())
    dg <- pipeline()$degs
    if (nrow(dg) < 2) {
      plot.new(); text(0.5, 0.5, "Not enough DEGs for heat-map.")
    } else {
      pheatmap(pipeline()$norm[rownames(dg), , drop = FALSE],
               scale = "row", show_rownames = FALSE)
    }
  })
  
  ## Downloads
  output$download_degs <- downloadHandler(
    filename = function() paste0("DEGs_FDR", FDR_CUTOFF, "_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(pipeline()$degs, file, row.names = TRUE)
  )
  output$download_norm <- downloadHandler(
    filename = function() paste0("norm_counts_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(pipeline()$norm, file, row.names = TRUE)
  )
  
  ## Summary & report (with citation)
  output$analysis_summary <- renderPrint({
    req(pipeline())
    cat("R-ODAF DESeq2 Pipeline summary\n\n")
    cat("Design column:    ", input$design_col, "\n")
    cat("Conditions:       ", input$condition1, "vs", input$condition2, "\n")
    cat("Min coverage:     ", MIN_COVERAGE, "\n")
    cat("Min CPM thresh:   ", MIN_CPM, "\n")
    cat("FDR cutoff:       ", FDR_CUTOFF, "\n")
    cat("Outlier var %:    ", VAR_THRESHOLD, "\n\n")
    cat("Total genes tested:", nrow(pipeline()$all_res), "\n")
    cat("Final DEGs:       ", nrow(pipeline()$degs), "\n\n")
    cat("Please cite:\n")
    cat("Verheijen, M. C., Meier, M. J., Asensio, J. O., Gant, T. W., Tong, W.,\n")
    cat("Yauk, C. L., & Caiment, F. (2022). R-ODAF: Omics data analysis framework\n")
    cat("for regulatory application. Regulatory Toxicology and Pharmacology,\n")
    cat("131, 105143.\n")
  })
  output$download_report <- downloadHandler(
    filename = function() paste0("R-ODAF_summary_", Sys.Date(), ".txt"),
    content  = function(file) {
      sink(file)
      output$analysis_summary()
      sink()
    }
  )
}

shinyApp(ui, server)
