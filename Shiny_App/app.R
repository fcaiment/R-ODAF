##############################################################################
#  shinyapp_with_homepage.R – R-ODAF DESeq2 pipeline with homepage + tutorial
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
library(stringr)        # for str_wrap()
library(grid)           # for unit() in panel spacing

## For HGNC mapping (human)
library(AnnotationDbi)
library(org.Hs.eg.db)

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
        arrange(desc(.[, 1]))
      dmat      <- dist(vals[, 1]) %>% as.matrix()
      diag_vals <- dmat[row(dmat) == col(dmat) + 1]
      flags     <- ((pca_obj$variance[pc] * diag_vals) / max_pc) > var_thresh
      if (any(flags)) out <- c(out, rownames(vals)[which(flags)])
    }
  }
  unique(out)
}

## Map Ensembl IDs to gene symbols (human)
map_to_hgnc <- function(gene_ids) {
  if (length(gene_ids) == 0) return(character(0))
  # strip version numbers: ENSG000001.12 -> ENSG000001
  ids_clean <- sub("\\.\\d+$", "", gene_ids)
  
  res <- suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      keys    = unique(ids_clean),
      keytype = "ENSEMBL",
      columns = c("SYMBOL")
    )
  )
  sym_map <- setNames(res$SYMBOL, res$ENSEMBL)
  out     <- sym_map[ids_clean]
  names(out) <- gene_ids
  out
}

##############################  UI  ##########################################
ui <- navbarPage(
  title = "R-ODAF",
  id    = "mainNav",
  
  ## ───────────────────────────── Home tab ────────────────────────────────
  tabPanel(
    "Home",
    fluidPage(
      tags$head(
        tags$style("body { overflow-y: scroll; }")
      ),
      fluidRow(
        column(
          width = 8,
          br(),
          h2("Welcome to the R-ODAF DESeq2 pipeline"),
          p("A Shiny interface for standardised RNA-seq differential expression based on the R-ODAF framework."),
          p("Use the Analysis tab to upload count and metadata tables, run the full pipeline, and download publication-ready results."),
          br(),
          actionButton("goToAnalysis", "Start analysis", class = "btn btn-primary")
        )
      ),
      hr(),
      fluidRow(
        column(
          width = 8,
          h4("What you can do"),
          tags$ul(
            tags$li("Perform QC with coverage checks and PCA-based outlier detection."),
            tags$li("Run DESeq2 with fixed R-ODAF filters and FDR control."),
            tags$li("Explore interactive plots and export DEG tables and normalized counts.")
          ),
          br(),
          h4("Inputs"),
          tags$ul(
            tags$li("Counts: genes × samples (CSV/TSV, row names = genes, columns = samples)."),
            tags$li("Metadata: samples × attributes (CSV/TSV, row names = sample IDs)."),
            tags$li("Specify design column, control label, and case label in the Analysis tab.")
          ),
          br(),
          h4("Data privacy"),
          p("All analyses run locally in your R session. No data are uploaded to external servers."),
          br(), br()
        )
      )
    )
  ),
  
  ## ─────────────────────────── Analysis tab ──────────────────────────────
  tabPanel(
    "Analysis",
    fluidPage(
      titlePanel("R-ODAF DESeq2 Pipeline"),
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
            tabPanel(
              "Data preview",
              h4("Counts preview"),   DTOutput("counts_preview"),
              h4("Metadata preview"), DTOutput("meta_preview")
            ),
            
            tabPanel(
              "PCA & outliers",
              h4("PCA before outlier removal"),
              plotlyOutput("pca_before"),
              fluidRow(
                column(
                  width = 4,
                  selectInput(
                    "pca_before_format", "Download format",
                    choices  = c("PNG" = "png",
                                 "TIFF" = "tiff",
                                 "SVG"  = "svg",
                                 "JPEG" = "jpeg"),
                    selected = "png"
                  ),
                  downloadButton("download_pca_before", "Download PCA (before)")
                )
              ),
              verbatimTextOutput("outlier_log"),
              h4("PCA after outlier removal"),
              plotlyOutput("pca_after"),
              fluidRow(
                column(
                  width = 4,
                  selectInput(
                    "pca_after_format", "Download format",
                    choices  = c("PNG" = "png",
                                 "TIFF" = "tiff",
                                 "SVG"  = "svg",
                                 "JPEG" = "jpeg"),
                    selected = "png"
                  ),
                  downloadButton("download_pca_after", "Download PCA (after)")
                )
              )
            ),
            
            tabPanel(
              "DE results",
              verbatimTextOutput("deg_summary"),
              DTOutput("deg_table"),
              br(),
              fluidRow(
                column(
                  width = 6,
                  h5("DEGs"),
                  radioButtons(
                    "deg_format", "Format",
                    choices  = c("CSV" = "csv", "Tab-separated TXT" = "txt"),
                    selected = "csv",
                    inline   = TRUE
                  ),
                  downloadButton("download_degs_any", "Download DEGs")
                ),
                column(
                  width = 6,
                  h5("Normalized counts"),
                  radioButtons(
                    "norm_format", "Format",
                    choices  = c("CSV" = "csv", "Tab-separated TXT" = "txt"),
                    selected = "csv",
                    inline   = TRUE
                  ),
                  downloadButton("download_norm_any", "Download Norm Counts")
                )
              ),
              br(),
              h4("Heat-map: all samples"),
              plotOutput("heatmap_all"),
              br(),
              fluidRow(
                column(
                  width = 4,
                  selectInput(
                    "heat_all_format", "Download format",
                    choices  = c("PNG" = "png",
                                 "TIFF" = "tiff",
                                 "SVG"  = "svg",
                                 "JPEG" = "jpeg"),
                    selected = "png"
                  ),
                  downloadButton("download_heat_all", "Download heatmap (all)")
                )
              ),
              br(),
              h4("Heat-map: final DEGs"),
              plotOutput("heatmap_degs"),
              br(),
              fluidRow(
                column(
                  width = 4,
                  selectInput(
                    "heat_degs_format", "Download format",
                    choices  = c("PNG" = "png",
                                 "TIFF" = "tiff",
                                 "SVG"  = "svg",
                                 "JPEG" = "jpeg"),
                    selected = "png"
                  ),
                  downloadButton("download_heat_degs", "Download heatmap (DEGs)")
                )
              )
            ),
            
            tabPanel(
              "Top gene expression",
              sliderInput(
                "nTopGenes", "Total DEGs to consider",
                min = 5, max = 100, step = 5, value = 20
              ),
              uiOutput("topRangeUI"),
              radioButtons(
                "top_label_type", "Gene labels",
                choices  = c("Ensembl ID" = "ens", "Gene symbol" = "sym"),
                selected = "ens",
                inline   = TRUE
              ),
              plotlyOutput("topGenePlots"),
              br(),
              fluidRow(
                column(
                  width = 4,
                  selectInput(
                    "top_expr_format", "Download format",
                    choices  = c("PNG" = "png",
                                 "TIFF" = "tiff",
                                 "SVG"  = "svg",
                                 "JPEG" = "jpeg"),
                    selected = "png"
                  ),
                  downloadButton("download_top_expr", "Download plot")
                )
              )
            ),
            
            tabPanel(
              "Avg expression",
              sliderInput(
                "nAvgGenes", "Total DEGs to consider",
                min = 5, max = 100, step = 5, value = 20
              ),
              uiOutput("avgRangeUI"),
              radioButtons(
                "avg_label_type", "Gene labels",
                choices  = c("Ensembl ID" = "ens", "Gene symbol" = "sym"),
                selected = "ens",
                inline   = TRUE
              ),
              plotlyOutput("avgExprPlot"),
              br(),
              fluidRow(
                column(
                  width = 4,
                  selectInput(
                    "avg_expr_format", "Download format",
                    choices  = c("PNG" = "png",
                                 "TIFF" = "tiff",
                                 "SVG"  = "svg",
                                 "JPEG" = "jpeg"),
                    selected = "png"
                  ),
                  downloadButton("download_avg_expr", "Download plot")
                )
              )
            ),
            
            tabPanel(
              "Summary",
              h4("Analysis summary"),
              verbatimTextOutput("analysis_summary"),
              br(),
              downloadButton("download_report",      "Download summary (.txt)"),
              downloadButton("download_sessioninfo", "Download session info (.txt)")
            )
          )
        )
      )
    )
  ),
  
  ## ─────────────────────────── Tutorials tab ─────────────────────────────
  tabPanel(
    "Tutorials",
    fluidPage(
      br(),
      if (file.exists("www/tutorials.md")) {
        includeMarkdown("www/tutorials.md")
      } else {
        tagList(
          h3("Tutorials"),
          p("The tutorial file 'www/tutorials.md' was not found."),
          p("Place your markdown tutorial in a folder called 'www' next to this app script,"),
          p("and ensure the images are in 'www/tutorials/'.")
        )
      }
    )
  ),
  
  ## ───────────────────────── Documentation tab ───────────────────────────
  tabPanel(
    "Documentation",
    fluidPage(
      br(),
      h3("Documentation"),
      p("This section summarises the main assumptions and steps of the pipeline."),
      h4("Input format"),
      tags$ul(
        tags$li("Counts: genes × samples, CSV/TSV; row names = gene IDs."),
        tags$li("Metadata: samples × attributes, CSV/TSV; row names = sample IDs."),
        tags$li("Sample IDs must match between counts columns and metadata row names.")
      ),
      h4("Fixed filter settings"),
      tags$ul(
        tags$li(paste("Min coverage per sample:", format(MIN_COVERAGE, scientific = FALSE))),
        tags$li(paste("Min CPM:", MIN_CPM)),
        tags$li(paste("FDR cutoff:", FDR_CUTOFF)),
        tags$li(paste("Outlier variance cutoff (%):", VAR_THRESHOLD))
      )
    )
  ),
  
  ## ─────────────────────────── Citation tab ──────────────────────────────
  tabPanel(
    "Citation",
    fluidPage(
      br(),
      h3("How to cite"),
      p("If you use this app in a publication, please cite:"),
      tags$pre(
        "Lodhi, S., Verheijen, M., de Kok, T. M., Caiment, F., Jennen, D. (2025).\n",
        "R-ODAF-Shiny: An Interactive and Reproducible Framework for RNA-seq\n",
        "Transcriptomics Analysis. (manuscript ready for submission).\n\n",
        "For the original R-ODAF framework, refer to:\n\n",
        "Verheijen, M. C., Meier, M. J., Asensio, J. O., Gant, T. W., Tong, W.,\n",
        "Yauk, C. L., & Caiment, F. (2022). R-ODAF: Omics data analysis framework\n",
        "for regulatory application. Regulatory Toxicology and Pharmacology,\n",
        "131, 105143.\n",
        style = "white-space: pre-wrap;"
      )
    )
  ),
  
  ## ─────────────────────────── License tab ───────────────────────────────
  tabPanel(
    "License",
    fluidPage(
      br(),
      h3("License"),
      tags$pre(
        "MIT License\n\n",
        "Copyright (c) 2025 Saad Lodhi\n\n",
        "Permission is hereby granted, free of charge, to any person obtaining a copy\n",
        "of this software and associated documentation files (the \"Software\"), to deal\n",
        "in the Software without restriction, including without limitation the rights\n",
        "to use, copy, modify, merge, publish, distribute, sublicense, and/or sell\n",
        "copies of the Software, and to permit persons to whom the Software is\n",
        "furnished to do so, subject to the following conditions:\n\n",
        "The above copyright notice and this permission notice shall be included in all\n",
        "copies or substantial portions of the Software.\n\n",
        "THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR\n",
        "IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,\n",
        "FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE\n",
        "AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER\n",
        "LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,\n",
        "OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE\n",
        "SOFTWARE.",
        style = "white-space: pre-wrap;"
      )
    )
  ),
  
  ## ─────────────────────── Contact / Help tab ────────────────────────────
  tabPanel(
    "Contact / Help",
    fluidPage(
      br(),
      h3("Contact / Help"),
      p("For questions, suggestions, or bug reports:"),
      tags$ul(
        tags$li(
          "GitHub issues: ",
          tags$a(
            href   = "https://github.com/your-org/your-repo/issues",
            "open an issue on GitHub",
            target = "_blank"
          )
        ),
        tags$li(
          "Email: ",
          tags$a(href = "email@maastrichtuniversity.nl", "email@maastrichtuniversity.nl")
        )
      )
    )
  )
)

##############################  server  ######################################
server <- function(input, output, session) {
  WINDOW_SIZE <- 5
  
  observeEvent(input$goToAnalysis, {
    updateNavbarPage(session, "mainNav", selected = "Analysis")
  })
  
  ## dynamic range sliders
  output$topRangeUI <- renderUI({
    sliderInput(
      "topRange", "Gene range",
      min = 1, max = input$nTopGenes,
      value = c(1, min(WINDOW_SIZE, input$nTopGenes)), step = 1
    )
  })
  output$avgRangeUI <- renderUI({
    sliderInput(
      "avgRange", "Gene range",
      min = 1, max = input$nAvgGenes,
      value = c(1, min(WINDOW_SIZE, input$nAvgGenes)), step = 1
    )
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
      cnts0 <- counts_data(); meta0 <- meta_data()
      common <- intersect(colnames(cnts0), rownames(meta0))
      validate(
        need(
          length(common) >= 2,
          "No overlapping sample IDs between counts (columns) and metadata (row names)."
        )
      )
      cnts <- cnts0[, common, drop = FALSE]
      meta <- meta0[common, , drop = FALSE]
      cnts[is.na(cnts)] <- 0
      
      validate(
        need(
          all(c(input$condition1, input$condition2) %in% meta[[input$design_col]]),
          "Both condition labels must be present in the selected design column of the metadata."
        )
      )
      
      incProgress(0.10, detail = "Filtering low-coverage samples...")
      keep_samp <- colSums(cnts) > MIN_COVERAGE
      cnts <- cnts[, keep_samp, drop = FALSE]
      meta <- meta[keep_samp, , drop = FALSE]
      
      incProgress(0.15, detail = "Fitting initial DESeq2...")
      design_f <- as.formula(paste0("~", input$design_col))
      dds_init <- DESeqDataSetFromMatrix(round(cnts), meta, design_f)
      dds_init <- DESeq(dds_init, quiet = TRUE)
      
      incProgress(0.20, detail = "PCA outlier check...")
      dds_cur <- dds_init; cnts_cur <- cnts; meta_cur <- meta
      removed <- character(); conds <- c(input$condition1, input$condition2)
      repeat {
        outs <- check_outliers_vst(
          dds_cur, as.data.frame(colData(dds_cur)),
          input$design_col, conds, VAR_THRESHOLD
        )
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
      
      incProgress(0.30, detail = "Final DESeq2 & R-ODAF filters...")
      res_all <- results(
        dds_cur[keep_genes, ],
        contrast = c(input$design_col, input$condition2, input$condition1),
        alpha = FDR_CUTOFF,
        independentFiltering = FALSE,
        cooksCutoff = FALSE,
        pAdjustMethod = "fdr"
      )
      res_df <- as.data.frame(res_all)
      sig_df <- subset(res_df, padj < FDR_CUTOFF)
      
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
            if (max(v) == 0) FALSE else (max(v) / sum(v)) >= 1.4 * (length(cols))^(-0.66)
          })
          Filt[g, "spike"] <- as.integer(sum(bad_spike) <= 1)
        }
        keep_final <- rownames(Filt)[rowSums(Filt) == 3]
        degs_final <- sig_df[keep_final, , drop = FALSE]
      } else {
        degs_final <- sig_df
      }
      
      ## HGNC annotation
      all_sym <- map_to_hgnc(rownames(norm_counts))
      norm_hgnc <- data.frame(
        GeneID     = rownames(norm_counts),
        GeneSymbol = unname(all_sym),
        norm_counts,
        row.names   = NULL,
        check.names = FALSE
      )
      
      if (nrow(res_df) > 0) {
        res_sym <- map_to_hgnc(rownames(res_df))
        all_res_hgnc <- data.frame(
          GeneID     = rownames(res_df),
          GeneSymbol = unname(res_sym),
          res_df,
          row.names   = NULL,
          check.names = FALSE
        )
      } else {
        all_res_hgnc <- data.frame()
      }
      
      if (nrow(degs_final) > 0) {
        deg_sym <- map_to_hgnc(rownames(degs_final))
        degs_hgnc <- data.frame(
          GeneID     = rownames(degs_final),
          GeneSymbol = unname(deg_sym),
          degs_final,
          row.names   = NULL,
          check.names = FALSE
        )
      } else {
        degs_hgnc <- data.frame()
      }
      
      incProgress(1, detail = "Done.")
      list(
        dds_init      = dds_init,
        dds_final     = dds_cur,
        norm          = norm_counts,
        norm_hgnc     = norm_hgnc,
        all_res       = res_df,
        all_res_hgnc  = all_res_hgnc,
        degs          = degs_final,
        degs_hgnc     = degs_hgnc,
        removed       = unique(removed)
      )
    })
  })
  
  ## force evaluation when button is clicked so progress bar appears immediately
  observeEvent(input$runAnalysis, {
    pipeline()
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
  
  ## helper: mapping from Ensembl -> symbol
  symbol_map <- reactive({
    req(pipeline())
    nh <- pipeline()$norm_hgnc
    setNames(nh$GeneSymbol, nh$GeneID)
  })
  
  ## Top-gene expression (with label choice; x tick labels removed to avoid overlap)
  topExprDF <- reactive({
    req(topGenesWin(), pipeline())
    mat       <- pipeline()$norm[topGenesWin(), , drop = FALSE]
    gene_ids  <- rownames(mat)
    sym_map   <- symbol_map()
    
    labels <- gene_ids
    if (!is.null(input$top_label_type) && input$top_label_type == "sym") {
      tmp <- sym_map[gene_ids]
      labels <- ifelse(is.na(tmp) | tmp == "", gene_ids, tmp)
    }
    
    df <- as.data.frame(mat)
    df$gene <- labels
    long <- pivot_longer(df, -gene, names_to = "sample", values_to = "expression")
    meta <- as.data.frame(colData(pipeline()$dds_final)); meta$sample <- rownames(meta)
    left_join(long, meta, by = "sample")
  })
  
  output$topGenePlots <- renderPlotly({
    req(topExprDF())
    p <- ggplot(
      topExprDF(),
      aes_string(
        x = input$design_col,
        y = "expression",
        color = input$design_col,
        text  = "sample"
      )
    ) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.6) +
      facet_wrap(
        ~ gene,
        ncol    = WINDOW_SIZE,
        scales  = "free_y",
        labeller = labeller(gene = function(x) str_wrap(x, width = 10))
      ) +
      scale_x_discrete(labels = NULL) +   # remove repeated group labels to avoid overlap
      labs(x = input$design_col, y = "Normalized count") +
      theme_minimal() +
      theme(
        strip.text      = element_text(size = 8),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        panel.spacing.x = unit(1.1, "lines")
      )
    ggplotly(p, tooltip = c("text", "gene", "expression"))
  })
  
  ## Avg expression (with label choice)
  avgExprDF <- reactive({
    req(topAvgGenesWin(), pipeline())
    mat      <- pipeline()$norm[topAvgGenesWin(), , drop = FALSE]
    gene_ids <- rownames(mat)
    sym_map  <- symbol_map()
    
    labels <- gene_ids
    if (!is.null(input$avg_label_type) && input$avg_label_type == "sym") {
      tmp <- sym_map[gene_ids]
      labels <- ifelse(is.na(tmp) | tmp == "", gene_ids, tmp)
    }
    label_vec <- setNames(labels, gene_ids)
    
    meta <- as.data.frame(colData(pipeline()$dds_final))
    do.call(
      rbind,
      lapply(gene_ids, function(g) {
        mu <- tapply(mat[g, ], meta[[input$design_col]], mean)
        data.frame(
          gene            = label_vec[g],
          condition       = names(mu),
          mean_expression = as.numeric(mu),
          stringsAsFactors = FALSE
        )
      })
    )
  })
  
  output$avgExprPlot <- renderPlotly({
    req(avgExprDF())
    p <- ggplot(
      avgExprDF(),
      aes_string(x = "gene", y = "mean_expression", fill = "condition")
    ) +
      geom_col(position = position_dodge()) +
      labs(x = "Gene", y = "Mean normalized count") +
      theme_minimal() +
      theme(
        axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
        axis.title.x = element_text(margin = margin(t = 10))
      )
    ggplotly(p, tooltip = c("gene", "condition", "mean_expression"))
  })
  
  ## PCA & outliers
  output$pca_before <- renderPlotly({
    req(pipeline())
    vsd <- vst(pipeline()$dds_init, blind = TRUE) |> assay()
    md  <- as.data.frame(colData(pipeline()$dds_init))
    plotly_build(
      biplot(
        pca(vsd, metadata = md, removeVar = 0.1),
        x     = "PC1",
        y     = "PC2",
        colby = input$design_col
      ) +
        coord_fixed()
    )
  })
  output$outlier_log <- renderPrint({
    req(pipeline())
    if (length(pipeline()$removed) == 0) {
      cat("No outliers removed.")
    } else {
      cat("Removed samples:", paste(pipeline()$removed, collapse = ", "))
    }
  })
  output$pca_after <- renderPlotly({
    req(pipeline())
    vsd2 <- vst(pipeline()$dds_final, blind = TRUE) |> assay()
    md2  <- as.data.frame(colData(pipeline()$dds_final))
    plotly_build(
      biplot(
        pca(vsd2, metadata = md2, removeVar = 0.1),
        x     = "PC1",
        y     = "PC2",
        colby = input$design_col
      ) +
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
    datatable(pipeline()$degs_hgnc, options = list(scrollX = TRUE, pageLength = 5))
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
      pheatmap(
        pipeline()$norm[rownames(dg), , drop = FALSE],
        scale         = "row",
        show_rownames = FALSE
      )
    }
  })
  
  ## Downloads with format choice
  output$download_degs_any <- downloadHandler(
    filename = function() {
      ext <- if (is.null(input$deg_format) || input$deg_format == "csv") ".csv" else ".txt"
      paste0("DEGs_FDR", FDR_CUTOFF, "_", Sys.Date(), ext)
    },
    content = function(file) {
      if (is.null(input$deg_format) || input$deg_format == "csv") {
        write.csv(pipeline()$degs_hgnc, file, row.names = FALSE)
      } else {
        write.table(
          pipeline()$degs_hgnc,
          file      = file,
          sep       = "\t",
          quote     = FALSE,
          row.names = FALSE
        )
      }
    }
  )
  
  output$download_norm_any <- downloadHandler(
    filename = function() {
      ext <- if (is.null(input$norm_format) || input$norm_format == "csv") ".csv" else ".txt"
      paste0("norm_counts_", Sys.Date(), ext)
    },
    content = function(file) {
      if (is.null(input$norm_format) || input$norm_format == "csv") {
        write.csv(pipeline()$norm_hgnc, file, row.names = FALSE)
      } else {
        write.table(
          pipeline()$norm_hgnc,
          file      = file,
          sep       = "\t",
          quote     = FALSE,
          row.names = FALSE
        )
      }
    }
  )
  
  ## Summary & report
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
    cat("Lodhi, S., Verheijen, M., de Kok, T. M., Caiment, F., Jennen, D. (2025).\n")
    cat("R-ODAF-Shiny: An Interactive and Reproducible Framework for RNA-seq\n")
    cat("Transcriptomics Analysis. (manuscript ready for submission).\n\n")
    cat("For the original R-ODAF framework, refer to:\n")
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
  
  ## sessionInfo download
  output$download_sessioninfo <- downloadHandler(
    filename = function() paste0("R-ODAF_sessionInfo_", Sys.Date(), ".txt"),
    content  = function(file) {
      sink(file)
      utils::sessionInfo()
      sink()
    }
  )
  
  ## PCA downloads
  output$download_pca_before <- downloadHandler(
    filename = function() {
      paste0("PCA_before_", Sys.Date(), ".", input$pca_before_format)
    },
    content = function(file) {
      req(pipeline())
      vsd <- vst(pipeline()$dds_init, blind = TRUE) |> assay()
      md  <- as.data.frame(colData(pipeline()$dds_init))
      p <- biplot(
        pca(vsd, metadata = md, removeVar = 0.1),
        x     = "PC1",
        y     = "PC2",
        colby = input$design_col
      ) + coord_fixed()
      ggsave(file, plot = p, dpi = 300, width = 6, height = 5, units = "in")
    }
  )
  
  output$download_pca_after <- downloadHandler(
    filename = function() {
      paste0("PCA_after_", Sys.Date(), ".", input$pca_after_format)
    },
    content = function(file) {
      req(pipeline())
      vsd2 <- vst(pipeline()$dds_final, blind = TRUE) |> assay()
      md2  <- as.data.frame(colData(pipeline()$dds_final))
      p <- biplot(
        pca(vsd2, metadata = md2, removeVar = 0.1),
        x     = "PC1",
        y     = "PC2",
        colby = input$design_col
      ) + coord_fixed()
      ggsave(file, plot = p, dpi = 300, width = 6, height = 5, units = "in")
    }
  )
  
  ## Top gene expression download (same styling as on-screen)
  output$download_top_expr <- downloadHandler(
    filename = function() {
      paste0("TopGeneExpression_", Sys.Date(), ".", input$top_expr_format)
    },
    content = function(file) {
      req(topExprDF())
      p <- ggplot(
        topExprDF(),
        aes_string(
          x = input$design_col,
          y = "expression",
          color = input$design_col,
          text  = "sample"
        )
      ) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.6) +
        facet_wrap(
          ~ gene,
          ncol    = WINDOW_SIZE,
          scales  = "free_y",
          labeller = labeller(gene = function(x) str_wrap(x, width = 10))
        ) +
        scale_x_discrete(labels = NULL) +
        labs(x = input$design_col, y = "Normalized count") +
        theme_minimal() +
        theme(
          strip.text      = element_text(size = 8),
          axis.text.x     = element_blank(),
          axis.ticks.x    = element_blank(),
          panel.spacing.x = unit(1.1, "lines")
        )
      ggsave(file, plot = p, dpi = 300, width = 7, height = 5, units = "in")
    }
  )
  
  ## Avg expression download
  output$download_avg_expr <- downloadHandler(
    filename = function() {
      paste0("AvgExpression_", Sys.Date(), ".", input$avg_expr_format)
    },
    content = function(file) {
      req(avgExprDF())
      p <- ggplot(
        avgExprDF(),
        aes_string(x = "gene", y = "mean_expression", fill = "condition")
      ) +
        geom_col(position = position_dodge()) +
        labs(x = "Gene", y = "Mean normalized count") +
        theme_minimal() +
        theme(
          axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
          axis.title.x = element_text(margin = margin(t = 10))
        )
      ggsave(file, plot = p, dpi = 300, width = 7, height = 5, units = "in")
    }
  )
  
  ## Heatmap downloads
  output$download_heat_all <- downloadHandler(
    filename = function() {
      paste0("Heatmap_all_", Sys.Date(), ".", input$heat_all_format)
    },
    content = function(file) {
      req(pipeline())
      ext <- input$heat_all_format
      if (ext == "svg") {
        svg(file, width = 7, height = 5)
      } else if (ext == "tiff") {
        tiff(file, width = 7, height = 5, units = "in", res = 300, compression = "lzw")
      } else if (ext == "jpeg") {
        jpeg(file, width = 7, height = 5, units = "in", res = 300, quality = 100)
      } else {
        png(file, width = 7, height = 5, units = "in", res = 300)
      }
      corrplot(cor(pipeline()$norm), method = "color", tl.cex = 0.7)
      dev.off()
    }
  )
  
  output$download_heat_degs <- downloadHandler(
    filename = function() {
      paste0("Heatmap_DEGs_", Sys.Date(), ".", input$heat_degs_format)
    },
    content = function(file) {
      req(pipeline())
      ext <- input$heat_degs_format
      if (ext == "svg") {
        svg(file, width = 7, height = 5)
      } else if (ext == "tiff") {
        tiff(file, width = 7, height = 5, units = "in", res = 300, compression = "lzw")
      } else if (ext == "jpeg") {
        jpeg(file, width = 7, height = 5, units = "in", res = 300, quality = 100)
      } else {
        png(file, width = 7, height = 5, units = "in", res = 300)
      }
      dg <- pipeline()$degs
      if (nrow(dg) < 2) {
        plot.new(); text(0.5, 0.5, "Not enough DEGs for heat-map.")
      } else {
        pheatmap(
          pipeline()$norm[rownames(dg), , drop = FALSE],
          scale         = "row",
          show_rownames = FALSE
        )
      }
      dev.off()
    }
  )
}

shinyApp(ui, server)
