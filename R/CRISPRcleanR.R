#### interface to CRISPRcleanR-WebApp

# Whole CRISPRcleanR pipeline
ccr.AnalysisPipeline <- function(
  # Library releated parameters
  library_builtin = NULL,
  library_file = NULL,

  # Counts / FCs related parameters
  file_counts = NULL,

  # FASTQ / BAM options
  files_FASTQ_controls = NULL,
  files_FASTQ_samples = NULL,
  files_BAM_controls = NULL,
  files_BAM_samples = NULL,
  aligner = "Rsubreads",
  maxMismatches = 0,
  nTrim5 = "0",
  nTrim3 = "0",
  nBestLocations = 2,
  strand = "F",
  duplicatedSeq = "keep",
  nthreads = 1,
  indexMemory = 2000,
  fastqc_plots = FALSE,

  # Main analysis parameters
  EXPname = "",
  outdir = "./",
  ncontrols = 1,
  min_reads = 30,
  method = "ScalingByTotalReads",
  FDRth = 0.05,
  retrun_data = FALSE,

  # Correction parameters
  min.ngenes = 3,
  alpha = 0.01,
  nperm = 10000,
  p.method = "hybrid",
  min.width = 2,
  kmax = 25,
  nmin = 200,
  eta = 0.05,
  trim = 0.025,
  undo.splits = "none",
  undo.prune = 0.05,
  undo.SD = 3,

  # Run MAGeCK
  run_mageck = FALSE,
  path_to_mageck = "mageck",

  # Other options undocumented
  is_web = FALSE,
  nseed = 0xA5EED,
  verbose = -1,
  columns_map = c(
    id = "id",
    CHR = "chr",
    startp = "pos_start",
    endp = "pos_end",
    idx_start.x = "idx_start",
    idx_end = "idx_end",
    guideIdx = "idx",
    genes = "gene",
    avgFC = "avg_logFC",
    avg.logFC = "avg_logFC",
    adjFC = "avg_logFC_adj"
  )
) {
  # Create main run folder
  step_name <- "Initialize parameters"
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # Create the subfolder structures
  outdir_data <- "./data/"
  outdir_pdf <- "./pdf/"
  withr::with_dir(
    outdir,
    dir.create(outdir_data, recursive = TRUE, showWarnings = FALSE)
  )
  withr::with_dir(
    outdir,
    dir.create(outdir_pdf, recursive = TRUE, showWarnings = FALSE)
  )
  status <- 0

  # Copy the description of the output files
  withr::with_dir(
      outdir,
      file.copy(
      file.path(
        system.file("extdata", package = "CRISPRcleanR"),
        "OUTPUT_README"
      ),
      file.path(outdir_data, "README")
    )
  )

  # Create folder to store web app data
  if (is_web) {
    outdir_json <- "./json/"
    withr::with_dir(
      outdir,
      dir.create(outdir_json, recursive = TRUE, showWarnings = FALSE)
    )
    pipeline_meta_file <- "pipleline_web_v1.0.0.json"
  } else {
    outdir_json <- ""
    pipeline_meta_file <- "pipleline_main_v1.0.0.json"
  }

  # Convert path from rel. to abs. for library file
  if (!is.null(library_file)) {
    library_file <- sapply(
      library_file,
      tools::file_path_as_absolute
    )
  }

  # Convert path from rel. to abs. for count file
  if (!is.null(file_counts)) {
    file_counts <- sapply(
      file_counts,
      tools::file_path_as_absolute
    )
  }

  # Unlist list of files and set control and samples FASTQ files
  if (!is.null(files_FASTQ_controls)) {
    files_FASTQ_controls <- unlist(files_FASTQ_controls)
    if (is.null(names(files_FASTQ_controls))) {
      files_FASTQ_controls <- as.character(sapply(
        files_FASTQ_controls,
        tools::file_path_as_absolute
      ))
    } else {
      files_FASTQ_controls <- sapply(
        files_FASTQ_controls,
        tools::file_path_as_absolute
      )
    }

    files_FASTQ_samples <- unlist(files_FASTQ_samples)
    if (is.null(names(files_FASTQ_samples))) {
      files_FASTQ_samples <- as.character(sapply(
        files_FASTQ_samples,
        tools::file_path_as_absolute
      ))
    } else {
      files_FASTQ_samples <- sapply(
        files_FASTQ_samples,
        tools::file_path_as_absolute
      )
    }

    # Set control samples number based on the input files
    ncontrols <- length(files_FASTQ_controls)

    # Create BAM folder
    outdir_bam <- "./bam/"
    withr::with_dir(
      outdir,
      dir.create(outdir_bam, recursive = TRUE, showWarnings = FALSE)
    )
  } else {
    outdir_bam <- ""
  }

  # Unlist list of files and set control and samples BAM files
  if (!is.null(files_BAM_controls)) {
    files_BAM_controls <- unlist(files_BAM_controls)
    if (is.null(names(files_BAM_controls))) {
      files_BAM_controls <- as.character(sapply(
        files_BAM_controls,
        tools::file_path_as_absolute
      ))
    } else {
      files_BAM_controls <- sapply(
        files_BAM_controls,
        tools::file_path_as_absolute
      )
    }

    files_BAM_samples <- unlist(files_BAM_samples)
    if (is.null(names(files_BAM_samples))) {
      files_BAM_samples <- as.character(sapply(
        files_BAM_samples,
        tools::file_path_as_absolute
      ))
    } else {
      files_BAM_samples <- sapply(
        files_BAM_samples,
        tools::file_path_as_absolute
      )
    }

    # Set control samples number based on the input files
    ncontrols <- length(files_BAM_controls)
  }

  # Remove
  try({
    EXPname <- iconv(EXPname, to = "UTF-8", from = "ASCII//TRANSLIT")
    EXPname <- stringr::str_replace_all(EXPname, "[[:punct:]]", "_")
    EXPname <- stringr::str_replace_all(EXPname, "[^[:alnum:]]", "_")
  })

  # Print pipeline parameters
  try(ccr.PrintPipelineParams(
    # Main analysis parameters
    EXPname = EXPname,
    outdir = outdir,
    ncontrols = ncontrols,
    min_reads = min_reads,
    method = method,
    FDRth = FDRth,

    # Library releated parameters
    library_builtin = library_builtin,
    library_file = library_file,

    # Counts / FCs related parameters
    file_counts = file_counts,

    # FASTQ / BAM options
    files_FASTQ_controls = files_FASTQ_controls,
    files_FASTQ_samples = files_FASTQ_samples,
    files_BAM_controls = files_BAM_controls,
    files_BAM_samples = files_BAM_samples,
    aligner = aligner,
    maxMismatches = maxMismatches,
    nBestLocations = nBestLocations,
    nTrim5 = nTrim5,
    nTrim3 = nTrim3,
    strand = strand,
    duplicatedSeq = duplicatedSeq,
    nthreads = nthreads,
    indexMemory = indexMemory,
    fastqc_plots = fastqc_plots,

    # Correction parameters
    min.ngenes = min.ngenes,
    alpha = alpha,
    nperm = nperm,
    p.method = p.method,
    min.width = min.width,
    kmax = kmax,
    nmin = nmin,
    eta = eta,
    trim = trim,
    undo.splits = undo.splits,
    undo.prune = undo.prune,
    undo.SD = undo.SD,

    # Run MAGeCK
    run_mageck = run_mageck,
    path_to_mageck = path_to_mageck,

    # Other options udocumented
    is_web = is_web,
    nseed = nseed
  ))

  # Start pipeline with error handling
  status <- tryCatch(
    expr = {
      withr::with_dir(outdir, {
        # Load essential and signature genes for QC
        BAGEL_essential <- vector()
        BAGEL_nonEssential <- vector()
        data(BAGEL_essential, envir = environment())
        data(BAGEL_nonEssential, envir = environment())
        data(EssGenes.ribosomalProteins, envir = environment())
        data(EssGenes.DNA_REPLICATION_cons, envir = environment())
        data(EssGenes.KEGG_rna_polymerase, envir = environment())
        data(EssGenes.PROTEASOME_cons, envir = environment())
        data(EssGenes.SPLICEOSOME_cons, envir = environment())
        SIGNATURES <- list(
          Ribosomal_Proteins = EssGenes.ribosomalProteins,
          DNA_Replication = EssGenes.DNA_REPLICATION_cons,
          RNA_Polymerase = EssGenes.KEGG_rna_polymerase,
          Proteasome = EssGenes.PROTEASOME_cons,
          Spliceosome = EssGenes.SPLICEOSOME_cons,
          BAGEL_Essential = BAGEL_essential,
          BAGEL_Non_Essential = BAGEL_nonEssential
        )

        # Load pipeline steps from file
        pipleline_steps <- jsonlite::fromJSON(txt = file.path(
          system.file("extdata", package = "CRISPRcleanR"),
          pipeline_meta_file
        ))

        # Remove MAGeCK related steps if run_mageck is F
        if (!run_mageck) {
          pipleline_steps <- pipleline_steps[
            !names(pipleline_steps) %in% c(
              "mageck_uncorrected",
              "mageck_corrected",
              "imacpt_on_phenotype"
            )
          ]
        }

        # Run pipeline steps
        for (step_name in names(pipleline_steps)) {
          pipleline_steps[[step_name]][["output"]] <- ccr.ExecPipelineStep(
            # Pipeline step data
            step_name = step_name,
            step_desc = pipleline_steps[[step_name]][["desc"]],
            step_function = eval(parse(
              text = paste0(
                "CRISPRcleanR::",
                pipleline_steps[[step_name]][["FUN"]]
              )
            )),
            step_pdf = pipleline_steps[[step_name]][["pdf"]],
            step_objs = pipleline_steps[[step_name]][["objs"]],
            step_files = pipleline_steps[[step_name]][["files"]],

            # Folder strcture
            outdir = "./",
            outdir_data = outdir_data,
            outdir_pdf = outdir_pdf,
            outdir_json = outdir_json,
            outdir_bam = outdir_bam,
            filename = NULL,

            # Main analysis parameters
            EXPname = EXPname,
            CL = EXPname,
            TITLE = EXPname,
            ncontrols = ncontrols,
            min_reads = min_reads,
            method = method,
            FDRth = FDRth,
            th = FDRth,
            sigFDR = FDRth,

            # Visialization options
            display = TRUE,
            saveToFig = ifelse(
              step_name == "norm",
              TRUE,
              FALSE
            ),
            saveTO = outdir_pdf,
            plotFCprofile = TRUE,

            # Library releated parameters
            libraryAnnotation = pipleline_steps[["library"]][["output"]],
            library_builtin = library_builtin,
            library_file = library_file,

            # Counts / FCs related parameters
            file_counts = file_counts,
            counts = pipleline_steps[["counts"]][["output"]],
            Dframe = pipleline_steps[["counts"]][["output"]],
            normalised_counts = (
              pipleline_steps[["norm"]][["output"]][["norm_counts"]]
            ),
            foldchanges = if (step_name == "sort") {
              pipleline_steps[["norm"]][["output"]][["logFCs"]]
            } else {
              pipleline_steps[["correct_LFC"]][["output"]][["corrected_logFCs"]]
            },
            gwSortedFCs = pipleline_steps[["sort"]][["output"]],
            sgRNA_FCprofile = pipleline_steps[["mean_FCs_sgRNA"]][["output"]],
            correctedFCs_and_segments = (
              pipleline_steps[["correct_LFC"]][["output"]]
            ),
            FCsprofile = if (regexpr("_by_gene$", step_name) > -1) {
              pipleline_steps[["mean_FCs_gene"]][["output"]]
            } else {
              pipleline_steps[["mean_FCs_sgRNA"]][["output"]]
            },

            # FASTQ / BAM options
            files_FASTQ_controls = files_FASTQ_controls,
            files_FASTQ_samples = files_FASTQ_samples,
            files_BAM_controls = files_BAM_controls,
            files_BAM_samples = files_BAM_samples,
            aligner = aligner,
            maxMismatches = maxMismatches,
            nTrim5 = nTrim5,
            nTrim3 = nTrim3,
            nBestLocations = nBestLocations,
            strand = strand,
            duplicatedSeq = duplicatedSeq,
            nthreads = nthreads,
            indexMemory = indexMemory,
            fastqc_plots = fastqc_plots,

            # Correction parameters
            min.ngenes = min.ngenes,
            minTargetedGenes = min.ngenes,
            alpha = alpha,
            nperm = nperm,
            p.method = p.method,
            min.width = min.width,
            kmax = kmax,
            nmin = nmin,
            eta = eta,
            trim = trim,
            undo.splits = undo.splits,
            undo.prune = undo.prune,
            undo.SD = undo.SD,
            return.segments.unadj = TRUE,
            return.segments.adj = TRUE,

            # QC parameters
            positives = if (regexpr("_by_gene$", step_name) > -1) {
              BAGEL_essential
            } else {
              ccr.genes2sgRNAs(
                pipleline_steps[["library"]][["output"]], BAGEL_essential
              )
            },
            negatives = if (regexpr("_by_gene$", step_name) > -1) {
              BAGEL_nonEssential
            } else {
              ccr.genes2sgRNAs(
                pipleline_steps[["library"]][["output"]], BAGEL_nonEssential
              )
            },
            SIGNATURES = SIGNATURES,
            pIs = 6,
            nIs = 7,

            # Run MAGeCK
            run_mageck = run_mageck,
            mgckInputFile = if (step_name == "run_mageck_uncorrected") {
              file.path(outdir_data, "mageck_uncorrected_sgRNA_count.tsv")
            } else {
              file.path(outdir_data, "mageck_corrected_sgRNA_count.tsv")
            },
            normMethod = "none",
            expName = if (step_name %in% c(
              "mageck_uncorrected", "mageck_corrected"
            )) {
              if (step_name == "run_mageck_uncorrected") {
                "mageck_uncorrected"
              } else {
                "mageck_corrected"
              }
            } else {
              EXPname
            },
            outputPath = outdir_data,
            path_to_mageck = path_to_mageck,

            # Other options udocumented
            is_web = is_web,
            nseed = nseed,
            verbose = verbose,
            columns_map = columns_map
          )
        }
      })
    },

    # Error handling
    error = function(exec_error) {
      # Return the error
      # ERROR: [ERROR_NUMBER] | TYPE:[CLIENT/INTERNAL] | MSG: Msg custom
      if (step_name == "Initialize parameters") {
        exec_error_step <- 0
        exec_error_step_desc <- step_name
      } else {
        exec_error_step <- seq_along(pipleline_steps)[
          names(pipleline_steps) == step_name
        ]
        exec_error_step_desc <- pipleline_steps[[step_name]][["desc"]]
      }
      exec_error_type <- ifelse(
        exec_error_step <= 3,
        "CLIENT",
        "INTERNAL"
      )
      exec_error_message <- paste0(
        exec_error,
        " in step ",
        exec_error_step_desc
      )
      stop(paste0(
        "ERROR: ", exec_error_step, " | ",
        "TYPE: ", exec_error_type, " | ",
        "MSG: ", exec_error_message
      ))
    }
  )

  # Return data / status
  if (retrun_data) {
    return(pipleline_steps)
  } else {
    return(status)
  }
}

# Execute pipeline function step
ccr.ExecPipelineStep <- function(
  step_name,
  step_desc,
  step_function,
  step_pdf,
  step_objs,
  step_files,
  outdir_data,
  outdir_json,
  outdir_pdf,
  ...
) {
  res <- NULL

  # List all arguments
  arguments <- list(...)
  arguments[["outdir_data"]] <- outdir_data
  arguments[["outdir_json"]] <- outdir_json
  arguments[["outdir_pdf"]] <- outdir_pdf

  # Open PDF file to collect grafical function output
  if (!is.null(step_pdf)) {
    pdf(paste0(outdir_pdf, step_pdf, ".pdf"), width = 10, height = 10)
  }

  # Run with the required arguments
  res <- do.call(
    step_function,
    arguments[intersect(
      names(formals(step_function)),
      names(arguments)
    )]
  )

  # Export results as data.frame
  if (!is.null(step_files)) {
    if (length(step_files) > 1) {
      for (step_obj in step_objs) {
        ccr.dfToFile(
          res[[step_obj]],
          file.path(
            outdir_data,
            step_files[[which(step_objs == step_obj)]]
          ),
          "TSV"
        )
      }
    } else {
      if (step_name == "signatures_by_gene") {
        ccr.dfToFile(
          data.frame(
            gene_set = names(res),
            score = res,
            stringsAsFactors = FALSE
          ),
          file.path(outdir_data, step_files),
          "TSV"
        )
      } else {
        if (!step_name %in% c(
          "ROC_by_sgRNA", "ROC_by_gene",
          "PrRc_by_sgRNA", "PrRc_by_gene"
        )) {
          ccr.dfToFile(
            res,
            file.path(outdir_data, step_files),
            "TSV"
          )
        }
      }
    }
  }

  # Close PDF after function call
  if (!is.null(step_pdf)) {
    dev.off()
  }

  # Extra steps after normalization
  if (step_name == "norm") {
    # Move PDFs created by the norm step
    file.rename(
      paste0(arguments[["EXPname"]], "_fcs.pdf"),
      paste0(outdir_pdf, "fcs.pdf")
    )
    file.rename(
      paste0(arguments[["EXPname"]], "_normCounts.pdf"),
      paste0(outdir_pdf, "normCounts.pdf")
    )
  }

  # Extra steps after correction
  if (step_name == "correct_counts") {
    # Export file in MAGeckFormat
    ccr.PlainTsvFile(
      sgRNA_count_object = res,
      fprefix = "mageck_corrected",
      path = outdir_data
    )
  }

  # Export QC stats for the web
  if (arguments[["is_web"]]) {

    # Calculate norm summary stats for the web
    if (step_name == "norm") {
      ccr.dfToFile(list(
        raw = ccr.countToDist(
          arguments[["counts"]],
          arguments[["ncontrols"]]
        ),
        norm = ccr.countToDist(
          res[["norm_counts"]],
          arguments[["ncontrols"]]
        )),
        file.path(outdir_json, "normCounts"), "JSON")
      ccr.dfToFile(list(
        norm = ccr.countToDist(
          res[["logFCs"]],
          0
        )),
        file.path(outdir_json, "fcs"), "JSON")
    }

    # Calculate corrected summary stats for the web
    if (step_name == "correct_LFC") {
      ccr.ExportCorrectedFCsToWeb(
        correctedFCs = res,
        columns_map = arguments[["columns_map"]],
        outdir_json = outdir_json
      )
    }

    # Calculate corrected summary stats for the web
    if (step_name %in% c("ROC_by_sgRNA", "ROC_by_gene")) {
      metrics_df <- data.frame(
        AUC = ifelse(
          is.null(res[["AUC"]]),
          NA,
          res[["AUC"]]
        ),
        THR = ifelse(
          is.null(res[["sigthreshold"]]),
          NA,
          res[["sigthreshold"]]
        ),
        Recall = ifelse(
          is.null(res[["Recall"]]),
          NA,
          res[["Recall"]]
        )
      )
      colnames(metrics_df) <- c(
        "AUC",
        paste0(round(arguments[["FDRth"]] * 100, 2), "% FDR logFC threshold"),
        paste0("Recall at ", round(arguments[["FDRth"]] * 100, 2), "% FDR")
      )
      ccr.dfToFile(
        list(
          metrics = metrics_df,
          curve = data.frame(res[["curve"]])
        ),
        file.path(outdir_json, step_files),
        "JSON"
      )
    }

    if (step_name %in% c("PrRc_by_sgRNA", "PrRc_by_gene")) {
      metrics_df <- data.frame(
        AUC = ifelse(is.null(res[["AUC"]]), NA, res[["AUC"]]),
        FDR = ifelse(
          is.null(res[["sigthreshold"]]),
          NA,
          res[["sigthreshold"]]
        ),
        precision = ifelse(
          is.null(res[["Recall"]]),
          NA,
          head(res[["curve"]][
            abs(res[["curve"]][, "recall"] - res[["Recall"]]) == min(
              abs(res[["curve"]][, "recall"] - res[["Recall"]])
            ),
            "precision"],
            1
          )
        ),
        recall = ifelse(is.null(res[["Recall"]]), NA, res[["Recall"]])
      )
      colnames(metrics_df) <- c(
        "AUC",
        paste0(round(arguments[["FDRth"]] * 100, 2), "% FDR logFC threshold"),
        paste0("Precision at ", round(arguments[["FDRth"]] * 100, 2), "% FDR"),
        paste0("Recall at ", round(arguments[["FDRth"]] * 100, 2), "% FDR")
      )
      ccr.dfToFile(
        list(
          metrics = metrics_df,
          curve = data.frame(res[["curve"]])
        ),
        file.path(outdir_json, step_files),
        "JSON"
      )
    }

    if (step_name == "signatures_by_gene") {
      geneFCs <- sort(arguments[["FCsprofile"]], decreasing = FALSE)
      Recall_scores_df_for_web <- list(
        metrics = data.frame(threshold = arguments[["th"]]),
        curve = data.frame(
          gene = names(geneFCs),
          rank = seq_along(geneFCs),
          logFC = geneFCs,
          stringsAsFactors = FALSE
        ),
        gene_set_array = list()
      )
      for (gene_set in names(res)) {
        Recall_scores_df_for_web[["gene_set_array"]][[gene_set]] <- list(
          score = res[[gene_set]],
          genes = Recall_scores_df_for_web[["curve"]][
            Recall_scores_df_for_web[["curve"]][,
              "gene"
            ] %in% arguments[["SIGNATURES"]][[gene_set]],
          ]
        )
        rownames(
          Recall_scores_df_for_web[["gene_set_array"]][[gene_set]][["genes"]]
        ) <- NULL
      }
      rownames(Recall_scores_df_for_web[["curve"]]) <- NULL
      ccr.dfToFile(
        Recall_scores_df_for_web,
        file.path(outdir_json, step_files),
        "JSON"
      )
    }
  }

  # Return results
  return(res)
}

ccr.RemoveExtraFiles <- function(
  is_web = FALSE,
  file_counts = NULL,
  files_FASTQ_controls = NULL,
  files_FASTQ_samples = NULL,
  files_BAM_controls = NULL,
  files_BAM_samples = NULL,
  outdir_data = NULL
) {
    # Remove unecessary files
  for (file_current in c(
    file_counts,
    files_FASTQ_controls,
    files_FASTQ_samples,
    files_BAM_controls,
    files_BAM_samples,
    list.files(
      outdir_data,
      pattern = "mageck_corrected"
    )
  )) {
    file_to_remove <- file.path(
      outdir_data,
      basename(file_current)
    )
    if (
      is_web &
      file.exists(file_to_remove) &
      !basename(file_to_remove) %in% c(
        "raw_counts.tsv",
        "count_norm.tsv",
        "counts_corrected.tsv",
        "mageck_corrected_sgRNA_count.tsv",
        "mageck_corrected.gene_summary.txt",
        "mageck_corrected.sgrna_summary.txt"
      )
    ) {
      file.remove(file_to_remove)
    }
  }
}

# Export the LFCs and segments
ccr.ExportCorrectedFCsToWeb <- function(
  correctedFCs,
  columns_map,
  outdir_json
) {
  data_web <- list()
  for (data_type in names(correctedFCs)[
    names(correctedFCs) != "SORTED_sgRNAs"]
  ) {
    if (data_type == "corrected_logFCs") {
      data_type <- "sgRNA_array"

      #format LFC data
      data_web[[data_type]] <- correctedFCs[["corrected_logFCs"]][
        correctedFCs[["SORTED_sgRNAs"]],
      ]
      data_web[[data_type]][, "adjFC"] <- ifelse(
        abs(
          data_web[[data_type]][, "correctedFC"] -
          data_web[[data_type]][, "avgFC"]
        ) < 0.01,
        NA,
        data_web[[data_type]][, "correctedFC"])
      data_web[[data_type]][, "id"] <- rownames(data_web)
      rownames(data_web[[data_type]]) <- NULL
    } else {

      #format segment data
      data_web[[data_type]] <- correctedFCs[[data_type]]
      data_web[[data_type]][, "idx_start"] <- vapply(
        as.character(data_web[[data_type]][, "guideIdx"]),
        function(s) as.integer(strsplit(s, ",")[[1]][[1]]),
        integer(length = 1)
      )
      data_web[[data_type]][, "idx_end"] <- vapply(
        as.character(data_web[[data_type]][, "guideIdx"]),
        function(s) as.integer(strsplit(s, ",")[[1]][[2]]),
        integer(length = 1)
      )
      data_web[[data_type]] <- merge(
        data_web[[data_type]],
        data.frame(aggregate(
          idx_start ~ CHR,
          data_web[[data_type]],
          min
        )),
        by = "CHR"
      )
      data_web[[data_type]][, "idx_start.x"] <- (
        data_web[[data_type]][, "idx_start.x"] -
        data_web[[data_type]][, "idx_start.y"] + 1
      )
      data_web[[data_type]][, "idx_end"] <- (
        data_web[[data_type]][, "idx_end"] -
        data_web[[data_type]][, "idx_start.y"] + 1
      )
    }

    #selct columns to export
    data_web[[data_type]] <- data_web[[data_type]][,
      names(columns_map)[names(columns_map) %in%
      colnames(data_web[[data_type]])]
    ]
    colnames(data_web[[data_type]]) <- columns_map[
      colnames(data_web[[data_type]])
    ]

    #split data by CHR
    data_web[[data_type]] <- split(
      data_web[[data_type]],
      data_web[[data_type]][, "chr"]
    )
  }

  #export JSON files
  for (chrN in names(data_web[[1]])) {
    ccr.dfToFile(
      lapply(data_web, function(l) l[[chrN]]),
      file.path(outdir_json, chrN),
      "JSON"
    )
  }
}


# Get library annotation data
ccr.getLibrary <- function(
  library_builtin,
  library_file,
  verbose = FALSE
) {

  libraryAnnotation <- NULL
  if (verbose > -1) pb <- txtProgressBar(min = 0, max = 1, style = 3)
  if (
    all(is.null(library_builtin), is.null(library_file)) |
    all(!is.null(library_builtin), !is.null(library_file))
  ) {
    if (all(is.null(library_builtin), is.null(library_file))) {
      stop("No library is defined")
    }
    if (all(!is.null(library_builtin), !is.null(library_file))) {
      stop("Only one of the library_builtin - library_file should be defined")
    }
  } else {
    if (!is.null(library_file)) {
      if (class(library_file) == "character") {
        if (file.exists(library_file)) {
          field_sep <- tools::file_ext(library_file)
          if (field_sep %in% c("txt", "tsv")) {
            libraryAnnotation <- read.table(
              library_file,
              sep = "\t",
              header = TRUE,
              stringsAsFactors = FALSE
            )
          } else {
            if (field_sep == "csv") {
              libraryAnnotation <- read.table(
                library_file,
                sep = ",",
                header = TRUE,
                stringsAsFactors = FALSE
              )
            } else {
              stop("Unknown file format loading library")
            }
          }
        }
      }
      if (class(library_file) == "data.frame") {
        libraryAnnotation <- library_file
      }
    } else {
      if (
        file.exists(system.file(
          "data",
          paste0(library_builtin, ".RData"),
          package = "CRISPRcleanR")
        )) {
        load(system.file(
          "data",
          paste0(library_builtin, ".RData"),
          package = "CRISPRcleanR"),
          environment()
        )
        libraryAnnotation <- get(library_builtin, environment())
      }
    }
    if (is.null(libraryAnnotation)) stop("Missing annotation library")
  }
  if (verbose > -1) setTxtProgressBar(pb, 1)
  stopifnot(class(libraryAnnotation) == "data.frame")
  stopifnot(all(table(libraryAnnotation[, "CODE"]) == 1))
  stopifnot(all(
    c("CODE", "GENES", "CHRM", "STARTpos", "ENDpos") %in%
    colnames(libraryAnnotation)
  ))

  # Set row names to CODE (sgRNA IDs)
  rownames(libraryAnnotation) <- libraryAnnotation[, "CODE"]

  if (verbose > -1) close(pb)
  return(libraryAnnotation)
}


# Load count data
ccr.getCounts <- function(
  file_counts,
  files_FASTQ_controls,
  files_FASTQ_samples,
  files_BAM_controls,
  files_BAM_samples,
  libraryAnnotation,
  maxMismatches,
  nTrim5,
  nTrim3,
  nthreads,
  nBestLocations,
  duplicatedSeq,
  indexMemory,
  strand,
  EXPname,
  outdir_data,
  outdir_bam,
  aligner,
  fastqc_plots,
  verbose
) {
  counts <- NULL
  if (verbose > -1) pb <- txtProgressBar(min = 0, max = 1, style = 3)
  if (sum(
    !is.null(file_counts),
    !all(is.null(files_FASTQ_controls), is.null(files_FASTQ_samples)),
    !all(is.null(files_BAM_controls), is.null(files_BAM_samples))
  ) != 1) {
    stop("Only one type of data should be defined")
  } else {
    # Align FASTQ files and get counts
    if (!all(is.null(files_FASTQ_controls), is.null(files_FASTQ_samples))) {
      if (all(!is.null(files_FASTQ_controls), !is.null(files_FASTQ_samples))) {

        # Get counts from FASTQ files
        counts <- ccr.FASTQ2counts(
          c(files_FASTQ_controls, files_FASTQ_samples),
          libraryAnnotation,
          maxMismatches = maxMismatches,
          nTrim5 = nTrim5,
          nTrim3 = nTrim3,
          nthreads = nthreads,
          nBestLocations = nBestLocations,
          strand = strand,
          duplicatedSeq = duplicatedSeq,
          indexMemory = indexMemory,
          EXPname = EXPname,
          outdir = outdir_bam,
          aligner = aligner,
          fastqc_plots = fastqc_plots,
          export_counts = TRUE,
          overwrite = FALSE
        )
        file.copy(
          file.path(outdir_bam, paste0(EXPname, ".counts.tsv")),
          file.path(outdir_data, "raw_counts.tsv")
        )
        file.copy(
          file.path(outdir_bam, paste0(EXPname, ".alignments_stats.tsv")),
          file.path(outdir_data, "alignments_stats.tsv")
        )
      } else {
        stop("FASTQ files should be suppled for both controls and samples")
      }
    }

    # Load counts from bam files
    if (!all(is.null(files_BAM_controls), is.null(files_BAM_samples))) {
      if (all(!is.null(files_BAM_controls), !is.null(files_BAM_samples))) {

        # Get counts from BAM files
        counts <- ccr.BAM2counts(
          c(files_BAM_controls, files_BAM_samples),
          libraryAnnotation,
          maxMismatches = maxMismatches,
          strand = strand,
          EXPname = EXPname,
          outdir = outdir_data,
          export_counts = TRUE,
          overwrite = TRUE
        )
        file.copy(
          file.path(outdir_bam, paste0(EXPname, ".counts.tsv")),
          file.path(outdir_data, "raw_counts.tsv")
        )
      } else {
        stop("FASTQ files should be suppled for both controls and samples")
      }
    }

    # Load counts from text matrix
    if (!is.null(file_counts)) {
      if (class(file_counts) == "character") {
        if (file.exists(file_counts)) {
          if (tools::file_ext(file_counts) %in% c("gz", "zip")) {
            field_sep <- tools::file_ext(tools::file_path_sans_ext(file_counts))
          } else {
            field_sep <- tools::file_ext(file_counts)
          }
          if (field_sep %in% c("txt", "tsv")) {
            counts <- read.table(
              file_counts,
              sep = "\t",
              header = TRUE,
              stringsAsFactors = FALSE
            )
          } else {
            if (field_sep == "csv") {
              counts <- read.table(
                file_counts,
                sep = ",",
                header = TRUE,
                stringsAsFactors = FALSE
              )
            } else {
              stop("Unknown file format loading count data")
            }
          }
        } else {
          stop("Count file doesn't exist")
        }
      }
      if (class(file_counts) == "data.frame") {
        counts <- file_counts
      }
    }

    # Rise error if no data is supplied
    if (is.null(counts)) stop("At least one type of data should be defined")
  }
  if (verbose > -1) setTxtProgressBar(pb, 1)
  if (verbose > -1) close(pb)
  stopifnot(all(c("sgRNA", "gene") %in% colnames(counts)))
  stopifnot(all(table(counts[, "sgRNA"]) == 1))

  # Set counts row names
  rownames(counts) <- counts[, "sgRNA"]

  return(counts)
}

#Check consistency between library and count files
ccr.checkCounts <- function(
  counts,
  libraryAnnotation,
  ncontrols = 1,
  min_reads = 30
) {
  # Check in the counts has enough overlap with the library
  if (any(!rownames(libraryAnnotation) %in% counts[, "sgRNA"])) {
    if ((
      sum(rownames(libraryAnnotation) %in% counts[, "sgRNA"]) /
      nrow(libraryAnnotation)
    ) < 0.8) {
      stop("The Count file include less than 80% of the library sgRNAs")
    } else {
      warning(paste0(
        "Count file contains ",
        sum(!rownames(libraryAnnotation) %in% counts[, "sgRNA"]),
        " sgRNAs without data"
      ))
    }
  }

  # Check in the library has enough overlap with the counts
  if (any(!counts[, "sgRNA"] %in% rownames(libraryAnnotation))) {
    if ((
      sum(!counts[, "sgRNA"] %in% rownames(libraryAnnotation)) /
      nrow(libraryAnnotation)
    ) > 0.2) {
      stop(paste0(
        "More than 20% of the sgRNAs",
        " in the count file are missing in the library"
      ))
    } else {
      warning(paste0(
        "Count miss ",
        sum(!counts[, "sgRNA"] %in% rownames(libraryAnnotation)),
        " sgRNAs from annotation library"
      ))
    }
  }

  # Check in the overlapping probes have enough reads
  if ((sum(rowMeans(counts[
    counts[, "sgRNA"] %in% libraryAnnotation[, "CODE"],
    seq(from = 3, to = 2 + ncontrols),
    drop = FALSE],
    na.rm = TRUE
  ) >= min_reads) / nrow(libraryAnnotation)) < 0.8) {
    stop(paste0(
      "Less than 80% of the library sgRNAs",
      " have more than ", min_reads, " reads"
    ))
  }

  # Return TRUE if there aren't errors
  return(TRUE)
}


# Print the run parameters
ccr.PrintPipelineParams <- function(
  ...
) {

  # Get the list of arguments
  arguments <- list(...)

  # Print Parameter list
  print("############################################")
  print("Input parameters")
  print("############################################")
  for (par_name in names(arguments)) {
    if (!par_name %in% c(
      "is_web", "columns_map",
      "GDSC.geneLevCNA", "CCLE.gisticCNA", "RNAseq.fpkms"
    )) {
      if (!is.null(arguments[[par_name]])) {
        if (par_name %in% c(
          "file_counts", "outdir",
          "files_FASTQ_controls", "files_FASTQ_samples",
          "files_BAM_controls", "files_BAM_controls"
        ) & arguments[["is_web"]] == TRUE &
          is.character(arguments[[par_name]])
        ) {
          arguments[[par_name]] <- basename(arguments[[par_name]])
        }
        print(paste0(
          "Parameter ", par_name, ": ",
          arguments[[par_name]],
          " (class: ", class(arguments[[par_name]]),
          "-", typeof(arguments[[par_name]]), ")"
        ))
      }
    }
  }
}

#### END interface to CRISPRcleanR-WebApp


#### sgRNAcount generation from low-level sequencing Files

#### Create index file for the library
ccr.CreateLibraryIndex <- function(
  libraryAnnotation,
  duplicatedSeq = "keep",
  EXPname = "",
  indexMemory = 2000,
  overwrite = FALSE
) {

  # Deal with duplicated sequences
  if (duplicatedSeq %in% c("exclude", "keep")) {
    if (duplicatedSeq == "exclude") {
      libraryToUse <- libraryAnnotation[
        !libraryAnnotation[["seq"]] %in%
        names(table(libraryAnnotation[["seq"]]))[
          table(libraryAnnotation[["seq"]]) > 1
        ],
        c("CODE", "seq")
      ]
    } else {
      libraryToUse <- libraryAnnotation[
        !duplicated(libraryAnnotation[["seq"]]),
        c("CODE", "seq")
      ]
    }
  } else {
    stop("duplicatedSeq param should be exclude or keep.")
  }

  sgRNAs <- Biostrings::DNAStringSet(libraryToUse[["seq"]])
  names(sgRNAs) <- libraryToUse[["CODE"]]
  libraryName <- paste0("Library_", EXPname)
  libraryFile <- paste0(libraryName, ".fa")
  if (file.exists(libraryFile) & !overwrite) {
    warning(paste0("Library ", libraryFile, " already exist"))
  } else {
    Biostrings::writeXStringSet(sgRNAs, file = libraryFile)
    write.table(
        libraryAnnotation[, c("CODE", "seq", "GENES")],
        file = paste0(libraryName, ".txt"),
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )
  }
  if (file.exists(paste0(libraryName, ".reads")) & !overwrite) {
    warning(paste0("Library index ", libraryName, " already exist"))
  } else {
    Rsubread::buildindex(
      libraryName,
      libraryFile,
      TH_subread = nrow(libraryToUse) + 1,
      gappedIndex = FALSE,
      memory = indexMemory,
      indexSplit = TRUE
    )
  }
  return(libraryName)
}

# Extract counts data from BAM files
ccr.BAM2counts <- function(
  BAMfileList,
  libraryAnnotation,
  maxMismatches = 0,
  strand = "F",
  EXPname = "",
  outdir = "./",
  export_counts = TRUE,
  overwrite = TRUE
) {
  counts <- data.frame(sgRNA = character(0), stringsAsFactors = FALSE)

  # Add seqlen to library for mismatch check
  libraryAnnotation[["seqlen"]] <- nchar(libraryAnnotation[["seq"]])

  # Create BAM sample names if missing
  if (is.null(names(BAMfileList))) {
    names(BAMfileList) <- tools::file_path_sans_ext(basename(BAMfileList))
  }

  # Stop if duplicates sample names
  if (any(table(names(BAMfileList)) > 1)) {
    stop("Duplicated BAM sample names")
  }

  # Define strand to use
  if (strand %in% c("F", "R", "*")) {
    if (strand == "F") strand <- "-"
    if (strand == "R") strand <- "+"
    if (strand == "*") strand <- c("-", "+")
  } else {
    stop("Strand should be one of F, R, *")
  }

  # Loop trough BAM samples files
  minWidth <- (min(libraryAnnotation[["seqlen"]]) - maxMismatches)
  for (BAMsample in BAMfileList) {
    if (file.exists(BAMsample)) {
      BAMsampleName <- names(BAMfileList)[BAMfileList == BAMsample]

      # Read alignment and get mapping quality / cigar
      myAlignedReads <- GenomicAlignments::readGAlignments(BAMsample)

      # Read counts by sgRNA ID
      myAlignedReads <- as.data.frame(table(
        seqnames(myAlignedReads[
          strand(myAlignedReads) %in% strand &
          width(myAlignedReads) >= minWidth
        ])
      ))
      colnames(myAlignedReads) <- c("sgRNA", BAMsampleName)
      myAlignedReads <- merge(
        libraryAnnotation[, c("CODE", "GENES")],
        myAlignedReads,
        by.x = "CODE",
        by.y = "sgRNA",
        all.y = TRUE
      )
      colnames(myAlignedReads) <- c("sgRNA", "gene", BAMsampleName)

      # Add sample to count matrix
      counts <- merge(counts, myAlignedReads, all = TRUE)
    } else {
      dev.off()
      stop(paste0("Missing BAM file: ", BAMsample))
    }
  }

  # Fix missing values
  counts[is.na(counts)] <- 0

  # Set counts row names
  rownames(counts) <- counts[, "sgRNA"]

  # Export count table
  if (export_counts) {
    if (
      !file.exists(file.path(outdir, paste0(EXPname, ".counts.tsv"))) |
      overwrite
    ) {
      write.table(
        counts,
        file = file.path(outdir, paste0(EXPname, ".counts.tsv")),
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE
      )
    }
  }

  return(counts)
}

# Extract counts data from FASTQ files with MAGeCK
ccr.MAGeCK2counts <- function(
  FASTQfileList,
  libraryAnnotation,
  maxMismatches = 0,
  nTrim5 = 0,
  strand = "F",
  fastqc_plots = TRUE,
  EXPname = "",
  outdir = "./",
  aligner = "Rsubreads",
  path_to_mageck = "mageck",
  overwrite = FALSE
) {
    textbunch <- paste0(
      path_to_mageck, " count ",
      "--list-seq ", paste0("Library_", EXPname), ".txt ",
      "--fastq ", paste0(FASTQfileList, collapse = " "), " ",
      "--norm-method none",
      "--sample-label ", paste0(names(FASTQfileList), collapse = ","), " ",
      "-n ", file.path(outdir, EXPname), " ",
      "--trim-5 ", nTrim5,
      if (maxMismatches > 0) "--count-n ",
      if (strand == "R") "--reverse-complement ",
      if (fastqc_plots) "--pdf-report "
    )
    output_text <- system(textbunch, intern = TRUE, wait = TRUE)
    print(output_text)

    counts <- read.delim(
      file.path(outdir, paste0(EXPname, ".counts.tsv"))
    )

    return(counts)
}

# Extract counts data from FASTQ files
ccr.FASTQ2counts <- function(
  FASTQfileList,
  libraryAnnotation,
  maxMismatches = 0,
  nTrim5 = "0",
  nTrim3 = "0",
  nthreads = 1,
  nBestLocations = 2,
  strand = "F",
  indexMemory = 2000,
  duplicatedSeq = "keep",
  EXPname = "",
  outdir = "./",
  aligner = "Rsubreads",
  fastqc_plots = TRUE,
  export_counts = TRUE,
  overwrite = FALSE
) {
  # Create FASTQ sample names if missing
  if (is.null(names(FASTQfileList))) {
    names(FASTQfileList) <- ifelse(
      tools::file_ext(basename(FASTQfileList)) == "gz",
      tools::file_path_sans_ext(
        tools::file_path_sans_ext(basename(FASTQfileList))
      ),
      tools::file_path_sans_ext(basename(FASTQfileList))
    )
  }

  # Stop if duplicates sample names
  if (any(table(names(FASTQfileList)) > 1)) {
    stop("Duplicated FASTQ sample names")
  }

  # Create library index
  libraryIndex <- ccr.CreateLibraryIndex(
    libraryAnnotation = libraryAnnotation,
    duplicatedSeq = duplicatedSeq,
    EXPname = EXPname,
    indexMemory = indexMemory,
    overwrite = overwrite
  )

  # Get count using Rsubreads alingments
  if (aligner == "Rsubreads") {
    alignment_stats <- data.frame(metrics = character())
    BAMfileList <- vector()

    # Loop though fastq files for alignment and QC
    for (FASTQsample in FASTQfileList) {
      if (file.exists(FASTQsample)) {
        FASTQsampleName <- names(FASTQfileList)[
          FASTQfileList == FASTQsample
        ]

        # Run QC
        if (file.exists(file.path(
            outdir,
            paste0(FASTQsampleName, ".html")
          )) & !overwrite
        ) {
          warning(paste(
            "FASTQ QC report for sample",
            FASTQsampleName,
            "already exists"
          ))
        } else {
          if (fastqc_plots) {
            qcRes <- Rqc::rqc(
              path = dirname(FASTQsample),
              sample = FALSE,
              pattern = paste0("^", basename(FASTQsample), "$"),
              n = 500000,
              outdir = outdir,
              file = FASTQsampleName,
              workers = 1,
              openBrowser = FALSE
            )
          }
        }

        if (file.exists(file.path(
          outdir,
          paste0(FASTQsampleName, ".bam")
          )) & !overwrite
        ) {
          warning(paste(
            "BAM file for sample",
            FASTQsampleName,
            "already exists"
          ))
        } else {
          myMapped <- Rsubread::align(
            # index for reference sequences
            index = libraryIndex,

            # input reads and output
            readfile1 = FASTQsample,
            type = "DNA",
            input_format = "gzFASTQ",
            output_format = "BAM",
            output_file = file.path(outdir, paste0(FASTQsampleName, ".bam")),

            # offset value added to Phred quality scores of read bases
            phredOffset = 33,

            # thresholds for mapping
            nsubreads = 10,
            TH1 = 1,
            maxMismatches = maxMismatches,

            # distance and orientation of paired end reads
            minFragLength = 0,
            maxFragLength = 100,
            PE_orientation = "rr",

            # unique mapping and multi-mapping
            unique = FALSE,
            nBestLocations = nBestLocations,

            # indel detection
            indels = as.integer(maxMismatches > 0),
            complexIndels = FALSE,

            # read trimming
            nTrim5 = as.integer(nTrim5),
            nTrim3 = as.integer(nTrim3),

            # read order
            keepReadOrder = FALSE,
            sortReadsByCoordinates = TRUE,

            # number of CPU threads
            nthreads = nthreads,

            # dynamic programming
            DP_GapOpenPenalty = -1,
            DP_GapExtPenalty = 0,
            DP_MismatchPenalty = 0,
            DP_MatchScore = 2,

            # detect structural variants
            detectSV = FALSE,

            # gene annotation
            useAnnotation = FALSE,
            chrAliases = NULL
          )

          # Export alignemnts stats
          myMapped[, "metrics"] <- rownames(myMapped)
          write.table(
            myMapped,
            file = file.path(
              outdir,
              paste0(FASTQsampleName, "_stats.txt")
            ),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE
          )
        }

        # Collect alignemnts stats
        if (file.exists(file.path(
          outdir,
          paste0(FASTQsampleName, "_stats.txt")
        ))) {
          alignment_stats <- merge(
            alignment_stats,
            read.delim(
              file.path(
                outdir,
                paste0(FASTQsampleName, "_stats.txt")
              ),
              header = TRUE,
              sep = "\t",
              stringsAsFactors = FALSE
            ),
            all = TRUE
          )
        }

        # Add BAM to BAM list
        BAMfileList <- c(
          BAMfileList,
          file.path(outdir, paste0(FASTQsampleName, ".bam"))
        )
        names(BAMfileList)[length(BAMfileList)] <- FASTQsampleName
      } else {
        stop(paste0("Missing FASTQ file: ", FASTQsample))
      }
    }

    # Export alignments stats
    write.table(
      alignment_stats,
      file = file.path(
        outdir,
        paste0(EXPname, ".alignments_stats.tsv")
      ),
      sep = "\t",
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE
    )

    # Read counts from BAM files
    counts <- ccr.BAM2counts(
      BAMfileList = BAMfileList,
      libraryAnnotation = libraryAnnotation,
      strand = strand,
      EXPname = EXPname,
      outdir = outdir,
      export_counts = export_counts,
      overwrite = overwrite
    )
  }

  # Get count using MAGeCK
  if (aligner == "mageck")  {
    counts <- ccr.MAGeCK2counts(
      FASTQfileList,
      libraryAnnotation,
      maxMismatches = maxMismatches,
      nTrim5 = nTrim5,
      strand = strand,
      fastqc_plots = fastqc_plots,
      EXPname = "",
      outdir = "./",
      path_to_mageck = "mageck",
      overwrite = FALSE
    )
  }
  return(counts)
}

#### END sgRNAcount generation from low-level sequencing Files

#### Analysis
ccr.NormfoldChanges <- function(
  filename,
  Dframe = NULL,
  display = TRUE,
  saveToFig = FALSE,
  outdir = "./",
  min_reads = 30,
  EXPname = "",
  libraryAnnotation,
  ncontrols = 1,
  method = "ScalingByTotalReads"
) {
  if (length(Dframe) == 0) {
    counts <- read.table(
      filename,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE
    )
  }else{
    counts <- Dframe
  }

  counts <- counts[is.element(counts[["sgRNA"]], rownames(libraryAnnotation)), ]

  if (saveToFig) {
    display <- TRUE
    pdf(
      paste0(outdir, EXPname, "_normCounts.pdf"),
      width = 10,
      height = 10
    )
  }

  withr::with_par(
    list(mfrow = c(1, 2)), {
      if (display) {
        ccr.boxplot(
          counts[, 3:ncol(counts)],
          main = paste(EXPname, "Raw sgRNA counts"),
          names = c(
            paste("CTRL", seq_len(ncontrols)),
            paste("library r", seq_len(ncol(counts) - 2 - ncontrols))
          )
        )
      }

      numd <- counts[
        is.element(
          counts[["sgRNA"]],
          rownames(libraryAnnotation)
        ),
        seq(from = 3, to = ncol(counts))
      ]

      IDX <- which(rowMeans(
        numd[, seq_len(ncontrols), drop = FALSE],
      ) >= min_reads)
      numd <- numd[IDX, ]
      counts <- counts[IDX, ]

      if (method == "MedRatios") {
        numd <- numd + 0.5
        pseudo_ref_sample <- apply(
          numd,
          MARGIN = 1,
          function(x) {
            prod(x) ^ (1 / length(x))
          }
        )
        pseudo_ref_mat <- matrix(
          rep(pseudo_ref_sample, ncol(numd)),
          length(pseudo_ref_sample), ncol(numd)
        )
        numd <- numd / pseudo_ref_mat
        normFact <- t(matrix(
          rep(colSums(numd), nrow(numd)),
          ncol(counts) - 2,
          nrow(numd)
        ))
      } else {
        if (method == "ScalingByTotalReads") {
            normFact <- t(matrix(
              rep(colSums(numd), nrow(numd)),
              ncol(counts) - 2,
              nrow(numd)
            ))
        } else {
          if (!is.element(
            method,
            libraryAnnotation[["GENES"]]
          )) {
            print("Invalid normalisation method")
            return()
          } else {
            gidx <- rownames(libraryAnnotation)[
              which(libraryAnnotation[["GENES"]] == method)
            ]
            gidx <- intersect(gidx, counts[["sgRNA"]])
            gidx <- match(gidx, counts[["sgRNA"]])
          }
          normFact <- t(matrix(
            rep(colSums(numd[gidx, ]), nrow(numd)),
            ncol(counts) - 2,
            nrow(numd)
          ))
        }
      }
      numd <- numd / normFact * 10000000
      normed <- cbind(counts[, seq_len(2)], numd)

      if (display) {
        ccr.boxplot(
          numd,
          main = paste(EXPname, "normalised sgRNA counts"),
          names = c(
            paste("CTRL", seq_len(ncontrols)),
            paste("library r", seq_len(ncol(counts) - 2 - ncontrols))
          )
        )
      }
    },
    no.readonly = TRUE
  )

  if (saveToFig) {
    dev.off()
  }

  nsamples <- ncol(counts) - (2 + ncontrols)

  for (i in seq_len(nsamples)) {
    c_foldchanges <- log2(
      (normed[, 2 + ncontrols + i] + 0.5) /
      (rowMeans(normed[,
        seq(from = 3, to = 2 + ncontrols),
        drop = FALSE
      ]) + 0.5)
    )

    c_foldchanges <- matrix(
      c_foldchanges,
      length(c_foldchanges),
      1,
      dimnames = list(normed[["sgRNA"]], colnames(normed)[2 + ncontrols + i])
    )

    if (i == 1) {
      foldchanges <- c_foldchanges
    } else {
      foldchanges <- cbind(foldchanges, c_foldchanges)
    }
  }

  if (saveToFig) {
    pdf(
      paste0(outdir, EXPname, "_fcs.pdf"),
      width = 8,
      height = 10
    )
  }

  if (display) {
    withr::with_par(
      list(mfrow = c(1, 1)), {
        ccr.boxplot(
          foldchanges,
          main = paste(EXPname, " sgRNA log fold changes"),
          names = paste("library r", seq_len(ncol(counts) - 2 - ncontrols))
        )
      },
      no.readonly = TRUE
    )
    if (saveToFig) {
      dev.off()
    }
  }

  foldchanges <- cbind(normed[, 1:2], foldchanges)

  save(
    normed,
    file = paste0(outdir, EXPname, "_normCounts.RData")
  )
  save(
    foldchanges,
    file = paste0(outdir, EXPname, "_foldChanges.RData")
  )

  return(list(norm_counts = normed, logFCs = foldchanges))
}

ccr.logFCs2chromPos <- function(
  foldchanges,
  libraryAnnotation
) {

  sgRNAsIds <- foldchanges[["sgRNA"]]

  genes <- as.character(libraryAnnotation[sgRNAsIds, "GENES"])

  chrN <- as.character(libraryAnnotation[sgRNAsIds, "CHRM"])

  if (any(chrN == "X")) chrN[which(chrN == "X")] <- "23"
  if (any(chrN == "Y")) chrN[which(chrN == "Y")] <- "24"
  chrN <- as.numeric(chrN)

  startp <- as.numeric(as.character(libraryAnnotation[sgRNAsIds, "STARTpos"]))
  endp <- as.numeric(as.character(libraryAnnotation[sgRNAsIds, "ENDpos"]))

  if (ncol(foldchanges) > 3) {
    converted <- data.frame(
      chrN,
      startp,
      endp,
      genes,
      rowMeans(foldchanges[, 3:ncol(foldchanges)]),
      stringsAsFactors = FALSE
    )
  }else{
    converted <- data.frame(
      chrN,
      startp,
      endp,
      genes,
      foldchanges[, 3:ncol(foldchanges)],
      stringsAsFactors = FALSE
    )
  }

  rownames(converted) <- foldchanges[["sgRNA"]]
  converted <- converted[order(converted[["chrN"]], converted[["startp"]]), ]
  colnames(converted)[[5]] <- "avgFC"
  colnames(converted)[[1]] <- "CHR"
  BP <- converted[["startp"]] + (
    converted[["endp"]] - converted[["startp"]]
  ) / 2
  converted <- cbind(converted, BP)

  return(converted)
}

ccr.cleanChrm <- function(
  gwSortedFCs,
  CHR,
  display = TRUE,
  label = "",
  saveTO = NULL,
  min.ngenes = 3,
  ignoredGenes = NULL,
  capped = FALSE,
  corrMet = "mean",
  alpha = 0.01,
  nperm = 10000,
  p.method = "hybrid",
  min.width = 2,
  kmax = 25,
  nmin = 200,
  eta = 0.05,
  trim = 0.025,
  undo.splits = "none",
  undo.prune = 0.05,
  undo.SD = 3,

  # Start update for web version
  return.segments.unadj = TRUE,
  return.segments.adj = FALSE,
  nseed = 0xA5EED,
  verbose = 1
  # End update for web version

) {

  gwSortedFCs <- as.data.frame(gwSortedFCs)

  ID <- which(gwSortedFCs[["CHR"]] == CHR)
  gwSortedFCs <- gwSortedFCs[ID, ]

  # Start update for web version
  if (length(nseed) > 0) {
    set.seed(nseed)
  }
  # End update for web version

  my.CNA.object <- CNA(
    cbind(gwSortedFCs[["avgFC"]]),
    gwSortedFCs[["CHR"]],
    gwSortedFCs[["BP"]],
    data.type = "logratio",
    sampleid = paste(label, "Chr", CHR, "sgRNA FCs")
  )

  my.smoothed.CNA.object <- smooth.CNA(my.CNA.object)
  my.segment.smoothed.CNA.object <- segment(
    my.smoothed.CNA.object,
    verbose = verbose,
    alpha = alpha,
    nperm = nperm,
    p.method = p.method,
    min.width = min.width,
    kmax = kmax,
    nmin = nmin,
    eta = eta,
    trim = trim,
    undo.splits = undo.splits,
    undo.prune = undo.prune,
    undo.SD = undo.SD
  )


  # Start update for web version
  regions.unadj <- my.segment.smoothed.CNA.object[["output"]]
  nsegments <- nrow(regions.unadj)
  # End update for web version

  nGeneInSeg <- vector()
  guides <- vector()
  newFC <- gwSortedFCs[["avgFC"]]
  correction <- rep(0, length(newFC))

  for (i in seq_len(nsegments)) {
    idxs <- seq(
      my.segment.smoothed.CNA.object[["segRows"]][i, 1],
      my.segment.smoothed.CNA.object[["segRows"]][i, 2]
    )
    includedGenes <- unique(gwSortedFCs[idxs, "genes"])
    includedGuides <- ID[range(idxs)]

    # Start update for web version
    regions.unadj[i, "loc.start"] <- gwSortedFCs[min(idxs), "startp"]
    regions.unadj[i, "loc.end"] <- gwSortedFCs[max(idxs), "endp"]
    # End update for web version

    if (length(ignoredGenes) > 0) {
      nGeneInSeg[i] <- length(setdiff(includedGenes, ignoredGenes))
    } else {
      nGeneInSeg[i] <- length(includedGenes)
    }

    if (nGeneInSeg[i] >= min.ngenes) {
      orSign <- sign(newFC[idxs])

      if (corrMet == "mean") {
        newFC[idxs] <- newFC[idxs] - mean(newFC[idxs])
      } else {
        newFC[idxs] <- newFC[idxs] - median(newFC[idxs])
      }

      if (capped) {
        newFC[idxs[which(sign(newFC[idxs]) != orSign)]] <- 0
      }

      correction[idxs] <- -sign(mean(newFC[idxs]))
    }

    guides[i] <- paste(includedGuides, collapse = ", ")
  }

  # Start update for web version
  regions.unadj <- cbind(regions.unadj[, 2:ncol(regions.unadj)], guides)
  colnames(regions.unadj) <- c(
    "CHR", "startp", "endp", "n.sgRNAs", "avg.logFC", "guideIdx"
  )
  # End update for web version

  norm.my.CNA.object <- CNA(
    cbind(newFC),
    gwSortedFCs[["CHR"]],
    gwSortedFCs[["BP"]],
    data.type = "logratio",
    sampleid = paste(label, "Chr", CHR, "sgRNA FCs - post CRISPRcleanR")
  )

  norm.my.smoothed.CNA.obj <- smooth.CNA(norm.my.CNA.object)
  norm.my.seg.smoothed.CNA.obj <- segment(
    norm.my.smoothed.CNA.obj,
    verbose = verbose
  )

  if (!display) {
    saveTO <- NULL
  }

  if (length(saveTO)) {
    display <- TRUE
    path <- paste0(saveTO, label, "/")
    if (!file.exists(path)) {
      dir.create(path)
    }
    pdf(paste0(path, CHR, ".pdf"), width = 7.5, height = 7.5)
  }

  if (display) {
    withr::with_par(
      list(mfrow = c(2, 1)), {
        plot(
          my.segment.smoothed.CNA.object,
          main = paste(label, "Chr", CHR, "sgRNA FCs"),
          segcol = "black",
          xmaploc = FALSE,
          xlab = "",
          ylab = "logFC"
        )
        plot(
          norm.my.seg.smoothed.CNA.obj,
          main = paste(label, "Chr", CHR, "sgRNA FCs, post CRISPRcleanR"),
          segcol = "black",
          xmaploc = FALSE,
          xlab = "sgRNA index",
          ylab = "logFC"
        )
      },
      no.readonly = TRUE
    )
  }

  if (length(saveTO)) {
    dev.off()
  }

  correctedFC <- newFC
  gwSortedFCs <- cbind(gwSortedFCs, correction, correctedFC)

  # Start update for web version
  regions.adj <- norm.my.seg.smoothed.CNA.obj[["output"]]
  nsegments <- nrow(regions.adj)
  guides <- vector()
  for (i in seq_len(nsegments)) {
    idxs <- seq(
      norm.my.seg.smoothed.CNA.obj[["segRows"]][i, 1],
      norm.my.seg.smoothed.CNA.obj[["segRows"]][i, 2]
    )
    includedGenes <- unique(gwSortedFCs[idxs, "genes"])
    includedGuides <- ID[range(idxs)]
    regions.adj[i, "loc.start"] <- gwSortedFCs[min(idxs), "startp"]
    regions.adj[i, "loc.end"] <- gwSortedFCs[max(idxs), "endp"]
    guides[i] <- paste(includedGuides, collapse = ", ")
  }
  regions.adj <- cbind(regions.adj[, 2:ncol(regions.adj)], guides)
  colnames(regions.adj) <- c(
    "CHR", "startp", "endp", "n.sgRNAs", "avg.logFC", "guideIdx"
  )
  # End update for web version

  gwSortedFCs[, "guideIdx"] <- seq_len(nrow(gwSortedFCs))
  res <- list(correctedFCs = gwSortedFCs)
  if (return.segments.unadj) res[["regions"]] <- regions.unadj
  if (return.segments.adj) res[["regions.adj"]] <- regions.adj

  return(res)
}

ccr.GWclean <- function(
  gwSortedFCs,
  label = "",
  display = TRUE,
  saveTO = NULL,
  ignoredGenes = NULL,
  min.ngenes = 3,
  alpha = 0.01,
  nperm = 10000,
  p.method = "hybrid",
  min.width = 2,
  kmax = 25,
  nmin = 200,
  eta = 0.05,
  trim = 0.025,
  undo.splits = "none",
  undo.prune = 0.05,
  undo.SD = 3,
  return.segments.unadj = TRUE,
  return.segments.adj = FALSE,
  nseed = 0xA5EED,
  verbose = 1
) {

  gwSortedFCs <- as.data.frame(gwSortedFCs)
  CHRs <- as.character(sort(unique(gwSortedFCs[["CHR"]])))

  # Start update for web version
  corrected_logFCs <- NULL
  segments <- NULL
  segments_adj <- NULL
  if (verbose == 0) pb <- txtProgressBar(min = 0, max = length(CHRs), style = 3)
  # End update for web version

  for (i in seq_along(CHRs)) {
    CHR <- CHRs[i]
    res <- ccr.cleanChrm(
      gwSortedFCs = gwSortedFCs,
      CHR = CHR,
      display = display,
      label = label,
      saveTO = saveTO,
      ignoredGenes = ignoredGenes,
      min.ngenes = min.ngenes,
      alpha = alpha,
      nperm = nperm,
      p.method = p.method,
      min.width = min.width,
      kmax = kmax,
      nmin = nmin,
      eta = eta,
      trim = trim,
      undo.splits = undo.splits,
      undo.prune = undo.prune,
      undo.SD = undo.SD,

      # Start update for web version
      return.segments.unadj = return.segments.unadj,
      return.segments.adj = return.segments.adj,
      nseed = nseed,
      verbose = verbose
    )

    corrected_logFCs <- rbind(corrected_logFCs, res[["correctedFCs"]])
    if (return.segments.unadj) segments <- rbind(segments, res[["regions"]])
    if (return.segments.adj) segments_adj <- rbind(
      segments_adj,
      res[["regions.adj"]]
    )
    if (verbose == 0) setTxtProgressBar(pb, i)

  }
  if (verbose == 0) close(pb)

  ret <- list(
    corrected_logFCs = corrected_logFCs,
    SORTED_sgRNAs = rownames(gwSortedFCs)
  )
  if (return.segments.unadj) ret[["segments"]] <- segments
  if (return.segments.adj) ret[["segments_adj"]] <- segments_adj

  # End update for web version

  return(ret)
}

ccr.correctCounts <- function(
  CL,
  normalised_counts,
  correctedFCs_and_segments,
  libraryAnnotation,
  minTargetedGenes = 3,
  OutDir = "./",
  ncontrols = 1,
  verbose = 1
) {

  normalised_counts <- normalised_counts[
    which(!is.na(normalised_counts[, 1])),
  ]
  rownames(normalised_counts) <- normalised_counts[["sgRNA"]]
  numdata <- normalised_counts[, 3:ncol(normalised_counts)]
  rownames(numdata) <- normalised_counts[["sgRNA"]]

  segments <- correctedFCs_and_segments[["segments"]]

  segment_guides <- list()

  nsegments <- nrow(segments)

  correctedCounts <- numdata

  # Start update for web version
  if (verbose == 0) pb <- txtProgressBar(min = 0, max = nsegments, style = 3)
  # End update for web version

  for (i in seq_len(nsegments)) {
    if (verbose > 0) print(paste0(
      CL, ": correcting counts on segment ",
      i, " (of ", nsegments, ")"
    ))
    current_idxs <- as.character(segments[["guideIdx"]][i])
    current_idxs <- as.numeric(unlist(str_split(current_idxs, ", ")[[1]]))
    segment_guides[[i]] <- correctedFCs_and_segments[["SORTED_sgRNAs"]][
      current_idxs[[1]]:current_idxs[[2]]
    ]

    ntarg <- length(
      unique(as.character(libraryAnnotation[segment_guides[[i]], "GENES"]))
    )

    guides <- segment_guides[[i]]
    if (ntarg >= minTargetedGenes) {

      if (ncontrols == 1) {
        c <- numdata[unlist(guides), 1]
      }else{
        c <- rowMeans(numdata[unlist(guides), seq_len(ncontrols)])
      }

      N <- correctedFCs_and_segments[["corrected_logFCs"]][
        guides,
        "correctedFC"
      ]
      reverted <- c * 2 ^ N

      nreps <- ncol(numdata) - ncontrols

      if (nreps > 1) {
        proportions <- (
          numdata[guides, (ncontrols + 1):ncol(numdata)] + 1
        ) / (
          matrix(
            rep(
              rowSums(numdata[guides, (ncontrols + 1):ncol(numdata)] + 1),
              nreps
            ),
            length(guides),
            nreps
          )
        )
        revertedCounts <- (reverted * nreps) * proportions
      }else{
        revertedCounts <- reverted
      }

      correctedCounts[
        guides,
        (ncontrols + 1):ncol(numdata)
      ] <- revertedCounts
    }

    # Start update for web version
    if (verbose == 0) setTxtProgressBar(pb, i)
  }
  if (verbose == 0) close(pb)
  # End update for web version

  adjusted <- as.data.frame(cbind(
    normalised_counts[["sgRNA"]],
    normalised_counts[["gene"]],
    correctedCounts
  ), stringAsFactors = FALSE)
  colnames(adjusted) <- colnames(normalised_counts)
  rownames(adjusted) <- NULL
  correctedCounts <- adjusted

  save(
    correctedCounts,
    file = paste0(OutDir, CL, "_correctedCounts.RData")
  )
  return(correctedCounts)
}

#### Utils
ccr.genes2sgRNAs <- function(
  libraryAnnotation,
  genes
) {
  notIncludedGenes <- genes[
    which(!is.element(genes, libraryAnnotation[["GENES"]]))
  ]
  if (length(notIncludedGenes) > 0) {
    warning(paste(
      "No sgRNAs targeting the following genes in this library:",
      paste(notIncludedGenes, collapse = ", ")
    ))
  }
  sgRNAs <- rownames(libraryAnnotation)[
    which(is.element(libraryAnnotation[["GENES"]], genes))
  ]
  return(sgRNAs)
}

ccr.geneMeanFCs <- function(
  sgRNA_FCprofile,
  libraryAnnotation
) {
  FCsprofile <- sgRNA_FCprofile
  FCsprofile <- aggregate(
    sgRNA_FCprofile ~ libraryAnnotation[names(sgRNA_FCprofile), "GENES"],
    FUN = mean
  )
  nn <- as.character(
    FCsprofile[["libraryAnnotation[names(sgRNA_FCprofile), \"GENES\"]"]]
  )
  FCsprofile <- FCsprofile[["sgRNA_FCprofile"]]
  names(FCsprofile) <- as.character(nn)

  return(FCsprofile)
}

ccr.sgRNAmeanFCs <- function(
  foldchanges
) {
  FCs <- foldchanges[["correctedFC"]]
  names(FCs) <- rownames(foldchanges)

  return(FCs)
}

ccr.get.gdsc1000.AMPgenes <- function(
  cellLine,
  minCN = 8,
  exact = FALSE,
  GDSC.geneLevCNA = NULL,
  GDSC.CL_annotation = NULL
) {

  if (is.null(GDSC.geneLevCNA)) {
    data(GDSC.geneLevCNA, envir = environment())
  }

  if (is.null(GDSC.CL_annotation)) {
    data(GDSC.CL_annotation, envir = environment())
  }

  if (
    !is.element(cellLine, GDSC.CL_annotation[["CL.name"]]) &
    !is.element(cellLine, GDSC.CL_annotation[["COSMIC.ID"]])
  ) {
    stop(paste0(
      "Cell line not found:",
      " Please provide a valid cell line cosmic identifier or name",
      " (see available cell lines here: http://www.cancerrxgene.org",
      "/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS1E.xlsx)"),
      call. = TRUE,
      domain = NULL
    )
  }

  if (is.element(cellLine, GDSC.CL_annotation[["COSMIC.ID"]])) {
    cid <- cellLine
  } else {
    cid <- as.character(
      GDSC.CL_annotation[["COSMIC.ID"]][
        GDSC.CL_annotation[["CL.name"]] == cellLine
      ]
    )
  }

  genes <- rownames(GDSC.geneLevCNA)

  if (!exact) {
    cnRange <- c(minCN:16)
  } else {
    cnRange <- minCN
  }

  if (length(intersect(cid, colnames(GDSC.geneLevCNA))) > 0) {
    amplifiedGenes <- NULL

    for (ar in cnRange) {
      amplifiedGenes <- c(
        amplifiedGenes,
        genes[unique(c(grep(paste0(",", ar, ","), GDSC.geneLevCNA[, cid])))]
      )
    }

    amplifiedGenes <- sort(amplifiedGenes)

    GENES <- amplifiedGenes
    CNAval <- GDSC.geneLevCNA[GENES, cid]

    CNAval <- unlist(str_split(CNAval, ","))

    if (length(GENES) > 0) {
      CN <- CNAval[seq(2, length(GENES) * 4, 4)]

      amplifiedGenes <- cbind(GENES, CN)
      amplifiedGenes <- as.data.frame(amplifiedGenes, stringsAsFactors = FALSE)

      colnames(amplifiedGenes)[[1]] <- "Gene"
      amplifiedGenes[, 2] <- as.numeric(amplifiedGenes[, 2])
      return(amplifiedGenes)
    }else{
      print(paste0("No amplified genes for this cell line",
        " according to the selected threshold"
      ))
      return(NULL)
    }
  }else{
    print(paste0(
      "There is no data available for this cell line",
      " in the built in CNA data frame. Provide gene level CN values in input"
    ))
    return(NULL)
  }
}

ccr.get.nonExpGenes <- function(
  cellLine,
  th = 0.05,
  amplified = FALSE,
  minCN = 8,
  RNAseq.fpkms = NULL,
  GDSC.CL_annotation = NULL
) {

  if (is.null(GDSC.CL_annotation)) {
    data(GDSC.CL_annotation, envir = environment())
  }

  if (is.null(RNAseq.fpkms)) {
    data(RNAseq.fpkms, envir = environment())
  }

  if (
    !is.element(cellLine, GDSC.CL_annotation[["CL.name"]]) &
    !is.element(cellLine, GDSC.CL_annotation[["COSMIC.ID"]])
  ) {
    stop(paste0(
      "Cell line not found:",
      " Please provide a valid cell line cosmic identifier or name",
      " (see available cell lines here: http://www.cancerrxgene.org/",
      "gdsc1000/GDSC1000_WebResources//Data/suppData/TableS1E.xlsx)"),
    call. = TRUE,
    domain = NULL)
  }

  if (is.element(cellLine, GDSC.CL_annotation[["COSMIC.ID"]])) {
    cid <- cellLine
  } else {
    cid <- as.character(
      GDSC.CL_annotation[["COSMIC.ID"]][
        GDSC.CL_annotation[["CL.name"]] == cellLine
      ]
    )
  }

  if (amplified) {
    amplifiedGenes <- ccr.get.gdsc1000.AMPgenes(
      cellLine = cellLine,
      minCN = minCN
    )
    amplifiedGenes <- unique(amplifiedGenes[["Gene"]])
  }

  if (length(intersect(cid, colnames(RNAseq.fpkms))) > 0) {
    expProf <- RNAseq.fpkms[, cid]
    notExpGenes <- names(which(expProf < th))

    if (amplified) {
      notExpGenes <- intersect(notExpGenes, amplifiedGenes)
    }
    return(notExpGenes)
  }else{
    print("No RNAseq data available for this cell line")
    return(NULL)
  }
}

ccr.get.CCLEgisticSets <- function(
  cellLine,
  CCLE.gisticCNA = NULL,
  GDSC.CL_annotation = NULL
) {

  if (is.null(GDSC.CL_annotation)) {
    data(GDSC.CL_annotation, envir = environment())
  }

  if (is.null(CCLE.gisticCNA)) {
    data(CCLE.gisticCNA, envir = environment())
  }

  if (
    !is.element(cellLine, GDSC.CL_annotation[["CL.name"]]) &
    !is.element(cellLine, GDSC.CL_annotation[["COSMIC.ID"]])
  ) {
    stop(paste0(
      "Cell line not found:",
      " Please provide a valid cell line cosmic identifier or name ",
      "(see available cell lines here: http://www.cancerrxgene.org/",
      "gdsc1000/GDSC1000_WebResources//Data/suppData/TableS1E.xlsx)"),
      call. = TRUE,
      domain = NULL
    )
  }

  if (is.element(cellLine, GDSC.CL_annotation[["COSMIC.ID"]])) {
    cid <- cellLine
  } else {
    cid <- as.character(
      GDSC.CL_annotation[["COSMIC.ID"]][
        GDSC.CL_annotation[["CL.name"]] == cellLine
      ]
    )
  }

  if (length(intersect(cid, colnames(CCLE.gisticCNA))) > 0) {
    gisticProf <- CCLE.gisticCNA[, cid]
    tmpGe <- rownames(CCLE.gisticCNA)[which(gisticProf == -2)]
    gistic_min2 <- tmpGe

    tmpGe <- rownames(CCLE.gisticCNA)[which(gisticProf == -1)]
    gistic_min1 <- tmpGe

    tmpGe <- rownames(CCLE.gisticCNA)[which(gisticProf == 0)]
    gistic_wt <- tmpGe

    tmpGe <- rownames(CCLE.gisticCNA)[which(gisticProf == 1)]
    gistic_plus1 <- tmpGe

    tmpGe <- rownames(CCLE.gisticCNA)[which(gisticProf == 2)]
    gistic_plus2 <- tmpGe

    return(list(
      gm2 = gistic_min2,
      gm1 = gistic_min1,
      gz = gistic_wt,
      gp1 = gistic_plus1,
      gp2 = gistic_plus2
    ))
  }else{
    print("No gistic CNA scores available for this cell line")
    return(NULL)
  }
}

ccr.PlainTsvFile <- function(
  sgRNA_count_object,
  fprefix = "",
  path = "./"
) {

    currentFileContent <- sgRNA_count_object

    fname <- paste0(path, fprefix, "_sgRNA_count.tsv")

    write.table(
      currentFileContent,
      quote = FALSE,
      row.names = FALSE,
      sep = "\t",
      file = fname
    )

    return(fname)
}

ccr.ExecuteMageck <- function(
  mgckInputFile,
  expName = "expName",
  normMethod = "none",
  outputPath = "./",
  ncontrols = 1,
  path_to_mageck = "mageck",
  verbose = 1
) {
  fc <- read.table(mgckInputFile, sep = "\t", header = TRUE)

  # Start update for web version
  Cnames <- colnames(fc)[3:(2 + ncontrols)]
  Tnames <- colnames(fc)[(3 + ncontrols):ncol(fc)]
  textbunch <- paste0(path_to_mageck, " test -k")
  textbunch <- paste0(
    textbunch, " ",
    mgckInputFile,
    " -c \"", paste0(Cnames, collapse = ","),
    "\" -t \"", paste0(Tnames, collapse = ","),
    "\" -n ", outputPath, expName,
    " --norm-method ", normMethod
  )

  system(textbunch, intern = TRUE, wait = TRUE)
  # End update for web version

  geneSummaryFN <- paste0(outputPath, expName, ".gene_summary.txt")

  return(geneSummaryFN)
}

# End update for web version

#### Assessment and visualisation
ccr.ROC_Curve <- function(
  FCsprofile,
  positives,
  negatives,
  display = TRUE,
  FDRth = NULL,
  expName = NULL
) {

  FCsprofile <- FCsprofile[
    intersect(c(positives, negatives),
    names(FCsprofile))
  ]

  predictions <- FCsprofile
  observations <- is.element(names(FCsprofile), positives) + 0
  names(observations) <- names(predictions)

  RES <- roc(observations, predictions, direction = ">", quiet = TRUE)

  if (display) {
    plot(
      RES,
      col = "blue",
      lwd = 3,
      xlab = "TNR",
      ylab = "Recall",
      main = expName
    )
  }

  SENS <- NULL
  threshold <- NULL
  COORS <- coords(
    RES,
    "all",
    ret = c("threshold", "ppv", "sensitivity", "specificity"),
    transpose = TRUE
  )
  if (length(FDRth) > 0) {

    FDR5percTh <- max(
      COORS["threshold", which(COORS["ppv", ] >= (1 - FDRth))]
    )

    threshold <- COORS[
      "threshold",
      min(which(COORS["threshold", ] <= FDR5percTh))
    ]

    SENS <- COORS[
      "sensitivity",
      min(which(COORS["threshold", ] <= FDR5percTh))
    ]
    SPEC <- COORS[
      "specificity",
      min(which(COORS["threshold", ] <= FDR5percTh))
    ]
    if (display) {
      abline(h = SENS, lty = 2)
    }
  }

  if (display) {
    if (length(SENS) == 0) {
      legend(
        "bottomright",
        paste("AUC = ", format(RES[["auc"]], digits = 3)),
        bty = "n"
      )
    }else{
      legend(
        "bottomright",
        c(
          paste0(
            "- - Recall ", 100 * FDRth,
            "%FDR = ", format(SENS, digits = 3)
          ),
          paste("AUC = ", format(RES[["auc"]], digits = 3))
        ),
        bty = "n"
      )
    }
  }

  COORS <- t(COORS[c("specificity", "sensitivity", "threshold"), ])
  RES <- list(
    AUC = RES[["auc"]],
    Recall = SENS,
    sigthreshold = threshold,
    curve = COORS
  )

  ### threshold, and recall at fixed FDR to be returned
  return(RES)
}

ccr.PrRc_Curve <- function(
  FCsprofile,
  positives,
  negatives,
  display = TRUE,
  FDRth = NULL,
  expName = NULL
) {

  FCsprofile <- FCsprofile[
    intersect(c(positives, negatives),
    names(FCsprofile))
  ]

  predictions <- -FCsprofile
  observations <- is.element(names(FCsprofile), positives) + 0
  names(observations) <- names(predictions)

  prc <- pr.curve(
    scores.class0 = predictions,
    weights.class0 = observations,
    curve = TRUE,
    sorted = TRUE
  )

  PRECISION <- prc[["curve"]][, 2]
  RECALL <- prc[["curve"]][, 1]

  if (display) {
    plot(
      RECALL,
      PRECISION,
      col = "blue",
      lwd = 3,
      xlab = "Recall",
      ylab = "Precision",
      type = "l",
      xlim = c(0, 1),
      ylim = c(0, 1),
      main = expName
    )
  }

  SENS <- NULL
  threshold <- NULL
  if (length(FDRth) > 0) {

    res <- ccr.ROC_Curve(
      FCsprofile,
      positives,
      negatives,
      display = FALSE,
      FDRth = FDRth
    )

    SENS <- res[["Recall"]]
    threshold <- res[["sigthreshold"]]

    if (display) {
        abline(h = 1 - FDRth, lty = 1)
        abline(v = SENS, lty = 2)
    }
  }

  RND <- sum(observations) / length(observations)
  if (display) {
    if (length(SENS) == 0) {
      legend(
        "bottomleft",
        paste("AUC = ",
        format(prc[["auc.integral"]], digits = 3)),
        bty = "n"
      )
    } else {
        legend(
          "bottomleft",
          c(
            paste0(
              "- - Recall ", 100 * FDRth,
              "%FDR = ", format(SENS, digits = 3)
            ),
            paste("AUC = ", format(prc[["auc.integral"]], digits = 3))
          ),
          bty = "n")
    }
    abline(h = RND, col = "lightgray")
  }
  curve <- prc[["curve"]]
  colnames(curve) <- c("recall", "precision", "threshold")
  RES <- list(
    AUC = prc[["auc.integral"]],
    Recall = SENS,
    sigthreshold = threshold,
    curve = curve,
    RND = RND
  )
  # ### threshold, and recall at fixed FDR to be returned
  return(RES)
}

ccr.randomised_ROC <- function(
  FCs,
  PERCrandn,
  ntrials,
  positives,
  negatives,
  LibraryAnnotation
) {

  posGuides <- ccr.genes2sgRNAs(
    libraryAnnotation = LibraryAnnotation,
    genes = positives
  )
  negGuides <- ccr.genes2sgRNAs(
    libraryAnnotation = LibraryAnnotation,
    genes = negatives
  )

  guidesToShuffle <- intersect(union(posGuides, negGuides), names(FCs))

  idGuidesToShuffle <- match(guidesToShuffle, names(FCs))

  nguides <- round(length(guidesToShuffle) * PERCrandn / 100)

  rndSENSITIVITY <- NULL
  rndSPECIFICITY <- NULL
  rndAUROC <- NULL

  rndRECALL <- NULL
  rndPRECISION <- NULL
  rndAUPRC <- NULL

  rndRECALL_at_fixedFDR <- NULL

  for (i in seq_len(ntrials)) {
    print(i)

    RNDFCs <- FCs
    mid <- sample(length(guidesToShuffle), nguides)
    toShuffle <- idGuidesToShuffle[mid]

    otherGuides <- setdiff(seq_along(FCs), toShuffle)

    toReplace <- otherGuides[sample(length(otherGuides), nguides)]

    names(RNDFCs)[
      c(toShuffle, toReplace)
    ] <- names(FCs)[
      c(toReplace, toShuffle)
    ]

    RNDgeneFCs <- ccr.geneMeanFCs(RNDFCs, LibraryAnnotation)

    RESroc <- ccr.ROC_Curve(
      RNDgeneFCs,
      positives,
      negatives,
      FDRth = 0.05,
      display = FALSE
    )
    rndSENSITIVITY <- rbind(rndSENSITIVITY, RESroc[["curve"]][, "sensitivity"])
    rndSPECIFICITY <- rbind(rndSPECIFICITY, RESroc[["curve"]][, "specificity"])
    rndAUROC <- c(rndAUROC, RESroc[["AUC"]])

    RESprrc <- ccr.PrRc_Curve(
      RNDgeneFCs,
      positives,
      negatives,
      FDRth = 0.05,
      display = FALSE
    )
    rndRECALL <- rbind(rndRECALL, RESprrc[["curve"]][, "recall"])
    rndPRECISION <- rbind(rndPRECISION, RESprrc[["curve"]][, "precision"])
    rndAUPRC <- c(rndAUPRC, RESprrc[["AUC"]])
    rndRECALL_at_fixedFDR <- c(rndRECALL_at_fixedFDR, RESprrc[["Recall"]])
  }

  return(list(
    rndSENSITIVITY = rndSENSITIVITY,
    rndSPECIFICITY = rndSPECIFICITY,
    rndAUROC = rndAUROC,
    rndRECALL = rndRECALL,
    rndPRECISION = rndPRECISION,
    rndAUPRC = rndAUPRC,
    rndRECALL_at_fixedFDR = rndRECALL_at_fixedFDR
  ))
}

ccr.VisDepAndSig <- function(
  FCsprofile,
  SIGNATURES,
  TITLE = "",
  pIs = NULL,
  nIs = NULL,
  th = 0.05,
  plotFCprofile = TRUE
) {

  sigNames <- names(SIGNATURES)
  withr::with_par(
    list(mar = c(5, 5, 8, 0)), {
      nsig <- length(SIGNATURES)

      if (length(pIs) > 0 & length(nIs) > 0) {
        RES <- ccr.fixedFDRthreshold(
          FCsprofile = FCsprofile,
          TruePositives = SIGNATURES[[pIs]],
          TrueNegatives = SIGNATURES[[nIs]],
          th = th
        )
        FDR5percRANK <- RES[["RANK"]]
      } else {
        FDR5percRANK <- NULL
      }

      layout(t(matrix(
        c(rep(1, 5), seq_len(nsig + 1), rep(1, 5), seq_len(nsig + 1)),
        nsig + 6,
        2
      )))

      nelements <- length(FCsprofile)

      if (plotFCprofile) {
        plot(
          sort(FCsprofile, decreasing = TRUE),
          rev(seq_along(FCsprofile)),
          ylim = c(nelements, 1),
          pch = 16,
          frame.plot = FALSE,
          yaxt = "n",
          xlim = c(min(FCsprofile), max(FCsprofile) + 1),
          ylab = "Depletion Rank",
          xlab = "Log FC",
          log = "y",
          cex = 1.2,
          cex.lab = 1.5,
          cex.axis = 1.2,
          main = TITLE,
          cex.main = 1.5
        )
      } else {
        plot(
          0,
          0,
          ylim = c(nelements, 1),
          pch = 16,
          frame.plot = FALSE,
          yaxt = "n",
          xlim = c(min(FCsprofile), max(FCsprofile) + 1),
          ylab = "Depletion Rank",
          xlab = "Log FC",
          log = "y",
          cex = 1.2,
          cex.lab = 1.5,
          cex.axis = 1.2,
          main = TITLE,
          cex.main = 1.5
        )
      }
      lines(x = c(0, 0), y = c(1, (length(FCsprofile)) + 10000), lty = 2)

      axis(
        side = 2,
        c(1, 10, 100, 1000, 10000, nelements),
        labels = c(1, 10, 100, 1000, 10000, nelements),
        las = 2,
        cex.axis = 1
      )

      withr::with_par(
        list(xpd = TRUE), {
          lines(
            c(min(FCsprofile), max(FCsprofile) + 1),
            c(FDR5percRANK, FDR5percRANK),
            col = "red",
            lwd = 3,
            lty = 2
          )

          text(
            x = max(FCsprofile),
            y = FDR5percRANK,
            "5%FDR",
            pos = 3,
            col = "red",
            cex = 1
          )
        },
        no.readonly = TRUE
      )

      TPR <- vector()

      withr::with_par(
        list(xpd = TRUE, mar = c(5, 0, 8, 0)), {
          for (i in seq_along(SIGNATURES)) {
            plot(
              1,
              1,
              xlim = c(0, 1),
              ylim = c(length(FCsprofile), 1),
              xaxt = "n",
              yaxt = "n",
              ylab = "",
              xlab = "",
              col = NA,
              log = "y",
              frame.plot = FALSE
            )

            hitPositions <- match(SIGNATURES[[i]], names(sort(FCsprofile)))
            hitPositions <- hitPositions[!is.na(hitPositions)]

            abline(
              h = hitPositions[hitPositions <= FDR5percRANK],
              lwd = 2,
              col = "blue"
            )
            abline(
              h = hitPositions[hitPositions > FDR5percRANK],
              lwd = 2,
              col = "gray"
            )

            TPR[i] <- (
              length(which(hitPositions <= FDR5percRANK)) /
              length(hitPositions)
            )

            lines(
              c(0, 1),
              c(FDR5percRANK, FDR5percRANK),
              col = "red",
              lwd = 3,
              lty = 2
            )

            text(
              0.4,
              1,
              labels = sigNames[i],
              pos = 4,
              offset = 0,
              srt = 80,
              cex = 1
            )

            if (i < (length(SIGNATURES) - 1)) {
              text(
                0.6,
                nelements + 50000,
                paste0(round(100 * TPR[i]), "%"),
                cex = 1.2,
                col = "blue"
              )
            }
          }
        },
        no.readonly = TRUE
      )
      names(TPR) <- sigNames
    },
    no.readonly = TRUE
  )
  return(TPR)
}

ccr.perf_statTests <- function(
  cellLine,
  libraryAnnotation,
  correctedFCs,
  outDir = "./",
  GDSC.geneLevCNA = NULL,
  CCLE.gisticCNA = NULL,
  RNAseq.fpkms = NULL,
  GDSC.CL_annotation = NULL,

  # Start update for web version
  verbose = 1
  # End update for web version
) {

  nnames <- c(
    "Dep (PN)",
    "Dep (Gistic)",
    "notExp",
    "Amp (Gistic +1)",
    "Amp (Gistic +2)",
    "Amp (Gistic +1) notExp",
    "Amp (Gistic +2) notExp",
    "Amp (PNs >= 2)",
    "Amp (PNs >= 4)",
    "Amp (PNs >= 8)",
    "Amp (PNs >= 10)",
    "Amp (PNs >= 2) notExp",
    "Amp (PNs >= 4) notExp",
    "Amp (PNs >= 8) notExp",
    "Amp (PNs >= 10) notExp",
    "Essential",
    "BAGEL Essential",
    "BAGEL Ess.Only",
    "BAGEL nonEssential",
    "whole-library"
  )

  guideSets <- ccr.get.guideSets(
    cellLine,
    GDSC.geneLevCNA,
    CCLE.gisticCNA,
    RNAseq.fpkms,
    libraryAnnotation = libraryAnnotation,
    GDSC.CL_annotation = GDSC.CL_annotation
  )
  PVALS <- vector()
  PVALSn <- vector()

  EFFsizes <- vector()
  EFFsizesN <- vector()

  SIGNS <- vector()
  SIGNSn <- vector()

  # Start update for web version
  if (verbose == 0) {
    pb <- txtProgressBar(min = 0, max = length(guideSets), style = 3)
  }

  for (i in seq_along(guideSets)) {
    if (verbose > 0) {
      print(paste("Testing sgRNAs targeting:", nnames[i], "genes"))
    }
    # End update for web version

    currentSet <- guideSets[[i]]

    if (length(currentSet) > 1) {
      tt <- t.test(
        correctedFCs[currentSet, "avgFC"],
        correctedFCs[setdiff(rownames(correctedFCs), currentSet), "avgFC"]
      )

      ttn <- t.test(
        correctedFCs[currentSet, "correctedFC"],
        correctedFCs[setdiff(rownames(correctedFCs), currentSet), "correctedFC"]
      )

      EFFsizes[i] <- ccr.cohens_d(
        x = correctedFCs[
          currentSet,
          "avgFC"
        ],
        y = correctedFCs[
          setdiff(rownames(correctedFCs), currentSet),
          "avgFC"
        ]
      )
      EFFsizesN[i] <- ccr.cohens_d(
        x = correctedFCs[
          currentSet,
          "correctedFC"
        ],
        y = correctedFCs[
          setdiff(rownames(correctedFCs), currentSet),
          "correctedFC"
        ]
      )

      PVALS[i] <- tt[["p.value"]]
      PVALSn[i] <- ttn[["p.value"]]

      SIGNS[i] <- sign(tt[["estimate"]][[2]] - tt[["estimate"]][[1]])
      SIGNSn[i] <- sign(ttn[["estimate"]][[2]] - ttn[["estimate"]][[1]])

    }else{
      PVALS[i] <- NA
      PVALSn[i] <- NA

      EFFsizes[i] <- NA
      EFFsizesN[i] <- NA

      SIGNS[i] <- NA
      SIGNSn[i] <- NA
    }

    # Start update for web version
    if (verbose == 0) {
      setTxtProgressBar(pb, i)
    }
  }
  if (verbose == 0) {
    close(pb)
  }
  # End update for web version

  pdf(paste0(outDir, cellLine, "_bp.pdf"), width = 15, height = 10)
  AT <- c(
    1, 2, 3, 4.5, 5.5, 6.5, 7.5, 9, 10,
    11, 12, 13, 14, 15, 16, 17.5, 18.5, 19.5, 20.5, 22
  )
  withr::with_par(
    list(mfrow = c(2, 1), mar = c(6, 4, 4, 1), xpd = NA), {
      plotpar <- boxplot(
        correctedFCs[guideSets[[1]], "avgFC"],
        correctedFCs[guideSets[[2]], "avgFC"],
        correctedFCs[guideSets[[3]], "avgFC"],
        correctedFCs[guideSets[[4]], "avgFC"],
        correctedFCs[guideSets[[5]], "avgFC"],
        correctedFCs[guideSets[[6]], "avgFC"],
        correctedFCs[guideSets[[7]], "avgFC"],
        correctedFCs[guideSets[[8]], "avgFC"],
        correctedFCs[guideSets[[9]], "avgFC"],
        correctedFCs[guideSets[[10]], "avgFC"],
        correctedFCs[guideSets[[11]], "avgFC"],
        correctedFCs[guideSets[[12]], "avgFC"],
        correctedFCs[guideSets[[13]], "avgFC"],
        correctedFCs[guideSets[[14]], "avgFC"],
        correctedFCs[guideSets[[15]], "avgFC"],
        correctedFCs[guideSets[[16]], "avgFC"],
        correctedFCs[guideSets[[17]], "avgFC"],
        correctedFCs[guideSets[[18]], "avgFC"],
        correctedFCs[guideSets[[19]], "avgFC"],
        correctedFCs[["avgFC"]],
        at = AT,
        names = nnames,
        las = 2,
        outline = FALSE,
        frame.plot = FALSE,
        main = cellLine,
        ylab = "pre-CRISPRcleanR sgRNA log2 FC",
        col = c(
          "#B0B0B0", "#898989", "#626262",
          c("red2", "red4"),
          "cornflowerblue", "blue",
          c("#FF6EB4", "#FF4978", "#FF243C", "#FF0000"),
          "#8EE5EE", "#5E98CD", "#2F4CAC", "#00008B",
          "darkgreen", "green", "darkcyan", "bisque4", "white"
        )
      )

      sigVec <- rep("", 19)
      sigVec[which(PVALS < 0.05)] <- "."
      sigVec[which(PVALS < 0.05 & EFFsizes > 0.5)] <- "*"
      sigVec[which(PVALS < 0.05 & EFFsizes > 1)] <- "**"
      sigVec[which(PVALS < 0.05 & EFFsizes > 2)] <- "***"

      posY <- (
        max(correctedFCs[["avgFC"]]) - min(correctedFCs[["avgFC"]])
      ) / 30 + max(correctedFCs[["avgFC"]])

      for (i in 1:19) {
        text(AT[i], posY, sigVec[i], cex = 1.5)
      }
      dd <- mean(correctedFCs[["avgFC"]])
      lines(x = c(0.2, 22.8), y = c(dd, dd), lty = 2, col = "red", lwd = 2)
      withr::with_par(
        list(mar = c(2, 4, 8, 1), xpd = NA), {
          boxplot(
            correctedFCs[guideSets[[1]], "correctedFC"],
            correctedFCs[guideSets[[2]], "correctedFC"],
            correctedFCs[guideSets[[3]], "correctedFC"],
            correctedFCs[guideSets[[4]], "correctedFC"],
            correctedFCs[guideSets[[5]], "correctedFC"],
            correctedFCs[guideSets[[6]], "correctedFC"],
            correctedFCs[guideSets[[7]], "correctedFC"],
            correctedFCs[guideSets[[8]], "correctedFC"],
            correctedFCs[guideSets[[9]], "correctedFC"],
            correctedFCs[guideSets[[10]], "correctedFC"],
            correctedFCs[guideSets[[11]], "correctedFC"],
            correctedFCs[guideSets[[12]], "correctedFC"],
            correctedFCs[guideSets[[13]], "correctedFC"],
            correctedFCs[guideSets[[14]], "correctedFC"],
            correctedFCs[guideSets[[15]], "correctedFC"],
            correctedFCs[guideSets[[16]], "correctedFC"],
            correctedFCs[guideSets[[17]], "correctedFC"],
            correctedFCs[guideSets[[18]], "correctedFC"],
            correctedFCs[guideSets[[19]], "correctedFC"],
            correctedFCs[["correctedFC"]],
            at = AT,
            outline = FALSE,

            # Start update for web version
            names = rep("", 20),
            frame.plot = FALSE,
            main = "",
            ylab = "post-CRISPRcleanR sgRNA log2 FC",
            # End update for web version

            col = c(
              "#B0B0B0", "#898989", "#626262",
              c("red2", "red4"),
              "cornflowerblue", "blue",
              c("#FF6EB4", "#FF4978", "#FF243C", "#FF0000"),
              "#8EE5EE", "#5E98CD", "#2F4CAC", "#00008B",
              "darkgreen", "green", "darkcyan", "bisque4", "white"
            )
          )

          sigVec1 <- rep("", 19)
          sigVec1[which(PVALSn < 0.05)] <- "."
          sigVec1[which(PVALSn < 0.05 & EFFsizesN > 0.5)] <- "*"
          sigVec1[which(PVALSn < 0.05 & EFFsizesN > 1)] <- "**"
          sigVec1[which(PVALSn < 0.05 & EFFsizesN > 2)] <- "***"

          posY <- (
            max(correctedFCs[["correctedFC"]]) -
            min(correctedFCs[["correctedFC"]])
          ) / 30 + max(correctedFCs[["correctedFC"]])

          for (i in 1:19) {
            text(AT[i], posY, sigVec1[i], cex = 1.5)
          }
          dd <- mean(correctedFCs[["correctedFC"]])
          lines(x = c(0.2, 22.8), y = c(dd, dd), lty = 2, col = "red", lwd = 2)

          legend(
            "bottomleft",
            legend = c(
              ".     p < 0.05",
              "*     p < 0.05, Effect size > 0.5",
              "**    p < 0.05, Effect size > 1",
              "***   p < 0.05, Effect size > 2"
            )
          )
        },
        no.readonly = TRUE
      )
    },
    no.readonly = TRUE
  )
  dev.off()

  EFFsizes <- rbind(EFFsizes, EFFsizesN)
  rownames(EFFsizes) <- c("pre-CRISPRcleanR", "post-CRISPRcleanR")
  colnames(EFFsizes) <- nnames[1:19]

  PVALS <- rbind(PVALS, PVALSn)
  rownames(PVALS) <- c("pre-CRISPRcleanR", "post-CRISPRcleanR")
  colnames(PVALS) <- nnames[1:19]

  SIGNS <- rbind(SIGNS, SIGNSn)
  rownames(SIGNS) <- c("pre-CRISPRcleanR", "post-CRISPRcleanR")
  colnames(SIGNS) <- nnames[1:19]

  return(list(PVALS = PVALS, SIGNS = SIGNS, EFFsizes = EFFsizes))
}

ccr.multDensPlot <- function(
  TOPLOT,
  COLS,
  XLIMS,
  TITLE,
  LEGentries,
  XLAB = ""
) {
  YM <- vector()
  for (i in seq_along(TOPLOT)) {
    YM[i] <- max(TOPLOT[[i]][["y"]], na.rm = TRUE)
  }

  Ymax <- max(YM, na.rm = TRUE)

  plot(
    0,
    0,
    col = NA,
    ylab = "density",
    xlab = XLAB,
    xlim = XLIMS,
    ylim = c(0, Ymax),
    type = "l",
    main = TITLE
  )

  for (i in seq_along(TOPLOT)) {
    cord.x <- c(TOPLOT[[i]][["x"]])
    cord.y <- c(TOPLOT[[i]][["y"]])
    rgbc <- col2rgb(COLS[i])
    currCol <- rgb(
      rgbc[[1]],
      rgbc[[2]],
      rgbc[[3]],
      alpha = 100,
      maxColorValue = 255
    )
    polygon(cord.x, cord.y, col = currCol, border = NA)
    lines(TOPLOT[[i]], col = COLS[i], lwd = 3)
    if (i == 1) {
      legend("topleft", legend = LEGentries, col = COLS, lwd = 3, bty = "n")
    }
  }
}

ccr.perf_distributions <- function(
  cellLine,
  correctedFCs,
  GDSC.geneLevCNA = NULL,
  CCLE.gisticCNA = NULL,
  RNAseq.fpkms = NULL,
  minCNs = c(8, 10),
  libraryAnnotation = NULL,
  GDSC.CL_annotation = NULL
) {
  guideSets <- ccr.get.guideSets(
    cellLine,
    GDSC.geneLevCNA,
    CCLE.gisticCNA,
    RNAseq.fpkms,
    libraryAnnotation = libraryAnnotation,
    GDSC.CL_annotation = GDSC.CL_annotation
  )

  names(guideSets)[[16]] <- "MSigDB CFEs"

  gs <- guideSets
  Xmin <- min(
    min(correctedFCs[, "avgFC"]),
    min(correctedFCs[, "correctedFC"])
  ) - 1
  Xmax <- max(
    max(correctedFCs[, "avgFC"]),
    max(correctedFCs[, "correctedFC"])
  ) + 1

  whatToPlot <- c("avgFC", "correctedFC")

  geneSetToPlot <- list(
    AMPLIFIEDpn = c(paste0("Amp (PNs >= ", minCNs, ")")),
    AMPLIFIEDgistic = c("Amp (Gistic +1)", "Amp (Gistic +2)"),
    AMPLIFIEDnotExp = c(
      "Amp (Gistic +2) notExp",
      paste0("Amp (PNs >= ", minCNs, ") notExp")
    ),
    REFERENCE = c("MSigDB CFEs", "BAGEL Essential", "BAGEL nonEssential")
  )

  COLS_list <- list(
    AMPLIFIEDpn = c("#FF243C", "#FF0000", "darkgrey"),
    AMPLIFIEDgistic = c("red2", "red4", "darkgrey"),
    AMPLIFIEDnotExp = c("blue", "#2F4CAC", "#00008B", "darkgray"),
    REFERENCE = c("darkgreen", "green", "bisque4")
  )

  withr::with_options(
    list(warn = -1), {
      for (i in seq_len(4)) {
        withr::with_par(
          list(mfrow = c(2, 1)), {
            for (j in seq_len(2)) {

              plotType <- geneSetToPlot[[i]]

              densities <- list()
              for (k in seq_along(plotType)) {
                if (length(gs[[plotType[k]]] > 0)) {
                  densities[[k]] <- density(
                    correctedFCs[gs[[plotType[k]]], whatToPlot[j]],
                    na.rm = TRUE
                  )
                }else{
                  densities[[k]] <- list(x = NA, y = NA)
                }
              }

              densOthers <- density(
                correctedFCs[
                  setdiff(rownames(correctedFCs), unlist(gs[plotType])),
                  whatToPlot[j]
                ],
                na.rm = TRUE
              )

              COLS <- COLS_list[[i]]
              toPlot <- append(densities, list(densOthers))

              if (j == 1) {
                xlab <- ""
                title <- "pre CRISPRcleanR"
              }else{
                xlab <- "sgRNA log FC"
                title <- "post CRISPRcleanR"
              }

              withr::with_par(
                list(mar = c(4, 4, 2, 1)), {
                  withr::with_options(
                    list(warn = -1), {
                      ccr.multDensPlot(
                        TOPLOT = toPlot,
                        COLS = COLS,
                        XLIMS = c(Xmin, Xmax),
                        TITLE = title,
                        LEGentries = c(plotType, "others"),
                        XLAB = xlab
                      )
                    }
                  )
                }
              )
            }
          },
          no.readonly = TRUE
        )
      }
    }
  )
}

ccr.RecallCurves <- function(
  cellLine,
  correctedFCs,
  GDSC.geneLevCNA = NULL,
  RNAseq.fpkms = NULL,
  minCN = 8,
  libraryAnnotation = NULL,
  GeneLev = FALSE,
  GDSC.CL_annotation = NULL,
  verbose = 1
) {

  # Start update for web version
  guideSets <- ccr.get.guideSets(
    cellLine,
    GDSC.geneLevCNA = GDSC.geneLevCNA,
    CCLE.gisticCNA = NULL,
    RNAseq.fpkms = RNAseq.fpkms,
    libraryAnnotation = libraryAnnotation,
    GDSC.CL_annotation = GDSC.CL_annotation
  )

  toPlot <- c("avgFC", "correctedFC")
  if (verbose == 0) {
    pb <- txtProgressBar(min = 0, max = 4, style = 3)
  }
  # End update for web version

  if (GeneLev) {
    guideSets <- lapply(guideSets, function(x) {
      libraryAnnotation[x, "GENES"]
    })
  }

  COLS <- c("bisque4", "green", "red", "blue")

  withr::with_par(
    list(mfrow = c(2, 1), mar = c(4, 4, 4, 1)), {
      AUCcombo <- NULL
      for (j in seq_along(toPlot)) {
        if (j == 1) {
          MAIN <- paste(cellLine, "pre-CRISPRcleanR")
        }else{
          MAIN <- paste(cellLine, "post-CRISPRcleanR")
        }

        if (!GeneLev) {
          O1 <- order(correctedFCs[, toPlot[j]])
          predictions1 <- rownames(correctedFCs)[O1]
          XLAB <- "sgRNA logFC percentile"
        }else{
          sgProf <- correctedFCs[, toPlot[j]]
          names(sgProf) <- rownames(correctedFCs)
          gProf <- ccr.geneMeanFCs(sgProf, libraryAnnotation)
          O1 <- order(gProf)
          predictions1 <- names(gProf)[O1]
          XLAB <- "gene Avg logFC percentile"
        }

        toTest <- c(
          "BAGEL nonEssential",
          "BAGEL Essential",
          paste0("Amp (PNs >= ", minCN, ")"),
          paste0("Amp (PNs >= ", minCN, ") notExp")
        )

        AUCs <- vector()
        for (i in seq_along(toTest)) {
          currentSet <- guideSets[[toTest[i]]]
          currentSet <- intersect(currentSet, predictions1)
          pp1 <- cumsum(
            is.element(predictions1, currentSet)
          ) / length(currentSet)

          if (i == 1) {
            plot(
              100 * (seq_along(predictions1)) / length(predictions1),
              100 * pp1,
              type = "l",
              xlim = c(0, 100),
              ylim = c(0, 100),
              ylab = "% Recall",
              xlab = XLAB,
              col = COLS[i],
              main = MAIN,
              lwd = 2
            )
          }else{
            lines(
              100 * (seq_along(predictions1)) / length(predictions1),
              100 * pp1,
              type = "l",
              xlim = c(0, 100),
              ylim = c(0, 100),
              col = COLS[i],
              lwd = 2
            )
          }

          AUCs[i] <- trapz(seq_along(predictions1) / length(predictions1), pp1)

          # Start update for web version

          if (verbose == 0) {
            setTxtProgressBar(pb, i)
          }
        }

        AUCcombo <- rbind(AUCcombo, AUCs)
        if (j == 1) {
          legend("bottomright", toTest, cex = 0.8, lwd = 2, col = COLS)
        }
      }
    },
    no.readonly = TRUE
  )

  if (verbose == 0) {
    close(pb)
  }
  # End update for web version

  rownames(AUCcombo) <- c("pre-CRISPRcleanR", "post-CRISPRcleanR")
  colnames(AUCcombo) <- toTest

  return(AUCcombo)
}

ccr.impactOnPhenotype <- function(
  MO_uncorrectedFile,
  MO_correctedFile,
  sigFDR = 0.05,
  expName = "expName",
  display = TRUE
) {

  pre <- read.table(
    MO_uncorrectedFile,
    sep = "\t",
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
  )
  pre <- pre[order(rownames(pre)), ]

  post <- read.table(
    MO_correctedFile,
    sep = "\t",
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
  )
  post <- post[order(rownames(post)), ]

  preD <- which(pre[["neg.fdr"]] < sigFDR & pre[["pos.fdr"]] >= sigFDR)
  preE <- which(pre[["neg.fdr"]] >= sigFDR & pre[["pos.fdr"]] < sigFDR)
  preNULL <- setdiff(seq_len(nrow(pre)), c(preD, preE))

  postD <- which(post[["neg.fdr"]] < sigFDR & post[["pos.fdr"]] >= sigFDR)
  postE <- which(post[["neg.fdr"]] >= sigFDR & post[["pos.fdr"]] < sigFDR)
  postNULL <- setdiff(seq_len(nrow(post)), c(postD, postE))

  aDD <- length(intersect(preD, postD))
  aDN <- length(intersect(preD, postNULL))
  aDE <- length(intersect(preD, postE))

  aND <- length(intersect(preNULL, postD))
  aNN <- length(intersect(preNULL, postNULL))
  aNE <- length(intersect(preNULL, postE))

  aED <- length(intersect(preE, postD))
  aEN <- length(intersect(preE, postNULL))
  aEE <- length(intersect(preE, postE))

  cm <- matrix(
    c(aDD, aDN, aDE, aND, aNN, aNE, aED, aEN, aEE),
    3,
    3,
    dimnames = list(c("cD", "cN", "cE"), c("uD", "uN", "uE"))
  )
  cm[is.na(cm)] <- 0

  IMPACTEDg <- 100 * sum(triu(cm, 1) + tril(cm, -1)) / sum(c(cm))
  IMPACTED_phenGenes <- 100 * (
    cm[2, 1] + cm[2, 3] + cm[3, 1] + cm[3, 2]) / sum(c(cm[, c(1, 3)])
  )

  IMPACTED_Depletions <- 100 * (cm[2, 1] + cm[2, 3]) / sum(cm[, 1])
  IMPACTED_Enrichments <- 100 * (cm[2, 3] + cm[1, 3]) / sum(cm[, 3])

  DISTORTEDg <- 100 * (cm[1, 3] + cm[3, 1]) / sum(c(cm))
  DISTORTED_phenGenes <- 100 * (cm[1, 3] + cm[3, 1]) / sum(c(cm[, c(1, 3)]))

  DISTORTED_Depletions <- 100 * cm[3, 1] / sum(cm[, 1])
  DISTORTED_Enrichments <- 100 * cm[1, 3] / sum(cm[, 3])

  geneCounts <- cm

  colnames(cm) <- paste(
    colSums(cm),
    c("loss of fitness", "no phenotype", "gain of fitness"),
    sep = "\n"
  )
  cm <- cm / t(matrix(rep(colSums(cm), nrow(cm)), 3, 3))

  if (display) {
    withr::with_par(
      list(mar = c(5, 4, 4, 10), xpd = TRUE), {
        barplot(
          100 * cm,
          col = c("red", "gray", "blue"),
          border = FALSE,
          main = expName,
          ylab = "%",
          xlab = "original counts"
        )
        legend(
          "right",
          c("loss of fitness", "no phenotype", "gain of fitness"),
          inset = c(-.5, 0),
          title = "Corrected counts",
          fill = c("red", "gray", "blue"),
          border = NA
        )
      },
      no.readonly = TRUE
    )

    withr::with_par(
      list(mfrow = c(2, 2), mar = c(0, 0, 2, 0), xpd = TRUE), {
        pie(
          c(IMPACTEDg, 100 - IMPACTEDg),
          col = c("blue", "white"),
          border = "gray",
          labels = c(paste0(format(IMPACTEDg, digits = 4), "%"), ""),
          main = "Overall impact"
        )
        pie(
          c(DISTORTEDg, 100 - DISTORTEDg),
          col = c("blue", "white"),
          border = "gray",
          labels = c(paste0(format(DISTORTEDg, digits = 4), "%"), ""),
          main = "Overall distortion"
        )
        pie(
          c(IMPACTED_phenGenes, 100 - IMPACTED_phenGenes),
          col = c("darkgreen", "white"),
          border = "gray",
          labels = c(paste0(
            format(IMPACTED_phenGenes, digits = 4), "%"
          ), ""),
          main = "Impact (G/L fitness genes)"
        )
        pie(
          c(DISTORTED_phenGenes, 100 - DISTORTED_phenGenes),
          col = c("darkgreen", "white"),
          border = "gray",
          labels = c(paste0(
            format(DISTORTED_phenGenes, digits = 4), "%"
          ), ""),
          main = "Distortion (G/L fitness genes)"
        )
      },
      no.readonly = TRUE
    )
  }

  dimnames(geneCounts) <- list(
    `corrected counts` = c("dep.", "null", "enr."),
    `original counts` = c("dep.", "null", "enr.")
  )

  id <- intersect(preD, postE)
  to_bind <- cbind(
    pre[id, c("neg.fdr", "pos.fdr")],
    post[id, c("neg.fdr", "pos.fdr")]
  )

  id <- intersect(preE, postD)
  to_bind <- rbind(
    to_bind,
    cbind(
      pre[id, c("neg.fdr", "pos.fdr")],
      post[id, c("neg.fdr", "pos.fdr")]
    )
  )
  colnames(to_bind) <- paste0(
    c("", "", "ccr.", "ccr."),
    colnames(to_bind)
  )

  id <- intersect(preD, postD)
  to_bind_c <- cbind(
    pre[id, c("neg.fdr", "pos.fdr")],
    post[id, c("neg.fdr", "pos.fdr")]
  )

  id <- intersect(preE, postE)
  to_bind_c <- rbind(
    to_bind_c,
    cbind(
      pre[id, c("neg.fdr", "pos.fdr")],
      post[id, c("neg.fdr", "pos.fdr")]
    )
  )
  colnames(to_bind_c) <- paste0(
    c("", "", "ccr.", "ccr."),
    colnames(to_bind_c)
  )

  id <- intersect(preD, postNULL)
  to_bind_A <- cbind(
    pre[id, c("neg.fdr", "pos.fdr")],
    post[id, c("neg.fdr", "pos.fdr")]
  )

  id <- intersect(preE, postNULL)
  to_bind_A <- rbind(
    to_bind_A,
    cbind(
      pre[id, c("neg.fdr", "pos.fdr")],
      post[id, c("neg.fdr", "pos.fdr")]
    )
  )
  colnames(to_bind) <- paste0(
    c("", "", "ccr.", "ccr."),
    colnames(to_bind)
  )

  return(list(
    `GW_impact %` = IMPACTEDg,
    `Phenotype_G_impact %` = IMPACTED_phenGenes,
    `Depleted_G_impact %` = IMPACTED_Depletions,
    `Enriched_G_impact %` = IMPACTED_Enrichments,
    `GW_distortion %` = DISTORTEDg,
    `Phenotype_G_distortion %` = DISTORTED_phenGenes,
    `Depleted_G_distortion %` = DISTORTED_Depletions,
    `Enriched_G_distortion %` = DISTORTED_Enrichments,
    geneCounts = geneCounts,
    distortion = to_bind,
    impact = to_bind_A
  ))
}

ccr.geneSummary <- function(
  sgRNA_FCprofile,
  libraryAnnotation,
  FDRth = 0.05
) {

    geneLevFCs <- ccr.geneMeanFCs(sgRNA_FCprofile, libraryAnnotation)

    data(BAGEL_essential)
    data(BAGEL_nonEssential)

    ROCresults <- ccr.PrRc_Curve(
      FCsprofile = geneLevFCs,
      positives = BAGEL_essential,
      negative = BAGEL_nonEssential,
      display = FALSE,
      FDRth = FDRth
    )

    SigDepletedVector <- (geneLevFCs < ROCresults[["sigthreshold"]])

    oo <- order(geneLevFCs)
    geneLevFCs <- geneLevFCs[oo]
    SigDepletedVector <- SigDepletedVector[oo]

    res <- data.frame(
      stringsAsFactors = FALSE,
      row.names = names(geneLevFCs),
      logFC = geneLevFCs,
      "SigDep" = SigDepletedVector
    )

    return(res)
}

## other exported non documented functions

#### Assessment and visualisation

### Utils

## not exported functions
ccr.boxplot <- function(
  toPlot,
  main,
  names
) {

  boxplot(toPlot, main = main, names = names, las = 2)
  MEANS <- apply(toPlot, MARGIN = 2, FUN = "mean")
  SD <- apply(toPlot, MARGIN = 2, FUN = "sd")
  for (i in seq_along(MEANS)) {
    lines(
      x = c(i - 0.3, i + 0.3),
      y = c(MEANS[i], MEANS[i]) + SD[i],
      col = "blue",
      lwd = 2
    )
    lines(
      x = c(i - 0.2, i + 0.2),
      y = c(MEANS[i], MEANS[i]),
      col = "red",
      lwd = 5
    )
    lines(
      x = c(i - 0.3, i + 0.3),
      y = c(MEANS[i], MEANS[i]) - SD[i],
      col = "blue",
      lwd = 2
    )
  }
}

ccr.cohens_d <- function(
  x,
  y
) {
  lx <- length(x) - 1
  ly <- length(y) - 1

  md  <- abs(mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE))
  csd <- lx * var(x, na.rm = TRUE) + ly * var(y, na.rm = TRUE)

  csd <- csd / (lx + ly)
  csd <- sqrt(csd)                     ## common sd computation

  cd  <- md / csd                        ## cohen"s d
  return(cd)
}

ccr.get.guideSets <- function(
  cellLine,
  GDSC.geneLevCNA = NULL,
  CCLE.gisticCNA = NULL,
  RNAseq.fpkms = NULL,
  libraryAnnotation = NULL,
  GDSC.CL_annotation = NULL
) {
  if (is.null(GDSC.geneLevCNA)) {
    data(GDSC.geneLevCNA, envir = environment())
  }
  if (is.null(CCLE.gisticCNA)) {
    data(CCLE.gisticCNA, envir = environment())
  }
  if (is.null(RNAseq.fpkms)) {
    data(RNAseq.fpkms, envir = environment())
  }

  GDSC.cna <- ccr.get.gdsc1000.AMPgenes(
    cellLine = cellLine,
    GDSC.geneLevCNA = GDSC.geneLevCNA,
    minCN = 0,
    GDSC.CL_annotation = GDSC.CL_annotation
  )
  CCLE.cna <- ccr.get.CCLEgisticSets(
    cellLine = cellLine,
    CCLE.gisticCNA = CCLE.gisticCNA,
    GDSC.CL_annotation = GDSC.CL_annotation
  )
  notExp <- ccr.get.nonExpGenes(
    cellLine = cellLine,
    RNAseq.fpkms = RNAseq.fpkms,
    GDSC.CL_annotation = GDSC.CL_annotation
  )

  data(BAGEL_essential, envir = environment())
  data(BAGEL_nonEssential, envir = environment())

  data(EssGenes.DNA_REPLICATION_cons, envir = environment())
  data(EssGenes.KEGG_rna_polymerase, envir = environment())
  data(EssGenes.PROTEASOME_cons, envir = environment())
  data(EssGenes.ribosomalProteins, envir = environment())
  data(EssGenes.SPLICEOSOME_cons, envir = environment())

  CFEgenes <- unique(c(
    EssGenes.DNA_REPLICATION_cons,
    EssGenes.KEGG_rna_polymerase,
    EssGenes.PROTEASOME_cons,
    EssGenes.ribosomalProteins,
    EssGenes.SPLICEOSOME_cons
  ))

  BAGEL_essOnly <- setdiff(BAGEL_essential, CFEgenes)

  withr::with_options(
    list(warn = -1), {
      BAGEL_essential <- ccr.genes2sgRNAs(
        BAGEL_essential,
        libraryAnnotation = libraryAnnotation
      )
      BAGEL_nonEssential <- ccr.genes2sgRNAs(
        BAGEL_nonEssential,
        libraryAnnotation = libraryAnnotation
      )
      EssGenes.DNA_REPLICATION_cons <- ccr.genes2sgRNAs(
        EssGenes.DNA_REPLICATION_cons,
        libraryAnnotation = libraryAnnotation
      )
      EssGenes.KEGG_rna_polymerase <- ccr.genes2sgRNAs(
        EssGenes.KEGG_rna_polymerase,
        libraryAnnotation = libraryAnnotation
      )
      EssGenes.PROTEASOME_cons <- ccr.genes2sgRNAs(
        EssGenes.PROTEASOME_cons,
        libraryAnnotation = libraryAnnotation
      )
      EssGenes.ribosomalProteins <- ccr.genes2sgRNAs(
        EssGenes.ribosomalProteins,
        libraryAnnotation = libraryAnnotation
      )
      EssGenes.SPLICEOSOME_cons <- ccr.genes2sgRNAs(
        EssGenes.SPLICEOSOME_cons,
        libraryAnnotation = libraryAnnotation
      )
      CFEgenes <- ccr.genes2sgRNAs(
        CFEgenes,
        libraryAnnotation = libraryAnnotation
      )
      BAGEL_essOnly <- ccr.genes2sgRNAs(
        BAGEL_essOnly,
        libraryAnnotation = libraryAnnotation
      )
      DelGDSC <- ccr.genes2sgRNAs(
        libraryAnnotation,
        GDSC.cna[["Gene"]][GDSC.cna[["CN"]] == 0]
      )
      DelCCLE <- ccr.genes2sgRNAs(
        libraryAnnotation,
        CCLE.cna[["gm2"]]
      )
      notExpG <- ccr.genes2sgRNAs(
        libraryAnnotation,
        notExp
      )
      gistic1 <- ccr.genes2sgRNAs(
        libraryAnnotation,
        CCLE.cna[["gp1"]]
      )
      gistic2 <- ccr.genes2sgRNAs(
        libraryAnnotation,
        CCLE.cna[["gp2"]]
      )
      gistic1ne <- ccr.genes2sgRNAs(
        libraryAnnotation,
        intersect(CCLE.cna[["gp1"]], notExp)
      )
      gistic2ne <- ccr.genes2sgRNAs(
        libraryAnnotation,
        intersect(CCLE.cna[["gp2"]], notExp)
      )
      pn2 <- ccr.genes2sgRNAs(
        libraryAnnotation,
        GDSC.cna[["Gene"]][GDSC.cna[["CN"]] >= 2]
      )
      pn4 <- ccr.genes2sgRNAs(
        libraryAnnotation,
        GDSC.cna[["Gene"]][GDSC.cna[["CN"]] >= 4]
      )
      pn8 <- ccr.genes2sgRNAs(
        libraryAnnotation,
        GDSC.cna[["Gene"]][GDSC.cna[["CN"]] >= 8]
      )
      pn10 <- ccr.genes2sgRNAs(
        libraryAnnotation,
        GDSC.cna[["Gene"]][GDSC.cna[["CN"]] >= 10]
      )
      pn2ne <- ccr.genes2sgRNAs(
        libraryAnnotation,
        intersect(GDSC.cna[["Gene"]][GDSC.cna[["CN"]] >= 2], notExp))
      pn4ne <- ccr.genes2sgRNAs(
        libraryAnnotation,
        intersect(GDSC.cna[["Gene"]][GDSC.cna[["CN"]] >= 4], notExp))
      pn8ne <- ccr.genes2sgRNAs(
        libraryAnnotation,
        intersect(GDSC.cna[["Gene"]][GDSC.cna[["CN"]] >= 8], notExp)
      )
      pn10ne <- ccr.genes2sgRNAs(
        libraryAnnotation,
        intersect(GDSC.cna[["Gene"]][GDSC.cna[["CN"]] >= 10], notExp)
      )
    }
  )

  guideSets <- list(
    DelGDSC,
    DelCCLE,
    notExpG,
    gistic1,
    gistic2,
    gistic1ne,
    gistic2ne,
    pn2,
    pn4,
    pn8,
    pn10,
    pn2ne,
    pn4ne,
    pn8ne,
    pn10ne,
    CFEgenes,
    BAGEL_essential,
    BAGEL_essOnly,
    BAGEL_nonEssential
  )

  nnames <- c(
    "Dep (PN)",
    "Dep (Gistic)",
    "notExp",
    "Amp (Gistic +1)",
    "Amp (Gistic +2)",
    "Amp (Gistic +1) notExp",
    "Amp (Gistic +2) notExp",
    "Amp (PNs >= 2)",
    "Amp (PNs >= 4)",
    "Amp (PNs >= 8)",
    "Amp (PNs >= 10)",
    "Amp (PNs >= 2) notExp",
    "Amp (PNs >= 4) notExp",
    "Amp (PNs >= 8) notExp",
    "Amp (PNs >= 10) notExp",
    "Essential",
    "BAGEL Essential",
    "BAGEL Ess.Only",
    "BAGEL nonEssential"
  )
  names(guideSets) <- nnames
  return(guideSets)
}

ccr.fixedFDRthreshold <- function(
  FCsprofile,
  TruePositives,
  TrueNegatives,
  th
) {
  presentGenes <- intersect(c(TruePositives, TrueNegatives), names(FCsprofile))
  predictions <- FCsprofile[presentGenes]
  observations <- is.element(presentGenes, TruePositives) + 0
  names(observations) <- presentGenes
  RES <- roc(observations, predictions, direction = ">", quiet = TRUE)
  COORS <- coords(RES, "all", ret = c("threshold", "ppv"), transpose = TRUE)
  FDRpercTh <- max(COORS["threshold", which(COORS["ppv", ] >= (1 - th))])
  FDRpercRANK <- max(which(sort(FCsprofile) <= FDRpercTh))
  return(list(FCth = FDRpercTh, RANK = FDRpercRANK))
}

ccr.sd_guideFCs <- function(
  FCs,
  distorted_genes,
  libraryAnnotation
) {
  SDS <- aggregate(FCs, by = list(libraryAnnotation[names(FCs), "GENES"]), "sd")
  XSDS <- SDS[["x"]]
  names(XSDS) <- SDS[["Group.1"]]
  boxplot(XSDS[setdiff(names(XSDS), distorted_genes)], XSDS[distorted_genes])
  print(
    t.test(XSDS[setdiff(names(XSDS), distorted_genes)], XSDS[distorted_genes])
  )
}

# Start update for web version

ccr.dfToFile <- function(
  df,
  outfile,
  data_format
) {
  if (!file.exists(dirname(outfile))) {
    dir.create(dirname(outfile), recursive = TRUE)
  }
  if (!toupper(data_format) %in% c("TXT", "TSV", "CSV", "JSON")) {
    warning("Export format not recongnized! Exports switched to TSV.")
    data_format <- "TSV"
  } else {
    if (toupper(data_format) == "JSON") {
      if (!any(rownames(installed.packages()) == "jsonlite")) {
        warning("Missing package [jsonlite]! No Export availalbe for web app.")
        data_format <- "MISSING"
      }
    }
  }
  if (toupper(data_format) == "JSON") {
    write(
      jsonlite::toJSON(df, pretty = TRUE),
      paste0(outfile, ".json")
    )
  }
  if (toupper(data_format) %in% c("TXT", "TSV", "CSV")) {
    write.table(
      df,
      quote = toupper(data_format) == "CSV",
      col.names = TRUE,
      row.names = any(is.na(as.integer(rownames(df)))),
      sep = ifelse(toupper(data_format) == "CSV", ",", "\t"),
      file = paste0(outfile, ".", tolower(data_format))
    )
  }
  return(data_format)
}

ccr.countToDist <- function(
  counts,
  ncontrols
) {
    counts_dist <- NULL
    sample_n <- 0
    for (sample in colnames(counts)[- c(1:2)]) {
      sample_n <- sample_n + 1
      sample_dist_params <- c(
        mean =  mean(counts[, sample], na.rm = TRUE),
        sd = sd(counts[, sample], na.rm = TRUE),
        quantile(counts[, sample], c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
      )
      tmp <- counts[
        counts[, sample] > sample_dist_params[["75%"]] + (
          sample_dist_params[["75%"]] -
          sample_dist_params[["25%"]]
        ) * 1.5 |
        counts[, sample] < sample_dist_params[["25%"]] - (
          sample_dist_params[["75%"]] -
          sample_dist_params[["25%"]]
        ) * 1.5,
        c(colnames(counts)[1:2], sample)
      ]
      colnames(tmp) <- c("sgRNA", "gene", "value")
      counts_dist[[sample_n]] <- list(
          info = data.frame(
          name = sample,
          type = ifelse(sample_n <= ncontrols, "control", "sample"),
          order = sample_n,
          stringsAsFactors = FALSE
        ),
        dist = data.frame(
          W1 = max(
            sample_dist_params[["0%"]], sample_dist_params[["25%"]] - (
              sample_dist_params[["75%"]] -
              sample_dist_params[["25%"]]
            ) * 1.5
          ),
          Q1 = sample_dist_params[["25%"]],
          median = sample_dist_params[["50%"]],
          Q3 = sample_dist_params[["75%"]],
        W2 = min(
          sample_dist_params[["100%"]],
          sample_dist_params[["75%"]] + (
            sample_dist_params[["75%"]] -
            sample_dist_params[["25%"]]
          ) * 1.5
        ),
        sd1 = sample_dist_params[["mean"]] - sample_dist_params[["sd"]],
        mean = sample_dist_params[["mean"]],
        sd2 = sample_dist_params[["mean"]] + sample_dist_params[["sd"]],
        stringsAsFactors = FALSE
      ),
      outliers = tmp
    )
  }
  return(counts_dist)
}

# End update for web version
