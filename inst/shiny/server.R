################################################################################
################################################################################
library("shiny")
library("shinyFiles")
library("DT")
library("data.table")
library("ggplot2")
theme_set(theme_bw())
library("magrittr")
################################################################################
shinyServer(function(input, output, session){
  ########################################
  # Hard-coded params
  ########################################
  # The user-exposed options for searching. Might be OS brittle.
  volumes <- c(
    Downloads = "~/Downloads/",
    Home = "~/",
    Root = "/")
  # Just an example of a restricted path. Not immediately useful.
  restrictedPaths = system.file(package='base')
  ########################################
  # Inputs and reactives and logic and stuff
  ########################################
  
  # input folder selection for OG
  shinyDirChoose(input, 'dirOG', roots=volumes, session=session, restrictions = restrictedPaths)
  inputDirPath = reactive({
    validate(need(input$dirOG, message = "..."))
    inputDirPath = parseDirPath(volumes, input$dirOG)
    if(!file.exists(inputDirPath)){
      inputDirPath <- ""
    }
    return(inputDirPath)
  })
  infoFilePath = reactive({
    infoFilePath = file.path(inputDirPath(), "info.txt")
    if(!file.exists(infoFilePath)){
      infoFilePath <- "NO `info.txt` file found!"
    }
    return(infoFilePath)
  })
  
  # input folder selection for learnErrors, which uses FT
  shinyDirChoose(input, 'dirFT', roots=volumes, session=session, restrictions = restrictedPaths)
  inputDirPathFT = reactive({
    validate(need(input$dirFT, message = "..."))
    inputDirPathFT = parseDirPath(volumes, input$dirFT)
    if(!file.exists(inputDirPathFT)){
      inputDirPathFT <- ""
    }
    return(inputDirPathFT)
  })
  infoFilePathFT = reactive({
    infoFilePathFT = file.path(inputDirPathFT(), "info.txt")
    if(!file.exists(infoFilePathFT)){
      infoFilePathFT <- "NO `info.txt` file found!"
    }
    return(infoFilePathFT)
  })
  fastqFilesTabFT = reactive({
    fread(input = infoFilePathFT())
  })
  # The OG fastq files tab
  fastqFilesTab = reactive({
    fread(input = infoFilePath())
  })
  # The samples
  includeSamples = reactive({
    validate(need(fastqFilesTab(), message = "..."))
    validate(need(input$NSamples, message = "Need number of files for downsampled quality evaluation."))
    fastqFilesTab <- copy(fastqFilesTab())
    # Define a subset of sample names
    NSamples = input$NSamples
    # The minimum of either the prescribed number of samples and the number available
    NSamples <- c(NSamples, fastqFilesTab[, uniqueN(Sample)]) %>% min(na.rm = TRUE)
    includeSamples = fastqFilesTab[, unique(Sample)] %>% 
      sample(size = NSamples, replace = FALSE) %>% 
      sort
    return(includeSamples)
  })
  
  # input folder selection for Run DADA, which uses FT
  shinyDirChoose(input, 'dirRD', roots=volumes, session=session, restrictions = restrictedPaths)
  inputDirPathRD = reactive({
    validate(need(input$dirRD, message = "..."))
    inputDirPathRD = parseDirPath(volumes, input$dirRD)
    if(!file.exists(inputDirPathRD)){
      inputDirPathRD <- ""
    }
    return(inputDirPathRD)
  })
  infoFilePathRD = reactive({
    infoFilePathRD = file.path(inputDirPathRD(), "info.txt")
    if(!file.exists(infoFilePathRD)){
      infoFilePathRD <- "NO `info.txt` file found!"
    }
    return(infoFilePathRD)
  })
  
  ########################################
  # Quality Sub-Sampling 
  # Tabulate q-values into a long granular table for each file.
  ########################################
  TabulateQuality = reactive({
    # How many reads, at most, to read from each file.
    # validate(need(input$nreads, message = "..."))
    input$nreads %>% need(message = "...") %>% validate
    nReads = 10^(input$nreads)
    # Don't move forward if no sequence table yet.
    validate(need(fastqFilesTab(), message = "Select Info File."))
    fastqFilesTab <- copy(fastqFilesTab())
    includeSamples = includeSamples()
    # Progress bar...
    # progress <- shiny::Progress$new(session, min=1, max=15)
    # on.exit(progress$close())
    withProgress(
      message = 'Sampling quality profiles from sequence files...',
      detail = 'This may take a while...',
      value = 0, 
      expr = {
        incProgUnit = 1 / (fastqFilesTab[, 1, by = c("File", "Direction")] %>% nrow)
        qtab <- fastqFilesTab[
          (Sample %chin% includeSamples),
          {
            message("Tabulating from:\n", File)
            incProgress(
              amount = incProgUnit,
              message = Direction[1],
              detail = File[1])
            tabulate_quality(
              fastqFile = file.path(inputDirPath(), File),
              nReads = nReads)
          }, by = c("File", "Direction")]
      })
    return(qtab)
  })
  maxLength = reactive({
    qtab = TabulateQuality() %>% copy
    return(
      qtab[, max(Cycle)]
    )
  })
  ########################################
  # Summarize Sub-Sampled Quality, reactively
  ########################################
  QualSummReact = reactive({
    validate(need(input$qquantile, message = "..."))
    # Quantile to show
    desiredQuantile = input$qquantile
    qtab = TabulateQuality() %>% copy
    # `CycleCounts` is collapsed for distributional summary at each cycle
    CycleCounts = qtab[, .(Count = sum(Count, na.rm = TRUE)),
                       by = c("Direction", "Cycle", "Score")]
    setorderv(CycleCounts, c("Cycle", "Score"), order = c(1, 1))
    # Define quantile
    CycleCounts[, Quantile := cumsum(Count)/sum(Count, na.rm = TRUE),
                by = c("Direction", "Cycle")]
    # Define proportion
    CycleCounts[, Proportion := (Count)/sum(Count, na.rm = TRUE),
                by = c("Direction", "Cycle")]
    # Collect summary statistic
    cycsum = CycleCounts[,
                         .(
                           mean = sum(Score * Count, na.rm = TRUE) / sum(Count, na.rm = TRUE),
                           median = Score[(Quantile >= 0.50)][1],
                           QuantN = Score[(Quantile >= desiredQuantile)][1]
                         ),
                         by = c("Direction", "Cycle")]
    nameQuantN = paste0("Quantile_", desiredQuantile)
    setnames(cycsum, "QuantN", nameQuantN)
    suppressWarnings({
      CycleStats <- melt.data.table(cycsum,
                                    id.vars = c("Direction", "Cycle"),
                                    variable.name = "Statistic",
                                    value.name = "Quality")
    })
    fitTab = CycleStats[, .(fit = list(fit = loess(Quality ~ Cycle, data = .SD))),
                        by = c("Direction", "Statistic")]
    # Define smoothed values at each entry
    # Form: A[B, bb:=i.b, on='a']
    # CycleStats[fitTab, Smooth := predict(i.fit[[1]], newdata = Cycle), on = c("Direction", "Statistic")]
    CycleStats[fitTab, Fit := i.fit, on = c("Direction", "Statistic")]
    CycleStats[, Smooth := predict(Fit[[1]], newdata = Cycle), by = c("Direction", "Statistic", "Cycle")]
    CycleStats[, Fit := NULL]
    return(
      list(
        CycleStats = CycleStats,
        CycleCounts = CycleCounts
      )
    )
  })
  # Compute default/predicted trim values
  suggestedTrimTable = reactive({
    # Min-quality at the indicated quantile (smooth)
    # for computing the default trimming parameters
    # e.g. 15L
    minQual = input$minQual
    LeftTrimDefault = 10L
    CycleStats = QualSummReact()$CycleStats %>% copy
    # Use this to define default Right Trim
    RightTrimTable = CycleStats[, .(Cycle = min(max(Cycle, na.rm = TRUE), Cycle[(Smooth <= minQual)], na.rm = TRUE)), by = "Direction"]
    RightTrimTable[, Side := "Right"]
    # Define default Left Trim
    LeftTrimTable = copy(RightTrimTable)[, c("Side", "Cycle") := list("Left", LeftTrimDefault)]
    TrimTable = list(RightTrimTable, LeftTrimTable) %>% rbindlist
    return(TrimTable)
  })
  output$uiForward <- renderUI({
    ForwardRight = suggestedTrimTable()[(Direction == "F" & Side == "Right")]$Cycle
    return(
      sliderInput("Forward", "Forward", min = 1, max = maxLength(), step = 5, value = c(10, ForwardRight))
    )
  })
  output$uiReverse <- renderUI({
    ReverseRight = suggestedTrimTable()[(Direction == "R" & Side == "Right")]$Cycle
    return(
      sliderInput("Reverse", "Reverse", min = 1, max = maxLength(), step = 5, value = c(10, ReverseRight))
    )
  })
  # Populate the trim table to be used in the chart and in the executed filtertrimming
  # from the user-input, the dual-slider for each direction.
  TrimTable = reactive({
    data.table(
      Direction = c("F", "F", "R", "R"), 
      Side = c("Left", "Right", "Left", "Right"),
      Cycle = c(input$Forward[1],
                input$Forward[2],
                input$Reverse[1],
                input$Reverse[2])
    )
  })
  ########################################
  # FilterTrim Quality Graphic
  ########################################
  # Define the main graphic as a ggplot2 object.
  pQbyC = reactive({
    qualitySummaryList = QualSummReact()
    TrimTable = copy(TrimTable())
    # Modify TrimTable manually/hardcoded as-needed
    # to avoid the weird "cliff" in this particular dataset
    # TrimTable[(Direction == "R1" & Side == "Right"), Cycle := 252]
    pQbyC = plot_quality_by_cycle(CycleStats = qualitySummaryList$CycleStats,
                                  CycleCounts = qualitySummaryList$CycleCounts,
                                  TrimTable = TrimTable)
    return(pQbyC)
  })
  ########################################
  # FilterTrim Execution
  ########################################
  # The filtertrim directory
  ftDir = reactive({
    # Create the filtertrim directory if it doesn't exist yet
    ftDir = "FT" %>% file.path(inputDirPath(), .)
    if(!(ftDir %>% dir.exists)){
      ftDir %>% dir.create()
    }
    return(ftDir)
  })
  
  ftFilesTable = reactive({
    ftDir = ftDir()
    # Define the file paths of the filter-trimmed files
    filterTrimFilesTable = fastqFilesTab() %>% copy
    filterTrimFilesTable[, FileOG := File %>% file.path(inputDirPath(), .) %>% normalizePath]
    filterTrimFilesTable[, FileFT := paste0(Sample, "_", Direction, "-filtrim.fastq.gz") %>% 
                           file.path(ftDir, .)]
    message("FilterTrim Paths:\n\n",
            filterTrimFilesTable$FileFT %>% paste0(collapse = "\n"),
            "\n\n")
    TrimTable = TrimTable() %>% copy
    # The must be ordered such that forward read is always first among pairs of files
    setorderv(filterTrimFilesTable, c("Sample", "Direction"))
    # Before returning, write the info.txt table to the filtertrim dir
    infoTabFT = copy(filterTrimFilesTable)[, .(Sample, Direction, FileFT)]
    infoTabFT[, File := basename(FileFT)]
    infoTabFT[, FileFT := NULL]
    write_table_tab(x = infoTabFT, 
                    file = file.path(ftDir, "info.txt")
    )
    return(filterTrimFilesTable)
  })
  
  prepFT = reactive({
    # should return dummy text, as non-execution workaround
    # to make FT dir and info file
    tab = ftFilesTable()
    return("filtertrim prepared.")
  })
  
  ExecuteFilterTrim = reactive({
    # This stops action unless filtertrim action button has been pressed.
    input$actionb_filtertrim %>% need(message = "...") %>% validate
    message("Filter Trim execution, iteration number:\n",
            input$actionb_filtertrim)
    # The output dir, created if need-be
    ftDir = ftDir()
    # Define the file paths of the filter-trimmed files
    filterTrimFilesTable = ftFilesTable() %>% copy
    # The table that defines the trimming params
    TrimTable = TrimTable() %>% copy
    # The must be ordered such that forward read is always first among pairs of files
    setorderv(filterTrimFilesTable, c("Sample", "Direction"))
    # The must be ordered such that forward read is always first among pairs of files
    setorderv(TrimTable, c("Direction", "Side"))
    (trimLeft <- TrimTable[(Side == "Left")]$Cycle)
    message("trimLeft:\n", trimLeft %>% paste0(collapse = ", "))
    (truncLen <- TrimTable[(Side == "Right")]$Cycle - trimLeft)
    message("truncLen:\n", truncLen %>% paste0(collapse = ", "))
    # write trimming parameters to the top-level folder
    write_table_tab(x = TrimTable, 
                    file = file.path(inputDirPath(), "trimtable.txt")
    )
    # Define the arguments list that will be passed to filterAndTrim
    ftArgsList = list(
      # File I/O
      fwd = filterTrimFilesTable[(Direction == "F")]$FileOG, 
      filt = filterTrimFilesTable[(Direction == "F")]$FileFT,
      rev = filterTrimFilesTable[(Direction == "R")]$FileOG, 
      filt.rev = filterTrimFilesTable[(Direction == "R")]$FileFT,
      multithread = input$multithread,
      # Filtering and trimming params
      maxLen = Inf,
      minLen = 0,
      trimLeft = trimLeft,
      truncLen = truncLen,
      maxEE = input$maxEE,
      minQ = 0,
      truncQ = input$truncQ,
      maxN = 0,
      n = input$n,
      rm.phix = TRUE,
      compress = TRUE,
      verbose = TRUE)
    # Save these arguments prior to running filterAndTrim
    jsonlite::write_json(x = ftArgsList,
                         path = file.path(inputDirPath(),
                                          "filterAndTrim.json"))
    withProgress(message = 'Filtering and trimming in progress...',
                 detail = 'This may take a while...',
                 value = 0, 
                 expr = {
                   # Execute the filter-trimming
                   incProgress(1/5)
                   do.call(what = "filterAndTrim",
                           args = ftArgsList)
                 })
    
    return("FilterTrim Executed!!!")
  })
  
  ########################################
  # Learn Errors Execution
  ########################################
  ExecuteLearnErrors = eventReactive(eventExpr = input$actionb_learnErrors,
                                     valueExpr = {
    # The FT dir, from which to read data
    ftDir = inputDirPathFT()
    
    fastqFilesTab = fastqFilesTabFT() %>% copy
    fastqFilesTab[, FileFull := ftDir %>% normalizePath %>% file.path(., File)]
    fastqFilesTab[, PassFilter := all(file.exists(FileFull)), by = "Sample"]
    fastqFilesTab[, Size := file.size(FileFull)]
    
    fastqFilesTab$FileFull %>% 
      head %>% 
      paste0(collapse = "\n") %>% 
      message("First five full paths for Learn Errors:\n", .)
    
    # Get the user-spec params
    learnSize = input$LE_learnSize
    minSize = 10^(input$LE_minSize)
    nReads = 10^(input$LE_nreads)
    multithread = input$LE_multithread
    
    # Define the samples to learn from (random)
    samplesToLearnFrom = fastqFilesTab[(Size > minSize & PassFilter), unique(Sample)]
    stopifnot(length(samplesToLearnFrom) > 0)
    learnSamples = sample(x = samplesToLearnFrom,
                          replace = FALSE,
                          size = min(learnSize, length(samplesToLearnFrom)))
    time0 = Sys.time()
    setorderv(fastqFilesTab, c("Sample", "Direction"))
    
    filesLearnForward = fastqFilesTab[(Sample %chin% learnSamples & Direction == "F")]$FileFull
    withProgress(
      message = 'Learning (forward) error rates...',
      detail = 'This can take a while...',
      value = 0, 
      expr = {
        incProgress(1/length(filesLearnForward))
        errF <- learnErrors(
          fls = filesLearnForward,
          nreads = nReads,
          multithread = multithread)
      }
    )
    # Learn reverse error rates
    filesLearnReverse = fastqFilesTab[(Sample %chin% learnSamples & Direction == "R")]$FileFull
    withProgress(
      message = 'Learning reverse error rates...',
      detail = 'This can take a while...',
      value = 0, 
      expr = {
        incProgress(1/length(filesLearnReverse))
        errR <- learnErrors(
          fls = filesLearnReverse,
          nreads = nReads,
          multithread = multithread)
      }
    )
    timeDADA2Learn = (time0 - Sys.time())
    message("\n\n Saving error matrices to FT directory...\n\n")
    # Save list of error matrices
    errs = list(forward = errF,
                reverse = errR)
    saveRDS(object = errs,
            file = file.path(ftDir, "errs.RDS"))
    # Return a final message
    msgLearnTime = paste("Learned errors from:\n", ftDir, "\n",
                         "time elapsed:\n", round(timeDADA2Learn, 3))
    message(msgLearnTime)
    return(msgLearnTime)
  })
  
  # Reactive holding the errors
  errs = reactive({
    ftDir = inputDirPathFT()
    errorsFile = file.path(ftDir, "errs.RDS")
    # Validate that file exists in order to show plot.
    file.exists(errorsFile) %>% 
      need(message = paste(errorsFile, " not found...")) %>% 
      validate
    # message the errors file path
    message("Learned Errors being stored at:\n",
            errorsFile)
    # A reactive object holding the learned errors.
    errs <- reactiveFileReader(intervalMillis = 1000,
                               session = session,
                               filePath = errorsFile,
                               readFunc = readRDS)
    # Expect to at least have $forward errors (even if not paired)
    errs()$forward %>% need(message = "$forward missing or malformed...") %>% validate
    return(errs())
  })
  
  pErrors = reactive({
    return(
      list(
        ForwardErrors = plotErrors(errs()$forward),
        ReverseErrors = plotErrors(errs()$reverse)
      )
    )
  })
  
  ########################################
  # DADA Execution
  ## DADA2::dada()
  ########################################
  testRunDada = eventReactive(
    eventExpr = input$actionb_rundada, 
    valueExpr = {
      message("test run dada2..., the button is working.")
    }
  )
  # Interpret relevant inputs to update info table.
  fastqFilesTabRD = reactive({
    fastqFilesTab = fread(input = infoFilePathRD())
    # Minimum file size
    minSize = 10^(input$RD_minSize)
    # The FT dir, from which to read data
    ftDir = inputDirPathRD()
    fastqFilesTab = fastqFilesTab %>% copy
    fastqFilesTab[, FileFull := ftDir %>% normalizePath %>% file.path(., File)]
    fastqFilesTab[, PassFilter := all(file.exists(FileFull)), by = "Sample"]
    fastqFilesTab[, Size := file.size(FileFull)]
    fastqFilesTab[, PassFilter := PassFilter & Size > minSize, by = "FileFull"]
    # Filter the files that don't pass the requirements of existence, and minimum size
    fastqFilesTabDADA2 = copy(fastqFilesTab)[(PassFilter)]
    
    samplesLost = nrow(fastqFilesTab) - nrow(fastqFilesTabDADA2)
    samplesLostPerc = 100 * (samplesLost / nrow(fastqFilesTab)) %>% round(digits = 1)
    # Status messages
    message("\n\nDADA2 start:\n Of the ", nrow(fastqFilesTab), " input files\n",
            samplesLost, " (", samplesLostPerc,
            "%) were lost due to existence-check and size filter.\n\n")
    fastqFilesTabDADA2$FileFull %>% 
      head %>% 
      paste0(collapse = "\n") %>% 
      message("First five full paths for Running DADA:\n", .)
    # Send it along to outputs
    return(fastqFilesTabDADA2)
  })
  # RUN DADA2
  ExecuteDADA2 = eventReactive(
    eventExpr = input$actionb_rundada,
    valueExpr = {
      message(".\n.\nRunning DADA2...\n.\n.")
      ftDir = inputDirPathRD()
      # Fail early if errors file not present
      errorsFile = file.path(ftDir, "errs.RDS")
      file.exists(errorsFile) %>% 
        need(message = "No learned errors file: `errs.RDS`!\nSee `Learn Errors` tab") %>% 
        validate
      errs = readRDS(errorsFile)
      
      # Run-related user-spec params
      nReads = 10^(input$RD_nreads)
      multithread = input$RD_multithread
      
      fastqFilesTabDADA2 = fastqFilesTabRD() %>% copy
      
      stopifnot(nrow(fastqFilesTabDADA2) > 0)
      
      # Set the dada2 cache directory, helpful for debugging
      dadaCacheDir = file.path(ftDir, "dadaCache")
      if(!dir.exists(dadaCacheDir)){
        dir.create(dadaCacheDir)
      }
      # Names of dada cache files
      fastqFilesTabDADA2[, FileDadaCache := file.path(dadaCacheDir,
                                                      paste0(Sample, "-", Direction, "-dada2.RDS"))]
      # RUN DADA via wrapper function.
      time0 = Sys.time()
      # Run on each sample read-pair
      setorderv(fastqFilesTabDADA2, c("Sample", "Direction"))
      setkeyv(fastqFilesTabDADA2, "Sample")
      # The multi-threaded one-sample-at-a-time approach:
      mergeTab = fastqFilesTabDADA2[(PassFilter),
                                    wrap_dada2_workflow(
                                      seqFiles = FileFull,
                                      dadaOutFiles = FileDadaCache,
                                      err = errs,
                                      # No need for selfConsist if we are confident
                                      # in convergence and consistency of our error model
                                      selfConsist = FALSE,
                                      # performance params
                                      multithread = multithread,
                                      nReads = nReads),
                                    by = c("Sample")]
      
      timeDADA2 = (time0 - Sys.time())
      
      saveRDS(mergeTab,
              file = file.path(ftDir, "dadaTab.RDS"))
      
      ####################
      # Remove Chimeras
      ####################
      seqtab = dcast.data.table(
        data = mergeTab,
        formula = Sample ~ sequence, 
        fun.aggregate = sum, 
        fill = 0, 
        value.var = "abundance")
      seqmat <- as(seqtab[, -1L, with = FALSE], "matrix")
      rownames(seqmat) <- seqtab$Sample
      # Chimera filter on the whole table
      seqmat <- removeBimeraDenovo(seqmat, multithread = NCores)
      dim(seqmat)
      saveRDS(seqmat, "seqmat.RDS")
      # Filter chimeras from mergeTab
      setkey(mergeTab, sequence)
      mergeTabNoChimera = mergeTab[colnames(seqmat)]
      Nseq = mergeTab[, uniqueN(sequence)]
      NseqNoChimera = mergeTabNoChimera[, uniqueN(sequence)]
      message("\nNumber of sequences determined to be chimeras (de novo):\n",
              (Nseq - NseqNoChimera), " out of ", Nseq, " denoised sequences (", 
              round(100*(Nseq - NseqNoChimera)/Nseq, digits = 1), "%)\n\n"
      )
      saveRDS(seqmat, 
              file = file.path(ftDir, "seqmat.RDS"))
      saveRDS(mergeTabNoChimera, 
              file = file.path(ftDir, "dadaTabBimeraFilt.RDS"))
      return(
        paste0("Time to execute DADA Workflow:\n", timeDADA2)
      )
    }
  )
  # Reactive holding the dada results
  DADAS = reactive({
    dadasFile = file.path(inputDirPathRD(), "dadaTabBimeraFilt.RDS")
    # Validate that file exists in order to show plot.
    file.exists(dadasFile) %>% need(message = "...") %>% validate
    # message the errors file path
    message("Learned Errors being stored at:\n",
            dadasFile)
    # A reactive object holding the learned errors.
    dadas <- reactiveFileReader(intervalMillis = 1000,
                               session = session,
                               filePath = dadasFile,
                               readFunc = readRDS)
    return(dadas())
  })
  
  pDADA = reactive({
    # plot_dada_summary(DADAS())
    dummy = DADAS()
    p = ggplot(data = data.frame(x = 1:10, y = 1:10), aes(x, y)) + geom_point()
    return(p)
  })
  
  ########################################
  # The output definitions
  ########################################
  # Filter Trim
  output$showInputPath <- renderText({inputDirPath()})
  output$showInfoPath <- renderText({infoFilePath() %>% basename})
  output$maxLength <- renderText({maxLength()})
  output$infoTable <- DT::renderDataTable({
    fastqFilesTab() 
  }, options = list(lengthChange = FALSE))
  output$includeSamples <- renderText({
    paste0(includeSamples(), collapse = ", ")
  })
  output$QbyC <- renderPlot({pQbyC()})
  output$PrepFT <- renderText({prepFT()})
  output$ExecFT <- renderText({ExecuteFilterTrim()})
  # Learn Errors
  output$showInputPathFT <- renderText({inputDirPathFT()})
  output$showInfoPathFT <- renderText({infoFilePathFT() %>% basename})
  output$infoTableFT <- DT::renderDataTable({
    fastqFilesTabFT() 
  }, options = list(lengthChange = FALSE))
  output$ExecLE <- renderText({ExecuteLearnErrors()})
  output$ForwardErrors <- renderPlot({pErrors()$ForwardErrors})
  output$ReverseErrors <- renderPlot({pErrors()$ReverseErrors})
  # Run DADA
  output$showInputPathRD <- renderText({inputDirPathRD()})
  output$showInfoPathRD <- renderText({infoFilePathRD() %>% basename})
  output$ExecRD <- renderText({ExecuteDADA2()})
  output$dadaMain = renderPlot({pDADA()})
  output$infoTableRD <- DT::renderDataTable({
    fastqFilesTab = fastqFilesTabRD() %>% copy
    fastqFilesTab[, c("FileFull", "PassFilter") := NULL]
    return(fastqFilesTab)
  }, options = list(lengthChange = FALSE))
})
