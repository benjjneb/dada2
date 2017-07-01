library("shiny")
library("shinyFiles")

################################################################################
# UI FilterTrim
################################################################################

# Sidebar Panel FilterTrim
sbp_filtertrim = sidebarPanel(
  width = 4,
  fluidRow(
    column(width = 6, h4("Data Selection")),
    column(width=6,
           shinyDirButton('dirOG', 'Select Data Folder', 'Please select a FastQ Seq Folder')
    )
  ),
  fluidRow(
    column(width=9,
           h5("Input Directory"),
           verbatimTextOutput('showInputPath', placeholder = TRUE)
    ),
    column(width=3,
           h5("Info File"),
           verbatimTextOutput('showInfoPath', placeholder = TRUE)
    )
  ),
  fluidRow(
    h4("Quality Downsample"),
    column(width = 5,
           sliderInput("nreads", "Reads / Sample (log10)", value = 3, min = 2, max = 5, step = 1)
    ),
    column(width = 5,
           sliderInput("NSamples", "Number of Samples", value = 10L, min = 1L, max = 50L, step = 1L, round = TRUE)
    ),
    column(width = 2,
           h5("Max Length"),
           verbatimTextOutput('maxLength', placeholder = TRUE)
    )
  ),
  fluidRow(
    h4("Trimming Prediction"),
    # input$minQual # 15
    column(width = 6,
           sliderInput("minQual", "Min Smoothed Quality", value = 15, min = 0, max = 40, step = 1)
    ),
    column(width = 6,
           sliderInput("qquantile", "Quality Quantile", value = 0.2, min = 0.0, max = 1.0, step = 0.05)
    )
  ),
  fluidRow(
    h4("Trimming Parameters"),
    column(width = 6,
           uiOutput("uiForward")
    ),
    column(width = 6,
           uiOutput("uiReverse")
    ),
    h4("Parameters: dada2::filterAndTrim()"),
    column(width = 3,
           numericInput("maxEE", "maxEE", value = 2L, min = 0L, step = 1L)
    ),
    column(width = 3,
           sliderInput("multithread", "multithread", value = NCores - 1L, min = 1L, max = NCores, step = 1L, round = TRUE)
    ),
    column(width = 3,
           sliderInput("n", "n [log10]", value = 5, min = 2, max = 8, step = 1)
    ),
    column(width = 3,
           numericInput("truncQ", "truncQ", value = 2L, min = 0L, step = 1L)
    )
  ),
  h4("Fastq Info Table"),
  DT::dataTableOutput('infoTable')
)

uiFilterTrim = fluidPage(
  h4("Filter and trim fastq sequences"),
  sbp_filtertrim,
  column(
    width=8,
    mainPanel(
      actionButton("actionb_filtertrim", "Execute FilterTrim", icon("filter")),
      h5("Samples included in quality profiling:"),
      textOutput("includeSamples"),
      h5("Quality by Cycle"),
      plotOutput("QbyC"),
      textOutput("PrepFT"),
      textOutput("ExecFT")
    )
  ),
  fluidRow(column(width = 12,
                  "On-page documentation here.",
                  "Directory should have... oligo-removed, sample-demultiplexed fastq sequences."))
)

################################################################################
# UI Learn Errors
################################################################################
# Sidebar Panel Learn Errors
sbp_learnerrors = sidebarPanel(
  width = 4,
  fluidRow(
    column(width = 6, h4("Data Selection")),
    column(width=6,
           shinyDirButton('dirFT', 'Select Filtertrim Folder',
                          'Please select the filtered and trimmed fastq folder.')
    )
  ),
  fluidRow(
    column(width=9,
           h5("Input Directory"),
           verbatimTextOutput('showInputPathFT', placeholder = TRUE)
    ),
    column(width=3,
           h5("Info File"),
           verbatimTextOutput('showInfoPathFT', placeholder = TRUE)
    )
  ),
  fluidRow(
    h4("Parameters: dada2::learnErrors()"),
    column(width = 3,
           sliderInput("LE_learnSize", "Number of Samples", value = 6, min = 1, max = 30, step = 1)
    ),
    column(width = 3,
           sliderInput("LE_minSize", "Min. File Size [log10]", value = 4, min = 3, max = 6, step = 1)
    ),
    column(width = 3,
           sliderInput("LE_nreads", "nreads [log10]", value = 5, min = 2, max = 8, step = 1)
    ),
    column(width = 3,
           sliderInput("LE_multithread", "multithread", value = NCores - 1L, min = 1L, max = NCores, step = 1L, round = TRUE)
    )
  ),
  h4("Fastq Info Table"),
  DT::dataTableOutput('infoTableFT')
)

uiLearnErrors = fluidPage(
  h4("Learn error matrix from a subset of samples in each direction"),
  sbp_learnerrors,
  column(
    width=8,
    mainPanel(
      actionButton("actionb_learnErrors", "Learn Errors",
                   icon("leanpub",
                        lib = "font-awesome")),
      h4("Forward Error Rates"),
      plotOutput("ForwardErrors"),
      h4("Reverse Error Rates"),
      plotOutput("ReverseErrors"),
      fluidRow(
        column(width = 12,
               textOutput("ExecLE"))
      )
    )
  ),
  fluidRow(column(width = 12, "On-page documentation here."))
)

################################################################################
# UI Run DADA2
################################################################################
# Sidebar Panel Run DADA2 (RD)
sbp_rundada = sidebarPanel(
  width = 4,
  fluidRow(
    column(width = 6, h4("Data Selection")),
    column(width=6,
           shinyDirButton('dirRD', 'Select Filtertrim Folder',
                          'Please select the filtered and trimmed fastq folder.')
    )
  ),
  fluidRow(
    column(width=9,
           h5("Input Directory"),
           verbatimTextOutput('showInputPathRD', placeholder = TRUE)
    ),
    column(width=3,
           h5("Info File"),
           verbatimTextOutput('showInfoPathRD', placeholder = TRUE)
    )
  ),
  fluidRow(
    h4("Parameters: dada2::dada()"),
    column(width = 3,
           sliderInput("RD_minSize", "Min. File Size [log10]", value = 4, min = 3, max = 6, step = 1)
    ),
    column(width = 3,
           sliderInput("RD_nreads", "nreads [log10]", value = 5, min = 2, max = 8, step = 1)
    ),
    column(width = 3,
           sliderInput("RD_multithread", "multithread", value = NCores - 1L, min = 1L, max = NCores, step = 1L, round = TRUE)
    )
  ),
  h4("Fastq Info Table"),
  DT::dataTableOutput('infoTableRD')
)

uiRunDADA = fluidPage(
  h4("Learn error matrix from a subset of samples in each direction"),
  sbp_rundada,
  column(
    width=8,
    mainPanel(
      fluidRow(
        column(
          width = 12,
          actionButton("actionb_rundada", "Run DADA",
                       icon("play",
                            lib = "font-awesome"))
        ),
        column(width = 12, textOutput("ExecRD")),
        plotOutput("dadaMain")
      )
      # (actionButton("actionb_chimerafilter", "Filter Chimeras",
      #               icon("group",
      #                    lib = "font-awesome")) %>% column(width = 12) %>% fluidRow),
      # h4("Chimera Summary here..."),
      # actionButton("actionb_shiftmeras", "Aggregate Shiftmeras",
      #              icon("align-center",
      #                   lib = "font-awesome")),
      # h4("Shiftmera Summary here...")
    )
  ),
  fluidRow(column(width = 12, "On-page documentation here."))
)

################################################################################
# Define general header tag list 
# List of tags to display as a common header above all tabPanels.
################################################################################
headerTagList = list(
  tags$style(type="text/css", ".phyloseq-print { font-size: 10px; }"),
  tags$base(target="_blank")
)
################################################################################
# Define the full user-interface, `ui`
################################################################################
ui = navbarPage(
  # The content/UI
  tabPanel("Filter Trim", uiFilterTrim),
  tabPanel("Learn Errors", uiLearnErrors),
  tabPanel("Run DADA", uiRunDADA),
  ## Meta
  # Top left corner title text
  title = "Shiny DADA2" %>% 
    a(href="http://benjjneb.github.io/dada2/",
      style="color:#F0F0F0") %>% 
    h4,
  # HTML headers for each panel
  header = headerTagList,
  collapsible = TRUE,
  windowTitle = "Shiny DADA2"
)
################################################################################
# There can only be one of these calls in the app code:
################################################################################
shinyUI(ui)
