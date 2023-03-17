library(shiny)
library(dplyr)
library(ape)
library(readr)
library(stringr)
library(tidyr)
library(purrr)
library(treeio)
library(ggtree)
library(ggtreeExtra)
library(ggstar)
library(ggnewscale)
library(ggplot2)
library(ggstance)
library(Biostrings)
library(ggmsa)
library(seqmagick)
library(XVector)
library(gtable)
library(grid)
library(DT)
library(cowplot)
library(ggplotify)
library(shinydashboard)
library(shinydashboardPlus)
library(fontawesome)
library(UpSetR)
library(ggvenn)
library(shinyjs)
library(Polychrome)
library(shinycssloaders)
library(svglite)
library(topGO)
library(Rgraphviz)
library(scales)
library(methods)

#####
# Ui
#####

#set shinycssloaders options
options(spinner.type=1)
options(spinner.color='#428bca')
options(spinner.size=2)

header <- shinydashboardPlus::dashboardHeader(
  #title = "A2TEA",
  ##adding logo to App
  #title = span(img(src = "a2tea_hexsticker.png", height = 35), "test"),
  #https://stackoverflow.com/questions/31440564/adding-a-company-logo-to-shinydashboard-header
  title =
          tags$a(
            href = 'https://github.com/tgstoecker/A2TEA.WebApp',
            target = "_blank",
            class = "dropdown",
            tags$img(
                src='A2TEA.WebApp/a2tea_hexsticker.png',
                height='60',
                style = "margin-top: -10px;")
               ),
#                 style = "padding-top:10px; padding-bottom:10px;")),
  #tags$li(a(href = 'https://www.google.com',
  #          tags$img(src = 'a2tea_hexsticker.png',
  #          title = "Company Home", height = "30px"),
  #          style = "padding-top:10px; padding-bottom:10px;"),
  #          class = "dropdown"),
  tags$li(
    class = "dropdown",
    # the following CSS keeps the header nice looking and on one line
    # the CSS at the beginning of the sidebar definitons is also required for this
    # sizes of individual elements and distan1002ces are important however the most important paramter for this to work is "height"
    # this is set to 60px in all elements and ensures nice design
    ## I'll probably have to adapt this a bit if/when I have a nice logo for A2TEA
    ## e.g. a hex design - https://github.com/GuangchuangYu/hexSticker
    #keep header fixed while scrolling
    tags$script(HTML("$('body').addClass('fixed');")),
    tags$style(".main-header {max-height: 60px;}"),
    tags$style(".main-header .logo {height: 60px;font-size:26px;font-weight:bold;padding:10px;}"),
    tags$style(".sidebar-toggle {height: 60px; padding-top: 12px !important;font-size: 26px;}"),
    tags$style(".navbar {min-height:60px !important}"),
    tags$style(".form-control {padding: 16px 10px;}"),
    tags$style(".btn-box-tool {font-size: 20px;}")),
  # add a github button in the header - to the repo or my profile
  tags$li(a(onclick = "onclick = window.open('https://github.com/tgstoecker')",
            href = NULL,
            icon("github"),
            title = "GitHub",
            style = "cursor: pointer; font-size: 32px; height: 60px;"),
          class = "dropdown"),
  leftUi = tagList(
    tags$li(tags$a("A2TEA",
              onclick = "onclick = window.open('https://github.com/tgstoecker')",
              href = NULL,
              title = "A2TEA",
              style = "cursor: pointer; font-size: 32px; height: 60px; padding-top: 24px; font-weight: 900;"),
            class = "dropdown"),
    #bookmark buton
    actionButton(inputId = "bookmark_click",
                 label = "Bookmark selection/s",
                 icon = shiny::icon("book-bookmark"),
                 style = "color: #ffffff; background-color: #428bca; border-color: #ffffff",
                 class = "bookmark_button"
    ),
    tags$style(".bookmark_button {font-size: 16px; padding: 10px; border-radius: 15px; margin-left: 30px;}"),
    #add custom notification placement - center of the screen
    tags$head(
      tags$style(
        HTML(".shiny-notification {
                       position:fixed;
                       top: calc(50%);
                       left: calc(50%);
                       font-size: 20px;}"
        )
      )
    )
  )
)


sidebar <- shinydashboardPlus::dashboardSidebar(
  # call useShinyjs() in order to use it!
  #useShinyjs(),
  #  actionButton(inputId = "toggle", "Hide/show plot"),
  # as explained in the header creation - necessary CSS to guarantee complete vertical alignment
  tags$style(".left-side, .main-sidebar {padding-top: 60px}"),
  # this code hides the upload and select panels when the the sidebar is toggled off
  # it keeps the sidebar tabs though ;D
  # this is achieved by specifically referencing (and thus hiding) ".form-group"
  # https://www.w3schools.com/css/css_display_visibility.asp
  ### horrible - removed most sliders when toggled to close; tags$style(".sidebar-collapse .form-group {display:none;}"),
  tags$style(".sidebar-collapse .logo {display:none;}"),
  #we can add custom css tags to finely control to ONLY hide sidebar select panels/ the upload field when toggling to closed mode
  #https://stackoverflow.com/questions/67942124/get-reactive-sliders-to-use-css-style-tags-in-r-shiny-instead-of-ignoring-them
  tags$style(".sidebar-collapse .hypothesis_select {display:none;}"),
  tags$style(".sidebar-collapse .data_upload {display:none;}"),
  # Add input of A2TE.RData object
  br(),
  # demo dataset
  tags$div(id = "demo_upload", class="data_upload",
           actionButton(inputId = "demo_A2TEA_choice",
                        label = "Try a demo A2TEA.RData file")
  ),
  #user upload
  tags$div(id = "data_upload", class="data_upload",
           fileInput(inputId = "A2TEA",
                     label = "Upload A2TEA.RData file:",
                     multiple = FALSE,
                     accept = ".RData")
  ),
  # Add selectors for Hypothesis & HOG
  uiOutput('select_Hypothesis'),
  ###  uiOutput('select_HOG'),
  sidebarMenu(id = "a2tea_sidebar_menu",
              menuItem("General", tabName = "general", icon = shiny::icon("house")),
              menuItem("TEA analysis", tabName = "tea", icon = shiny::icon("tree")), #badgeLabel = "new", badgeColor = "green")
              menuItem("Set analyses", tabName = "set_analyses", icon = shiny::icon("not-equal")),
              menuItem("GO term analyses", tabName = "go_term_analyses", icon = shiny::icon("bullseye")),
              menuItem("Bookmarks", tabName = "bookmarks", icon = shiny::icon("book-bookmark"))
  )
)

#create a controlbar - for the bookmarks
controlbar = dashboardControlbar(
  #very useful tip on modifying the controlbar - icon change not considered inside the function
  #https://stackoverflow.com/questions/71058072/remove-toggle-gear-icon-controlbar
  # change icon
  tags$script(
    HTML(
      'var e = document.querySelector("body > div.wrapper > header > nav > div:nth-child(4) > ul > li:last-child > a > i");
           e.setAttribute("class", "dropdown fa fa-book-bookmark");'
    )
  ),
  #add class and modify size via the parent element - icon itself breaks layout
  tags$script(
    HTML(
      'var e = document.querySelector("body > div.wrapper > header > nav > div:nth-child(4) > ul > li:last-child > a");
           e.setAttribute("class", "dropdown cb-changes");
           e.setAttribute("style", "font-size: 30px; padding-right: 25px");'
    )
  ),
  id = "controlbar",
  width = 650,
  br(),
  controlbarMenu(
    id = "cb_menu",
    controlbarItem(
      title = "Gene level bookmarks",
      box(width = 12,
          background ="yellow",
          gradient = TRUE,
          DT::DTOutput("cb_bookmark_genes_df")
      )
    ),
    controlbarItem(
      title = "OG level bookmarks",
      box(width = 12,
          background ="yellow",
          gradient = TRUE,
          DT::DTOutput("cb_bookmark_ogs_df")
      )
    )
  )
)

#footer definition - layout inspired by GeneTonic's approach
A2TEA_footer <- fluidRow(
  column(width = 1,
         background ="yellow",
         gradient = TRUE,
         align = "right",
         tags$a(
           href = "https://github.com/tgstoecker/A2TEA.WebApp",
           target = "_blank",
           tags$img(src = "A2TEA.WebApp/a2tea_hexsticker.png",
                    height = "70px",
                    style = "margin: -10px")
         )
  ),
  column(width = 11,
         background ="yellow",
         gradient = TRUE,
         align = "center",
         "A2TEA is developed by ",
         tags$a(href = "https://github.com/tgstoecker", "Tyll Stöcker,",
                target = "_blank"),
         " Carolin Uebermuth, Florian Boecker & Heiko Schoof in the Crop Bioinformatics group of the ",
         tags$a(href = "https://www.inres.uni-bonn.de/en", "INRES",
                target = "_blank"),
         "- Institute of Crop Science and Resource Conservation", br(),
         "License: ", tags$a(href = "https://opensource.org/licenses/MIT", "MIT",
                             target = "_blank",),
         "- The ",
         tags$a(href = "https://github.com/tgstoecker", "A2TEA workflow",
                target = "_blank"),
         " & ",
         tags$a(href = "https://github.com/tgstoecker", "A2TEA WebApp",
                target = "_blank"),
         "are developed & available on ",
         tags$a(href = "https://github.com/tgstoecker", "GitHub",
                target = "_blank")
  )
)

#CSS tags under dashboardBody()!
footer = dashboardFooter(
  left = A2TEA_footer,
  right = NULL
)

body <- dashboardBody(
  #actually belongs to footer but tags are not accepted there...
  tags$style(".main-footer {background: #F39C12; color: #FFF}"),
  #global change of hyperlink color - https://stackoverflow.com/questions/53399403/r-shiny-dashboard-change-color-for-all-hyperlinks
  tags$head(tags$style(HTML("a {color: #428bca}"))),
  #infoboxes additional styling
  tags$style(".info-box-icon {border-radius: 15px;}"),
  tags$style(".info-box, .bg-aqua {background: #428bca !important; border-radius: 15px;}"),
  tags$style(".info-box-text {text-transform: none !important; font-size: 16px}"),
  tags$style(".info-box-number {font-weight: 700; font-size: 18px;}"),
  #notification styling
  tags$style(".shiny-notification {color: #428bca; background-color: #F39C12; border: 5px solid #428bca}"),
  tags$style(".shiny-notification-warning {color: #5342ca; background-color: #FF4800; border: 5px solid #428bca}"),
  # call useShinyjs() in order to use it!
  #useShinyjs(),
  # tag "tea"
  tags$style(".progress {border-radius: 15px;}"),
  tags$style(".dataTable {color: #428bca}"),
  tags$style(".selectize-input, .option {color: #428bca; border-radius: 15px;}"),
  tags$style(".selectize-dropdown, .option {color: #428bca; border-radius: 15px;}"),
  tags$style(".btn {color: #428bca; border-radius: 15px;}"),
  tags$style(".form-control {color: #428bca; border-radius: 15px;}"),
  #for custom links buttons
  tags$style(".btn-custom-link {color: #fff !important; background: #428bca; border-radius: 15px;}"),
  tags$style(".dataTables_wrapper .dataTables_filter input {border-radius: 15px;}"),
  tags$style("div.dt-buttons {margin-left: 20px;}"),
  tags$style("button.dt-button {border-radius: 15px;}"),
  # tag "bookmark"
  tags$style(".box-footer {color: #428bca}"),
  tags$style(".bookmark_reset {margin-top: 57px}"),
  tags$style(".rdata_subset_export {color: white; background-color: #428bca}"),

  tabItems(
    # first tab: general information; species plot, upsetR/venn diagramme
    tabItem(tabName = "general",
            h2("General information on the experiment"),

            # row1
            fluidRow(
              # Column 1
              column(width = 7,
                     shinyjs::hidden(
                       div(id = "id_diff_exp_table_outer",
                           box(
                             style = "height: 600px; overflow-y: auto; overflow-x: auto",
                             background ="yellow",
                             gradient = TRUE,
                             width = NULL,
                             collapsed = TRUE,
                             collapsible = TRUE,
                             title = "Differential expression results",
                             radioButtons(
                               inputId = "ensembl_db_choice",
                               choices = c("Plants",
                                           "Bacteria",
                                           "Protists",
                                           "Metazoa",
                                           "Fungi",
                                           "Vertebrates"),
                               label = "Ensembl DB choice",
                               selected = "Plants",
                               inline = TRUE,
                               width = "100%"
                             ),
                             #hidden button that we use for complete datatable download
                             downloadButton("download_deg_full",
                                            "", style = "visibility: hidden;"),
                             shinycssloaders::withSpinner(
                               DT::DTOutput("DEG_table")
                             )
                           )
                       ))
              ),
              # Column 2
              column(width = 5,
                     # in order to use shinyjs on box elements - wrap them in div() expressions
                     # give the divs their own ids! in order to use shinyjs on the box
                     shinyjs::hidden(
                       div(id = "id_speciesPlot_outer",
                           box(
                             style = "height: 600px; overflow-y: auto; overflow-x: auto",
                             id = "id_speciesPlot",
                             background ="yellow",
                             gradient = TRUE,
                             collapsed = TRUE,
                             collapsible = TRUE,
                             sidebar = boxSidebar(
                               icon = fontawesome::fa("gears"),
                               width = 25,
                               startOpen = FALSE,
                               id = "mycardsidebar",
                               selectInput(inputId = "species_tree_choice",
                                           label = "All species or hypothesis subset:",
                                           choices = c("all species", "hypothesis species"),
                                           selected = "hypothesis species"),
                               numericInput(inputId = "species_tree_height", "height", value = 600),
                               sliderInput(inputId = "species_tree_hexpand",
                                           "horizontal axis expansion",
                                           min = 0.1,
                                           max = 4.0,
                                           value = 0.7),
                               selectInput(inputId = "species_plot_export_choice",
                                           label = "Export format",
                                           choices = c("pdf", "jpeg", "tiff", "png", "bmp", "svg"),
                                           selected = "pdf"),
                               downloadButton(outputId = "species_tree_download", label = "Download Plot")
                             ),
                             width = NULL,
                             title = "Phylogenetic Tree - Orthofinder/STAG/STRIDE",
                             shinycssloaders::withSpinner(
                               plotOutput("speciesPlot")
                             )
                           )
                       )
                     )
              )
            ),

            # row 2
            fluidRow(
              column(width = 7,
                     shinyjs::hidden(
                       div(id = "id_func_anno_table_outer",
                           box(
                             style = "height: 600px; overflow-y: auto; overflow-x: auto",
                             background ="yellow",
                             gradient = TRUE,
                             width = NULL,
                             collapsed = TRUE,
                             collapsible = TRUE,
                             title = "Functional annnotation",
                             #hidden button that we use for complete datatable download
                             downloadButton("download_func_anno_full",
                                            "", style =
                                              "visibility: hidden;"),
                             shinycssloaders::withSpinner(
                               DTOutput("func_anno_table")
                             )
                           )
                       ))
              ),
              column(width = 5,
                     shinyjs::hidden(
                       div(id = "id_hypothesis_conservation_outer",
                           box(
                             style = "height: 600px; overflow-y: auto; overflow-x: auto",
                             background ="yellow",
                             gradient = TRUE,
                             width = NULL,
                             collapsed = TRUE,
                             collapsible = TRUE,
                             title = "Orthologous Groups - Conservation",
                             shinycssloaders::withSpinner(
                               plotOutput("hypothesis_conservation_plot")
                             ),
                             sidebar = boxSidebar(
                               icon = fontawesome::fa("gears"),
                               width = 25,
                               startOpen = FALSE,
                               id = "conservation_plot_sidebar",
                               uiOutput('select_hypothesis_conservation_plot'),
                               numericInput(inputId = "hypothesis_conservation_plot_height", "Plot Height", value = 600),
                               uiOutput('hypothesis_conservation_venn_hexpand_slider'),
                               selectInput(inputId = "conservation_plot_export_choice",
                                           label = "Export format",
                                           choices = c("pdf", "jpeg", "tiff", "png", "bmp", "svg"),
                                           selected = "pdf"),
                               downloadButton(outputId = "conservation_plot_download", label = "Download Plot")
                             )
                           )
                       )
                     )
              )
            ),
            #some space
            br(),
            br(),
            #3rd - + rows
            div(id = "id_hypothesis_infos_general",
                uiOutput("general_ui_infoboxes")
            )
    ),

    # second tab: Trait-specific evolutionary adaptation analysis
    tabItem(tabName = "tea",
            h2("Trait-specific evolutionary adaption analysis"),

            # first row
            fluidRow(
              column(width = 12,
                     shinyjs::hidden(
                       div(id = "hypothesis_HOG_level_table_outer",
                           box(style = "height: 600px; overflow-y: auto; overflow-x: auto",
                               background ="yellow",
                               gradient = TRUE,
                               width = NULL,
                               collapsed = TRUE,
                               collapsible = TRUE,
                               title = "OG table",
                               #hidden button that we use for complete datatable download
                               downloadButton("download_hog_full",
                                              "", style = "visibility: hidden;"),
                               shinycssloaders::withSpinner(
                                 DTOutput("hypothesis_HOG_level_table")
                               )
                           )
                       )
                     )
              )
            ),
            shinyjs::hidden(
              div(id = "id_tea_vis_tool_outer",
                  box(width = 12,
                      gradient = TRUE,
                      background = "yellow",
                      collapsed = FALSE,
                      collapsible = TRUE,
                      title = "TEA visualization tool",
                      #(visual) options inside the sidebar
                      sidebar = boxSidebar(
                        width = 25,
                        startOpen = FALSE,
                        icon = fontawesome::fa("gears"),
                        id = "tea_tree_mycardsidebar",
                        sliderInput(
                          inputId = "ortho_tree_hexpand",
                          "Horizontal Axis Expansion",
                          min = 0.1,
                          max = 5.0,
                          value = 0.1,
                          step = 0.1
                        ),
                        #choice where to put the legend
                        selectInput(
                          inputId = "tree_plot_legend_pos",
                          label = "Legend position",
                          choices = c("left", "top", "right", "bottom", "none"),
                          selected = "right"
                        ),
                        # define size of tree lines
                        sliderInput(
                          min = 0.1,
                          max = 2.5,
                          value = 1.0,
                          inputId = "tree_size",
                          label = "Tree Lines - Size"
                        ),
                        # setting plot height
                        sliderInput(
                          inputId = "tree_height",
                          "Plot height (px)",
                          value = 600,
                          min = 500,
                          max = 5000,
                          step = 50
                        ),
                        #sidebar options for adding visual changes to log2 FC layer
                        uiOutput("tea_log2FC_n_break_options"),
                        uiOutput("tea_log2FC_axis_text_options"),
                        uiOutput("tea_log2FC_vjust_options"),
                        uiOutput("tea_log2FC_linetype_options"),
                      ),
                      fluidRow(
                        # row 2 - column 1 (only 1)
                        column(
                          width = 4,
                          boxPad(
                            fluidRow(
                              column(width = 4,
                                     uiOutput('select_HOG')
                              ),
                              column(width = 8,
                                     checkboxInput(inputId = "bookmark_viewed_OG_choice", label = "Add OG to bookmark-choices?", value = FALSE, width = NULL)
                              )
                            ),
                            uiOutput('add_OGs_control'),
                            selectInput(
                              inputId = "tree_layout",
                              label = "Tree layout",
                              # removed choices 'dendrogram', 'inward_circular', 'equal_angle', 'daylight', 'ape'
                              # either not supported for the layers I want to display or not useful for TEA analysis
                              choices = c('rectangular', 'slanted', 'ellipse', 'roundrect', 'circular', 'fan', 'radial')
                            ),
                            # turn off/on branch_length option
                            selectInput(
                              inputId = "branch_length",
                              label = "Branch length",
                              choices = c('branch.length', 'none')
                            ),
                            splitLayout(
                              checkboxInput(inputId = "log2FC_layer",
                                            label = "Add log2FC layer",
                                            value = FALSE
                              ),
                              checkboxInput(inputId = "tea_log2FC_add_sign_stars_choice",
                                            label = "Mark significant DEGs",
                                            value = FALSE
                              ),
                            ),
                            splitLayout(
                              uiOutput("tea_log2FC_offset_options"),
                              uiOutput("tea_log2FC_pwidth_options"),
                            ),
                            checkboxInput(inputId = "msa_layer",
                                          label = "Add MSA layer - only recommended for NON-circular views",
                                          value = FALSE),
                            uiOutput("tea_msa_width_options"),
                            splitLayout(
                              uiOutput("tea_msa_offset_options"),
                              uiOutput("tea_msa_pwidth_options"),
                            ),
                            selectInput(inputId = "tree_plot_export_choice",
                                        label = "Export format",
                                        choices = c("pdf", "jpeg", "tiff", "png", "bmp", "svg"),
                                        selected = "pdf"),
                            downloadButton(outputId = "tree_plot_download", label = "Download Plot"),
                          )
                        ),
                        column(width = 8,
                               shinyjs::hidden(
                                 div(id = "id_treePlot_outer",
                                     # https://stackoverflow.com/questions/49129194/r-shiny-dynamic-box-height
                                     # height: XX, where XX can be any CSS unit (I personally prefer viewpoint ratios, e.g. height: 33vh will make a box as high as one third of a screen, regardless of screen resolution), sets the height of your box;
                                     # overflow-y: auto adds a vertical scrollbar if required.
                                     # With this approach, when a user resizes the screen all boxes maintain their equal heights.
                                     #                    box(style = "height: 1000px; overflow-y: auto; overflow-x: auto",
                                     box(style = "height: 600px; overflow-y: auto; overflow-x: auto",
                                         background ="yellow",
                                         gradient = TRUE,
                                         width = NULL,
                                         title = "Tree Plot - Orthologous Group",
                                         shinycssloaders::withSpinner(
                                           plotOutput('treePlot'
                                           )
                                         )
                                     )
                                 )
                               )
                        )
                      )
                  ))),
            # row 3
            fluidRow(
              # row 3 - column 1
              column(width = 12,
                     shinyjs::hidden(
                       div(id = "id_treePlot_table_outer",
                           # https://stackoverflow.com/questions/49129194/r-shiny-dynamic-box-height
                           # height: XX, where XX can be any CSS unit (I personally prefer viewpoint ratios, e.g. height: 33vh will make a box as high as one third of a screen, regardless of screen resolution), sets the height of your box;
                           # overflow-y: auto adds a vertical scrollbar if required.
                           # With this approach, when a user resizes the screen all boxes maintain their equal heights.
                           box(style = "height: 600px; overflow-y: auto; overflow-x: auto",
                               background ="yellow",
                               gradient = TRUE,
                               width = NULL,
                               collapsed = TRUE, collapsible = TRUE,
                               title = "Expanded OG + extended BLAST hits",
                               #hidden button that we use for complete datatable download
                               downloadButton("download_blast_full",
                                              "", style = "visibility: hidden;"),
                               shinycssloaders::withSpinner(
                                 DTOutput("ortho_tree_table")
                               )
                           )
                       )
                     )
              )
            ),
            # row 4 - column 1
            shinyjs::hidden(
              div(id = "msa_vis_tool_outer",
                  box(width = 12, gradient = TRUE, background = "yellow", title = "MSA", collapsible = TRUE, collapsed = TRUE,
                      fluidRow(
                        column(width = 8,
                               shinyjs::hidden(
                                 div(id = "id_msa_solo_plot_outer",
                                     # https://stackoverflow.com/questions/49129194/r-shiny-dynamic-box-height
                                     # height: XX, where XX can be any CSS unit (I personally prefer viewpoint ratios, e.g. height: 33vh will make a box as high as one third of a screen, regardless of screen resolution), sets the height of your box;
                                     # overflow-y: auto adds a vertical scrollbar if required.
                                     # With this approach, when a user resizes the screen all boxes maintain their equal heights.
                                     #                    box(style = "height: 600px; overflow-y: auto; overflow-x: auto",
                                     box(style = "height: 800px; overflow-y: auto; overflow-x: auto",
                                         background ="yellow",
                                         gradient = TRUE,
                                         width = NULL,
                                         title = "MSA plot - choice of displayed Ortogroup/s in TEA analysis tool",
                                         conditionalPanel(
                                           condition = "input.msa_solo_button > 0",
                                           style = "display: none;",
                                           withSpinner(plotOutput("msa_solo_plot"))
                                         )
                                     )
                                 )
                               )
                        ),
                        column(width = 4,
                               shinyjs::hidden(
                                 div(id = "id_msa_solo_plot_options",
                                     boxPad(
                                       closable = TRUE,
                                       collapsible = TRUE,
                                       width = 75,
                                       startOpen = FALSE,
                                       id = "msa_sidebar",
                                       splitLayout(
                                         numericInput(inputId = "msa_plot_height", "height", value = 800),
                                         numericInput(inputId = "msa_plot_width", "width", value = 15000)
                                       ),
                                       splitLayout(
                                         checkboxInput(inputId = "msa_plot_add_conservation_choice",
                                                       label = "Add conservation plot",
                                                       value = FALSE
                                         ),
                                         selectInput(inputId = "msa_solo_color_scheme",
                                                     label = "Color scheme",
                                                     choices = c('Clustal', 'Chemistry_AA', 'Shapely_AA', 'Zappo_AA',
                                                                 'Taylor_AA', 'LETTER', 'CN6', 'Chemistry_NT', 'Shapely_NT', 'Zappo_NT', 'Taylor_NT'),
                                                     selected = "Chemistry_AA")
                                       ),
                                       sliderInput(inputId = "msa_solo_id_labels_margin", label = "ID labels margin", min = -100, max = 100, value = 0, step = 10),
                                       sliderInput(inputId = "msa_solo_id_labels_size", label = "ID labels size", min = 5, max = 150, value = 20, step = 5),
                                       sliderInput(inputId = "msa_solo_pos_labels_size", label = "Position labels size", min = 5, max = 150, value = 20, step = 5),
                                       selectInput(inputId = "msa_solo_plot_export_choice",
                                                   label = "Export format",
                                                   choices = c("pdf", "svg"),
                                                   selected = "pdf"),
                                       actionButton("msa_solo_button", "Create plot"),
                                       downloadButton(outputId = "msa_solo_plot_download", label = "Download Plot")
                                     )
                                 )
                               )
                        )
                      )
                  )
              )
            )
    ),
    tabItem(tabName = "set_analyses",
            h2("Set analyses"),

            fluidRow(
              shinyjs::hidden(
                div(id = "id_set_analyses_stats_options",
                    box(width = 9,
                        title = "Enrichment analysis",
                        background ="yellow",
                        gradient = TRUE,
                        column(width = 3,
                               sliderInput(inputId = "n_deg_set_choice",
                                           label = "# DEGs for set choices",
                                           min = 1,
                                           max = 10,
                                           value = 1,
                                           step = 1
                               ),
                               selectInput(inputId = "background_set_choice",
                                           label = "Choose background",
                                           choices = c("all OGs",
                                                       "all OGs with # DEGs from ANY species",
                                                       "conserved OGs",
                                                       "conserved OGs; # DEGs from ANY species",
                                                       "conserved OGs; # DEGs from ALL species",
                                                       #hypothesis specific but not strictly conserved
                                                       "hypothesis; all OGs; EXPANDED",
                                                       "hypothesis; all OGs; # DEGs from ANY species",
                                                       "hypothesis; all OGs; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; all OGs; # DEGs from ALL EXPANDED species",
                                                       "hypothesis; all OGs; EXPANDED; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; all OGs; EXPANDED; # DEGs from ALL EXPANDED species",
                                                       #hypothesis specific AND conserved
                                                       "hypothesis; conserved OGs; EXPANDED",
                                                       "hypothesis; conserved OGs; # DEGs from ANY species",
                                                       "hypothesis; conserved OGs; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; conserved OGs; # DEGs from ALL EXPANDED species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL EXPANDED species"
                                           ),
                                           selected = "all OGs"),

                               selectInput(inputId = "background_subset_choice",
                                           label = "Choose subset in background",
                                           choices = c("all OGs",
                                                       "all OGs with # DEGs from ANY species",
                                                       "conserved OGs",
                                                       "conserved OGs; # DEGs from ANY species",
                                                       "conserved OGs; # DEGs from ALL species",
                                                       #hypothesis specific but not strictly conserved
                                                       "hypothesis; all OGs; EXPANDED",
                                                       "hypothesis; all OGs; # DEGs from ANY species",
                                                       "hypothesis; all OGs; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; all OGs; # DEGs from ALL EXPANDED species",
                                                       "hypothesis; all OGs; EXPANDED; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; all OGs; EXPANDED; # DEGs from ALL EXPANDED species",
                                                       #hypothesis specific AND conserved
                                                       "hypothesis; conserved OGs; EXPANDED",
                                                       "hypothesis; conserved OGs; # DEGs from ANY species",
                                                       "hypothesis; conserved OGs; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; conserved OGs; # DEGs from ALL EXPANDED species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL EXPANDED species"
                                           ),
                                           selected = "all OGs with # DEGs from ANY species"),

                               selectInput(inputId = "interest_set_choice",
                                           label = "Choose set of interest",
                                           choices = c("all OGs",
                                                       "all OGs with # DEGs from ANY species",
                                                       "conserved OGs",
                                                       "conserved OGs; # DEGs from ANY species",
                                                       "conserved OGs; # DEGs from ALL species",
                                                       #hypothesis specific but not strictly conserved
                                                       "hypothesis; all OGs; EXPANDED",
                                                       "hypothesis; all OGs; # DEGs from ANY species",
                                                       "hypothesis; all OGs; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; all OGs; # DEGs from ALL EXPANDED species",
                                                       "hypothesis; all OGs; EXPANDED; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; all OGs; EXPANDED; # DEGs from ALL EXPANDED species",
                                                       #hypothesis specific AND conserved
                                                       "hypothesis; conserved OGs; EXPANDED",
                                                       "hypothesis; conserved OGs; # DEGs from ANY species",
                                                       "hypothesis; conserved OGs; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; conserved OGs; # DEGs from ALL EXPANDED species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL EXPANDED species"
                                           ),
                                           selected = "conserved OGs"),

                               selectInput(inputId = "interest_subset_choice",
                                           label = "Choose subset of interest",
                                           choices = c("all OGs",
                                                       "all OGs with # DEGs from ANY species",
                                                       "conserved OGs",
                                                       "conserved OGs; # DEGs from ANY species",
                                                       "conserved OGs; # DEGs from ALL species",
                                                       #hypothesis specific but not strictly conserved
                                                       "hypothesis; all OGs; EXPANDED",
                                                       "hypothesis; all OGs; # DEGs from ANY species",
                                                       "hypothesis; all OGs; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; all OGs; # DEGs from ALL EXPANDED species",
                                                       "hypothesis; all OGs; EXPANDED; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; all OGs; EXPANDED; # DEGs from ALL EXPANDED species",
                                                       #hypothesis specific AND conserved
                                                       "hypothesis; conserved OGs; EXPANDED",
                                                       "hypothesis; conserved OGs; # DEGs from ANY species",
                                                       "hypothesis; conserved OGs; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; conserved OGs; # DEGs from ALL EXPANDED species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY EXPANDED species",
                                                       "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL EXPANDED species"
                                           ),
                                           selected = "conserved OGs; # DEGs from ANY species")
                        ),
                        column(width = 4,
                               shinycssloaders::withSpinner(
                                 DTOutput("set_analysis_df")
                               )
                        ),
                        column(width = 5,
                               textOutput("phyper_header"),
                               tags$head(tags$style("#phyper_header{
                                     color: 'black';
                                     font-size: 20px;
                                     font-style: bold;
                                     }"
                               )
                               ),
                               verbatimTextOutput("phyper_code"),
                               br(),
                               br(),
                               shinycssloaders::withSpinner(
                                 verbatimTextOutput("phyper_choices")
                               ),
                               br(),
                               br(),
                               splitLayout(cellWidths = c("25%", "75%"),
                                           actionButton(inputId = "calc_set_phyper", label = "Calculate"),
                                           shinycssloaders::withSpinner(
                                             textOutput("phyper_value")
                                           ),
                                           tags$head(tags$style("#phyper_value{
                                       color: #428bca;
                                       font-size: 20px;
                                       font-style: italic;
                                       }"
                                           )
                                           )
                               )
                        )
                    )
                ))
            ),

            fluidRow(
              shinyjs::hidden(
                div(id = "id_set_analyses_circ_and_og_size",
                    column(width = 6,
                           box(width = 12,
                               title = "Circular set plot",
                               background ="yellow",
                               gradient = TRUE,
                               column(width = 2,
                                      numericInput(inputId = "set_circ_plot_height", "height", value = 450),
                                      selectInput(inputId = "set_circ_plot_export_choice",
                                                  label = "Export format",
                                                  choices = c("pdf", "jpeg", "tiff", "png", "bmp", "svg"),
                                                  selected = "pdf"),
                                      br(),
                                      br(),
                                      downloadButton(outputId = "set_circ_plot_download", label = "Download")
                               ),
                               column(width = 10,
                                      boxPad(style = "height: 500px; overflow-y: auto; overflow-x: auto",
                                             #title = "Circular set plot",
                                             shinycssloaders::withSpinner(
                                               plotOutput("set_circ_plot")
                                             )
                                      )
                               ),
                               #(visual) options inside the sidebar
                               sidebar = boxSidebar(
                                 width = 25,
                                 startOpen = FALSE,
                                 icon = fontawesome::fa("gears"),
                                 id = "set_circ_sidebar",
                                 sliderInput(inputId = "interest_ratio_centerdist_choice", label = "Interest ratio - dist. from center", min = 0.1, max = 4, value = 1.5, step = 0.1),
                                 sliderInput(inputId = "interest_ratio_clockpos_choice", label = "Interest ratio - clock position", min = 0.1, max = 1, value = 0.1, step = 0.1),
                                 selectInput(inputId = "interest_ratio_color", label = "Interest ratio - text color", choices = c("black", "white"), selected = "black"),
                                 sliderInput(inputId = "interest_ratio_size", label = "Interest ratio - text size", min = 1, max = 15, value = 8, step = 1),
                                 sliderInput(inputId = "background_ratio_centerdist_choice", label = "Background ratio - dist. from center", min = 0.1, max = 4, value = 3, step = 0.1),
                                 sliderInput(inputId = "background_ratio_clockpos_choice", label = "Background ratio - clock position", min = 0.1, max = 1, value = 0.1, step = 0.1),
                                 selectInput(inputId = "background_ratio_color", label = "Background ratio - text color", choices = c("black", "white"), selected = "black"),
                                 sliderInput(inputId = "background_ratio_size", label = "Background ratio - text size", min = 1, max = 15, value = 8, step = 1),
                                 textInput(inputId = "set_circ_custom_x_axis_label", label = "Custom x-axis label", value = "Ratio: Interest (inside) vs. background (outside)")
                               )
                           )
                    ),
                    column(width = 6,
                           box(width = 12,
                               title = "OG sizes distribution",
                               background ="yellow",
                               gradient = TRUE,
                               column(width = 3,
                                      splitLayout(
                                        numericInput(inputId = "set_og_size_plot_height", "height", value = 450),
                                        ""
                                      ),
                                      selectInput(inputId = "set_og_size_background_choice",
                                                  label = "Choose the background set",
                                                  choices = c("all OGs",
                                                              "all OGs with # DEGs from ANY species",
                                                              "conserved OGs",
                                                              "conserved OGs; # DEGs from ANY species",
                                                              "conserved OGs; # DEGs from ALL species",
                                                              #hypothesis specific but not strictly conserved
                                                              "hypothesis; all OGs; EXPANDED",
                                                              "hypothesis; all OGs; # DEGs from ANY species",
                                                              "hypothesis; all OGs; # DEGs from ANY EXPANDED species",
                                                              "hypothesis; all OGs; # DEGs from ALL EXPANDED species",
                                                              "hypothesis; all OGs; EXPANDED; # DEGs from ANY EXPANDED species",
                                                              "hypothesis; all OGs; EXPANDED; # DEGs from ALL EXPANDED species",
                                                              #hypothesis specific AND conserved
                                                              "hypothesis; conserved OGs; EXPANDED",
                                                              "hypothesis; conserved OGs; # DEGs from ANY species",
                                                              "hypothesis; conserved OGs; # DEGs from ANY EXPANDED species",
                                                              "hypothesis; conserved OGs; # DEGs from ALL EXPANDED species",
                                                              "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY species",
                                                              "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL species",
                                                              "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY EXPANDED species",
                                                              "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL EXPANDED species"
                                                  ),
                                                  selected = "hypothesis; conserved OGs; EXPANDED"),

                                      selectInput(inputId = "set_og_size_interest_choice",
                                                  label = "Choose subset of interest",
                                                  choices = c("all OGs",
                                                              "all OGs with # DEGs from ANY species",
                                                              "conserved OGs",
                                                              "conserved OGs; # DEGs from ANY species",
                                                              "conserved OGs; # DEGs from ALL species",
                                                              #hypothesis specific but not strictly conserved
                                                              "hypothesis; all OGs; EXPANDED",
                                                              "hypothesis; all OGs; # DEGs from ANY species",
                                                              "hypothesis; all OGs; # DEGs from ANY EXPANDED species",
                                                              "hypothesis; all OGs; # DEGs from ALL EXPANDED species",
                                                              "hypothesis; all OGs; EXPANDED; # DEGs from ANY EXPANDED species",
                                                              "hypothesis; all OGs; EXPANDED; # DEGs from ALL EXPANDED species",
                                                              #hypothesis specific AND conserved
                                                              "hypothesis; conserved OGs; EXPANDED",
                                                              "hypothesis; conserved OGs; # DEGs from ANY species",
                                                              "hypothesis; conserved OGs; # DEGs from ANY EXPANDED species",
                                                              "hypothesis; conserved OGs; # DEGs from ALL EXPANDED species",
                                                              "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY species",
                                                              "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL species",
                                                              "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY EXPANDED species",
                                                              "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL EXPANDED species"
                                                  ),
                                                  selected = "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL species"),
                                      selectInput(inputId = "set_og_size_plot_export_choice",
                                                  label = "Export format",
                                                  choices = c("pdf", "jpeg", "tiff", "png", "bmp", "svg"),
                                                  selected = "pdf"),
                                      downloadButton(outputId = "set_og_size_plot_download", label = "Download")
                               ),
                               column(width = 9,
                                      boxPad(style = "height: 500px; overflow-y: auto; overflow-x: auto",
                                             shinycssloaders::withSpinner(
                                               plotOutput("set_og_size_plot")
                                             )
                                      )
                               ),
                               #(visual) options inside the sidebar
                               sidebar = boxSidebar(
                                 width = 25,
                                 startOpen = FALSE,
                                 icon = fontawesome::fa("gears"),
                                 id = "set_og_size_sidebar",
                                 sliderInput(inputId = "set_og_size_plot_alpha_choice",
                                             label = "Fill alpha value",
                                             min = 0.1,
                                             max = 1,
                                             value = 0.9,
                                             step = 0.1)
                               )
                           )
                    )
                ))
            )
    ),
    tabItem(tabName = "go_term_analyses",
            h2("Gene ontology term - enrichment analysis"),

            fluidRow(
              column(width = 4,
                     shinyjs::hidden(
                       div(id = "id_go_term_set_choices_outer",
                           box(
                             background ="yellow",
                             gradient = TRUE,
                             width = 12,
                             title = "Go term set choices",
                             selectInput(inputId = "go_ontology_choice",
                                         label = "Ontology",
                                         choices = c('BP', 'MF', 'CC'),
                                         selected = "BP"),
                             #hypothesis choice + "conserved_vs_all"
                             #needs server side rendering since I don't know how many hypotheses are in the experiment...
                             selectInput(inputId = "go_expansion_choice",
                                         label = "Expansion",
                                         choices = c('yes', 'no'),
                                         selected = "yes"),
                             selectInput(inputId = "go_deg_choice",
                                         label = "DEG criterium",
                                         choices = c("none",
                                                     "at least # DEGs in any hypothesis species",
                                                     "at least # DEGs in ANY expanded species",
                                                     "at least # DEGs in ALL expanded species"),
                                         selected = "none"),
                             sliderInput(inputId = "go_n_deg_choice",
                                         label = "Minimum # DEGs for DEG criterium",
                                         min = 1,
                                         max = 10,
                                         value = 1,
                                         step = 1),
                             selectInput(inputId = "go_algo_choice",
                                         label = "Algorithm",
                                         choices = c('classic', 'elim', 'weight', 'weight01', 'parentchild'),
                                         selected = "weight01"),
                             sliderInput(inputId = "go_top_n_nodes_choice",
                                         label = "Top # of Nodes",
                                         min = 1,
                                         max = 100,
                                         step = 1,
                                         value = 33),
                             sliderInput(inputId = "go_sig_filter_choice",
                                         label = "Min. # of sig. GO term annotations",
                                         min = 1,
                                         max = 5,
                                         step = 1,
                                         value = 1),
                             actionButton("go_start_button", "Compute")
                           )
                       )
                     )
              ),
              column(width = 8,
                     shinyjs::hidden(
                       div(id = "id_go_term_output_table_outer",
                           box(style = "height: 700px; overflow-y: auto; overflow-x: auto",
                               background ="yellow",
                               gradient = TRUE,
                               width = NULL,
                               title = "Go term set choices",
                               shinycssloaders::withSpinner(
                                 DTOutput("go_term_output_table")
                               )
                           )
                       )
                     )
              )
            ),
            #Datatable of significant OGs
            fluidRow(
              column(
                width = 4,
                shinyjs::hidden(
                  div(id = "id_go_sig_OGs_df_outer",
                      box(style = "height: 450px; overflow-y: auto; overflow-x: auto",
                          width = 12,
                          background ="yellow",
                          gradient = TRUE,
                          collapsible = TRUE,
                          collapsed = TRUE,
                          title = "Sign. OGs/genes of OGs for chosen GO term",
                          splitLayout(
                            uiOutput("select_GO_term"),
                            radioButtons(inputId = "sig_go_terms_OGs_or_genes_choice",
                                         inline = TRUE,
                                         label = "OG or genes level",
                                         choices = c("OG level", "Gene level"),
                                         selected = "OG level")
                          ),
                          shinycssloaders::withSpinner(
                            DTOutput("sig_OGs_for_term_df")
                          )
                      )
                  )
                )
              ),
              column(
                width = 8,
                shinyjs::hidden(
                  div(id = "id_go_enrich_outer",
                      box(
                        width = 12,
                        background ="yellow",
                        gradient = TRUE,
                        collapsible = TRUE,
                        collapsed = TRUE,
                        title = "GO enrichment plots",
                        column(
                          width = 9,
                          div(id = "id_go_enrich_plots",
                              boxPad(style = "height: 450px; overflow-y: auto; overflow-x: auto",
                                     shinycssloaders::withSpinner(
                                       plotOutput("go_enrich_plot"),
                                     )
                              )
                          )
                        ),
                        column(
                          width = 3,
                          shinyjs::hidden(
                            div(id = "id_go_enrich_options",
                                boxPad(
                                  background ="yellow",
                                  gradient = TRUE,
                                  title = "Go term enrichment plots",
                                  startOpen = TRUE,
                                  numericInput(inputId = "go_enrich_plot_height", "height", value = 450),
                                  sliderInput(inputId = "go_enrich_N_GOs_choice",
                                              label = "# of top GO terms to display",
                                              min = 1,
                                              max = 30,
                                              step = 1,
                                              value = 10),
                                  radioButtons(inputId = "enrich_plot_choice",
                                               label = "Enrichment plot choice",
                                               choices = c("dot", "bar"),
                                               selected = "dot"),
                                  selectInput(inputId = "go_enrich_plot_export_choice",
                                              label = "Export format",
                                              choices = c("pdf", "jpeg", "tiff", "png", "bmp", "svg"),
                                              selected = "pdf"),
                                  splitLayout(
                                    actionButton(inputId = "go_enrich_button", "Create plot"),
                                    downloadButton(outputId = "go_enrich_plot_download", label = "Download")
                                  )
                                )
                            )
                          )
                        )
                      )
                  )
                )
              )
            ),
            #Graph vis. of GO terms
            shinyjs::hidden(
              div(id = "go_graph_tool_outer",
                  box(style = "height: 800px; overflow-y: auto; overflow-x: auto",
                      width = 12,
                      gradient = TRUE,
                      background = "yellow",
                      title = "Go graph",
                      collapsible = TRUE,
                      collapsed = TRUE,
                      fluidRow(
                        column(
                          width = 3,
                          shinyjs::hidden(
                            div(id = "id_go_graph_options",
                                boxPad(
                                  #closable = TRUE,
                                  #collapsible = TRUE,
                                  width = 75,
                                  #startOpen = TRUE,
                                  icon = fontawesome::fa("gears"),
                                  id = "id_go_graph_sidebar",
                                  splitLayout(
                                    numericInput(inputId = "go_graph_height", "height", value = 1000),
                                    numericInput(inputId = "go_graph_width", "width", value = 1200)
                                  ),
                                  sliderInput(inputId = "go_graph_top_n_nodes_choice",
                                              label = "Top # of nodes to display",
                                              min = 1,
                                              max = 100,
                                              step = 1,
                                              value = 1),
                                  splitLayout(
                                    sliderInput(inputId = "go_graph_n_max_chars_choice",
                                                label = "Max name length for GO terms",
                                                min = 1,
                                                max = 100,
                                                step = 1,
                                                value = 20),
                                    sliderInput(inputId = "go_graph_cex_choice",
                                                label = "Text size modifier (cex)",
                                                min = 0.1,
                                                max = 2,
                                                step = 0.1,
                                                value = 0.3)
                                  ),
                                  selectInput(inputId = "go_graph_export_choice",
                                              label = "Export format",
                                              choices = c("pdf", "jpeg", "tiff", "png", "bmp", "svg"),
                                              selected = "pdf"),
                                  splitLayout(
                                    actionButton(inputId = "go_graph_button",
                                                 label = "Create graph"),
                                    downloadButton(outputId = "go_graph_download", label = "Download Plot")
                                  )
                                ),
                            )
                          )
                        ),
                        column(width = 9,
                               shinyjs::hidden(
                                 div(id = "id_go_graph_plot",
                                     boxPad(
                                       style = "height: 800px; overflow-y: auto; overflow-x: auto",
                                       background ="yellow",
                                       gradient = TRUE,
                                       width = NULL,
                                       title = "Graph of GO terms",
                                       shinycssloaders::withSpinner(
                                         plotOutput("go_graph")
                                       ),
                                     )
                                 )
                               )
                        )
                      )
                  )
              )
            )
    ),
    tabItem(tabName = "bookmarks",
            h2("Bookmarks"),
            uiOutput("ui_bookmarks")
    )
  )
)

ui <- shinydashboardPlus::dashboardPage(
  #activate Shinyjs
  useShinyjs(),
  skin = "yellow",
  header = header,
  sidebar = sidebar,
  controlbar = controlbar,
  body = body,
  footer = footer
)
