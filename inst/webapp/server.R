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
# Server
#####
server <- function(input, output, session) {
  
  #############
  #create session specific dir for intermediary files etc. - deleted on close
  #############
  
  #create session directory - to be deleted after close
  currwd <- getwd()
  session_dir <- fn_dir(currwd)
  
  
  ###########
  # data upload and handling
  ###########    
  
  # change maximum file upload size = default 100Mb
  options(shiny.maxRequestSize=1000*1024^2)
  
  ## event reactive expression for loading the data once uploaded by the user
  # hypotheses.tsv (simple df with)
  hypotheses_tsv <- eventReactive(input$A2TEA, {
    t = new.env()
    load(input$A2TEA$datapath, envir = t)
    get("hypotheses", envir=t)
  })
  # core A2TEA object
  HYPOTHESES.a2tea <- eventReactive(input$A2TEA, {
    t = new.env()
    load(input$A2TEA$datapath, envir = t)
    get("HYPOTHESES.a2tea", envir=t)
  })
  # Gene level diff. exp. table
  HOG_DE.a2tea <- eventReactive(input$A2TEA, {
    t = new.env()
    load(input$A2TEA$datapath, envir = t)
    get("HOG_DE.a2tea", envir=t)
  })
  # List of HOGs
  HOG_level_list <- eventReactive(input$A2TEA, {
    t = new.env()
    load(input$A2TEA$datapath, envir = t)
    get("HOG_level_list", envir=t)
  })
  # List of HOGs
  SFA <- eventReactive(input$A2TEA, {
    t = new.env()
    load(input$A2TEA$datapath, envir = t)
    get("SFA", envir=t)
  }) 
  #all species tree
  all_speciesTree <- eventReactive(input$A2TEA, {
    t = new.env()
    load(input$A2TEA$datapath, envir = t)
    get("all_speciesTree", envir=t)
  })
  #pep fasta sequences
  A2TEA.fa.seqs <- eventReactive(input$A2TEA, {
    t = new.env()
    load(input$A2TEA$datapath, envir = t)
    get("A2TEA.fa.seqs", envir=t)
  })
  
  ## rendering of UI elements for the sidebar
  # choice of hypothesis to be displayed
  output$select_Hypothesis = renderUI({
    #we can add custom css tags to finely control to ONLY hide sidebar select panels/ the upload field when toggling to closed mode
    #https://stackoverflow.com/questions/67942124/get-reactive-sliders-to-use-css-style-tags-in-r-shiny-instead-of-ignoring-them
    tags$div(id = "hypothesis_select", class="hypothesis_select",
             selectInput(inputId = 'select_Hypothesis_server', 
                         label = 'Select Hypothesis:', 
                         choices = hypotheses_tsv()$name
             )
    )
  })
  
  #we define all cases in which we want or animation of "in-sliding" panels to occur
  #on read-in of data, changes between different sidebar menu options and when changing between hypotheses
  animation_req <- reactive({
    paste(input$A2TEA, input$a2tea_sidebar_menu, input$select_Hypothesis_server)
  })
  
  ## observe upload of data and only then show boxes in which the plots are displayed
  #also animate when changing between sidebar menu options ;D
  # by delaying the showing we circumvent an early quick popup when changing between sidebar pages
  # the logic here is quite neat: once we change to another menu via the sidebar only the new view is opened and all others closed;
  # this way on the next change the view is empty and the user sees the animation without prev. results being visible for  short period of time ;D
  # General page
  observeEvent(animation_req(), {
    req(input$A2TEA)
    shinyjs::hideElement("id_speciesPlot_outer", anim = FALSE)
    shinyjs::hideElement("id_hypothesis_conservation_outer", anim = FALSE)
    shinyjs::hideElement("id_diff_exp_table_outer", anim = FALSE)
    shinyjs::hideElement("id_func_anno_table_outer", anim = FALSE)
    shinyjs::hideElement("id_hypothesis_infos_general", anim = FALSE)
    if (input$a2tea_sidebar_menu == "general") {
      delay(ms = 100,
            c(shinyjs::showElement("id_speciesPlot_outer", anim = TRUE, animType = "slide", time = 1),
              shinyjs::showElement("id_hypothesis_conservation_outer", anim = TRUE, animType = "slide", time = 1),
              shinyjs::showElement("id_diff_exp_table_outer", anim = TRUE, animType = "slide", time = 1),
              shinyjs::showElement("id_func_anno_table_outer", anim = TRUE, animType = "slide", time = 1),
              shinyjs::showElement("id_hypothesis_infos_general", anim = TRUE, animType = "slide", time = 1)
            )
      )
    }
    # TEA analysis tab
    shinyjs::hideElement("hypothesis_HOG_level_table_outer", anim = FALSE)
    shinyjs::hideElement("id_treePlot_table_outer", anim = FALSE)
    shinyjs::hideElement("id_tea_vis_tool_outer", anim = FALSE)
    shinyjs::hideElement("id_treePlot_outer", anim = FALSE)
    shinyjs::hideElement("msa_vis_tool_outer", anim = FALSE)
    shinyjs::hideElement("id_msa_solo_plot_outer", anim = FALSE)
    shinyjs::hideElement("id_msa_solo_plot_options", anim = FALSE)
    if (input$a2tea_sidebar_menu == "tea") {
      delay(ms = 100,
            expr = c(shinyjs::showElement("hypothesis_HOG_level_table_outer", anim = TRUE, animType = "slide", time = 1),
                     shinyjs::showElement("id_treePlot_table_outer", anim = TRUE, animType = "slide", time = 1),
                     shinyjs::showElement("id_tea_vis_tool_outer", anim = TRUE, animType  = "slide", time = 1),
                     shinyjs::showElement("id_treePlot_outer", anim = TRUE, animType = "slide", time = 1),
                     shinyjs::showElement("msa_vis_tool_outer", anim = TRUE, animType = "slide", time = 1),
                     shinyjs::showElement("id_msa_solo_plot_outer", anim = TRUE, animType = "slide", time = 1),
                     shinyjs::showElement("id_msa_solo_plot_options", anim = TRUE, animType = "slide", time = 1)
            )
      )
    }
    #set analyses
    shinyjs::hideElement("id_set_analyses_stats_options", anim = FALSE)
    shinyjs::hideElement("id_set_analyses_circ_and_og_size", anim = FALSE)      
    if (input$a2tea_sidebar_menu == "set_analyses") {
      delay(ms = 100,
            expr = c(
              shinyjs::showElement("id_set_analyses_stats_options", anim = TRUE, animType = "slide", time = 1),
              shinyjs::showElement("id_set_analyses_circ_and_og_size", anim = TRUE, animType = "slide", time = 1)
            )
      )
    }
    #go term analyses
    shinyjs::hideElement("id_go_term_set_choices_outer", anim = FALSE)
    shinyjs::hideElement("id_go_term_output_table_outer", anim = FALSE)
    shinyjs::hideElement("id_go_sig_OGs_df_outer", anim = FALSE)
    shinyjs::hideElement("id_go_enrich_outer", anim = FALSE)
    shinyjs::hideElement("id_go_enrich_plots", anim = FALSE)
    shinyjs::hideElement("id_go_enrich_options", anim = FALSE)
    shinyjs::hideElement("go_graph_tool_outer", anim = FALSE)
    shinyjs::hideElement("id_go_graph_options", anim = FALSE)
    shinyjs::hideElement("id_go_graph_plot", anim = FALSE)
    if (input$a2tea_sidebar_menu == "go_term_analyses") {
      delay(ms = 100,
            expr = c(
              shinyjs::showElement("id_go_term_set_choices_outer", anim = TRUE, animType = "slide", time = 1),
              shinyjs::showElement("id_go_term_output_table_outer", anim = TRUE, animType = "slide", time = 1),
              shinyjs::showElement("id_go_sig_OGs_df_outer", anim = TRUE, animType = "slide", time = 1),
              shinyjs::showElement("id_go_enrich_outer", anim = TRUE, animType = "slide", time = 1),
              shinyjs::showElement("id_go_enrich_plots", anim = TRUE, animType = "slide", time = 1),
              shinyjs::showElement("id_go_enrich_options", anim = TRUE, animType = "slide", time = 1),
              shinyjs::showElement("go_graph_tool_outer", anim = TRUE, animType = "slide", time = 1),
              shinyjs::showElement("id_go_graph_options", anim = TRUE, animType = "slide", time = 1),
              shinyjs::showElement("id_go_graph_plot", anim = TRUE, animType = "slide", time = 1)
            )
      )
    }
    #bookmarks
    shinyjs::hideElement("id_bookmarks_outer", anim = FALSE)
    if (input$a2tea_sidebar_menu == "bookmarks") {
      delay(ms = 100,
            expr = c(
              shinyjs::showElement("id_bookmarks_outer", anim = TRUE, animType = "slide", time = 1)
            )
      )
    }
  })
  
  
  
  # a dropdown list of the expanded HOGs of the current hypothesis
  # -> observe choice of HOG; e.g. for building trees
  # first create reactive list of all exp. HOGs 
  hypothesisExpOGs <- reactive({
    # require the input$select_Hypothesis_server to be present, since otherwise for a brief second an error is thrown
    req(input$select_Hypothesis_server)
    vars <- all.vars(
      parse(
        text = names(HYPOTHESES.a2tea()[[hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server]]]@expanded_OGs)
      )
    )    
    vars <- as.list(vars)
    return(vars)    
  })
  
  
  
  output$select_HOG = renderUI({     
    selectInput('select_HOG_server', 'Select OG', hypothesisExpOGs())
  })
  
  
  
  
  ########################
  # General (start) page
  ########################
  
  
  #reactive DEG table
  res_deg_df  <- reactive({
    res_deg_df <- HOG_DE.a2tea() %>%
      #get rid of NAs since this will clash with datatable JS
      mutate(log2FoldChange = case_when(
        is.na(log2FoldChange) ~ 0,
        TRUE ~ log2FoldChange
      ) 
      ) %>%
      #round all numerical columns that or not p-value/padj
      mutate(across(3:6, round, 3)) %>%
      #rename HOG to OG
      dplyr::rename(OG = HOG)
    
    res_deg_df$NCBI <- createNCBILink(res_deg_df$gene)
    res_deg_df$Ensembl <- createEnsemblPlantsLink(val = res_deg_df$gene, 
                                                  ensembl_db = input$ensembl_db_choice)  
    return(res_deg_df)
  })
  
  
  # render table displaying functional annotation information
  # advanced filtering is enabled
  output$DEG_table <- DT::renderDT(server = TRUE, {
    req(res_deg_df())   
    
    DT::datatable(
      res_deg_df(),
      
      filter = list(position = 'top', clear = FALSE),
      extensions = "Buttons",
      options = list(
        dom = 'lBfrtip',
        lengthMenu = list(c(10, 25, 50, 100, 200), c('10', '25', '50', '100', '200')),
        pageLength = 10,
        buttons = list(
          'copy', 'print',
          list(extend = "csv",
               fieldSeparator = "\t",
               extension = ".tsv",
               text = "Download Current Page", filename = "deg_results",
               exportOptions = list(
                 modifier = list(page = "current")
               )
          )
        ),
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      #since rownames = FALSE the column position is n-1
      #we actually need only columns 2 and 10 to be selectable
      selection=list(mode="multiple",
                     target="cell",
                     #allow only to select cell from gene and OG columns
                     selectable = rbind(cbind(1:nrow(res_deg_df()), rep(1, nrow(res_deg_df()))),
                                        cbind(1:nrow(res_deg_df()), rep(9, nrow(res_deg_df()))))
      ),
      rownames=FALSE,
      #no escape - so hyperlinks are displayed
      escape = FALSE
    ) %>%
      formatStyle(
        "log2FoldChange",
        background = styleColorBar_divergent(
          res_deg_df()$log2FoldChange,
          scales::alpha("#428bca", 0.4),
          scales::alpha("#F39C12", 0.4)
        ),
        backgroundSize = "100% 90%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      )
  })  
  
  #compute reactive for transformed functional annotation table
  SFA_general <- eventReactive(SFA(),{
    
    SFA_general <- tibble(
      "species" = character(),
      "OG" = character(),
      "gene" = character(),
      "Blast-Hit-Accession" = character(),
      "Human-Readable-Description" = character(),
      "Gene-Ontology-Term" = character()
    )
    
    for (i in 1:length(SFA())) {
      
      current <- SFA()[[i]] %>%
        dplyr::select(-`AHRD-Quality-Code`) %>%
        dplyr::rename(OG = HOG, gene = `Protein-Accession`) %>%
        mutate(species = names(SFA())[i]) %>%
        dplyr::relocate(species, .before = OG) %>%
        mutate(OG = case_when(
          is.na(OG) ~ "singleton",
          TRUE ~ OG
        )
        )
      
      SFA_general <- dplyr::bind_rows(SFA_general, current)
    }
    return(SFA_general)
    
  })
  
  SFA_general_OG_level <- eventReactive(SFA_general(), {
    SFA_general_OG_level <- SFA_general() %>%
      #remove all singleton rows
      filter(OG != "singleton") %>%
      dplyr::select(OG, `Human-Readable-Description`, `Gene-Ontology-Term`) %>%
      group_by(OG) %>%
      #remove redundancy while collapsing with ", "
      summarise(`Gene-Ontology-Term` = paste(unique(`Gene-Ontology-Term`), collapse=", "),
                `Human-Readable-Description` = paste(unique(`Human-Readable-Description`), collapse=", ")) %>%
      #we should remove all rows that have NO GO terms associated
      #and remove all NAs from the Gene-Ontology-Term column if they have at least 1 GO term
      #remove inline NAs in Gene-Ontology-Term column
      mutate(
        `Gene-Ontology-Term` = str_remove_all(`Gene-Ontology-Term`, "NA, |, NA") 
      )
    
    return(SFA_general_OG_level)
  })
  
  
  # render table displaying functional annotation information
  # advanced filtering is enabled
  output$func_anno_table <- renderDT(server = TRUE, {
    datatable(
      SFA_general(), 
      filter = list(position = 'top', clear = FALSE),
      
      extensions = "Buttons",
      options = list(
        dom = 'lBfrtip',
        lengthMenu = list(c(10, 25, 50, 100, 200), c('10', '25', '50', '100', '200')),
        pageLength = 10,
        buttons = list(
          'copy', 'print',
          list(extend = "csv",
               fieldSeparator = "\t",
               extension = ".tsv",
               text = "Download Current Page", filename = "functional_annotation",
               exportOptions = list(
                 modifier = list(page = "current")
               )
          )
        ),
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      selection=list(mode="multiple", 
                     target="cell",
                     #allow only to select cell from gene and OG columns
                     selectable = rbind(cbind(1:nrow(SFA_general()), rep(1, nrow(SFA_general()))),
                                        cbind(1:nrow(SFA_general()), rep(2, nrow(SFA_general()))))
      ),
      rownames = FALSE
    )
  })
  
  # creating the species tree plot
  # make renderUI element reactive
  # access to element based on name in first element selectInput('')
  output$speciesPlot <- renderPlot(
    height = function() input$species_tree_height,
    {
      req(input$select_Hypothesis_server)
      req(input$species_tree_choice)
      
      if (input$species_tree_choice == "all species") {
        
        # assigning hypothesis specific species tree  
        p <- all_speciesTree()
        p <- ggtree::ggtree(p, aes(color = "#F39C12", show.legend = FALSE)) +
          ggtree::hexpand(input$species_tree_hexpand, direction = 1) +
          ggtree::hexpand(input$species_tree_hexpand, direction = -1) +
          ggtree::geom_tiplab(color = "#428bca", show.legend = FALSE) +
          theme(plot.margin=margin(6, 120, 6, 6),
                legend.position = "none")
        
      } else if (input$species_tree_choice == "hypothesis species") {
        
        # assigning hypothesis specific species tree  
        p <- HYPOTHESES.a2tea()[[hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server]]]@species_tree
        p <- ggtree::ggtree(p, aes(color = "#F39C12", show.legend = FALSE)) +
          ggtree::hexpand(input$species_tree_hexpand, direction = 1) +
          ggtree::hexpand(input$species_tree_hexpand, direction = -1) +
          ggtree::geom_tiplab(color = "#428bca", show.legend = FALSE) +
          theme(plot.margin=margin(6, 120, 6, 6),
                legend.position = "none")
      }
      
      
      ggsave(filename = paste0(session_dir, "/", "species_tree.", input$species_plot_export_choice), 
             plot = as.grob(p),
             device = input$species_plot_export_choice,
             limitsize = FALSE)
      
      return(p)
    })
  
  #download option species tree
  output$species_tree_download <- downloadHandler(
    filename = function() {
      paste0(session_dir, "/", "species_tree.", input$species_plot_export_choice)
    },
    content = function(file) {
      file.copy(paste0(session_dir, "/", "species_tree.", input$species_plot_export_choice), file, overwrite=TRUE)
    }
  )
  
  
  ## Creating upsetR and venn plots
  # based on hypothesis/ when it is changed, create specific list with incl. genotypes/species
  genotypeList <- eventReactive(input$select_Hypothesis_server, {
    
    workset <- HOG_level_list()[[hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server]]] %>% 
      dplyr::select(HOG, ends_with("_total")) #%>%
    #for testing purposes only - else remove!
    workset <- workset %>% setNames(names(workset) %>% stringr::str_replace("_total",""))
    
    # a vector with all species names for currently selected hypothesis
    species_vec <- names(workset %>% dplyr::select(-HOG))
    
    # using map_at from purrr
    # we perform the operation for all columns in species_vec
    # filter/only keep HOG "rows" per species if species possesses at least 1 gene inside the HOG
    # with purrr::list_modify("HOG" = NULL) the untouched pure HOG column (now a list) is removed
    genotypeList <- workset %>% purrr::map_at(species_vec, 
                                              ~filter(workset, .!= 0) %>% 
                                                dplyr::select(HOG) %>% 
                                                as_vector() %>%
                                                unname()
    ) %>% purrr::list_modify("HOG" = NULL)
    
    # return genotypeList
    genotypeList
  })
  
  
  #check how many species in hypothesis
  hypothesis_conservation_plot_choices <- reactive({
    req(input$select_Hypothesis_server)
    req(genotypeList)
    
    if (length(genotypeList()) > 4) {
      options = "UpSetR"
    } else {
      options = c("UpSetR", "Venn")
    }
    return(options)
  })
  
  output$select_hypothesis_conservation_plot = renderUI({     
    selectInput(inputId = 'hypothesis_OGs_conservation_vis_choice', 
                label = 'Venn (<= 4 species) / UpSetR', 
                choices = hypothesis_conservation_plot_choices())
  })
  
  
  #create horizontal expansion slider - only if Venn plot was chosen
  output$hypothesis_conservation_venn_hexpand_slider = renderUI({ 
    req(input$hypothesis_OGs_conservation_vis_choice == "Venn")
    
    sliderInput(inputId = "hypothesis_conservation_venn_hexpand", 
                label = "Horizontal axis expansion - labels", 
                min = 0.0, 
                max = 2.0,
                step = 0.1,
                value = 0.4)
  })
  
  #create conservation plot
  output$hypothesis_conservation_plot <- renderPlot(
    
    height = function() input$hypothesis_conservation_plot_height,
    #    width = function() input$hypothesis_conservation_plot_width,  
    {
      req(input$select_Hypothesis_server)
      req(input$hypothesis_OGs_conservation_vis_choice)
      
      
      if (input$hypothesis_OGs_conservation_vis_choice == "UpSetR") {
        
        # create upset plot - needs ggplotify to write to object and export
        
        p <- as.ggplot(
          upset(fromList(genotypeList()), 
                order.by = "freq",
                empty.intersections = "on",
                point.size = 3.5, 
                line.size = 2, 
                mainbar.y.label = "HOG Intersections", 
                text.scale = c(2, 1.3, 1, 1, 2, 2),
                main.bar.color = "navy",
                matrix.color = "navy",
                sets.bar.color = "navy"
          )
        )
      }
      
      else if (input$hypothesis_OGs_conservation_vis_choice == "Venn") {
        
        req(input$hypothesis_conservation_venn_hexpand)
        
        # create Venn plot 
        p <- ggvenn(genotypeList(), 
                    stroke_linetype = 2, 
                    text_size = 5, 
                    set_name_color = "navy",
                    stroke_color = "yellow3",
                    stroke_size = 0.8, 
                    set_name_size = 6, 
                    show_percentage = TRUE, 
                    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")
        ) + 
          scale_x_continuous(expand = expansion(mult = input$hypothesis_conservation_venn_hexpand))
      }
      
      conservation_plot <- function() {
        plot(p)
      }
      
      #write plot to file for export
      ggsave(filename = paste0(session_dir, "/", "conservation_plot.", input$conservation_plot_export_choice), 
             plot = as.grob(function() conservation_plot()),
             device = input$conservation_plot_export_choice,
             limitsize = FALSE)
      
      return(conservation_plot())
      
    })
  
  #download option for conservation plot
  output$conservation_plot_download <- downloadHandler(
    filename = function() {
      paste0(session_dir, "/", "conservation_plot.", input$conservation_plot_export_choice)
    },
    content = function(file) {
      file.copy(paste0(session_dir, "/", "conservation_plot.", input$conservation_plot_export_choice), file, overwrite=TRUE)
    }
  )
  
  
  #quick stats reactives for start page
  species_num <- reactive({
    req(input$select_Hypothesis_server)
    
    exp_in_num <- hypotheses_tsv() %>%
      filter(name == input$select_Hypothesis_server) %>%
      pull(expanded_in) %>%
      str_split(., pattern = ";") %>%
      unlist() %>%
      length()
    comp_to_num <- hypotheses_tsv() %>%
      filter(name == input$select_Hypothesis_server) %>%
      pull(compared_to) %>%
      str_split(., pattern = ";") %>%
      unlist() %>%
      length()
    
    species_num <- exp_in_num + comp_to_num
    
    return(species_num)
  })
  
  #num of expanded OGs can easily be computed by length of input choices in tree analysis
  
  #number of DEGs for all species part of a hypothesis
  species_degs_num <- reactive({
    req(input$select_Hypothesis_server)
    
    exp_in_names <- hypotheses_tsv() %>%
      filter(name == input$select_Hypothesis_server) %>%
      pull(expanded_in) %>%
      str_split(., pattern = ";") %>%
      unlist()
    
    comp_to_names <- hypotheses_tsv() %>%
      filter(name == input$select_Hypothesis_server) %>%
      pull(compared_to) %>%
      str_split(., pattern = ";") %>%
      unlist()
    
    species_names <- c(exp_in_names, comp_to_names)
    
    species_degs_num <- HOG_DE.a2tea() %>%
      filter(species %in% species_names) %>%
      filter(significant == "yes") %>%
      nrow()
    
    return(species_degs_num)
  })
  
  output$general_ui_infoboxes <- renderUI({
    req(input$A2TEA)

    fluidRow(
      column(width = 6,
             valueBoxOutput("general_infobox", width = 12),
             valueBoxOutput("general_infobox2", width = 12),
             valueBoxOutput("general_infobox3", width = 12)
      ),
      column(width=6,        
             #hypotheses table
             box(
               style = "height: 600px; overflow-y: auto; overflow-x: auto",
               background ="yellow",
               gradient = TRUE,
               width = NULL,
               collapsed = TRUE,
               collapsible = TRUE,
               title = "Hypothesis details:",
               shinycssloaders::withSpinner(
                 DTOutput("hypotheses_DT_table")
               )
             )
      )
    ) 
  })
  
  output$general_infobox <- renderInfoBox({
    req(input$select_Hypothesis_server)
    req(species_num())
    
    infoBox(icon = icon("chart-column", lib = "font-awesome"),
            width = 6,
            title = "# Species in hypothesis",
            value = paste0(species_num()),
            fill = TRUE
    )
  })
  
  output$general_infobox2 <- renderInfoBox({
    req(input$select_Hypothesis_server)
    req(hypothesisExpOGs())
    
    infoBox(icon = icon("project-diagram"),
            width = 6,
            title = "# of Orthologous Groups passing expansion criteria",
            value = paste(length(hypothesisExpOGs())),
            fill = TRUE
    )
  })
  
  output$general_infobox3 <- renderInfoBox({
    req(input$select_Hypothesis_server)
    req(species_degs_num())
    
    infoBox(icon = icon("chart-pie"),
            width = 6,
            title = "# Differentially expressed genes/transcripts",
            value = paste(species_degs_num()),
            fill = TRUE
    )
  })
  
  #table of hypothesis options in A2TEA pipeline
  output$hypotheses_DT_table <- renderDT(server = TRUE, {

    hypotheses_choice_DT <- hypotheses_tsv() %>%
      filter(name == input$select_Hypothesis_server) %>%
      #rename column header to easier understandable/longer phrases
      dplyr::rename("Index number of hypothesis" = hypothesis) %>%
      dplyr::rename("Name of hypothesis" = name) %>%
      dplyr::rename("Species checked for gene family expansion events" = expanded_in) %>%
      dplyr::rename("Species which were compared to" = compared_to) %>%
      dplyr::rename("Min. # of checked/expanded species that need to fulfill expansion criteria" = Nmin_expanded_in) %>%
      dplyr::rename("Min. # of compared to species of that need to fulfill expansion criteria ('contraction')" = Nmin_compared_to) %>%
      dplyr::rename("Min. factor of expansion" = min_expansion_factor) %>%
      dplyr::rename("Min. # of additional genes in exp. species" = min_expansion_difference) %>%
      dplyr::rename("Min. # of genes of expanded species required to consider the OG for expansion" = Nmin_expanded_genes) %>%
      dplyr::rename("Was ploidy normalization performed?" = ploidy_normalization) %>%
      dplyr::rename("Does every member of the expanded species have to be present (not expanded!) in the OG?" = expanded_in_all_found) %>%
      dplyr::rename("Does every member of the compared to species have to be present (not expanded!) in the OG?" = compared_to_all_found) %>%
      base::t(.)
    
    datatable(
      hypotheses_choice_DT, 
      filter = list(position = 'none', clear = FALSE),
      
      extensions = "Buttons",
      options = list(
        dom = 'Bt',
        pageLength = 12,
        buttons = list(
          'copy', 'print',
          list(extend = "csv",
               fieldSeparator = "\t",
               extension = ".tsv",
               text = "Download", filename = paste0("hypothesis_", 
                                                    hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server],
                                                    "_information"),
               exportOptions = list(
                 modifier = list(page = "current")
               )
          )
        ),
        ordering = FALSE,
        searching = FALSE#,
        #columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      colnames = c("", ""),
      rownames= TRUE
    ) 
  })
  
  
  ########################
  # Tea analysis 
  ########################
  
  # read-in hypothesis specific HOG_level_list
  hypothesis_HOG_level_list <- eventReactive(input$select_Hypothesis_server, {
    HOG_level_list()[[hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server]]]  
  })
  
  
  # the complete list of extended blast hits for this HOG - as_vector(tibble)
  blast_hits_table <- eventReactive(input$select_HOG_server, {
    HYPOTHESES.a2tea()[[hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server]]]@expanded_OGs[[input$select_HOG_server]]@blast_table  
  })
  
  ## for adding OGs analysis
  # adding closest complete OGs
  num_max_additional_OGs <- eventReactive(input$select_HOG_server, {
    nmaOGs <- length(HYPOTHESES.a2tea()[[hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server]]]@expanded_OGs[[input$select_HOG_server]]@add_OG_analysis)
    nmaOGs
  })
  
  
  xx_changes <- reactive({
    paste(input$select_HOG_server, input$add_OGs)
  })
  
  xx_changes <- xx_changes %>% debounce(2000)
  
  #I need abit redundant code in the following (best solution as of yet)
  #when switching between OGs on the adding OGs analysis I have to account for the fact that I don't always,
  #have the same amount of addtional sets available; this can lead to the problem that if the slider for
  #for additional sets is on e.g. 5 and I switch to another exp. OG with only 2 additonal sets I get a small error
  #while I have an observeEvent that resets the slider on such changes this is not evaluated before these reactive
  #expressions; the workaround for now is that I have dedicated eventReactive for:
  #case 1: switch of OG - set addtional OGs/sets to 0 (= set 1 = only the expanded OG itself)
  #case 2: adding/removing additional OGs - set the value to whatever the user choses
  #for both "tree_choice_add_OGs" & "all_genes_HOG_add_OGs" this means that  we can have more or less the same
  #code except this minor detail
  
  #building correct tree
  tree_choice_add_OGs <- eventReactive(xx_changes(), {#list(input$select_HOG_server, input$add_OGs), {
    #actually requires user choices/the standard value to be known here (defined further down)
    req(input$select_HOG_server)
    req(input$add_OGs)
    #we add 1 since the #1 set is just the exp. OG itself and we have substracted in the slider input accordingly
    HYPOTHESES.a2tea()[[hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server]]]@expanded_OGs[[input$select_HOG_server]]@add_OG_analysis[[input$add_OGs+1]]@tree
  }, ignoreNULL = FALSE, ignoreInit = FALSE) 
  
  
  all_genes_HOG_add_OGs <- eventReactive(xx_changes(), {
    req(input$select_HOG_server)
    req(input$add_OGs)
    #we add 1 since the #1 set is just the exp. OG itself and we have substracted in the slider input accordingly
    agHaOGs <- HYPOTHESES.a2tea()[[hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server]]]@expanded_OGs[[input$select_HOG_server]]@add_OG_analysis[[input$add_OGs+1]]@genes
    as.vector(agHaOGs)
  })
  
  
  #for correct coloring of the tree in the adding OGs analysis we need to base it 
  #on the set of genes in the largest set (so exp. OG + all additional OGs)
  max_all_genes_HOG_add_OGs <- eventReactive(input$select_HOG_server, {
    choice_OG <- HYPOTHESES.a2tea()[[hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server]]]@expanded_OGs[[input$select_HOG_server]]@add_OG_analysis
    n_all_sets <- length(choice_OG)
    #take max set for choice & associate with (H)OG or singleton
    magHaOGs <- as.vector(choice_OG[[n_all_sets]]@genes$value)
    magHaOGs
  })
  
  
  # the complete MSA for this OG (incl. pot. the chosen amount of additional closest OGs) - as_vector(tibble)
  msa <- eventReactive(xx_changes(), {
    req(input$select_HOG_server)
    req(input$add_OGs)
    HYPOTHESES.a2tea()[[hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server]]]@expanded_OGs[[input$select_HOG_server]]@add_OG_analysis[[input$add_OGs+1]]@msa    
  }) 
  
  
  ##for MSA solo plot!
  # the complete MSA for this OG (incl. pot. the chosen amount of additional closest OGs) - as_vector(tibble)
  msa_solo <- eventReactive(input$msa_solo_button, {#xx_changes(), {
    req(input$select_HOG_server)
    req(input$add_OGs)
    HYPOTHESES.a2tea()[[hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server]]]@expanded_OGs[[input$select_HOG_server]]@add_OG_analysis[[input$add_OGs+1]]@msa    
  })
  
  
  #### rendering HOG list with interesting values tea, cafe, expansion, etc.
  # advanced filtering is enabled
  output$hypothesis_HOG_level_table <- renderDT(server = TRUE, {
    datatable(
      hypothesis_HOG_level_list() %>% dplyr::rename(OG = HOG),
      filter = list(position = 'top', clear = FALSE),
      
      extensions = "Buttons",
      options = list(
        dom = 'lBfrtip',
        lengthMenu = list(c(10, 25, 50, 100, 200), c('10', '25', '50', '100', '200')),
        pageLength = 10,
        buttons = list(
          'copy', 'print',
          list(extend = "csv",
               fieldSeparator = "\t",
               extension = ".tsv",
               text = "Download Current Page", filename = "og_table",
               exportOptions = list(
                 modifier = list(page = "current")
               )
          )
        ),
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      rownames= FALSE
    )
  })
  
  
  ### renderUI elements
  output$add_blast_hit_control <- renderUI({
    req(input$select_Hypothesis_server)
    req(input$select_HOG_server)
    
    sliderInput(inputId = "add_blast_hits",
                label = "Additional Blast Hits", 
                value = 0, 
                min = 0,
                step = 1,
                max = num_max_extended_blast_hits()
    )        
  })
  
  
  #adding additional closest OGs analysis
  output$add_OGs_control <- renderUI({
    req(input$select_Hypothesis_server)
    req(input$select_HOG_server)

    sliderInput(inputId = "add_OGs",
                label = "Additional closest OGs", 
                value = 0, 
                min = 0,
                step = 1,
                #we substract 1 since the #1 set is just the exp. OG itself
                max = num_max_additional_OGs()-1
    )        
  })
  
  output$tea_log2FC_n_break_options <- renderUI({
    req(input$select_Hypothesis_server)
    req(input$select_HOG_server)
    req(input$log2FC_layer == TRUE)
    
    sliderInput(inputId = "tea_log2FC_n_break_choice",
                label = "log2FC - Breaks in Axis", 
                value = 3, 
                min = 0,
                step = 1,
                max = 15
    )
  })
  
  output$tea_log2FC_offset_options <- renderUI({
    req(input$select_Hypothesis_server)
    req(input$select_HOG_server)
    req(input$log2FC_layer == TRUE)
    
    sliderInput(inputId = "tea_log2FC_offset_choice",
                label = "log2FC - Offset", 
                value = 0.8, 
                min = 0,
                step = 0.05,
                max = 5
    )
  })
  
  
  output$tea_log2FC_pwidth_options <- renderUI({
    req(input$select_Hypothesis_server)
    req(input$select_HOG_server)
    req(input$log2FC_layer == TRUE)
    
    sliderInput(inputId = "tea_log2FC_pwidth_choice",
                label = "log2FC - layer width", 
                value = 0.2,
                min = 0.1,
                step = 0.1,
                max = 5
    )
  })
  
  output$tea_log2FC_axis_text_options <- renderUI({
    req(input$select_Hypothesis_server)
    req(input$select_HOG_server)
    req(input$log2FC_layer == TRUE)
    
    sliderInput(inputId = "tea_log2FC_axis_text_choice",
                label = "log2FC - Size of Axis Text", 
                value = 4,
                min = 1,
                step = 1,
                max = 8
    )
  })
  
  output$tea_log2FC_vjust_options <- renderUI({
    req(input$select_Hypothesis_server)
    req(input$select_HOG_server)
    req(input$log2FC_layer == TRUE)
    
    sliderInput(inputId = "tea_log2FC_vjust_choice",
                label = "log2FC - Axis Text dist. adjust", 
                value = 1,
                min = 0.1,
                step = 0.1,
                max = 5
    )
  })    
  
  output$tea_log2FC_linetype_options <- renderUI({
    req(input$select_Hypothesis_server)
    req(input$select_HOG_server)
    req(input$log2FC_layer == TRUE)
    
    sliderInput(inputId = "tea_log2FC_linetype_choice",
                label = "log2FC - Linetype", 
                value = 1,
                min = 0,
                step = 1,
                max = 6
    )
  })
  
  
  #MSA alignment - max value (width of alignment)
  output$tea_msa_width_options <- renderUI({
    req(input$select_Hypothesis_server)
    req(input$select_HOG_server)
    req(input$msa_layer == TRUE)
    req(msa())
    
    sliderInput(inputId = "tea_plot_msa_range_pos", 
                label = "Range of MSA:", 
                value = c(1, max(tidy_msa(msa())["position"])), 
                min = 1, 
                max = max(tidy_msa(msa())["position"]))
    
  }) 
  
  #render options for MSA layer in TEA plots conditionally
  #msa visual choices
  output$tea_msa_offset_options <- renderUI({
    req(input$select_Hypothesis_server)
    req(input$select_HOG_server)
    req(input$msa_layer == TRUE)
    req(msa())
    
    sliderInput(inputId = "tea_msa_offset_choice",
                label = "MSA layer offset", 
                value = 0.8, 
                min = 0,
                step = 0.05,
                max = 10
    )
  }) 
  
  output$tea_msa_pwidth_options <- renderUI({
    req(input$select_Hypothesis_server)
    req(input$select_HOG_server)
    req(input$msa_layer == TRUE)
    req(msa())
    
    sliderInput(inputId = "tea_msa_pwidth_choice",
                label = "MSA layer width", 
                value = 0.5,
                min = 0.1,
                step = 0.1,
                max = 10
    )
  })
  
  
  ### The Tree/expansion + annotation plot ###
  output$treePlot <- renderPlot(

    height = function() input$tree_height,
    {
      req(input$select_Hypothesis_server)
      req(input$tree_plot_export_choice)
      
      info <- HOG_DE.a2tea()                           
      
      tree <- tree_choice_add_OGs() 
      
      # rlang package has the function is_empty()
      # with it we can test for "character(0)"
      el <- list()
      
      for (i in 1:length(tree$tip.label)) {
        if (is_empty(info$HOG[info$gene == tree$tip.label[i]]) == FALSE) {
          first <- info$HOG[info$gene == tree$tip.label[i]]
          el <- c(el, list(filler = tree$tip.label[i]))
          names(el)[i] <- first
        }
        else {
          el <- c(el, list(singleton = tree$tip.label[i]))      
        }
      }
      
      HOG_vec_table <- info %>%
        filter(gene %in% max_all_genes_HOG_add_OGs()) %>%
        arrange(match(gene, max_all_genes_HOG_add_OGs())) %>%
        #number singletons so as to not to get confused in the index numbers (each singleton gets it's own color)
        mutate(HOG = case_when(
          HOG == "singleton" ~ paste("singleton", sep = "-", cumsum(.$HOG=="singleton")),
          TRUE ~ as.character(HOG)
        )
        ) 
      
      HOG_vec <- HOG_vec_table %>%
        #reduce OG/singleton column vector to unique set in order of first appearance -> colors for factor
        distinct(HOG) %>%
        pull()
      
      #by subsetting based on the number of sets we actually plot the legend is reduced to only
      #those OGs that the user is currently looking at ;D
      #brackets important here..
      HOG_vec <- HOG_vec[1:(input$add_OGs+1)]
      
      #new for add OG analysis - here singletons are individually numbered and colored                      
      for (i in 1:length(el)) {
        if (names(el)[i] == "singleton") {
          curr_singleton_id <- el[[i]]
          correct_singleton_num <- HOG_vec_table %>% filter(gene == curr_singleton_id) %>% pull
          names(el)[i] <- correct_singleton_num
        }
      }
      
      # reduce the list to unique tags
      el_split <- sapply(unique(names(el)), function(x) unname(unlist(el[names(el)==x])), simplify=FALSE)
      
      # length() of el_split gives us the number of distinct HOGs for the current tree
      # we can use this for an ifelse check and subsequent choice of coloring
      # 8 or less HOGs then we use the OkabeIto scale with black
      # https://www.chronicle.com/blogs/profhacker/color-blind-accessible-figures
      #palette_OkabeIto_black <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
      #                            "#0072B2", "#D55E00", "#CC79A7", "#000000")
      # more than 8 - then additional colors are used;
      # generated by P36 = Polychrome::createPalette(36, c("#E69F00", "#56B4E9", "#009E73"))
      # after a couple of tries these are the colors (first 8 are OkabeIto scale):
      cols <- c(
        '#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7','#000000',
        '#D8EE00','#A70D2E','#EEE2BB','#DF9AFD','#163BAF','#E6DEFD','#715600','#164F5D',
        '#91F58A','#F30DFE','#00FAFA','#FF9F99','#FC7516','#FD0D8D','#FE83CD','#840086',
        '#688B26','#225AFC','#8178B0','#0DFDC3','#B4EDDD','#CDBD4B','#88867C','#2A9AFE',
        '#9D0060','#86382A','#FCBAD3','#E556DA'
      )
      # more than 36?
      # the we extend with colorRampPalette
      # colorRampPalette(cols)(number of HOGs higher than 36)
      if (length(el_split) > 36) {
        cols <- colorRampPalette(cols)(length(el_split))
      }                                
      
      # naming the colors in the cols vector based on the current set of genes 
      names(cols) <- levels(factor(HOG_vec, levels = HOG_vec))
      
      
      ## THE TREE
      # 1) basic tree + layout, branch.length, size of branches
      # 2) getting rid of unwanted 0s in the HOG -tree association
      # 3) finishing tree:coloring, scales, legend, margins
      
      # basic tree
      p_OTU <- ggtree(tr = tree,
                      layout = input$tree_layout,
                      branch.length = input$branch_length,
                      size = input$tree_size,
                      inherit.aes = TRUE)
      
      
      # the HOGs can't be assigned in some cases which leads to the introduction of "0"s
      # we need to get rid of them somehow, as this leads to 1 additional HOG (called 0) being introduced
      # this shifts the colouring and leads to all kinds of problems
      # the solution here is to assign to each zero occurence the HOG of the next element (that is not a zero)
      # with this we can emphasize:
      # the expansion aspect as the 0 element is colored like the next element (further away from root node);looks nice
      gOTU <- groupOTU(p_OTU, el_split, 'HOG', overlap = "origin")
      
      # the whole thing can be written as a neat function;
      # input "OTU" = gOTU$data$HOG
      
      zero2HOG <- function(OTU) {
        # get index positions of all 0 occurences
        zero_index_pos <- OTU %>% `==`(0) %>% which()
        # create empty vector; length is length of just determined # of zero positions
        all_zero_index_pos_shift <- vector(length = length(zero_index_pos))
        
        # we iterate over the vector of zero positions
        for (i in 1:length(zero_index_pos)){
          # assign current zero index position to a "fresh" vector
          zero_index_pos_shift <- zero_index_pos[i]
          # while this position is referring to a zero keep doing the while loop;
          # stop once we don't find a 0 at checked position (so either a HOG or singleton)
          while (OTU[zero_index_pos_shift] == 0) {
            # add 1; (shift position one to the right)
            zero_index_pos_shift <- zero_index_pos_shift + 1
          }
          # add non zero position to ith element in result vector
          all_zero_index_pos_shift[i] <- zero_index_pos_shift
        }
        return(all_zero_index_pos_shift)
      }
      
      # only call function if we find the zeros in the HOG/singleton vector
      if (0 %in% gOTU$data$HOG) {
        # calling function and checking output
        zero_alt_HOG_pos <- zero2HOG(gOTU$data$HOG)
        # we can use the shifted vector to get the HOG elements at these positions 
        zero_alt_HOG = gOTU$data$HOG[zero_alt_HOG_pos]
        ## replace all occurences of 0 with the HOG of the next, non-zero element
        gOTU$data$HOG <- replace(gOTU$data$HOG, gOTU$data$HOG %in% 0, zero_alt_HOG)
      }
      
      # Finishing up the tree
      p <- gOTU + 
        #     p <-  groupOTU(p_OTU, el_split, 'HOG', overlap = "origin") +
        aes(color=HOG, inherit.aes=FALSE) +
        geom_tiplab(aes(color=HOG), show.legend = TRUE) +
        hexpand(input$ortho_tree_hexpand, direction = 1) +
        hexpand(input$ortho_tree_hexpand, direction = -1) +
        scale_color_manual(values=cols[1:length(HOG_vec)], 
                           breaks = levels(factor(HOG_vec, levels = HOG_vec)),
                           # guides solves ugly legend design
                           # helpful link: https://aosmith.rbind.io/2020/07/09/ggplot2-override-aes/
                           guide = guide_legend(keywidth=0.5,
                                                keyheight=0.5,
                                                order=1,
                                                override.aes = list(size = 10))) +
        #                              theme(plot.margin=margin(6, 120, 6, 6),  
        theme(legend.position = input$tree_plot_legend_pos)
      
      ### Additional Options                  
      ## adding additional layers to the tree
      # logfold change + sign. stars
      tree_dge <- info %>%
        filter(gene %in% tree$tip.label) %>%
        replace_na(list(log2FoldChange = 0)) %>% dplyr::select(gene, log2FoldChange, significant)
      
      tree_dge_only_sig <- info %>%
        filter(gene %in% tree$tip.label) %>%
        replace_na(list(log2FoldChange = 0)) %>% dplyr::select(gene, log2FoldChange, significant) %>%
        filter(significant == "yes")
      
      #get max pos/neg log fold change - for automatic breaks in y-axis scale
      max_lfc <- tree_dge %>% 
        mutate(arr = abs(log2FoldChange)) %>% #get absolute log2FoldChange
        arrange(desc(arr)) %>% #arrange according to decreasing absolute log2FoldChanges
        slice_head(n = 1) %>%
        pull(arr) %>% floor() 
      
      
      if (input$log2FC_layer == TRUE & input$tea_log2FC_add_sign_stars_choice == FALSE) { 
        req(input$log2FC_layer == TRUE)
        req(input$tea_log2FC_add_sign_stars_choice == FALSE)
        req(input$tea_log2FC_n_break_choice)
        req(input$tea_log2FC_offset_choice)
        req(input$tea_log2FC_pwidth_choice)
        req(input$tea_log2FC_axis_text_choice)
        req(input$tea_log2FC_vjust_choice)
        req(input$tea_log2FC_linetype_choice)
 
        p <- p + 
          new_scale_fill() + 
          geom_fruit(
            data=tree_dge,
            geom = geom_bar, 
            mapping = aes(
              y=gene,
              x=log2FoldChange
            ),
            orientation="y",
            stat="identity",
            fill="navy",
            colour = "navy",
            alpha=.7,
            inherit.aes=FALSE,
            offset = input$tea_log2FC_offset_choice,
            pwidth = input$tea_log2FC_pwidth_choice,
            axis.params = 
              list(
                axis = "x", 
                text.size = !!input$tea_log2FC_axis_text_choice, 
                # nbreak - variable or function accepted here?
                nbreak = !!input$tea_log2FC_n_break_choice, #input$nbreak_choice, 
                line.size = 0.5, 
                vjust = !!input$tea_log2FC_vjust_choice,
                line.color = "black",
                inherit.aes = FALSE,
              ),
            grid.params =
              list(
                size = 0.5,
                linetype = !!input$tea_log2FC_linetype_choice
              )
          )           
      }
      
      #log2FC layer + significance stars
      if (input$log2FC_layer == TRUE & input$tea_log2FC_add_sign_stars_choice == TRUE) {    
        req(input$log2FC_layer == TRUE)
        req(input$tea_log2FC_add_sign_stars_choice == TRUE)
        req(input$tea_log2FC_n_break_choice)
        req(input$tea_log2FC_offset_choice)
        req(input$tea_log2FC_pwidth_choice)
        req(input$tea_log2FC_axis_text_choice)
        req(input$tea_log2FC_vjust_choice)
        req(input$tea_log2FC_linetype_choice)
        
        p <- p + 
          new_scale_fill() +
          geom_fruit_list(
            geom_fruit(
              data=tree_dge,
              geom = geom_bar, 
              mapping = aes(
                y=gene,
                x=log2FoldChange
              ),
              orientation="y",
              stat="identity",
              fill="navy",
              colour = "navy",
              alpha=.7,
              inherit.aes=FALSE,
              offset = input$tea_log2FC_offset_choice,
              pwidth = input$tea_log2FC_pwidth_choice,
              axis.params = 
                list(
                  axis = "x", 
                  text.size = !!input$tea_log2FC_axis_text_choice, 
                  # nbreak - variable or function accepted here?
                  nbreak = !!input$tea_log2FC_n_break_choice, #input$nbreak_choice, 
                  line.size = 0.5, 
                  vjust = !!input$tea_log2FC_vjust_choice,
                  line.color = "black",
                  inherit.aes = FALSE,
                ),
              grid.params =
                list(
                  size = 0.5, #color = "navy"
                  linetype = !!input$tea_log2FC_linetype_choice
                )
            ),
            new_scale_fill(), # To initialize fill scale.
            geom_fruit(
              data=tree_dge_only_sig,
              geom = geom_star,
              mapping = aes(y=gene, x=log2FoldChange, fill=significant),
              offset = input$tea_log2FC_offset_choice,
              size = 4,
              color = "orange",
              starstroke = 1
            )
          ) + 
          new_scale_fill() # To initialize fill scale.
      }
      
      
      #adding a MSA layer
      if (input$msa_layer == TRUE) {
        req(input$msa_layer == TRUE)
        req(msa())
        req(input$tea_plot_msa_range_pos)
        req(input$tea_msa_offset_choice)
        req(input$tea_msa_pwidth_choice)
        
        msa_aped <- ape::as.AAbin(msa())
        p <- msaplot(p,
                     msa_aped, 
                     window=c(input$tea_plot_msa_range_pos[1], input$tea_plot_msa_range_pos[2]),
                     offset = input$tea_msa_offset_choice,
                     width = input$tea_msa_pwidth_choice,
                     bg_line = TRUE
        ) 
        
      }                   
      
      #write plot to file for export
      ggsave(filename = paste0(session_dir, "/",
                               "hypothesis_",
                               hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server],
                               "_",
                               input$select_HOG_server, 
                               ".", 
                               input$tree_plot_export_choice), 
             plot = p, 
             device = input$tree_plot_export_choice,
             limitsize = FALSE)
      
      # return the final tree plot - with/without additional layers    
      return(p)         
      
    })
  
  #download option for tree plot
  output$tree_plot_download <- downloadHandler(
    filename = function() {
      paste0(session_dir, "/",
             "hypothesis_",
             hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server],
             "_",
             input$select_HOG_server, 
             ".", 
             input$tree_plot_export_choice)
    },
    content = function(file) {
      file.copy(paste0(session_dir, "/",
                       "hypothesis_",
                       hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server],
                       "_",
                       input$select_HOG_server, 
                       ".", 
                       input$tree_plot_export_choice), 
                file, 
                overwrite=TRUE)
    }
  )
  
  # render table displaying HOG & extended BLAST hits - BLAST scores
  # advanced filtering is enabled
  output$ortho_tree_table <- renderDT(server = TRUE, {
    datatable(blast_hits_table(), 
              filter = list(position = 'top', clear = FALSE),
              
              extensions = "Buttons",
              options = list(
                dom = 'lBfrtip',
                lengthMenu = list(c(10, 25, 50, 100, 200), c('10', '25', '50', '100', '200')),
                pageLength = 10,
                buttons = list(
                  'copy', 'print',
                  list(extend = "csv",
                       fieldSeparator = "\t",
                       extension = ".tsv",
                       text = "Download Current Page", filename = paste0(
                         input$select_HOG_server,
                         "_blast_hits"),
                       exportOptions = list(
                         modifier = list(page = "current")
                       )
                  )
                ),
                columnDefs = list(list(className = 'dt-center', targets = "_all"))
              ),
              selection=list(mode="multiple", 
                             target="cell",
                             #allow only to select cell from gene and OG columns
                             selectable = rbind(cbind(1:nrow(blast_hits_table()), rep(1, nrow(blast_hits_table()))),
                                                cbind(1:nrow(blast_hits_table()), rep(3, nrow(blast_hits_table()))))
              ),
              rownames= FALSE
    )
  })
  
  # render MSA plot for currently selected set of genes (expanded HOG + further BLAST hits)
  
  output$msa_solo_plot <- renderPlot(
    height = function() input$msa_plot_height,
    width = function() input$msa_plot_width,
    {
      
      req(input$msa_solo_button)
      
      msa_set_solo <- tidy_msa(msa_solo(), input$msa_solo_range_pos[1], input$msa_solo_range_pos[2])
      
      msa_plot <- ggplot() + 
        geom_msa(data = msa_set_solo,  
                 color = input$msa_solo_color_scheme, 
                 seq_name = TRUE) +
        theme(
          panel.border = element_blank(), 
          panel.background = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.title=element_blank(),
          axis.text.y = element_text(margin = margin(l = 0, r = input$msa_solo_id_labels_margin, t = 0, b = 0), size = input$msa_solo_id_labels_size),
          axis.text.x = element_text(size = input$msa_solo_pos_labels_size),
          axis.ticks=element_blank()
        )          
      
      #optionally add conservation plot below
      if (input$msa_plot_add_conservation_choice == TRUE) {
        req(input$msa_plot_add_conservation_choice)
        
        msa_plot <- msa_plot + geom_msaBar()
      }
      
      #write plot to file for export
      ggsave(filename = paste0(session_dir, "/",
                               "msa_solo_plot.", input$msa_solo_plot_export_choice), 
             plot = as.grob(msa_plot),
             device = input$msa_solo_plot_export_choice,
             limitsize = FALSE)
      
      return(msa_plot)
      
    }) %>% bindEvent(input$msa_solo_button)
  
  
  #download option for set circ plot
  output$msa_solo_plot_download <- downloadHandler(
    filename = function() {
      paste0(session_dir, "/",
             "msa_solo_plot.", input$msa_solo_plot_export_choice)
    },
    content = function(file) {
      file.copy(paste0(session_dir, "/",
                       "msa_solo_plot.", input$msa_solo_plot_export_choice), 
                file,
                overwrite=TRUE)
    }
  )                                    
  
  ########################
  # GO term analyses
  ########################
  
  #quick solution to generate an all OGs set (no singletons! for now)
  
  all_set <- eventReactive(SFA(), { 
    bind_rows(SFA()) %>%
      dplyr::select(HOG, `Gene-Ontology-Term`) %>%
      group_by(HOG) %>%
      #remove redundancy while collapsing with ", "
      summarise(`Gene-Ontology-Term` = paste(unique(`Gene-Ontology-Term`), collapse=", ")) %>%
      #we should remove all rows that have NO GO terms associated
      #and remove all NAs from the Gene-Ontology-Term column if they have at least 1 GO term
      #remove inline NAs in Gene-Ontology-Term column
      mutate(
        `Gene-Ontology-Term` = str_remove_all(`Gene-Ontology-Term`, "NA, |, NA") 
      ) %>% 
      #remove lines with only NA in Gene-Ontology-Term column
      filter(!str_detect(`Gene-Ontology-Term`, 'NA'))
  })
  
  
  A2TEA_GO_list <- eventReactive(input$go_start_button, {
    req(input$go_start_button)
    req(all_set())
    req(input$go_ontology_choice)
    req(input$go_expansion_choice)
    req(input$go_deg_choice)
    req(input$go_n_deg_choice)
    req(input$go_algo_choice)
    req(input$go_top_n_nodes_choice)
    req(input$go_sig_filter_choice)
    
    #define "universe" of analysis -> background set of all OGs + GO annotation
    
    universe <- all_set() %>%
      mutate(universe_list = setNames(.[["Gene-Ontology-Term"]], pull(.["HOG"]))) %>%
      pull(universe_list)
    
    conserved_OGs <- HOG_level_list()$all_species_overview %>%
      #at least 1 gene from every species thus "conserved"
      filter_at(vars(c(ends_with("_total"))), all_vars(. > 0)) %>%
      pull(HOG)
    
    #subsetting the all universe to only the conserved OGs
    universe <- universe[names(universe) %in% conserved_OGs]
    
    #whatever we compare to the universe it is always a subset of it
    #we can start with an table for the chosen hypothesis -> reduced to the conserved set of OGs!
    #this is then reduced to the set of interest based on the parameters
    if (!is.na(hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server])) {
      
      int_set_df <- HOG_level_list()[[hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server]]] %>% filter(HOG %in% conserved_OGs)
      
      #allow filtering of interesting set to only OGs that show expansion
      #requires hypothesis
      if (input$go_expansion_choice=="yes") {
        
        int_set_df <- int_set_df %>%
          filter(expansion == "yes")      
      }
      
      #getting set of expanded species as vector from hypotheses table for a particular hypothesis
      #"hypotheses" is part of A2TEA.RData object
      exp_species <- hypotheses_tsv() %>%
        filter(hypothesis == !!hypotheses_tsv()$hypothesis[hypotheses_tsv()$name==input$select_Hypothesis_server]) %>%
        pull(expanded_in) %>%
        str_split(., ";") %>%
        unlist()
      exp_species_sigDE <- paste0(exp_species, "_sigDE")
      exp_species_total <- paste0(exp_species, "_total")
      
      
      if (input$go_deg_choice=="at least # DEGs in any species" & !is.na(input$go_n_deg_choice)) {
        int_set_df <- int_set_df %>%
          filter_at(vars(c(ends_with("_sigDE"))), any_vars(. >= input$go_n_deg_choice))
      } 
      else if (input$go_deg_choice=="at least # DEGs in ANY expanded species" & !is.na(input$go_n_deg_choice)) {
        #need to have hypothesis information here to know expanded species
        int_set_df <- int_set_df %>%
          #reduce to expanded species for easy filtering
          dplyr::select(HOG, any_of(exp_species_total), any_of(exp_species_sigDE)) %>%
          filter_at(vars(c(ends_with("_sigDE"))), any_vars(. >= input$go_n_deg_choice))
      } 
      else if (input$go_deg_choice=="at least # DEGs in ALL expanded species" & !is.na(input$go_n_deg_choice)) {
        #need to have hypothesis information here to know expanded species
        int_set_df <- int_set_df %>%
          #reduce to expanded species for easy filtering
          dplyr::select(HOG, any_of(exp_species_total), any_of(exp_species_sigDE)) %>%
          filter_at(vars(c(ends_with("_sigDE"))), all_vars(. >= input$go_n_deg_choice))
      } 
      
    }
    
    #extract names of OGs of interesting final set from int_set_df
    int_set <- int_set_df %>% pull(HOG)  
    
    #important to create true gene 2 GO mapping; GO column char vector!
    universe <- universe %>%
      strsplit(split = ", ")
    
    inGenes <- factor(as.integer(names(universe) %in% int_set))
    names(inGenes) <- names(universe)
    
    GOdata <- new("topGOdata", 
                  ontology=input$go_ontology_choice, 
                  allGenes=inGenes,
                  annot=annFUN.gene2GO, 
                  gene2GO=universe)
    
    #Perform Fisher's exact test:
    
    resultFisher <- runTest(GOdata, 
                            algorithm=input$go_algo_choice, 
                            statistic="fisher")
    
    #We can list the top significant results found:
    allRes <- GenTable(GOdata, 
                       classicFisher = resultFisher, 
                       orderBy = "resultFisher", 
                       ranksOf = "classicFisher", 
                       topNodes = input$go_top_n_nodes_choice,
                       numChar=1000)
    allRes <- allRes[allRes$Significant>input$go_sig_filter_choice,]  
    
    #store results in output list - no need to get fancy with S4 here
    A2TEA_GO_list <- list("int_set_df" = int_set_df, 
                          "GOdata" = GOdata,
                          "resultFisher" = resultFisher, 
                          "allRes" = allRes)
    
    return(A2TEA_GO_list)
    
  })                     
  
  output$go_term_output_table <- renderDT(server = FALSE, {
    
    res_df <- A2TEA_GO_list()$allRes
    res_df$GO.ID <- createAmiGOLink(res_df$GO.ID)
    
    DT::datatable(
      res_df,
      filter = list(position = 'top', clear = FALSE),
      
      extensions = "Buttons",
      options = list(
        dom = 'lBfrtip',
        lengthMenu = list(c(10, 25, 50, 100, 200), c('10', '25', '50', '100', '200')),
        pageLength = 10,
        buttons = list(
          'copy', 'print',
          list(extend = "csv",
               fieldSeparator = "\t",
               extension = ".tsv",
               text = "Download Current Page", filename = "GO_term_enrichment_results",
               exportOptions = list(
                 modifier = list(page = "current")
               )
          ),
          list(extend = "csv",
               fieldSeparator = "\t",
               extension = ".tsv", 
               text = "Download Full Results", filename = "GO_term_enrichment_results",
               exportOptions = list(
                 modifier = list(page = "all")
               )
          )
        ),
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      rownames= FALSE,
      escape = FALSE
    )
  })
  
  
  ggdata <- reactive({
    
    req(A2TEA_GO_list())  
    
    goEnrichment <- A2TEA_GO_list()$allRes
    goEnrichment$classicFisher <- as.numeric(goEnrichment$classicFisher)
    
    ntop <- input$go_enrich_N_GOs_choice
    ggdata <- goEnrichment[1:ntop,]
    ggdata$Term <- factor(ggdata$Term, levels = rev(ggdata$Term)) # fixes order
    
    return(ggdata)
  })                              
  
  output$go_enrich_plot <- renderPlot(
    height = function() input$go_enrich_plot_height,
    {
      
      req(ggdata())  
      req(input$enrich_plot_choice)
      
      if (input$enrich_plot_choice == "bar") {
        
        p <- ggplot(ggdata(), aes(x=Term, y=Significant, fill = -log10(classicFisher))) + 
          geom_bar(stat = "identity") + 
          coord_flip() + 
          theme_bw(base_size = 24) +
          scale_fill_continuous(low = "#428bca", high = "#F39C12")  
        
        # Adjust theme components
        p <- p + theme(axis.text.x = element_text(colour = "black", vjust = 1, angle = 0, size = 12, face = 'bold'), 
                       axis.text.y = element_text(colour = "black", hjust = 1, angle = 0, size = 12, face = 'bold'),
                       axis.title = element_text(color = "black", margin = margin(10, 5, 0, 0), size = 12, face = 'bold'),
                       axis.title.x = element_text(size = 12, face = 'bold'),
                       axis.title.y = element_text(angle = 90, size = 12, face = 'bold'),#, face = 'bold'),
                       legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
                       legend.text = element_text(size = 14, face = "bold"), # Text size
                       title = element_text(size = 14, face = "bold"),
                       axis.line = element_line(colour = 'black'))
      }
      else if (input$enrich_plot_choice == "dot") {
        
        p <- ggplot(ggdata(),
                    aes(x = Term, y = -log10(classicFisher), size = Annotated, fill = -log10(classicFisher))) +
          expand_limits(y = 1) +
          geom_point(shape = 21) +
          scale_size(range = c(2.5,12.5)) +
          scale_fill_continuous(low = "#428bca", high = "#F39C12") +
          xlab('') + ylab('Enrichment score') +
          labs(
            subtitle = 'Top n terms ordered by Fisher p-value',
            caption = 'Cut-off lines drawn at equivalents of p=0.05, p=0.01, p=0.001') +
          geom_hline(yintercept = c(-log10(0.05), -log10(0.01), -log10(0.001)),
                     linetype = c("dotted", "longdash", "solid"),
                     colour = c("black", "black", "black"),
                     size = c(0.5, 1.5, 3)) +
          theme_bw(base_size = 24) +
          theme(
            legend.position = 'right',
            legend.background = element_rect(),
            plot.title = element_text(angle = 0, size = 16, face = 'bold', vjust = 1),
            plot.subtitle = element_text(angle = 0, size = 14, face = 'bold', vjust = 1),
            plot.caption = element_text(angle = 0, size = 12, face = 'bold', vjust = 1),
            axis.text.x = element_text(angle = 0, size = 12, face = 'bold', hjust = 1.10),
            axis.text.y = element_text(angle = 0, size = 12, face = 'bold', vjust = 0.5),
            axis.title = element_text(size = 12, face = 'bold'),
            axis.title.x = element_text(size = 12, face = 'bold'),
            axis.title.y = element_text(size = 12, face = 'bold'),
            axis.line = element_line(colour = 'black'),
            #Legend
            legend.key = element_blank(), # removes the border
            legend.key.size = unit(1, "cm"), # Sets overall area/size of the legend
            legend.text = element_text(size = 14, face = "bold"), # Text size
            title = element_text(size = 14, face = "bold")) +
          coord_flip()     
        
      }
      
      #write plot to file for export
      ggsave(filename = paste0(session_dir, "/",
                               "go_enrich_plot.", input$go_enrich_plot_export_choice), 
             plot = as.grob(p),
             device = input$go_enrich_plot_export_choice,
             limitsize = FALSE)
      
      return(p)
      
    })                               
  
  #download option for go graph
  output$go_enrich_plot_download <- downloadHandler(
    filename = function() {
      paste0(session_dir, "/",
             "go_enrich_plot.", input$go_enrich_plot_export_choice)
    },
    content = function(file) {
      file.copy(paste0(session_dir, "/",
                       "go_enrich_plot.", input$go_enrich_plot_export_choice), 
                file,
                overwrite=TRUE)
    }
  )
  
  # a dropdown list of the expanded HOGs of the current hypothesis
  # first create reactive list of all GO terms that were generated for the current GO analysis 
  sig_GO_terms <- reactive({
    
    req(A2TEA_GO_list())
    vars <- as.list(A2TEA_GO_list()$allRes$GO.ID)
    return(vars)    
  })
  
  output$select_GO_term = renderUI({     
    selectInput(inputId = 'select_go_term_server', label = 'Select GO term', choices = sig_GO_terms())
  })
  
  
  # Getting sig. OGs for an enriched GO term
  #with GO term
  #all_sig_OGs_or_genes_with_term_df <- eventReactive(input$go_sig_OGs_or_genes_df_button, {
  all_sig_OGs_or_genes_with_term_df <- reactive({
    
    req(A2TEA_GO_list())
    req(input$select_go_term_server)
    req(input$sig_go_terms_OGs_or_genes_choice)
    
    all_OGs_with_term <- genesInTerm(A2TEA_GO_list()$GOdata, input$select_go_term_server) %>%
      unname() %>%
      unlist()
    
    all_sig_OGs_or_genes_with_term_df <- A2TEA_GO_list()$int_set_df %>%
      filter(HOG %in% all_OGs_with_term)
    
    if (input$sig_go_terms_OGs_or_genes_choice == "Gene level") {
      all_sig_OGs_or_genes_with_term_df <- HOG_DE.a2tea() %>%
        filter(HOG %in% all_sig_OGs_or_genes_with_term_df$HOG)
    }
    
    all_sig_OGs_or_genes_with_term_df <- all_sig_OGs_or_genes_with_term_df %>%
      #add the OG term itself to the table
      mutate(GO.ID = input$select_go_term_server) %>%
      #relocate without second param moves column to first position
      relocate(GO.ID) %>%
      #turn into links
      mutate(GO.ID = createAmiGOLink(.$GO.ID)) %>%
      #rename HOG to OG
      dplyr::rename(OG = HOG)
    
    return(all_sig_OGs_or_genes_with_term_df)
  })
  
  output$sig_OGs_for_term_df <- renderDT(server = FALSE, {
    DT::datatable(
      all_sig_OGs_or_genes_with_term_df(),
      filter = list(position = 'top', clear = FALSE),
      
      extensions = "Buttons",
      options = list(
        dom = 'lBfrtip',
        lengthMenu = list(c(10, 25, 50), c('10', '25', '50')),
        pageLength = 10,
        buttons = list(
          'copy', 'print',
          list(extend = "csv",
               fieldSeparator = "\t",
               extension = ".tsv",
               text = "Download Current Page", filename = paste0("sigs_", input$sig_go_terms_OGs_or_genes_choice, "_", input$select_go_term_server),
               exportOptions = list(
                 modifier = list(page = "current")
               )
          ),
          list(extend = "csv",
               fieldSeparator = "\t",
               extension = ".tsv", 
               text = "Download Full Results", filename = paste0("sigs_", input$sig_go_terms_OGs_or_genes_choice, "_", input$select_go_term_server),
               exportOptions = list(
                 modifier = list(page = "all")
               )
          )
        ),
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      rownames= FALSE,
      escape = FALSE
    )
  })
  
  
  #plotting the GO graph
  output$go_graph <- renderPlot(
    height = function() input$go_graph_height,
    width = function() input$go_graph_width,
    {
      req(A2TEA_GO_list())
      #   req(GO_data_graph())
      req(input$go_graph_export_choice)
      
      #topGO can output a graph of the GO structure of terms
      graph <- showSigOfNodes(A2TEA_GO_list()$GOdata,
                              topGO::score(A2TEA_GO_list()$resultFisher),
                              firstSigNodes = input$go_graph_top_n_nodes_choice, 
                              useInfo = 'all',
                              .NO.CHAR = input$go_graph_n_max_chars_choice,
                              swPlot = FALSE)
      
      plot_graph <- function() {
        par(cex = input$go_graph_cex_choice)
        plot(graph$complete.dag)
      }
      plot_graph()
      
      #write plot to file for export
      ggsave(filename = paste0(session_dir, "/", "go_graph.", input$go_graph_export_choice), 
             plot = as.grob(function() plot_graph()),
             device = input$go_graph_export_choice,
             limitsize = FALSE)
      
    }) %>% bindEvent(input$go_graph_button)          
  
  
  #download option for go graph
  output$go_graph_download <- downloadHandler(
    filename = function() {
      paste0(session_dir, "/", "go_graph.", input$go_graph_export_choice)
    },
    content = function(file) {
      file.copy(paste0(session_dir, "/", "go_graph.", input$go_graph_export_choice), file, overwrite=TRUE)
    }
  )
  
  
  ##################
  # Set analyses
  ##################
  
  set_fun_changes <- reactive({
    paste(input$select_Hypothesis_server, input$n_deg_set_choice)
  })
  
  set_fun <- eventReactive(set_fun_changes(), {
    req(input$select_Hypothesis_server)

    #hypothesis choice by user
    hypothesis <- hypotheses_tsv() %>% 
      filter(name == input$select_Hypothesis_server) %>%
      pull(hypothesis)
    
    #need to define number of min DEGs for filter conditions
    n_deg = input$n_deg_set_choice
    
    #  #hypotheses object needed to define expanded species per hypothesis
    exp_species <- hypotheses_tsv() %>%
      filter(hypothesis == !!hypothesis) %>%
      pull(expanded_in) %>%
      str_split(., ";") %>%
      unlist()
    exp_species_sigDE <- paste0(exp_species, "_sigDE")

    ##################
    ##################
    
    
    #all OGs
    hyper_all_OGs <- HOG_level_list()$all_species_overview #%>%
    #    nrow()
    
    #all_OGs with at least N DEGs from ANY species
    hyper_all_OGs_N_DEGs_ANY_species <- HOG_level_list()$all_species_overview %>%
      #require at least N DEG from any species
      filter_at(vars(c(ends_with("_sigDE"))), any_vars(. >= n_deg))

    #all conserved OGs
    hyper_conserved_OGs <- HOG_level_list()$all_species_overview %>%
      #at least 1 gene from every species thus "conserved"
      filter_at(vars(c(ends_with("_total"))), all_vars(. > 0))

    #all conserved OGs with at least N DEGs from ANY species
    hyper_conserved_OGs_N_DEGs_ANY_species <- HOG_level_list()$all_species_overview %>%
      #at least 1 gene from every species thus "conserved"
      filter_at(vars(c(ends_with("_total"))), all_vars(. > 0)) %>%
      #require at least 1 DEG from any species
      filter_at(vars(c(ends_with("_sigDE"))), any_vars(. >= n_deg))

    #all conserved OGs with at least N DEGs from ALL species
    hyper_conserved_OGs_N_DEGs_ALL_species <- HOG_level_list()$all_species_overview %>%
      #at least 1 gene from every species thus "conserved"
      filter_at(vars(c(ends_with("_total"))), all_vars(. > 0)) %>%
      #require at least 1 DEG from any species
      filter_at(vars(c(ends_with("_sigDE"))), all_vars(. >= n_deg))

    ### Hypothesis specific sets for all (not exclusively conserved) OGs
    #hypothesis - all conserved OGs + hypothesis EXPANSION
    hyper_all_OGs_hypothesis_expansion <- HOG_level_list()[[hypothesis]] %>%
      #adding expansion filter
      filter(expansion == "yes")
    
    #hypothesis - all conserved OGs + hypothesis N DEGs from ANY species
    #at least # DEGs in any species
    hyper_all_OGs_hypothesis_N_DEGs_ANY_species <- HOG_level_list()[[hypothesis]] %>%
      #require at least 1 DEG in any hypothesis species
      filter_at(vars(c(ends_with("_sigDE"))), any_vars(. >= n_deg))
    
    #at least # DEGs in ANY expanded species
    hyper_all_OGs_hypothesis_N_DEGs_ANY_EXPANDED_species <- HOG_level_list()[[hypothesis]] %>%
      #require at least N DEG in any hypothesis from ANY EXPANDED species
      dplyr::select(HOG, all_of(ends_with("_total")), any_of(exp_species_sigDE)) %>%
      filter_at(vars(c(ends_with("_sigDE"))), any_vars(. >= n_deg))

    #at least # DEGs in ALL expanded species
    hyper_all_OGs_hypothesis_N_DEGs_ALL_EXPANDED_species <- HOG_level_list()[[hypothesis]] %>%
      #require at least 1 DEG in any hypothesis from ALL expanded species
      dplyr::select(HOG, all_of(ends_with("_total")), any_of(exp_species_sigDE)) %>%
      filter_at(vars(c(ends_with("_sigDE"))), all_vars(. >= n_deg))
    
    #hypothesis - all conserved OGs + hypothesis EXPANSION + hypothesis DEGs
    #at least # DEGs in any species
    hyper_all_OGs_hypothesis_expansion_N_DEGs_ANY_species <- HOG_level_list()[[hypothesis]] %>%
      #adding expansion filter
      filter(expansion == "yes") %>%
      #require at least 1 DEG in any hypothesis species
      filter_at(vars(c(ends_with("_sigDE"))), any_vars(. >= n_deg))

    #at least # DEGs in ANY expanded species
    hyper_all_OGs_hypothesis_expansion_N_DEGs_ANY_EXPANDED_species <- HOG_level_list()[[hypothesis]] %>%
      #adding expansion filter
      filter(expansion == "yes") %>%
      #require at least N DEG in any hypothesis EXPANDED species
      dplyr::select(HOG, all_of(ends_with("_total")), any_of(exp_species_sigDE)) %>%
      filter_at(vars(c(ends_with("_sigDE"))), any_vars(. >= n_deg))

    #at least # DEGs in ALL expanded species
    hyper_all_OGs_hypothesis_expansion_N_DEGs_ALL_EXPANDED_species <- HOG_level_list()[[hypothesis]] %>%
      #adding expansion filter
      filter(expansion == "yes") %>%
      #require at least 1 DEG in any hypothesis species
      dplyr::select(HOG, all_of(ends_with("_total")), any_of(exp_species_sigDE)) %>%
      filter_at(vars(c(ends_with("_sigDE"))), all_vars(. >= n_deg))
    
    ### Hypothesis specific sets for conserved OGs
    #hypothesis - all conserved OGs + hypothesis EXPANSION
    hyper_conserved_OGs_hypothesis_expansion <- HOG_level_list()[[hypothesis]] %>%
      #at least 1 gene from every species thus "conserved"
      filter_at(vars(c(ends_with("_total"))), all_vars(. > 0)) %>%
      #adding expansion filter
      filter(expansion == "yes")
    
    #hypothesis - all conserved OGs + hypothesis N DEGs from ANY species
    #at least # DEGs in any species
    hyper_conserved_OGs_hypothesis_N_DEGs_ANY_species <- HOG_level_list()[[hypothesis]] %>%
      #at least 1 gene from every species thus "conserved"
      filter_at(vars(c(ends_with("_total"))), all_vars(. > 0)) %>%
      #require at least 1 DEG in any hypothesis species
      filter_at(vars(c(ends_with("_sigDE"))), any_vars(. >= n_deg))

    #at least # DEGs in ANY expanded species
    hyper_conserved_OGs_hypothesis_N_DEGs_ANY_EXPANDED_species <- HOG_level_list()[[hypothesis]] %>%
      #at least 1 gene from every species thus "conserved"
      filter_at(vars(c(ends_with("_total"))), all_vars(. > 0)) %>%
      #require at least N DEG in any hypothesis from ANY EXPANDED species
      dplyr::select(HOG, all_of(ends_with("_total")), any_of(exp_species_sigDE)) %>%
      filter_at(vars(c(ends_with("_sigDE"))), any_vars(. >= n_deg))

    #at least # DEGs in ALL expanded species
    hyper_conserved_OGs_hypothesis_N_DEGs_ALL_EXPANDED_species <- HOG_level_list()[[hypothesis]] %>%
      #at least 1 gene from every species thus "conserved"
      filter_at(vars(c(ends_with("_total"))), all_vars(. > 0)) %>%
      #require at least 1 DEG in any hypothesis from ALL expanded species
      dplyr::select(HOG, all_of(ends_with("_total")), any_of(exp_species_sigDE)) %>%
      filter_at(vars(c(ends_with("_sigDE"))), all_vars(. >= n_deg))
    
    #hypothesis - all conserved OGs + hypothesis EXPANSION + hypothesis DEGs
    #at least # DEGs in any species
    hyper_conserved_OGs_hypothesis_expansion_N_DEGs_ANY_species <- HOG_level_list()[[hypothesis]] %>%
      #at least 1 gene from every species thus "conserved"
      filter_at(vars(c(ends_with("_total"))), all_vars(. > 0)) %>%
      #adding expansion filter
      filter(expansion == "yes") %>%
      #require at least 1 DEG in any hypothesis species
      filter_at(vars(c(ends_with("_sigDE"))), any_vars(. >= n_deg))

    #at least # DEGs in any species
    hyper_conserved_OGs_hypothesis_expansion_N_DEGs_ALL_species <- HOG_level_list()[[hypothesis]] %>%
      #at least 1 gene from every species thus "conserved"
      filter_at(vars(c(ends_with("_total"))), all_vars(. > 0)) %>%
      #adding expansion filter
      filter(expansion == "yes") %>%
      #require at least 1 DEG in any hypothesis species
      filter_at(vars(c(ends_with("_sigDE"))), all_vars(. >= n_deg))

    #at least # DEGs in ANY expanded species
    hyper_conserved_OGs_hypothesis_expansion_N_DEGs_ANY_EXPANDED_species <- HOG_level_list()[[hypothesis]] %>%
      #at least 1 gene from every species thus "conserved"
      filter_at(vars(c(ends_with("_total"))), all_vars(. > 0)) %>%
      #adding expansion filter
      filter(expansion == "yes") %>%
      #require at least N DEG in any hypothesis EXPANDED species
      dplyr::select(HOG, all_of(ends_with("_total")), any_of(exp_species_sigDE)) %>%
      filter_at(vars(c(ends_with("_sigDE"))), any_vars(. >= n_deg))

    #at least # DEGs in ALL expanded species
    hyper_conserved_OGs_hypothesis_expansion_N_DEGs_ALL_EXPANDED_species <- HOG_level_list()[[hypothesis]] %>%
      #at least 1 gene from every species thus "conserved"
      filter_at(vars(c(ends_with("_total"))), all_vars(. > 0)) %>%
      #adding expansion filter
      filter(expansion == "yes") %>%
      #require at least 1 DEG in any hypothesis species
      dplyr::select(HOG, all_of(ends_with("_total")), any_of(exp_species_sigDE)) %>%
      filter_at(vars(c(ends_with("_sigDE"))), all_vars(. >= n_deg))
    
    
    set_fun_nrows <- list(
      #general sets
      "all OGs" = nrow(hyper_all_OGs), 
      "all OGs with # DEGs from ANY species" = nrow(hyper_all_OGs_N_DEGs_ANY_species),
      "conserved OGs" = nrow(hyper_conserved_OGs), 
      "conserved OGs; # DEGs from ANY species" = nrow(hyper_conserved_OGs_N_DEGs_ANY_species),
      "conserved OGs; # DEGs from ALL species" = nrow(hyper_conserved_OGs_N_DEGs_ALL_species),
      #hypothesis specific but not strictly conserved
      "hypothesis; all OGs; EXPANDED" = nrow(hyper_all_OGs_hypothesis_expansion),
      "hypothesis; all OGs; # DEGs from ANY species" = nrow(hyper_all_OGs_hypothesis_N_DEGs_ANY_species),
      "hypothesis; all OGs; # DEGs from ANY EXPANDED species" = nrow(hyper_all_OGs_hypothesis_N_DEGs_ANY_EXPANDED_species),
      "hypothesis; all OGs; # DEGs from ALL EXPANDED species" = nrow(hyper_all_OGs_hypothesis_N_DEGs_ALL_EXPANDED_species),
      "hypothesis; all OGs; EXPANDED; # DEGs from ANY species" = nrow(hyper_all_OGs_hypothesis_expansion_N_DEGs_ANY_species),
      "hypothesis; all OGs; EXPANDED; # DEGs from ANY EXPANDED species" = nrow(hyper_all_OGs_hypothesis_expansion_N_DEGs_ANY_EXPANDED_species),
      "hypothesis; all OGs; EXPANDED; # DEGs from ALL EXPANDED species" = nrow(hyper_all_OGs_hypothesis_expansion_N_DEGs_ALL_EXPANDED_species),
      #hypothesis specific AND conserved
      "hypothesis; conserved OGs; EXPANDED" = nrow(hyper_conserved_OGs_hypothesis_expansion),
      "hypothesis; conserved OGs; # DEGs from ANY species" = nrow(hyper_conserved_OGs_hypothesis_N_DEGs_ANY_species),
      "hypothesis; conserved OGs; # DEGs from ANY EXPANDED species" = nrow(hyper_conserved_OGs_hypothesis_N_DEGs_ANY_EXPANDED_species),
      "hypothesis; conserved OGs; # DEGs from ALL EXPANDED species" = nrow(hyper_conserved_OGs_hypothesis_N_DEGs_ALL_EXPANDED_species),
      "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY species" = nrow(hyper_conserved_OGs_hypothesis_expansion_N_DEGs_ANY_species),
      "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL species" = nrow(hyper_conserved_OGs_hypothesis_expansion_N_DEGs_ALL_species),
      "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY EXPANDED species" = nrow(hyper_conserved_OGs_hypothesis_expansion_N_DEGs_ANY_EXPANDED_species),
      "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL EXPANDED species" = nrow(hyper_conserved_OGs_hypothesis_expansion_N_DEGs_ALL_EXPANDED_species)
    )
    
    set_fun_dfs <- list(
      #general sets
      "all OGs" = hyper_all_OGs, 
      "all OGs with # DEGs from ANY species" = hyper_all_OGs_N_DEGs_ANY_species,
      "conserved OGs" = hyper_conserved_OGs, 
      "conserved OGs; # DEGs from ANY species" = hyper_conserved_OGs_N_DEGs_ANY_species,
      "conserved OGs; # DEGs from ALL species" = hyper_conserved_OGs_N_DEGs_ALL_species,
      #hypothesis specific but not strictly conserved
      "hypothesis; all OGs; EXPANDED" = hyper_all_OGs_hypothesis_expansion,
      "hypothesis; all OGs; # DEGs from ANY species" = hyper_all_OGs_hypothesis_N_DEGs_ANY_species,
      "hypothesis; all OGs; # DEGs from ANY EXPANDED species" = hyper_all_OGs_hypothesis_N_DEGs_ANY_EXPANDED_species,
      "hypothesis; all OGs; # DEGs from ALL EXPANDED species" = hyper_all_OGs_hypothesis_N_DEGs_ALL_EXPANDED_species,
      "hypothesis; all OGs; EXPANDED; # DEGs from ANY species" = hyper_all_OGs_hypothesis_expansion_N_DEGs_ANY_species,
      "hypothesis; all OGs; EXPANDED; # DEGs from ANY EXPANDED species" = hyper_all_OGs_hypothesis_expansion_N_DEGs_ANY_EXPANDED_species,
      "hypothesis; all OGs; EXPANDED; # DEGs from ALL EXPANDED species" = hyper_all_OGs_hypothesis_expansion_N_DEGs_ALL_EXPANDED_species,
      #hypothesis specific AND conserved
      "hypothesis; conserved OGs; EXPANDED" = hyper_conserved_OGs_hypothesis_expansion,
      "hypothesis; conserved OGs; # DEGs from ANY species" = hyper_conserved_OGs_hypothesis_N_DEGs_ANY_species,
      "hypothesis; conserved OGs; # DEGs from ANY EXPANDED species" = hyper_conserved_OGs_hypothesis_N_DEGs_ANY_EXPANDED_species,
      "hypothesis; conserved OGs; # DEGs from ALL EXPANDED species" = hyper_conserved_OGs_hypothesis_N_DEGs_ALL_EXPANDED_species,
      "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY species" = hyper_conserved_OGs_hypothesis_expansion_N_DEGs_ANY_species,
      "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL species" = hyper_conserved_OGs_hypothesis_expansion_N_DEGs_ALL_species,
      "hypothesis; conserved OGs; EXPANDED; # DEGs from ANY EXPANDED species" = hyper_conserved_OGs_hypothesis_expansion_N_DEGs_ANY_EXPANDED_species,
      "hypothesis; conserved OGs; EXPANDED; # DEGs from ALL EXPANDED species" = hyper_conserved_OGs_hypothesis_expansion_N_DEGs_ALL_EXPANDED_species
    )
    
    #return a list with [1] being the num of rows and [2] being the  dfs themselves
    set_fun <- list(
      "set_fun_nrows" = set_fun_nrows,
      "set_fun_dfs" = set_fun_dfs
    )
    
    
    return(set_fun)
    
  })
  
  #success in sample
  output$set_analysis_df <- renderDT(server = FALSE, {
    req(input$n_deg_set_choice)
    req(input$background_set_choice)
    req(input$background_subset_choice)
    req(input$interest_set_choice)
    req(input$interest_subset_choice)
    
    DT::datatable(
      data.frame(
        Set = c("background set",
                "background subset",
                "set of interest",
                "subset of interest"),
        Size = c(set_fun()[[1]][[input$background_set_choice]],
                 set_fun()[[1]][[input$background_subset_choice]],
                 set_fun()[[1]][[input$interest_set_choice]],
                 set_fun()[[1]][[input$interest_subset_choice]])
      ),
      
      extensions = "Buttons",
      options = list(
        dom = 'Brti',
        buttons = list(
          'copy', 'print',
          list(extend = "csv",
               fieldSeparator = "\t",
               extension = ".tsv",
               text = "Download", filename = "set_comparison",
               exportOptions = list(
                 modifier = list(page = "current")
               )
          )
        ),
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      rownames= FALSE
    )
  })
  
  
  phyper_values <- eventReactive(input$calc_set_phyper, {
    req(input$n_deg_set_choice)
    req(input$background_set_choice)
    req(input$background_subset_choice)
    req(input$interest_set_choice)
    req(input$interest_subset_choice)
    
    
    q <- set_fun()[[1]][[input$interest_subset_choice]] #success in sample
    m <- set_fun()[[1]][[input$background_subset_choice]] #success in background
    n <- set_fun()[[1]][[input$background_set_choice]] - set_fun()[[1]][[input$background_subset_choice]] #failure in background
    k <- set_fun()[[1]][[input$interest_set_choice]] #sample size
    
    #perform hypergeometric test for enrichment
    p <- phyper((q-1), m, n, k, lower.tail=F)
    p
    
    phyper_values <- list(
      "q" = q,
      "m" = m,
      "n" = n,
      "k" = k,
      "p" = p
    )
    
  })
  
  output$phyper_header <- renderText({
    "Perform enrichment test (hypergeometric distribution):"
  }) 
  
  output$phyper_code <- renderText({
    
    HTML(paste("phyper('success in sample - 1',", 
               "       'success in background',", 
               "       'failure in background',", 
               "       'sample size',",
               "       lower.tail= FALSE)", 
               sep = "\n"))
  })
  
  output$phyper_choices <- renderText({
    req(input$n_deg_set_choice)
    req(input$background_set_choice)
    req(input$background_subset_choice)
    req(input$interest_set_choice)
    req(input$interest_subset_choice)
    
    HTML(paste(paste0("phyper(", set_fun()[[1]][[input$interest_subset_choice]], " - ", 1, ","), 
               paste0("       ", set_fun()[[1]][[input$background_subset_choice]], ","),
               paste0("       ", set_fun()[[1]][[input$background_set_choice]], " - ", set_fun()[[1]][[input$background_subset_choice]], ","),
               paste0("       ", set_fun()[[1]][[input$interest_set_choice]], ","),
               "       lower.tail= FALSE)", 
               sep = "\n"))
  }) 
  
  output$phyper_value <- renderText({ paste0("p = ", phyper_values()[["p"]])  })
  
  output$set_circ_plot <- renderPlot(
    height = function() input$set_circ_plot_height,
    #width = function() input$set_circ_plot_width,  
    {
      req(input$n_deg_set_choice)
      req(input$background_set_choice)
      req(input$background_subset_choice)
      req(input$interest_set_choice)
      req(input$interest_subset_choice)
      req(input$set_circ_custom_x_axis_label)
      req(input$interest_ratio_centerdist_choice)
      req(input$interest_ratio_clockpos_choice)
      req(input$interest_ratio_color)
      req(input$interest_ratio_size)
      req(input$background_ratio_centerdist_choice)
      req(input$background_ratio_clockpos_choice)
      req(input$background_ratio_color)
      req(input$background_ratio_size)
      
      
      #allow user to set x-axis title themselves?
      #with log x axis values (correspond to number of OGs)
      p <- ggplot() +
        geom_bar(aes(x=log10(set_fun()[[1]][[input$background_set_choice]])/2,
                     y=c(set_fun()[[1]][[input$background_subset_choice]], (set_fun()[[1]][[input$background_set_choice]]-set_fun()[[1]][[input$background_subset_choice]]))),
                 width = log10(set_fun()[[1]][[input$background_set_choice]]), 
                 position="fill", 
                 stat="identity", 
                 fill=c("#F39C12", "#428bca"), 
                 color="white", 
                 alpha=0.5)+
        geom_bar(aes(x=log10(set_fun()[[1]][[input$interest_set_choice]])/2, 
                     y=c(set_fun()[[1]][[input$interest_subset_choice]], set_fun()[[1]][[input$interest_set_choice]]-set_fun()[[1]][[input$interest_subset_choice]])), 
                 width = log10(set_fun()[[1]][[input$interest_set_choice]]), 
                 position="fill", 
                 stat="identity", 
                 fill=c("#F39C12", "#428bca"), 
                 color="white", 
                 alpha=0.6)+
        coord_polar(theta = "y", start=0)+
        scale_y_continuous(name = input$set_circ_custom_x_axis_label, breaks = NULL, expand = c(0,0))+
        scale_x_continuous(name ="Number of OGs, log10-scaling", 
                           breaks = seq(0,0,0), 
                           expand = c(0,0))+
        theme(
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title = element_text(size = 20),
          panel.background = element_rect(fill="transparent", colour=NA),
          panel.grid  = element_blank(),
          axis.title.y = element_text(vjust = -5),
          axis.title.x = element_text(vjust = 5))+
        geom_text(aes(x=input$interest_ratio_centerdist_choice, y=input$interest_ratio_clockpos_choice, 
                      label=paste(set_fun()[[1]][[input$interest_subset_choice]], "/", set_fun()[[1]][[input$interest_set_choice]], sep = " ")), 
                  color=input$interest_ratio_color, 
                  size=input$interest_ratio_size) +
        geom_text(aes(x=input$background_ratio_centerdist_choice, y=input$background_ratio_clockpos_choice, 
                      label=paste(set_fun()[[1]][[input$background_subset_choice]], "/", set_fun()[[1]][[input$background_set_choice]], sep = " ")), 
                  color=input$background_ratio_color, 
                  size=input$background_ratio_size)
      
      #write plot to file for export
      ggsave(filename = paste0(session_dir, "/",
                               "set_circ_plot.", input$set_circ_plot_export_choice), 
             plot = as.grob(p),
             device = input$set_circ_plot_export_choice,
             limitsize = FALSE)
      
      return(p)
      
    }) #%>% bindEvent(input$set_circ_button)   
  
  #download option for set circ plot
  output$set_circ_plot_download <- downloadHandler(
    filename = function() {
      paste0(session_dir, "/",
             "set_circ_plot.", input$set_circ_plot_export_choice)
    },
    content = function(file) {
      file.copy(paste0(session_dir, "/",
                       "set_circ_plot.", input$set_circ_plot_export_choice), 
                file,
                overwrite=TRUE)
    }
  )
  
  #the distribution plot - histogram of OG sizes
  output$set_og_size_plot <- renderPlot(
    height = function() input$set_og_size_plot_height,
    #width = function() input$set_og_size_plot_width,
    {
      
      req(input$n_deg_set_choice)
      req(input$set_og_size_background_choice)
      req(input$set_og_size_interest_choice)
      req(input$set_og_size_plot_alpha_choice)
      req(input$set_og_size_plot_export_choice)
      
      #two sets
      a <- set_fun()[[2]][[input$set_og_size_background_choice]] %>%
        rowwise() %>% 
        mutate(OG_size = sum(c_across(ends_with("total")))) %>%
        dplyr::select(OG_size)
      
      #get sizes of OGs
      b <- set_fun()[[2]][[input$set_og_size_interest_choice]] %>%
        rowwise() %>% 
        mutate(OG_size = sum(c_across(ends_with("total")))) %>%
        dplyr::select(OG_size)
      
      
      combo <- dplyr::bind_rows(list(a,b),
                                .id="Dataset")
      
      combo <- combo %>%
        mutate(Dataset = 
                 case_when(
                   Dataset == 1 ~ input$set_og_size_background_choice,
                   Dataset == 2 ~ input$set_og_size_interest_choice
                 )
        )
      
      p <- ggplot(combo, aes(OG_size, fill = Dataset)) + 
        geom_histogram(position = "identity", alpha = input$set_og_size_plot_alpha_choice, binwidth=1) +
        scale_fill_manual(values = c("#428bca", "#F39C12")) +
        theme_classic() +
        theme(
          legend.position="bottom"
        )  
      
      #write plot to file for export
      ggsave(filename = paste0(session_dir, "/",
                               "set_og_size_plot.", input$set_og_size_plot_export_choice), 
             plot = as.grob(p),
             device = input$set_og_size_plot_export_choice,
             limitsize = FALSE)
      
      return(p)
      
    })
  
  #download option for set circ plot
  output$set_og_size_plot_download <- downloadHandler(
    filename = function() {
      paste0(session_dir, "/",
             "set_og_size_plot.", input$set_og_size_plot_export_choice)
    },
    content = function(file) {
      file.copy(paste0(session_dir, "/",
                       "set_og_size_plot.", input$set_og_size_plot_export_choice), 
                file,
                overwrite=TRUE)
    }
  )
  
  
  ####################
  #bookmarks
  ####################
  
  #initialize reactive bookmarking
  reactive_bookmarks <- reactiveValues()
  
  reactive_bookmarks$genes <- c()
  reactive_bookmarks$ogs <- c()
  
  #general tab bookmark reactives
  reactive_bookmarks$general_deg_table_genes <- c()
  reactive_bookmarks$func_anno_table_genes <- c()
  reactive_bookmarks$general_deg_table_ogs <- c()
  reactive_bookmarks$func_anno_table_ogs <- c()
  
  #tea tab bookmark reactives
  reactive_bookmarks$tea_og_table_ogs <- c()
  reactive_bookmarks$tea_og_tree_view_ogs <- c()
  #reactive_bookmarks$tea_blast_table_qseqid_genes <- c()
  #reactive_bookmarks$tea_blast_table_sseqid_genes <- c()          
  reactive_bookmarks$stable_bookmarks_blast_table_genes_df <- data.frame(gene=character(), OG=character())
  
  #go term enrichment bookmark reactives
  reactive_bookmarks$sig_go_gene_level_df <- data.frame(gene=character(), OG=character())
  reactive_bookmarks$sig_go_og_level_df <- data.frame(OG=character())
  
  #initializing proxies to have access to the tables
  #we use this to deselect cells once the bookmark button is pushed
  #general tab proxies
  deg_table_proxy <- dataTableProxy(outputId = "DEG_table")
  func_anno_table_proxy <- dataTableProxy(outputId = "func_anno_table")        
  
  #tea tab proxies
  og_table_proxy <- dataTableProxy(outputId = "hypothesis_HOG_level_table")
  blast_table_proxy <- dataTableProxy(outputId = "ortho_tree_table")
  
  #set analysis tab proxies
  
  #go enrichment tab proxies
  sig_OGs_for_term_df_proxy <- dataTableProxy(outputId = "sig_OGs_for_term_df")   
  
  #bookmarking functionality via a reactive
  observeEvent(input$bookmark_click, {
    
    sel_genes <- NULL
    sel_ogs <- NULL
    
    sel_deg_table_genes <- NULL
    sel_func_anno_genes <- NULL
    
    sel_deg_table_ogs <- NULL
    sel_func_anno_ogs <- NULL
    
    sel_og_table_ogs <- NULL
    sel_og_tree_view <- NULL
    sel_blast_table_qseqid_genes <- NULL
    sel_blast_table_sseqid_genes <- NULL
    
    sel_go_og_level_ogs <- NULL
    sel_go_gene_level_genes <- NULL
    
    if (input$a2tea_sidebar_menu != "bookmarks") {
      
      if (input$a2tea_sidebar_menu == "general") {
        
        #deg table bookmark choices
        #checking whether no choice was made; if never any click then special case since variable does not exist!
        #the trick is not to use req() to check existence of input but rather isTruthy which does not enforce downstream existence even outside the if else!!
        if (isTruthy(input$DEG_table_cells_selected) != 0 && nrow(input$DEG_table_cells_selected) != 0) {
          
          sel_deg_table <- as.data.frame(input$DEG_table_cells_selected)     
          
          sel_deg_table_genes <- sel_deg_table %>%
            filter(.[[2]] == 1) %>%
            pull(1)
          
          sel_deg_table_ogs <- sel_deg_table %>%
            filter(.[[2]] == 9) %>%
            pull(1)
        }
        
        #functional annotation table bookmark click choices
        #if (nrow(input$func_anno_table_cells_selected) != 0) {
        if (isTruthy(input$func_anno_table_cells_selected) != 0 && nrow(input$func_anno_table_cells_selected) != 0) {
          sel_func_anno_table <- as.data.frame(input$func_anno_table_cells_selected)     
          
          sel_func_anno_ogs <- sel_func_anno_table %>%
            filter(.[[2]] == 1) %>%
            pull(1)
          
          sel_func_anno_genes <- sel_func_anno_table %>%
            filter(.[[2]] == 2) %>%
            pull(1)
        }
        
        ##merge of user choices for general tab
        sel_genes <- unique(c(sel_deg_table_genes, sel_func_anno_genes))
        sel_ogs <- unique(c(sel_deg_table_ogs, sel_func_anno_ogs))
        
      }
      
      else if (input$a2tea_sidebar_menu == "tea") {
        
        if (isTruthy(input$hypothesis_HOG_level_table_rows_selected) != 0 && length(input$hypothesis_HOG_level_table_rows_selected) != 0) {
          
          sel_og_table_ogs <- input$hypothesis_HOG_level_table_rows_selected 
          
        }
        
        #option for user to add currently vied OG easily to bookmarks
        if (input$bookmark_viewed_OG_choice == TRUE) {
          
          sel_og_tree_view <- input$select_HOG_server
          
          reactive_bookmarks$tea_og_tree_view_ogs <- c(reactive_bookmarks$tea_og_tree_view_ogs, sel_og_tree_view)
          
          #once this is done uncheck the checkbox
          updateCheckboxInput(session=session, inputId="bookmark_viewed_OG_choice", value = FALSE)
          
        }
        
        if (isTruthy(input$ortho_tree_table_cells_selected) != 0 && nrow(input$ortho_tree_table_cells_selected) != 0) {
          
          #sel_blast_table_genes <- input$ortho_tree_table_cell_selected
          sel_blast_table <- as.data.frame(input$ortho_tree_table_cells_selected)
          
          sel_blast_table_qseqid_genes <- sel_blast_table %>%
            filter(.[[2]] == 1) %>%
            pull(1)
          
          sel_blast_table_sseqid_genes <- sel_blast_table %>%
            filter(.[[2]] == 3) %>%
            pull(1)
          
          ##blast table - qseqid
          bookmarks_blast_table_qseqid_genes_df <- blast_hits_table()[sel_blast_table_qseqid_genes, ] %>%
            dplyr::select(qseqid_name) %>%
            dplyr::rename(gene = qseqid_name) %>%
            dplyr::left_join(., res_deg_df(), by = "gene") %>%
            dplyr::select(gene, OG)
          
          ##blast table - sseqid
          bookmarks_blast_table_sseqid_genes_df <- blast_hits_table()[sel_blast_table_sseqid_genes, ] %>%
            dplyr::select(sseqid_name) %>%
            dplyr::rename(gene = sseqid_name) %>%
            dplyr::left_join(., res_deg_df(), by = "gene") %>%
            dplyr::select(gene, OG)
          
          bookmarks_blast_table_genes_df <- bind_rows(bookmarks_blast_table_qseqid_genes_df, 
                                                      bookmarks_blast_table_sseqid_genes_df)
          
          
          
          reactive_bookmarks$stable_bookmarks_blast_table_genes_df <- bind_rows(reactive_bookmarks$stable_bookmarks_blast_table_genes_df,
                                                                                bookmarks_blast_table_genes_df)
          
        }
        
        sel_genes <- unique(c(sel_blast_table_qseqid_genes, sel_blast_table_sseqid_genes))
        sel_ogs <- unique(c(sel_og_table_ogs, sel_og_tree_view))
        
      }
      
      #else if (input$a2tea_sidebar_menu == "set_analyses") {
      #            
      #  XXXXXX <- NULL
      # 
      #}
      
      else if (input$a2tea_sidebar_menu == "go_term_analyses") {
        
        if (isTruthy(input$sig_OGs_for_term_df_rows_selected) != 0 && length(input$sig_OGs_for_term_df_rows_selected) != 0) {
          
          if (input$sig_go_terms_OGs_or_genes_choice == "OG level") {
            
            sel_go_og_level_ogs <- input$sig_OGs_for_term_df_rows_selected
            
            bookmarks_sig_go_og_level_df <- all_sig_OGs_or_genes_with_term_df()[sel_go_og_level_ogs, ] %>%
              dplyr::select(OG)
            
            reactive_bookmarks$sig_go_og_level_df <- bind_rows(reactive_bookmarks$sig_go_og_level_df, 
                                                               bookmarks_sig_go_og_level_df)
            
          }
          else if (input$sig_go_terms_OGs_or_genes_choice == "Gene level") {
            
            sel_go_gene_level_genes <- input$sig_OGs_for_term_df_rows_selected
            
            bookmarks_sig_go_gene_level_df <- all_sig_OGs_or_genes_with_term_df()[sel_go_gene_level_genes, ] %>%
              dplyr::select(gene, OG)
            
            reactive_bookmarks$sig_go_gene_level_df <- bind_rows(reactive_bookmarks$sig_go_gene_level_df, 
                                                                 bookmarks_sig_go_gene_level_df)
          }
        }
        
        sel_genes <- sel_go_gene_level_genes
        sel_ogs <- sel_go_og_level_ogs
        
      }    
      
      #check if any user selections were made
      if (length(sel_genes) == 0 && length(sel_ogs) == 0) {
        showNotification("No genes or orthologous groups selected!", type = "warning", duration = 4)
      }
      
      #check for newness of user selections
      if (length(sel_genes) > 0 && all(sel_genes %in% reactive_bookmarks$genes)) {
        showNotification("All selected genes already bookmarked!", type = "warning", duration = 4)
      }
      if (length(sel_ogs) > 0 && all(sel_ogs %in% reactive_bookmarks$ogs)) {
        showNotification("All selected orthologous groups already bookmarked!", type = "warning", duration = 4)
      }
      
      #if here sth. of interest was bookmarked -> give feedback to user
      if (length(sel_genes) > 0 && 
          any(!sel_genes %in% reactive_bookmarks$genes) |
          length(sel_ogs) > 0 && 
          any(!sel_ogs %in% reactive_bookmarks$ogs)
      ) {
        showNotification("Non-redundant selections were bookmarked", duration = 4)
      }    
      
      #general gene reactives   
      reactive_bookmarks$general_deg_table_genes <- unique(c(reactive_bookmarks$general_deg_table_genes, sel_deg_table_genes))
      reactive_bookmarks$func_anno_table_genes <- unique(c(reactive_bookmarks$func_anno_table_genes, sel_func_anno_genes))
      
      #general OG reactives
      reactive_bookmarks$general_deg_table_ogs <- unique(c(reactive_bookmarks$general_deg_table_ogs, sel_deg_table_ogs))
      reactive_bookmarks$func_anno_table_ogs <- unique(c(reactive_bookmarks$func_anno_table_ogs, sel_func_anno_ogs))         
      
      #tea gene reactives   
      #reactive_bookmarks$tea_blast_table_qseqid_genes <- sel_blast_table_qseqid_genes
      #reactive_bookmarks$tea_blast_table_sseqid_genes <- sel_blast_table_sseqid_genes
      
      #tea OG reactives
      reactive_bookmarks$tea_og_table_ogs <- unique(c(reactive_bookmarks$tea_og_table_ogs, sel_og_table_ogs))
      
      #GO term enrichment gene reactives
      #reactive_bookmarks$sig_go_og_level_df <- unique(c(reactive_bookmarks$sig_go_og_level_df, sel_go_og_level_ogs))
      
      #GO term enrichment OG reactives
      #reactive_bookmarks$sig_go_gene_level_df <- unique(c(reactive_bookmarks$sig_go_gene_level_df, sel_go_gene_level_genes))
      
      #summed reactives
      reactive_bookmarks$genes <- c(reactive_bookmarks$general_deg_table_genes,
                                    reactive_bookmarks$func_anno_table_genes,
                                    reactive_bookmarks$tea_blast_table_qseqid_genes,
                                    reactive_bookmarks$tea_blast_table_sseqid_genes)
      reactive_bookmarks$ogs <- c(reactive_bookmarks$general_deg_table_ogs,
                                  reactive_bookmarks$func_anno_table_ogs,
                                  reactive_bookmarks$tea_og_tree_view_ogs,
                                  reactive_bookmarks$tea_og_table_ogs)
      
      #using the proxy we can de-select once the bookmark button is clicked
      #general tab
      DT::selectCells(deg_table_proxy, input$DEG_table_cells_selected[!input$DEG_table_cells_selected %in% sel_deg_table_genes])
      DT::selectCells(deg_table_proxy, input$DEG_table_cells_selected[!input$DEG_table_cells_selected %in% sel_deg_table_ogs])
      DT::selectCells(func_anno_table_proxy, input$func_anno_table_cells_selected[!input$func_anno_table_cells_selected %in% sel_func_anno_genes])
      DT::selectCells(func_anno_table_proxy, input$func_anno_table_cells_selected[!input$func_anno_table_cells_selected %in% sel_func_anno_ogs])  
      #tea tab
      DT::selectRows(og_table_proxy, input$hypothesis_HOG_level_table_rows_selected[!input$hypothesis_HOG_level_table_rows_selected %in% sel_og_table_ogs])
      DT::selectCells(blast_table_proxy, input$ortho_tree_table_cell_selected[!input$ortho_tree_table_cell_selected %in% sel_func_anno_genes])
      #go term enrichment tab
      DT::selectRows(sig_OGs_for_term_df_proxy, input$sig_OGs_for_term_df_rows_selected[!input$sig_OGs_for_term_df_rows_selected %in% sel_go_og_level_ogs])
      DT::selectRows(sig_OGs_for_term_df_proxy, input$sig_OGs_for_term_df_rows_selected[!input$sig_OGs_for_term_df_rows_selected %in% sel_go_gene_level_genes])
    }
    else if (input$a2tea_sidebar_menu == "bookmarks") {
      showNotification("No bookmarking on the bookmarking tab ;D", type = "warning", duration = 4)
    }
    
  })
  
  
  #genes and OGs can appear in many tables and OGs can even be redundant...
  #we have to build reactive tables with redundancies removed
  #the DT rendering is completely seperate and only after this work is done!
  #this way we can also ensure correct counting of genes/OGs for the infoboxes
  reac_bookmark_genes_df <- reactive({
    
    reac_bookmark_genes_df <- data.frame(gene=character(), OG=character())       
    
    bookmarks_general_deg_table_genes_df <- res_deg_df()[reactive_bookmarks$general_deg_table_genes, ] %>%
      dplyr::select(gene, OG)  
    
    bookmarks_general_func_anno_genes_df <- SFA_general()[reactive_bookmarks$func_anno_table_genes, ] %>%
      dplyr::select(gene, OG)
    
    reac_bookmark_genes_df <- bind_rows(reac_bookmark_genes_df, bookmarks_general_deg_table_genes_df)
    reac_bookmark_genes_df <- bind_rows(reac_bookmark_genes_df, bookmarks_general_func_anno_genes_df)
    
    reac_bookmark_genes_df <- bind_rows(reac_bookmark_genes_df, reactive_bookmarks$stable_bookmarks_blast_table_genes_df)
    
    reac_bookmark_genes_df <- bind_rows(reac_bookmark_genes_df, reactive_bookmarks$sig_go_gene_level_df)
    
    #reac_bookmark_genes_df <- bind_rows(reac_bookmark_genes_df, bookmarks_blast_table_qseqid_genes_df)
    #reac_bookmark_genes_df <- bind_rows(reac_bookmark_genes_df, bookmarks_blast_table_sseqid_genes_df)
    
    #remove gene redundancy if same gene from different tables was bookmarked
    reac_bookmark_genes_df <- reac_bookmark_genes_df %>%
      group_by(gene) %>% 
      slice_head(n = 1) %>%
      ungroup()
    
    #add functional annotation info
    reac_bookmark_genes_df <- dplyr::left_join(reac_bookmark_genes_df, SFA_general(), by = "gene") %>%
      dplyr::select(gene, OG.x, `Human-Readable-Description`) %>%
      dplyr::rename(OG = OG.x)
    
    #perform final sorting by gene name
    reac_bookmark_genes_df <- reac_bookmark_genes_df %>%
      arrange(gene)
    
    return(reac_bookmark_genes_df)
    
  })
  
  reac_bookmark_ogs_df <- reactive({
    
    reac_bookmark_ogs_df <- data.frame(OG=character())       
    
    bookmarks_general_deg_table_ogs_df <- res_deg_df()[reactive_bookmarks$general_deg_table_ogs, ] %>%
      dplyr::select(OG)  
    
    bookmarks_general_func_anno_ogs_df <- SFA_general()[reactive_bookmarks$func_anno_table_ogs, ] %>%
      dplyr::select(OG)    
    
    bookmarks_tea_og_table_ogs_df <- hypothesis_HOG_level_list()[reactive_bookmarks$tea_og_table_ogs, ] %>%
      dplyr::rename(OG = HOG) %>%
      dplyr::select(OG)
    
    #with the tree view choice we actually get the OG/s itself returned
    bookmarks_tea_og_tree_view_ogs_df <- data.frame(OG = reactive_bookmarks$tea_og_tree_view_ogs)
    
    reac_bookmark_ogs_df <- bind_rows(reac_bookmark_ogs_df, bookmarks_general_deg_table_ogs_df)
    reac_bookmark_ogs_df <- bind_rows(reac_bookmark_ogs_df, bookmarks_general_func_anno_ogs_df)
    reac_bookmark_ogs_df <- bind_rows(reac_bookmark_ogs_df, bookmarks_tea_og_table_ogs_df)
    reac_bookmark_ogs_df <- bind_rows(reac_bookmark_ogs_df, reactive_bookmarks$sig_go_og_level_df)
    reac_bookmark_ogs_df <- bind_rows(reac_bookmark_ogs_df, bookmarks_tea_og_tree_view_ogs_df)
    
    #remove OG redundancy if same OG from different tables was bookmarked or
    #same OG for different genes are bookmarked in the same table
    reac_bookmark_ogs_df <- reac_bookmark_ogs_df %>%
      group_by(OG) %>% 
      slice_head(n = 1) %>%
      ungroup()
    
    #add functional annotation info
    reac_bookmark_ogs_df <- dplyr::left_join(reac_bookmark_ogs_df, SFA_general_OG_level(), by = "OG") %>%
      dplyr::select(OG, `Human-Readable-Description`)
    
    #remove singleton choices and perform final sorting by OG name
    reac_bookmark_ogs_df <- reac_bookmark_ogs_df %>%
      filter(OG != "singleton") %>%
      arrange(OG)  
    
    return(reac_bookmark_ogs_df)
  })
  
  observeEvent(input$reset_genes_choice, {
    reactive_bookmarks$genes <- c()
    reactive_bookmarks$general_deg_table_genes <- c()
    reactive_bookmarks$func_anno_table_genes <- c()
    reactive_bookmarks$stable_bookmarks_blast_table_genes_df <- data.frame(gene=character(), OG=character())
    reactive_bookmarks$sig_go_gene_level_df <- data.frame(gene=character(), OG=character())
    showNotification("Resetting list of bookmarked genes", duration = 4)
  })
  
  observeEvent(input$reset_ogs_choice, {
    reactive_bookmarks$ogs <- c()
    reactive_bookmarks$general_deg_table_ogs <- c()
    reactive_bookmarks$func_anno_table_ogs <- c()
    reactive_bookmarks$tea_og_tree_view_ogs <- c()
    reactive_bookmarks$tea_og_table_ogs <- c()
    reactive_bookmarks$sig_go_og_level_df <- data.frame(OG=character())
    showNotification("Resetting list of bookmarked orthologous groups", duration = 4)
  })          
  
  ##creating the bookmarking datatables
  #general gene DT table
  #order table by gene name - since ordering by addition is currently not possible
  output$bookmark_genes_df <- DT::renderDT(server = FALSE, {
    
    DT::datatable(
      reac_bookmark_genes_df(),
      
      extensions = "Buttons",
      options = list(
        dom = 'lBfrtip',
        lengthMenu = list(c(10, 25, -1), c('10', '25', 'All')),
        pageLength = 10,
        buttons = list(
          'copy', 'print',
          list(extend = "csv",
               fieldSeparator = "\t",
               extension = ".tsv",
               text = "Download Current Page", filename = "bookmarked_genes",
               exportOptions = list(
                 modifier = list(page = "current")
               )
          ),
          list(extend = "csv",
               fieldSeparator = "\t",
               extension = ".tsv", 
               text = "Download Full Results", filename = "bookmarked_genes",
               exportOptions = list(
                 modifier = list(page = "all")
               )
          )
        ),
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      rownames= FALSE
    )
  })
  
  #a controlbar version of the gene level bookmarks table
  output$cb_bookmark_genes_df <- DT::renderDT(server = FALSE, {
    DT::datatable(
      reac_bookmark_genes_df(),
      options = list(
        dom = 'ltp',
        lengthMenu = list(c(10, 25), c('10', '25')),
        pageLength = 10,
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      rownames= FALSE
    )
  })
  
  #general OG DT table
  output$bookmark_ogs_df <- DT::renderDT(server = FALSE, {
    
    DT::datatable(
      reac_bookmark_ogs_df(),
      
      extensions = "Buttons",
      options = list(
        dom = 'lBfrtip',
        lengthMenu = list(c(10, 25, -1), c('10', '25', 'All')),
        pageLength = 10,
        buttons = list(
          'copy', 'print',
          list(extend = "csv",
               fieldSeparator = "\t",
               extension = ".tsv",
               text = "Download Current Page", filename = "bookmarked_orthologous_groups",
               exportOptions = list(
                 modifier = list(page = "current")
               )
          ),
          list(extend = "csv",
               fieldSeparator = "\t",
               extension = ".tsv", 
               text = "Download Full Results", filename = "bookmarked_orthologous_groups",
               exportOptions = list(
                 modifier = list(page = "all")
               )
          )
        ),
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      rownames= FALSE
    ) 
  })
  
  #a controlbar version of og bookmarks table
  output$cb_bookmark_ogs_df <- DT::renderDT(server = FALSE, {
    DT::datatable(
      reac_bookmark_ogs_df(),
      options = list(
        dom = 'ltp',
        lengthMenu = list(c(10, 25), c('10', '25')),
        pageLength = 10,
        columnDefs = list(list(className = 'dt-center', targets = "_all"))
      ),
      rownames= FALSE
    )
  })
  
  #render the info boxes
  output$infobox_book_genes <- renderInfoBox({
    req(reac_bookmark_genes_df())
    infoBox(
      title = "Bookmarked genes",
      value = nrow(reac_bookmark_genes_df()),
      icon = shiny::icon("book-bookmark"),
      fill = TRUE,
      width = 12
    )
  })  
  
  output$infobox_book_OGs <- renderInfoBox({
    req(reac_bookmark_ogs_df())
    infoBox(
      title = "Bookmarked OGs",
      value = nrow(reac_bookmark_ogs_df()),
      icon = shiny::icon("book-bookmark"),
      fill = TRUE,
      width = 12
    )
  })
  
  
  output$ui_bookmarks <- renderUI({
    shinyjs::hidden(
      div(id = "id_bookmarks_outer",
          tagList(
            fluidRow(width = 12,
                     column(
                       width = 6,
                       infoBoxOutput("infobox_book_genes", width = 6),
                       actionButton("reset_genes_choice", 
                                    class = "bookmark_reset",
                                    label ="",
                                    icon = icon("trash")
                       )
                     ),
                     column(
                       width = 6,
                       infoBoxOutput("infobox_book_OGs", width = 6),
                       actionButton("reset_ogs_choice",
                                    class = "bookmark_reset",
                                    label ="",
                                    icon = icon("trash")    
                       )
                     )
            ),
            fluidRow(width = 12,
                     column(width = 6,
                            box(
                              style = "height: 500px; overflow-y: auto; overflow-x: auto", 
                              background ="yellow",
                              gradient = TRUE,
                              width = NULL,
                              collapsed = FALSE,
                              title = "Bookmarked genes",
                              DTOutput("bookmark_genes_df")
                            )
                     ),
                     column(width = 6,
                            box(
                              style = "height: 500px; overflow-y: auto; overflow-x: auto", 
                              background ="yellow",
                              gradient = TRUE,
                              width = NULL,
                              collapsed = FALSE,
                              title = "Bookmarked orthologous groups",
                              DTOutput("bookmark_ogs_df")
                            )
                     )
            ),
            fluidRow(width = 12,
                     column(width = 4
                     ),
                     column(width = 4,
                            downloadButton(outputId = "button_bookmarks_RData_export", class = "rdata_subset_export",
                                           label = "Generate .RData subset for all bookmarks and related data",
                                           icon = icon("save", lib = "font-awesome"))
                     ),
                     column(width = 4
                     )
            )
          )
      ))
  })
  
  ### subset creation based on chosen bookmarks and after button click          
  #############
  #exporting to subset .RData object containing all info bookmarked items
  #-> for both genes and OGs we get info on OG level since this way nothing is lost
  
  #1 - get all user chosen OGs - extract from bookmark tables (so also parent OGs from bookmarked genes)
  #2 - check if OGs are expanded in any hypothesis (bookmarked OGs may simply be interesting because of GO terms)
  #2.1 - we make it easy for ourselves by simply removing every entry in a copy of HYPOTHESES.a2tea
  #2.2 - since the addtional OGs are only referred to by their genes in additional sets we iterate over the new subset HYPOTHESES.a2tea and get all genes of bookmarked OGs and their addtional OGs
  #3 based on bookmarked OGs and the gene list gained in #2.2 we create a master list of OGs and genes necessary for the subset object
  #3.1 the problem is however that we at this point have (for our purposes) an incomplete set of OGs and an incomplete set of genes
  #3.2 by extracting both the rows of the bookmarked OGs as well as those for the genes of additional OGs and then removing duplicates we can generate a gene level table with all genes/OGs for the subset
  #3.3 this thus corresponds to the HOG_DE.a2tea table of the subset
  #4 - subset objects with list of all OGs based on pre. decisions          
  
  
  #establish # of hypotheses
  num_hypotheses <- reactive({
    length(HYPOTHESES.a2tea())
  })
  
  #master sets of bookmarked OGs
  bookmarked_OGs <- reactive({
    bookmarked_OGs <- unique(c(reac_bookmark_genes_df()$OG, reac_bookmark_ogs_df()$OG))
    return(bookmarked_OGs)
  })
  
  
  subset_HYPOTHESES.a2tea <- reactive({
    
    #create copy of HYPOTHESES.a2tea
    subset_HYPOTHESES.a2tea <- HYPOTHESES.a2tea()
    
    for (n in 1:num_hypotheses()) {
      
      cur_hypo_exp_OGs <- names(subset_HYPOTHESES.a2tea[[n]]@expanded_OGs)
      
      for (exp_OG in cur_hypo_exp_OGs) {
        
        if (exp_OG %in% bookmarked_OGs()) {
        } else {
          subset_HYPOTHESES.a2tea[[n]]@expanded_OGs[[exp_OG]] <- NULL
        }
      }
    }
    return(subset_HYPOTHESES.a2tea)
  })
  
  
  subset_HOG_DE.a2tea <- reactive({
    #always the last set is the biggest containing all genes of exp. OG as well as all addtional OGs
    #this implementation is really not smart - no size pre-allocation and we prob. can get to the largest set somehow and not merge everything...
    #good enough for now - optimize later
    
    subset_all_genes_add_OGs <- c()
    
    for (n in 1:num_hypotheses()) {
      
      cur_hypo_exp_OGs <- names(subset_HYPOTHESES.a2tea()[[n]]@expanded_OGs)
      
      for (exp_OG in cur_hypo_exp_OGs) {
        
        exp_OG_num_og_sets <- length(subset_HYPOTHESES.a2tea()[[n]]@expanded_OGs[[exp_OG]]@add_OG_analysis)
        
        for (set in 1:exp_OG_num_og_sets) {
          
          genes <- subset_HYPOTHESES.a2tea()[[n]]@expanded_OGs[[exp_OG]]@add_OG_analysis[[set]]@genes %>% pull(value)
          subset_all_genes_add_OGs <- unique(c(subset_all_genes_add_OGs, genes))
          
        }
      }
    }
    
    #subset_all_genes_add_OGs
    
    # get subsets for bookmarked OGs
    subset_HOG_DE.a2tea_bookmarked_OGs <- HOG_DE.a2tea() %>%
      filter(HOG %in% bookmarked_OGs())
    
    # get subsets for bookmarked OGs
    subset_HOG_DE.a2tea_genes_add_OGs <- HOG_DE.a2tea() %>%
      filter(gene %in% subset_all_genes_add_OGs)
    
    subset_HOG_DE.a2tea <- dplyr::bind_rows(
      subset_HOG_DE.a2tea_bookmarked_OGs,
      subset_HOG_DE.a2tea_genes_add_OGs
    ) %>%
      distinct()
    
    #this is the new HOG_DE.a2tea
    return(subset_HOG_DE.a2tea)
    
  })
  
  subset_genes <- reactive({
    #gene vector for subset filtering
    subset_genes <- subset_HOG_DE.a2tea()$gene
    return(subset_genes)
  })
  
  subset_OGs <- reactive({
    #og vector for subset filtering
    subset_OGs <- subset_HOG_DE.a2tea() %>% filter(HOG != "singleton") %>% pull(HOG) #remove singletons...
    return(subset_OGs)
  })
  
  subset_HOG_level_list <- reactive({
    
    #getting subset for HOG_level_list
    subset_HOG_level_list <- HOG_level_list()
    
    #+1 since we include the all_species table
    for (n in 1:(num_hypotheses()+1)) {
      
      subset_HOG_level_list[[n]] <- subset_HOG_level_list[[n]] %>%
        filter(HOG %in% subset_OGs())
    }
    
    return(subset_HOG_level_list)
    
  })
  
  #subset pep fasta sequences
  subset_A2TEA.fa.seqs <- reactive({
    
    #getting subset for A2TEA.fa.seqs
    subset_A2TEA.fa.seqs <- A2TEA.fa.seqs()[subset_genes()]
    
    return(subset_A2TEA.fa.seqs)
  })
  
  #getting subset for HOG_level_list
  subset_SFA <- reactive({
    
    subset_SFA <- SFA()
    num_species <- length(subset_SFA)
    
    for (n in 1:num_species) {
      
      subset_SFA[[n]] <- subset_SFA[[n]] %>%
        filter(HOG %in% subset_OGs())
    }
    return(subset_SFA)
  })
  
  #the download button
  output$button_bookmarks_RData_export <- downloadHandler(
    filename = function() {
      "subset_A2TEA_finished.RData"
    }, content = function(file) {
      #write.csv(bookmarks_subset(), file)
      
      #save searches the object named bookmarks_subset(), and there's no such object. 
      #solution:
      hypotheses <- hypotheses_tsv()
      HYPOTHESES.a2tea <- subset_HYPOTHESES.a2tea()
      A2TEA.fa.seqs <- subset_A2TEA.fa.seqs()
      HOG_DE.a2tea <- subset_HOG_DE.a2tea()
      HOG_level_list <- subset_HOG_level_list()
      SFA <- subset_SFA()
      all_speciesTree <- all_speciesTree()
      
      save(hypotheses, 
           HYPOTHESES.a2tea, 
           A2TEA.fa.seqs,
           HOG_DE.a2tea, 
           HOG_level_list,
           SFA,
           all_speciesTree,
           file = file,
           compress = "xz", 
           compression_level = 9)
    }
  )
  
  ####################
  #actions to perform on end of session
  ####################
  
  session$onSessionEnded(function() {
    #remove session dir with intermediary plots, etc.
    print(paste0("Removing files in session dir: ", session_dir))
    unlink(x = session_dir, recursive = TRUE)
  })

}
