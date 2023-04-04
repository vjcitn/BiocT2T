
#' visualize PCA colored by populations for a selected SuperPopulation
#' @importFrom GGally ggpairs
#' @param pcs matrix of principal components, assumed to have colnames PC[1,...] typically PC1-PC10
#' @param superpop character(1)
#' @param selectedPCs numeric(), values typically drawn from 1:10
#' @param dropNAs logical(1) defaults to TRUE
#' @note Population and SuperPopulation definitions are at `https://www.internationalgenome.org/data-portal/population`,
#' mapping from sample to population/superpopulation is available via `BiocHail::path_1kg_annotations`.
#' 
viewStrat1KG = function(pcs, selectedPCs=1:4, superpop="AFR", dropNAs=TRUE) {
## ----getanno------------------------------------------------------------------
  npc = colnames(pcs)
  touse = paste0("PC", selectedPCs)
  stopifnot(all(touse %in% npc))
  anp = BiocHail::path_1kg_annotations()
  md = read.delim(anp, sep="\t")
  validsp = unique(md$SuperPopulation)
  if (!(superpop %in% validsp)) {
   message("valid superpopulations are:")
   print(validsp)
   stop(sprintf("superpop %s not valid", superpop))
   }
  validP = unique(md$Population)
  # library(dplyr) # Available via Depends
  pcdf = data.frame(pcs)
  pcdf$Sample = rownames(pcdf)
  lj = left_join(pcdf, md, by="Sample")
  if (dropNAs) lj = na.omit(lj)
  lj = dplyr::select(lj, all_of(c(touse, "SuperPopulation", "Population")))
  lj = dplyr::filter(lj, SuperPopulation==superpop)
  lj = dplyr::select(lj, -SuperPopulation)
  list(lj=lj) #, pl=ggpairs(lj, aes(colour=Population)))
}

#' multipanel app to view stratification patterns in various populations
#' @import shiny
#' @importFrom GGally ggpairs
#' @export
stratapp = function() {
 data(pc_3kloci, package="BiocT2T")
 data(pc_190kloci, package="BiocT2T")
 data(pc_990kloci, package="BiocT2T")
 pclist = list(k3 = pc_3kloci, k190 = pc_190kloci,
      k990 = pc_990kloci)
 res = c("k3", "k190", "k990")
 names(res) =  c("3k", "190k", "990k")
 superpops = c("EUR", "EAS", "AMR", "SAS", "AFR")
 splpops = list(AFR = c("ACB", "GWD", "ESN", "MSL", "YRI", "LWK", "ASW"), 
    AMR = c("PUR", "CLM", "PEL", "MXL"), EAS = c("CHS", "CDX", 
    "KHV", "CHB", "JPT"), EUR = c("GBR", "FIN", "IBS", "CEU", 
    "TSI"), SAS = c("PJL", "BEB", "STU", "ITU", "GIH"))

 server = function(input, output) {
  getPCs = reactive({
    viewStrat1KG(pclist[[input$res]], input$pcs, superpop=input$superpop)
    })
  output$pairsplot = renderPlot({
   dat = getPCs() #viewStrat1KG(pclist[[input$res]], input$pcs, superpop=input$superpop)
   dpop = dat[[1]]$Population
   dat = dplyr::select(dat[[1]], -Population)
   pairs(dat, pch=19, cex=.4, col=factor(dpop))
   }, height=600L)
  output$ggally = renderPlot({
   dat = getPCs()[[1]]
   GGally::ggpairs(dat, aes(colour=Population))
   }, height=650L)
  output$about = DT::renderDataTable({
   data(igsr_pops)
   igsr_pops |> dplyr::select(-`Data collections`)
   })

  observe({
      if (input$stopBtn > 0) {
        isolate({
          stopApp(returnValue = 0)
        })
      }
    })
  output$desc = renderPrint({ packageDescription("BiocT2T") })

 }

 ui = fluidPage(
  sidebarLayout(
   sidebarPanel(
    helpText("BiocT2T Population Stratification viewer, based on randomly chosen SNP sets on chr17."),
    helpText("See the 'About' tab for information on subpopulation codes used on the GGally tab."),
    radioButtons("superpop", "Superpop", choices=superpops, selected="AFR"),
#    uiOutput("subpop"),
    selectInput("pcs", "PCs to use", choices=1:10, selected=4:7, multiple=TRUE),
    radioButtons("res", "Resolution (# SNPS for PCs)", choices=res, selected="k990"), 
    actionButton("stopBtn", "Stop app"),
    width=2),
   mainPanel(
    tabsetPanel(
     tabPanel("2d",
      plotOutput("pairsplot")
      ),
     tabPanel("GGally",
      helpText("rendering detailed plot takes a moment..."),
      plotOutput("ggally")
      ),
     tabPanel("about", DT::dataTableOutput("about"), verbatimTextOutput("desc"))
     )
    )
   )
  ) # end ui
 shinyApp(ui=ui, server=server)
}
