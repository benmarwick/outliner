
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("EBImage")


library("shiny")
library("EBImage") # >= 4.19.3
library(imager)
library(sf)

ui <- fluidPage(# Application title
  titlePanel("Image outline, chords, and landmarks"),
  
  # Sidebar with a select input for the image
  sidebarLayout(
    sidebarPanel(
      fileInput("image", "Select image"),
      sliderInput(
        "slider5",
        label = "Crop left",
        min = 0,
        max = 1,
        value = 0
      ),
      sliderInput(
        "slider6",
        label = "Crop right",
        min = 0,
        max = 1,
        value = 1
      ),
      sliderInput(
        "slider7",
        label = "Crop top",
        min = 0,
        max = 1,
        value = 0
      ),
      sliderInput(
        "slider8",
        label = "Crop bottom",
        min = 0,
        max = 1,
        value = 1
      ),
      sliderInput(
        "slider11",
        label = "Gamma correction",
        min = 0,
        max = 10,
        value = 1,
        step = 0.1
      ),
      sliderInput(
        "slider9",
        label = "Contrast",
        min = 0,
        max = 2,
        value = 1,
        step = 0.1
      ),
      sliderInput(
        "slider10",
        label = "Brightness",
        min = -5,
        max = 5,
        value = 0,
        step = 0.1)
      ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Interactive browser", displayOutput("widget")),
        tabPanel("B&W raster", displayOutput("bw_raster")),
        tabPanel("Set scale", plotOutput("scale_raster", click = "plot_click")),
        tabPanel("Contours & chords", plotOutput("contour"))
      ),
        fluidRow(    
          downloadButton("downloadDataOutline",
                         label = "Download outline coords"),
          downloadButton("downloadDataPoints",
                         label = "Download point coords"),
          column(
            4,
            scale_factor <- numericInput(
              "scale_factor_number", 
              "How many mm between the two points on the scale bar?", 
              min = 0, max = 100, step = 0.1, value = 0
            ),
            tableOutput("scale_ratio"),
            actionButton("save_scale_factor", 
                         "Save scale factor"),
            uiOutput("saved_value"),
            actionButton("reset_points",
                         "Clear clicked points"),
            tableOutput("click_points")
          )
        )
      )
    )
  )


server <- function(input, output) {
  img <- reactive({
    f <- input$image
    if (is.null(f))
      return(NULL)
    readImage(f$datapath)
  })
  


  # modify the image: cropping
prepped_image <- reactive({
  req(img())
  x <- img()
  colorMode(x) = Grayscale
  x_dim <- dim(x)
  # contrast
  x1 <- x * input$slider9
  # brightness
  x2 <- x1 + input$slider10
  # Gamma Correction
  x3 <- x2 ^ input$slider11
  # Cropping
  x3[(input$slider5 * x_dim[1]):(input$slider6 * x_dim[1]),
           (input$slider7 * x_dim[2]):(input$slider8 * x_dim[2]), ]
  })
  
  # display image unaltered
  output$widget <- renderDisplay({
    req(img())
    EBImage::display(img())
  })
  
  # display image after preparations
  output$bw_raster <- renderDisplay({
    req(prepped_image())
    EBImage::display(prepped_image())
  })
  
  # display image after preparations so we can get the scale factor
  output$scale_raster <- renderPlot({
    req(prepped_image())
    req(point_coords)
    x4 <- prepped_image()
    plot(x4)
    ## show the points at the locations we click
    points(
      point_coords$x,
      point_coords$y,
      pch = 3,
      col = "black",
      cex = 5
    )
    
  })
  
  # compute contour points
  contour_points <- reactive({
    req(prepped_image())
    x5 <- prepped_image()
    
    # Canny edge detection
    canny <- cannyEdges(as.cimg(imageData(x5[,,1])))
    
    # find contours
    oc = ocontour(canny)
    oc[[1]]
  })
  
  # compute max length, max width, and mid-point width
  length_and_widths <- reactive({
    req(contour_points())
    points_ex <- contour_points()
    
    # from https://stackoverflow.com/questions/55854273/how-to-find-the-longest-chord-perpendicular-to-the-maximum-chord-through-a-polyg/55874609#55874609
    
    ## this returns the i,j of the largest elements in matrix `m`
    findmax <- function(m){
      v = which.max(m) - 1
      c(v %% nrow(m)+1, v %/% nrow(m)+1)
    }
    
    ## Return an sf line through a point at an angle of a given length
    pline <- function(pt, angle, length){
      st_linestring(
        cbind(
          pt[1] + c(length,-length) * sin(angle),
          pt[2] + c(length,-length) * cos(angle)
        )
      )
    }
    
    ### find the max length chord across all pairs of vertex points
    max_chord <- function(polygon){
      ## get polygon coordinates
      xy = st_coordinates(polygon)[,1:2]
      
      ## compute the distance matrix and find largest element
      df_dist = as.matrix(dist(xy))
      maxij = findmax(df_dist)
      
      ## those elements define the largest chord
      chord = rbind(
        xy[maxij[1],],
        xy[maxij[2],]
      )
      chord
    }
    
    ## return the line that is the chord at angle perp.angle of length 
    # through any of the polygon vertices
    max_perp_chord <- function(polygon, perp.angle, length){
      ## get polygon vertices
      pts = st_coordinates(polygon)[,c(1,2)]
      ## return the perpendicular lines
      perplines = lapply(1:nrow(pts), function(i){
        ## through the i-th vertex
        xy = pts[i,,drop=FALSE]
        perpline = pline(xy, perp.angle, length)
        ## intersect it with the polygon
        inters = st_intersection(polygon, perpline)
        inters
      }
      )
      
      ## get the vector of intersection lengths, find the largest
      perplengths = unlist(lapply(perplines, st_length))
      longest = which.max(perplengths)
      ## return the longest line
      perplines[[longest]]
      
    }
    
    ## get the coords of the midpoint on the chord
    find_mid_point_on_chord <- function(chord){
      x_mid = (chord[2,1] + chord[1,1])/2
      y_mid = (chord[2,2] + chord[1,2])/2
      
      m <- matrix(c(x_mid, y_mid), ncol = 2)
      
      st_point(m)
      
    }
    
    ### find the max length chord across all pairs of vertex points
    max_chord <- function(polygon){
      ## get polygon coordinates
      xy = st_coordinates(polygon)[,1:2]
      
      ## compute the distance matrix and find largest element
      df_dist = as.matrix(dist(xy))
      maxij = findmax(df_dist)
      
      ## those elements define the largest chord
      chord = rbind(
        xy[maxij[1],],
        xy[maxij[2],]
      )
      chord
    }
    
    ### find the max chord that is perpendicular to the most max chord
    find_max_chord <- function(spolygon, chord=max_chord(spolygon)){
      
      ## Now compute the length and angle of the longest chord
      chord.length = sqrt(diff(chord[,1])^2 + diff(chord[,2])^2)
      chord.theta = atan2(diff(chord[,1]), diff(chord[,2]))
      
      ## The perpendicular is at this angle plus pi/2 radians
      perp = chord.theta + pi/2
      max_perp_chord(spolygon, perp, chord.length)
    }
    
    
    ### find the perpendicular chord at the mid point of the max chord
    max_perp_chord_midpoint <- function(polygon){
      
      chord <- max_chord(polygon)
      
      ## Now compute the length and angle of the longest chord
      chord.length = sqrt(diff(chord[,1])^2 + diff(chord[,2])^2)
      chord.theta = atan2(diff(chord[,1]), diff(chord[,2]))
      perp = chord.theta + pi/2
      
      # compute mid point on chord
      mid_point_on_chord <- find_mid_point_on_chord(chord)
      
      ## get polygon vertices
      pts = st_coordinates(polygon)[,c(1,2)]
      ## return the perpendicular lines
      perplines = lapply(1:nrow(pts), function(i){
        ## through the i-th vertex
        xy = pts[i,,drop=FALSE]
        perpline = pline(xy, perp, chord.length)
        ## intersect it with the polygon
        inters = st_intersection(polygon, perpline)
        inters
      }
      )
      ## drop empties
      idx <- sapply(perplines, function(i) any(is.na(st_dimension(i))))
      perplines <- perplines[!idx]
      
      ## get the line closest to the mid point of the max chord
      closest <- which.min(sapply(perplines, function(i) st_distance(mid_point_on_chord, i)))
      perplines[[closest]]
      
    }
    
    ## construct an sf polygon from points:
    polygon = st_polygon(list(as.matrix(rbind(points_ex, points_ex[1,]))))
    
    ## compute chords
    chord = max_chord(polygon)
    mid_point_on_chord <- find_mid_point_on_chord(chord)
    pchord = find_max_chord(polygon, chord)
    max_perp_chord_midpoint_line <- max_perp_chord_midpoint(polygon)

    
    return(list(chord = chord,
                mid_point_on_chord = mid_point_on_chord,
                pchord = pchord,
                max_perp_chord_midpoint_line = max_perp_chord_midpoint_line))
    
  })
    
# show contour on object 
output$contour <- renderPlot({
  req(contour_points())
  req(prepped_image())
  req(length_and_widths())
    # plot original image and contours
    EBImage::display(prepped_image(),
                     method = "raster")
    lines(contour_points(),
          type = 'l',
          asp = 1,
          col = "green")
    
    points(contour_points(),
           cex = 0.5,
           col = "red")
    
    length_and_widths <- length_and_widths()
    points(length_and_widths$chord, col="green", cex=2)
    points(length_and_widths$mid_point_on_chord, col="green", cex=2)
    lines(length_and_widths$chord, col="green", lwd=2)
    plot(length_and_widths$pchord, add=TRUE, col = "blue")
    points(length_and_widths$pchord, col = "blue")
    plot(length_and_widths$max_perp_chord_midpoint_line, add = T, col = "green")
    points(length_and_widths$max_perp_chord_midpoint_line, add = T, col = "green")
    
    # label the chords
    
    # max chord length
    text(label = round(dist(length_and_widths$chord),2), 
         pos = 4,
         x = length_and_widths$chord[1,1], 
         y = length_and_widths$chord[1,2])
    
    # max perp chord length
    text(label = paste(round(st_length(length_and_widths$pchord),2)), 
         pos = 4,
         x = length_and_widths$pchord[3], 
         y = length_and_widths$pchord[4])
    
    # length mid point chord
    text(label = round(st_length(length_and_widths$max_perp_chord_midpoint_line),2), 
         pos = 4,
         x = length_and_widths$max_perp_chord_midpoint_line[3], 
         y = length_and_widths$max_perp_chord_midpoint_line[4])
    
    # length mid point on the max chord
    text(label = round(st_distance(length_and_widths$pchord, 
                                   length_and_widths$mid_point_on_chord),2), 
         pos = 1,
         x = length_and_widths$mid_point_on_chord[1], 
         y = length_and_widths$mid_point_on_chord[2])
    
  })
    

  ##  Coords of the locations we click on
  point_coords <- reactiveValues()
  
  ## Don't fire off the plot click too often
  plot_click_slow <- debounce(reactive(input$plot_click), 100)
  
  # collect points on click
  observeEvent(plot_click_slow(), {
    point_coords$x <- c(point_coords$x,
                                 (plot_click_slow()$x))
    point_coords$y <- c(point_coords$y,
                                 (plot_click_slow()$y))
  })
  
  # clear points when we click the button
  observeEvent(input$reset_points,{
    point_coords$x <- NULL
    point_coords$y <- NULL
  })
  
  # show a table of the coords of the locations we click on
  output$click_points <- renderTable({
    all_the_points <- data.frame(point_coords$x,
                                 point_coords$y)
    
    # show the data frame on the app
    all_the_points
  })
  
  
  # helpers to get the scaling value calculated, displayed, and stored
  all_the_points <- reactive({
    # get distance from last point to 2nd last point
    data.frame(point_coords$x,
               point_coords$y)
    
  })
  
  last_two_points <- reactive({
    req(all_the_points())
    all_the_points()[c(nrow(all_the_points()), 
                       (nrow(all_the_points()) - 1 )), ]
  })
  
  dist_last_two_points <- reactive({
    req(last_two_points())
    dist(last_two_points())
  })
  
  scale_factor <- reactive({
    req(dist_last_two_points())
    dist_last_two_points() / input$scale_factor_number
    
  })
  
  output$scale_ratio <-  renderUI({ 
    req(dist_last_two_points())
    req(scale_factor())
    paste0("Scale factor: ", 
           round(dist_last_two_points(), 3),  'px / ',
           input$scale_factor_number,
           "mm = ",
           round(scale_factor(), 3))
  })
  
  # store scale factor when we click the button
  scale_factor_storage <- eventReactive(input$save_scale_factor, {
    req(scale_factor())
    scale_factor()
  })
  
  # print the scale factor so we can see it
  output$saved_value <- renderUI({ 
    req(scale_factor_storage())
    paste0("This value has been stored: ", 
           round(scale_factor_storage(), 3))
  })
  
# Allows content from the Shiny application to be made available to the user as file downloads
# files for contour points and landmark points
  
  # file name to store point coords
  file_name_points <- reactive({
    inFile <- input$image
    
    if (is.null(inFile))
      return(NULL)
    paste0(inFile$name, sprintf("-point-coords-%s.csv",
                                format(Sys.time(), "%Y-%b-%d-%Hh%Mm%Ss")))
  })
  
  output$downloadDataOutline <- downloadHandler(
    filename <- function(){
      paste0(stringi::stri_extract_first(str = input$image$name, regex = ".*(?=\\.)"),
             sprintf("-outline-coords-%s.csv", format(Sys.time(), "%Y-%b-%d-%Hh%Mm%Ss")))
    },
  content <- function(file) {
    req(contour_points())
  write.csv(contour_points(), file)
  },
  contentType = "text/csv")

# this will get the landmarks when we have figured out how to clear the scale points
output$downloadDataPoints <- downloadHandler(
  filename <- function(){
    paste0(stringi::stri_extract_first(str = input$image$name, regex = ".*(?=\\.)"),
           sprintf("-landmark-coords-%s.csv", format(Sys.time(), "%Y-%b-%d-%Hh%Mm%Ss")))
  },
  content <- function(file) {
    req(all_the_points())
    write.csv(all_the_points(), file)
  },
  contentType = "text/csv")


}

# Run the application
shinyApp(ui = ui, server = server)