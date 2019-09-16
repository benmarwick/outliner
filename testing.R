library(EBImage)  
library(imager)
library(jpeg)
library(maptools)
library(autothresholdr)
library(ijtiff)
library(ggplot2)
library(dplyr)

img_file <- "test-images/levpoint.jpg"
x <- readImage(img_file)
# x <- getFrame(x, 3)
x_dim <- dim(x) 
x2 <- x * 0.7   # increase contrast
EBImage::display(x2)
plot(x2)

x_dim

# Cropping
EBImage::display(x[(x_dim[1] * 0) : (x_dim[1] * 1),
                   (x_dim[2] * 0) : (x_dim[2] * 1) , ])

EBImage::display(x)
colorMode(x) = Grayscale
EBImage::display(x)

kern = makeBrush(10, shape='diamond')
x= EBImage::erode(x, kern)
EBImage::display(x_erode)

y <- x > otsu(x_erode, levels=256)
EBImage::display(y, all=TRUE)

y1 <- thresh(y, 1, 1, 0.1)
EBImage::display(y1, all=TRUE)


###
img_file <- "test-images/gp.png"
x <- readImage(img_file)
x_dim <- dim(x) 
colorMode(x) = Color

x_cropped <- 
  x[(x_dim[1] * 0) : (x_dim[1] * 0.75),
    (x_dim[2] * 0) : (x_dim[2] * 1) , ]

x_gam <- x_cropped ^ 2.5
x_con <- x_gam + 0.35
x_bri <- x_con * 1

EBImage::display(x_bri, all=TRUE)

path <- tempfile(pattern = "temp", fileext = ".tif")
write_tif(as_ijtiff_img(x_bri), path, overwrite = TRUE)
y_tif <- round(read_tif(path) * 185)
ijtiff::display(y_tif)
EBImage::display(as_EBImage(y_tif))

auto_thresh_value <- auto_thresh(y_tif, method = 'Huang2')[1]

# plot image histogram
tibble(x = as.vector(y_tif)) %>%
  ggplot() + 
  aes(x) + 
  stat_density(bw = 3) +
  geom_vline(xintercept = auto_thresh_value, 
             colour = "red") +
  theme_minimal()

threshed <- auto_thresh_apply_mask( y_tif, method = auto_thresh_value )
ijtiff::display(threshed)
threshed_eb <- as_EBImage(threshed)

threshed_eb1 <- threshed_eb * 2
EBImage::display(threshed_eb1)

oc = ocontour(threshed_eb1)
ocs <- vapply(oc, length, integer(1))
idx <- as.integer(which.max(ocs))
plot(oc[[1]], 
     type='l', 
     asp = 1)
points(oc[[1]], 
       col=2, 
       cex = 0.25)

oc = ocontour(threshed_eb)
plot(oc[[1]], 
     type='l', 
     asp = 1)
points(oc[[1]], 
       col=2, 
       cex = 0.25)



####


library(imager)
canny <- cannyEdges(y_tif)
canny <- cannyEdges(as.cimg(imageData(x_bri[,,1])))

path <- tempfile(pattern = "temp", fileext = ".jpeg")
writeImage(x_bri, path)
x_bri_i <- load.image(path)
canny <- cannyEdges(x_bri_i)

canny_inv <- !canny
plot(canny_inv)

oc = ocontour(canny)
EBImage::display(x_bri, method = "raster")
lines(oc[[1]], 
     type='l', 
     asp = 1,
     add = T)
points(oc[[1]], 
       col=2, 
       cex = 0.25)


########
img_file <- "gp.png" # "test-images/gp.png"
x <- readImage(img_file)
x_dim <- dim(x) 

EBImage::display(x, 
                 method = "raster")

x_cropped <- 
  x[(x_dim[1] * 0) : (x_dim[1] * 0.7),
    (x_dim[2] * 0) : (x_dim[2] * 1) , ]

EBImage::display(x_cropped, 
                 method = "raster")

x_gam <- x_cropped ^ 1

EBImage::display(x_gam, 
                 method = "raster")

x_con <- x_gam + 0

EBImage::display(x_con, 
                 method = "raster")

x_bri <- x_con * 1

EBImage::display(x_bri, 
                 method = "raster")

canny <- cannyEdges(as.cimg(imageData(x_bri[,,1])))

oc = ocontour(canny)

lines(oc[[1]], 
      type='l', 
      asp = 1,
      col = "green")
points(oc[[1]], 
       cex = 0.5, 
       col = "red")

#####





MBR <- function(points) {
  # Analyze the convex hull edges                       
  a <- alphahull::ashape(points, alpha=1000)                 # One way to get a convex hull...
  e <- a$edges[, 5:6] - a$edges[, 3:4]            # Edge directions
  norms <- apply(e, 1, function(x) sqrt(x %*% x)) # Edge lengths
  v <- diag(1/norms) %*% e                        # Unit edge directions
  w <- cbind(-v[,2], v[,1])                       # Normal directions to the edges
  
  # Find the MBR
  vertices <- (points) [a$alpha.extremes, 1:2]    # Convex hull vertices
  minmax <- function(x) c(min(x), max(x))         # Computes min and max
  x <- apply(vertices %*% t(v), 2, minmax)        # Extremes along edges
  y <- apply(vertices %*% t(w), 2, minmax)        # Extremes normal to edges
  areas <- (y[1,]-y[2,])*(x[1,]-x[2,])            # Areas
  k <- which.min(areas)                           # Index of the best edge (smallest area)
  
  # Form a rectangle from the extremes of the best edge
  cbind(x[c(1,2,2,1,1),k], y[c(1,1,2,2,1),k]) %*% rbind(v[k,], w[k,])
}


# test with an actual artefact outline
plot(points <- read.csv("~/Downloads/levpoint-outline-coords-2019-Apr-25-09h32m35s.csv")[,2:3], asp = 1)
plot(points <- read.csv("~/Downloads/gp-outline-coords-2019-Apr-22-10h56m19s.csv")[,2:3], asp = 1)
mbr <- MBR(as.matrix(points))

# Plot the hull, the MBR, and the points
limits <-
  alphahull::apply(mbr, 2, function(x)
    c(min(x), max(x))) # Plotting limits
plot(
  alphahull::ashape(points, alpha = 1000),
  col = "Gray",
  pch = 20,
  xlim = limits[, 1],
  ylim = limits[, 2],
  asp = 1
)                # The hull
lines(mbr, col = "Blue", lwd = 3)                         # The MBR
points(points, pch = 19, cex = 0.5)                                # The points

# draw max dimension points and line
suppressPackageStartupMessages(library(tidyverse))
df_dist = data.frame(as.matrix(dist(cbind(points$V1,points$V2))))
df_dist_x = df_dist %>% 
  mutate(row.1 = 1:nrow(df_dist)) %>% 
  mutate(Y = paste0("Y", row_number())) %>%
  gather(X,  distance, X1:nrow(.)) %>% 
  select(X, Y, distance) %>% 
  mutate_at(vars(X, Y), parse_number)

df_dist_x_max <- 
df_dist_x %>% 
  dplyr::filter(distance == max(distance)) 

points(points[df_dist_x_max$X[1],], col = "red", cex = 2)
points(points[df_dist_x_max$X[2],], col = "red", cex = 2)

segments(points[df_dist_x_max$X[1], 'V1'], 
         points[df_dist_x_max$X[1], 'V2'],
         points[df_dist_x_max$X[2], 'V1'], 
         points[df_dist_x_max$X[2], 'V2'],
         col = "green")

# draw 2nd max dimension
p1 <- c(points[df_dist_x_max$X[1], 'V1'], 
  points[df_dist_x_max$X[1], 'V2'])
p2 <- c(points[df_dist_x_max$X[2], 'V1'], 
  points[df_dist_x_max$X[2], 'V2'])

x_mid = (p1[1] + p2[1])/2
y_mid = (p1[2] + p2[2])/2
mid_point <- c(x_mid, y_mid)
points(mid_point[1], mid_point[2], col = "red", cex = 1)

# rotate polygon so that max chord is vertical
points[df_dist_x_max$X[1],]
points[df_dist_x_max$X[2],]


# transform the points and lines into spatial objects
library(sf)
points_sf <- st_as_sf(points, coords = c("V1", "V2"))
library(sp)
library(rgeos)

newline = matrix(c(points[df_dist_x_max$X[1], 'V1'], 
                   points[df_dist_x_max$X[1], 'V2'],
                   points[df_dist_x_max$X[2], 'V1'], 
                   points[df_dist_x_max$X[2], 'V2']), byrow = T, nrow = 2)

spline <- as(st_as_sfc(st_as_text(st_linestring(newline))), "Spatial") # there is probably a more straighforward solution...
position <- gProject(spline, as(points_sf, "Spatial"))
position <-  coordinates(gInterpolate(spline, position))
colnames(position) <- c("X2", "Y2")

segments <- 
  data.frame(st_coordinates(points_sf), position)  

segments$dist <- NULL
for(i in 1:nrow(segments)){
  segments$dist[i] <- 
proxy::dist(data.frame(segments$X[i], segments$Y[i]),  
            data.frame(segments$X2[i], segments$Y2[i]))
}

# max width perpendicular to length axis
max_segment <- segments[which.max(segments$dist), ]
max_segment <- segments[segments$Y == max_segment$Y, ]

# width at midpoint of length axis
l_mid_margin <- segments[which.min(abs(segments$Y  - mid_point[2])), ]
r_mid_margin <- segments[which.min(abs(segments$Y2 - mid_point[2])), ]

points(l_mid_margin$X, l_mid_margin$Y, col = "pink")
points(r_mid_margin$X, r_mid_margin$Y, col = "pink")

segments(l_mid_margin$X, l_mid_margin$Y,
         r_mid_margin$X, r_mid_margin$Y,
         col = "pink")
  

library(ggplot2)
ggplot() +
  geom_sf(data = points_sf) +
  geom_point(data = l_mid_margin, 
             aes(X,Y), 
             colour = "green") +
  geom_point(data = r_mid_margin, 
             aes(X,Y), 
             colour = "blue") +
  geom_segment(data = data.frame(rbind(l_mid_margin, 
                                       r_mid_margin)),  
               aes(X,  Y, xend = X2, yend = Y2), colour = "purple") +
  geom_segment(data = max_segment,  aes(X, 
                                        Y,
                                        xend = X2, 
                                        yend = Y2), colour = "pink") +
  geom_line(data = as.data.frame(newline), aes(V1,V2)) +
  coord_sf()

# small example ###################################
points_ex <- points # points[sample(nrow(points), 700, TRUE), ]
plot(points_ex, asp = 1, cex = 0.5)

# draw max dimension points and line
suppressPackageStartupMessages(library(tidyverse))
df_dist = data.frame(as.matrix(dist(cbind(points_ex$V1,points_ex$V2))))
df_dist_x = df_dist %>% 
  mutate(row.1 = 1:nrow(df_dist)) %>% 
  mutate(Y = paste0("Y", row_number())) %>%
  gather(X,  distance, X1:nrow(.)) %>% 
  select(X, Y, distance) %>% 
  mutate_at(vars(X, Y), parse_number)

df_dist_x_max <- 
  df_dist_x %>% 
  dplyr::filter(distance == max(distance)) 

points(points_ex[df_dist_x_max$X[1],], col = "red", cex = 2)
points(points_ex[df_dist_x_max$X[2],], col = "red", cex = 2.5)

segments(points_ex[df_dist_x_max$X[1], 'V1'], 
         points_ex[df_dist_x_max$X[1], 'V2'],
         points_ex[df_dist_x_max$X[2], 'V1'], 
         points_ex[df_dist_x_max$X[2], 'V2'],
         col = "green")

# draw 2nd max dimension
p1 <- c(points_ex[df_dist_x_max$X[1], 'V1'], 
        points_ex[df_dist_x_max$X[1], 'V2'])
p2 <- c(points_ex[df_dist_x_max$X[2], 'V1'], 
        points_ex[df_dist_x_max$X[2], 'V2'])

x_mid = (p1[1] + p2[1])/2
y_mid = (p1[2] + p2[2])/2
mid_point <- c(x_mid, y_mid)
points(mid_point[1], mid_point[2], col = "red", cex = 1)

# transform the points and lines into spatial objects
library(sf)
library(sp)
library(rgeos)

points_sf <- st_as_sf(points_ex, coords = c("V1", "V2"))

newline = matrix(c(points_ex[df_dist_x_max$X[1], 'V1'], 
                   points_ex[df_dist_x_max$X[1], 'V2'],
                   points_ex[df_dist_x_max$X[2], 'V1'], 
                   points_ex[df_dist_x_max$X[2], 'V2']), byrow = T, nrow = 2)

spline <- as(st_as_sfc(st_as_text(st_linestring(newline))), "Spatial") # there is probably a more straighforward solution...
position <- gProject(spline, as(points_sf, "Spatial"))
position <-  coordinates(gInterpolate(spline, position))
colnames(position) <- c("X2", "Y2")

segments <- 
  data.frame(st_coordinates(points_sf), position)  

segments$dist <- NULL
for(i in 1:nrow(segments)){
  segments$dist[i] <- 
    proxy::dist(data.frame(segments$X[i], segments$Y[i]),  
                data.frame(segments$X2[i], segments$Y2[i]))
}

# max width perpendicular to length axis
max_segment <- segments[which.max(segments$dist), ]
max_segment <- segments[segments$Y == max_segment$Y, ]

segments(max_segment$X[1], max_segment$Y[1],
         max_segment$X2[1], max_segment$Y2[1],
         col = "purple")

segments(max_segment$X[2], max_segment$Y[2],
         max_segment$X2[2], max_segment$Y2[2],
         col = "purple")

# width at midpoint of length axis
l_mid_margin <- segments[which.min(abs(segments$Y  - mid_point[2])), ]
r_mid_margin <- segments[which.min(abs(segments$Y2 - mid_point[2])), ]

points(l_mid_margin$X, l_mid_margin$Y, col = "pink")
points(r_mid_margin$X, r_mid_margin$Y, col = "pink")

segments(l_mid_margin$X, l_mid_margin$Y,
         r_mid_margin$X, r_mid_margin$Y,
         col = "pink")

##-0--------------------------------------------------

# from https://stackoverflow.com/questions/55854273/how-to-find-the-longest-chord-perpendicular-to-the-maximum-chord-through-a-polyg/55874609#55874609

plot(points_ex <- read.csv("~/Downloads/gp-outline-coords-2019-Apr-22-10h56m19s.csv")[,2:3], asp = 1)
plot(points_ex <- read.csv("~/Downloads/levpoint-outline-coords-2019-Apr-25-09h32m35s.csv")[,2:3], asp = 1)

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


### drawing ###

## construct an sf polygon from points:
polygon = st_polygon(list(as.matrix(rbind(points_ex, points_ex[1,]))))
# Get the max length chord between vertices:

chord = max_chord(polygon)

# Plot the polygon and the chord and the chord points:

plot(polygon)
points(chord, col="green",cex=2)
lines(chord,col="green",lwd=2)

# find find_mid_point_on_chord(chord)
mid_point_on_chord <- find_mid_point_on_chord(chord)
plot(mid_point_on_chord, col="green",cex=2, add = T)

## show the max perpendicular chord
pchord = find_max_chord(polygon, chord)
plot(pchord, add=TRUE, col = "blue")
points(pchord, col = "blue")

## show the perpendicular chord at the mid point of the max chord
max_perp_chord_midpoint_line <- max_perp_chord_midpoint(polygon)
plot(max_perp_chord_midpoint_line, add = T, col = "green")
points(max_perp_chord_midpoint_line, add = T, col = "green")

# label the chords

# max chord length
text(label = round(dist(chord),2), 
     pos = 3,
     x = chord[1,1], 
     y = chord[1,2])

# max perp chord length
text(label = paste(round(st_length(pchord),2)), 
     pos = 4,
     x = pchord[3], 
     y = pchord[4])

# length mid point chord
text(label = round(st_length(max_perp_chord_midpoint_line),2), 
     pos = 4,
     x = max_perp_chord_midpoint_line[3], 
     y = max_perp_chord_midpoint_line[4])

# length mid point on the max chord
text(label = round(st_distance(pchord, mid_point_on_chord),2), 
     pos = 1,
     x = mid_point_on_chord[1], 
     y = mid_point_on_chord[2])





# resize images and overwrite 
# webshot::resize(img_file,  "600x")


