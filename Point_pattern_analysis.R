#PERU SPPA 

library(sf)
library(spatstat)
library(spatstat.geom)
library(dplyr)

#load in latrine points
latrine_pts <- st_read("C:\\Users\\lvbj3q\\Desktop\\SPPA\\SPPA_latrinepoints_nodup_100625.shp")
st_geometry_type(latrine_pts)
summary(latrine_pts)
plot(st_geometry(latrine_pts))

#extract coordinates
latrine_coords <- st_coordinates(latrine_pts)
class(latrine_coords) #should be function
head(latrine_coords) #should be lat & long

#projected ppp object 
window <- owin(xrange = range(latrine_coords[,1]),
               yrange = range(latrine_coords[,2]))

latrine_ppp <- ppp(x = latrine_coords[,1],
                   y = latrine_coords[,2],
                   window = window)
summary(latrine_ppp)
class(latrine_ppp) #should be planar point pattern (ppp) object

plot.ppp(latrine_ppp)
plot.ppp(latrine_ppp, main = "latrine locations", axes = T)
grid()


###G-function: evaluation of nearest neighbor distances - proportion of points with a neighbor within distance r
G <- Gest(latrine_ppp)
summary(G)
#r is distance, theo is theoretical G(r)
plot(G, xlim = c(0, 435.3), main = "G-functions - latrine locations") #x-axis set to largest observed r 

#compare observed to theoretical G-function for a random distribution of points (CSR)
plot(G, cbind(rs, theo) ~ theo, main = "G-function - latrine locations", xlab = "Theoretical G(r)", ylab = "Observed G(r)")

###K-function: evaluation of clustering or dispersion at different spatial scales

#monte-carlo envelopes w/ 99 simulations
envelope_latrine <- envelope(latrine_ppp, Kest, nsim = 99)
plot(envelope_latrine)
#-> simulated K differs sig from observed & is greater, so points are clustered
#-> because k increases linearly with distance (r) points are clustered at multiple scales


###Kernel density estimation 

#bandwidth selection using cross-validation
print(bw.o <- bw.diggle(latrine_ppp)) #choice
print(bw.CvL(latrine_ppp))
print(bw.ppl(latrine_ppp))

#compute density surface 
d1 <- density.ppp(latrine_ppp, sigma = bw.o, kernel = "quartic")
d05 <- density.ppp(latrine_ppp, sigma = bw.o, kernel = "quartic", adjust = 0.5)
d4 <- density.ppp(latrine_ppp, sigma = bw.o, kernel = "quartic", adjust = 4)
d2 <- density.ppp(latrine_ppp, sigma = bw.o, kernel = "quartic", adjust = 2)


par(mfrow = c(2, 2))
plot.im(d05, main = paste("Bandwidth=", round(bw.o * 0.4, 4), " (optimum*.5)", #looks the best 
                          sep = "", collapse = ""))
contour(d05, add = TRUE)
plot(d1, main = paste("Bandwidth=", round(bw.o, 4), " (optimum)", sep = "",
                      collapse = ""))
contour(d1, add = TRUE)
plot(d2, main = paste("Bandwidth=", round(bw.o * 1.5, 4), " (optimum*1.5)",
                      sep = "", collapse = ""))
contour(d2, add = TRUE)
plot(d4, main = paste("Bandwidth=", round(bw.o * 2, 4), " (optimum*2)",
                      sep = "", collapse = ""))
contour(d4, add = TRUE)
par(mfrow = c(1, 1))

d <- d05


