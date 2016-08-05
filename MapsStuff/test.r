library(rgdal)
library(threejs)

import("RenderData.png")
# Download MODIS 16-day 1 degree Vegetation Index data manually from
# http://neo.sci.gsfc.nasa.gov/view.php?datasetId=MOD13A2_M_NDVI
# or use the following cached copy from May 25, 2014
cache <- tempfile()
# writeBin(
#   readBin(
#     url("http://illposed.net/nycr2015/MOD13A2_E_NDVI_2014-05-25_rgb_360x180.TIFF",
#         open="rb"),
#     what="raw", n=1e6),
#   con=cache)
x <- readGDAL(import("RenderData.png"))
# Obtain lat/long coordinates and model values as a data.frame
x <- as.data.frame(cbind(coordinates(x), x@data[,1]))
names(x) <- c("long","lat","value")

# Remove ocean areas and NA values
x <- x[x$value < 255,]
x <- na.exclude(x)

# Cut the values up into levels corresponding to the
# 99th, 95th, 90th, percentiles and then all the rest.
x$q <- as.numeric(
         cut(x$value,
             breaks=quantile(x$value, probs=c(0,0.90,0.95,0.99,1)),
                            include.lowest=TRUE))

# Colors for each level
col = c("#0055ff","#00aaff","#00ffaa","#aaff00")[x$q]

# bling out the data
globejs(lat=x$lat, long=x$long,
         val=x$q^3,    # Bar height
         color=col,
         pointsize=0.5,
         atmosphere=TRUE)
