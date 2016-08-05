library(threejs)

# Plot flights to frequent destinations from
# Callum Prentice's global flight data set,
# http://callumprentice.github.io/apps/flight_stream/index.html
data(flights)

# Approximate locations as factors
dest   <- factor(sprintf("%.2f:%.2f",flights[,3], flights[,4]))

# A table of destination frequencies
freq <- sort(table(dest), decreasing=TRUE)

# The most frequent destinations in these data, possibly hub airports?
frequent_destinations <- names(freq)[1:10]

# Subset the flight data by destination frequency
idx <- dest %in% frequent_destinations
frequent_flights <- flights[idx, ]

# Lat/long of frequent destinations
latlong <- unique(frequent_flights[,3:4])

# Plot frequent destinations as bars, and the flights to and from
# them as arcs. Adjust arc width and color by frequency.
earth <- system.file("images/world.jpg",  package="threejs")
globejs(img=earth, lat=latlong[,1], long=latlong[,2], arcs=frequent_flights,
        arcsHeight=0.3, arcsLwd=2, arcsColor="#ffff00", arcsOpacity=0.15,
        atmosphere=TRUE)
