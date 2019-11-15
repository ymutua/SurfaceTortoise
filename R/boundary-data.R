#' Boundary polygon for field 10 at Bjertorp farm in southwest Sweden
#'
#' @usage data(boundary)
#'
#' @format A SpatialPolygonsDataFrame with the boundary polygon. The polygon has
#' been shrunk (i.e. the field is somewhat larger than the polygon)
#' Projected coordinate system Sweref99TM (epsg: 3006).
#'
#' @keywords datasets
#'
#' @references Piikki, K., Söderström, M., & Stenberg, B. (2013). Sensor data
#' fusion for topsoil clay mapping. Geoderma, 199, 106-116.
#'
#' @example
#' data(boundary)
#' raster::plot(boundary)
"boundary"
