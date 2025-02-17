#' @export
qlayer <- function (mapping = NULL, data = NULL, geom = ggplot2::GeomPoint, 
                    stat = ggplot2::StatIdentity, 
    position = ggplot2::position_identity(), ..., na.rm = FALSE, show.legend = NA, 
    inherit.aes = TRUE) 
{
    ggplot2::layer(data = data, mapping = mapping, geom = geom, 
        stat = stat, position = position, show.legend = show.legend, 
        inherit.aes = inherit.aes, params = rlang::list2(na.rm = na.rm, 
            ...))
}
