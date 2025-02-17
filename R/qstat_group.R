#' @export
qstat_group <- qstat

#' @export
qstat_layer <- function (compute_layer, ...){
    ggplot2::ggproto("StatTemp", ggplot2::Stat, compute_layer = compute_layer, 
        ...)
}

#' @export
qstat_panel <- function (compute_panel, ...){
    ggplot2::ggproto("StatTemp", ggplot2::Stat, compute_panel = compute_panel, 
        ...)
}
