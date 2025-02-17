#' @export
qstat <- function (compute_group = ggplot2::Stat$compute_group, ...) {
    ggplot2::ggproto("StatTemp", ggplot2::Stat, compute_group = compute_group, 
        ...)
}
