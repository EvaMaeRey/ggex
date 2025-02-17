proto_update <- function (`_class`, `_inherit`, default_aes_update = NULL, ...) {
    if (!is.null(default_aes_update)) {
        
      new_default_aes <- combine_aes(`_inherit`$default_aes, default_aes_update)
      
    } else {
      
      new_default_aes <- `_inherit`$default_aes
      
    }
  
    ggplot2::ggproto(`_class` = `_class`, 
                     `_inherit` = `_inherit`,
                     default_aes = new_default_aes, ...)
}

#' @export
qproto_update <- function (`_inherit`, default_aes_update = NULL, ...){
    proto_update("protoTemp", `_inherit`, default_aes_update = default_aes_update, 
        ...)
}
