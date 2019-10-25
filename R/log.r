#' Set the level of logger
#' @export set.log.level
#' @param level log level
#' @examples
#' \dontrun{
#' set.log.level(futile.logger::DEBUG)
#' set.log.level(futile.logger::INFO)
#' set.log.level(futile.logger::WARN)
#' }
set.log.level<-function(level){
  futile.logger::flog.threshold(level)
}
