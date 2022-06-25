MXP_locate_regex <- function(string, pattern) {
  MXPgreg <- gregexpr(pattern, string)[[1]]
  if (MXPgreg[1] > 0) {
    MXPgreg_lengthchar <- attributes(MXPgreg)$match.length
    #
    loc_mat <- do.call(rbind, lapply(1:length(MXPgreg_lengthchar), function(i) {
      c(MXPgreg[i], (MXPgreg[i] + MXPgreg_lengthchar[i] - 1))
    }))
  } else {
    loc_mat <- matrix(c(NA, NA), nrow = 1)
  }
  #
  return(loc_mat)
}