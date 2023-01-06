MXP_locate_regex <- function(string, pattern, ignore.case = FALSE,
                             perl = FALSE, fixed = FALSE, useBytes = FALSE) {
  ##
  MXPgreg <- gregexpr(pattern, string, ignore.case, perl, fixed, useBytes)[[1]]
  if (MXPgreg[1] > 0) {
    MXPgreg_lengthchar <- attributes(MXPgreg)$match.length
    #
    loc_mat <- do.call(rbind, lapply(1:length(MXPgreg_lengthchar), function(i) {
      c(MXPgreg[i], (MXPgreg[i] + MXPgreg_lengthchar[i] - 1))
    }))
  } else {
    loc_mat <- NULL
  }
  #
  return(loc_mat)
}
