# -------------------------------------------------------------------------------------------------------------------- #
#
#          PLOT HELPERS
#
# -------------------------------------------------------------------------------------------------------------------- #


# Mix two colors according to mixing parameter
cmix = function(col1,col2,weight,format = 'rgb'){

  if(ifelse(is.matrix(col1),nrow(col1),length(col1)) != length(weight) &
     length(col1) != 1){
    stop('Length of col1 must be either 1 or matching the length of weight')
  }
  if(ifelse(is.matrix(col2),nrow(col2),length(col2)) != length(weight) &
     length(col2) != 1){
    stop('Length of col1 must be either 1 or matching the length of weight')
  }
  if(length(weight) == 1){
    if(ifelse(is.matrix(col1),nrow(col1),length(col1)) !=
       ifelse(is.matrix(col2),nrow(col2),length(col2))){
      stop('If length of weight = 1, number of colors in col1 and col2 must match')
    }
  }

  nrows = max(c(ifelse(is.matrix(col1),nrow(col1),length(col1)),
                ifelse(is.matrix(col2),nrow(col2),length(col2)),
                length(weight)))

  if(is.character(col1)){
    if(length(col1) == 1){
      col1 = grDevices::col2rgb(col1)
      col1 = matrix(c(col1),ncol=3,nrow=nrows,byrow=T)
    } else {
      col1 = t(sapply(col1,grDevices::col2rgb))
    }
  } else{
    col1 = matrix(c(col1),ncol=3,nrow=nrows,byrow=F)
  }
  if(is.character(col2)){
    if(length(col2) == 1){
      col2 = grDevices::col2rgb(col2)
      col2 = matrix(c(col2),ncol=3,nrow=nrows,byrow=T)
    } else {
      col2 = t(sapply(col2,grDevices::col2rgb))
    }
  } else{
    col2 = matrix(c(col2),ncol=3,nrow=nrows,byrow=F)
  }


  col = col1 * (1-weight) + col2 * weight

  if(format == 'rgb') return(col)
  if(format == 'hex') return(grDevices::rgb(data.frame(col),maxColorValue = 255))
  if(!format %in% c('rgb','hex')) stop('Choose either "rgb" or "hex" as format')

}


# round even
round_even = function(x){
  rx = round(x)
  test = rx %% 2
  ifelse(test == 0,rx,
         ifelse( x < 0,
                 ifelse(x <  rx,floor(x),ceiling(x)),
                 ifelse(x >= rx,ceiling(x),floor(x))))
}


# find complementary color
compCol = function(col){

  # ----- extract color

  if(is.character(col)){
    col = grDevices::col2rgb(col)
    } else if(is.numeric(col)){
    if(length(col) == 3){
      col = matrix(col,nrow=3)
      } else stop('if numeric, col must be length 3')
    } else stop('col must be character or numeric')

  # ----- transform to complementary color

  compCol    = grDevices::rgb2hsv(col)
  compCol[1] = ifelse(compCol[1] > .5, compCol[1] - .5, compCol[1] + .5)
  compCol    = grDevices::hsv(compCol[1],compCol[2],compCol[3])

  # return
  return(compCol)
  }

