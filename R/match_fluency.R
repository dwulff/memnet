#' Match group fluency data
#'
#' Prune fluency data of two groups so that the groups'
#' distributions of number of productions match.
#'
#' @param data a data frame containing (at least) fluency productions,
#'   subject id, and a grouping variable.
#' @param labels a character vector containing the labels of variables
#'   containing the fluency productions, subject id, and grouping.
#'
#' @return
#' Function returns two lists with matching numbers of fluency
#' productions.
#'
#' @export
match_fluency = function(data, labels, type = 'match', tail = TRUE){

  # collect first row of each subject
  tmp_data = ddply(data, labels[-1], function(x) x[1,])

  # split groups and retrieve ids
  groups = split(tmp_data, tmp_data[,labels[3]] < median(tmp_data[,labels[3]]))
  ids_g1 = groups[[1]][,labels[2]]
  ids_g2 = groups[[2]][,labels[2]]
  tmp_g1 = subset(data, data[,labels[2]] %in% ids_g1)
  tmp_g2 = subset(data, data[,labels[2]] %in% ids_g2)

  # get responses
  resp_g1 = split(tmp_g1[,labels[1]],tmp_g1[,labels[2]])
  resp_g2 = split(tmp_g2[,labels[1]],tmp_g2[,labels[2]])

  # get lengths
  len_g1 = sapply(resp_g1,length)
  len_g2 = sapply(resp_g2,length)

  # order
  resp_g1 = resp_g1[order(len_g1)]
  resp_g2 = resp_g2[order(len_g2)]
  len_g1  = len_g1[order(len_g1)]
  len_g2  = len_g2[order(len_g2)]

  # create new lists
  resp_g1_new = list()
  resp_g2_new = list()
  for(i in 1:length(len_g1)){
    if(type == 'match'){
      id = which.min(abs(len_g1[i] - len_g2))
      } else {
      id = 1
      }
    min_len = min(len_g1[i],len_g2[id])
    if(tail == TRUE){
      resp_g1_new[[i]] = resp_g1[[i]][ 1:min_len]
      resp_g2_new[[i]] = resp_g2[[id]][1:min_len]
      } else {
      resp_g1_new[[i]] = resp_g1[[i]][ (1+len_g1[i] -min_len):len_g1[i]]
      resp_g2_new[[i]] = resp_g2[[id]][(1+len_g2[id]-min_len):len_g2[id]]
      }
    resp_g2 = resp_g2[-id]
    len_g2  = len_g2[-id]
    }
  out = list(resp_g1_new, resp_g2_new)
  names(out) = c('above','below')
  return(out)
  }
#
# a = match_fluency(d, labels, tail = T)
#
# lengths(a[[1]])
# lengths(a[[2]])
#
#lengths(a[[1]]) - lengths(resp_g1)
#lengths(a[[2]]) - lengths(resp_g2)
#

