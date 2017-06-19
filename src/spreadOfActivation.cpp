#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector staSearch(NumericMatrix x){
    int from                = x( 0, 0)-1;             // from node -> row
    int to                  = x( 0, 1)-1;             // to node -> column
    int nrun                = x( 0, 2);               // how many runs
    int stepmax             = x( 0, 3);               // step limit
    NumericVector steplist  = NumericVector(nrun);    // data container
    int nnodes              = x.nrow();               // number of nodes
    NumericMatrix adj       = x(_,Range(4,nnodes+3)); // adjacanecy matrix
    NumericVector edges;
    std::vector<int> neigh;
    int steps;
    int pos;
    int found;
    int j;
    int run;
    int sel;
    for(run = 0; run < nrun; run++){
      steps = 0;              
      pos   = from;
      found = 0;
      while(found == 0 & steps < stepmax){
        steps++;
        edges = adj(pos,_);
        neigh.clear();
        for(j = 0; j < nnodes; j++){
          if(edges(j) == 1){
            neigh.push_back(j);
            }
          }
        sel = rand() % neigh.size();
        pos = neigh[sel];
        if(pos == to){
          found = 1;
          }
        }
      steplist(run) = steps;
      }
    return steplist;}

// [[Rcpp::export]]
NumericVector dynSearch(NumericMatrix x){
    int from                = x( 0, 0)-1;             // from node -> row
    int to                  = x( 0, 1)-1;             // to node -> column
    int nrun                = x( 0, 2);               // how many runs
    int stepmax             = x( 0, 3);               // step limit
    NumericVector steplist  = NumericVector(nrun);    // data container
    int nnodes              = x.nrow();               // number of nodes
    NumericMatrix adj       = x(_,Range(4,nnodes+3)); // adjacanecy matrix
    NumericVector nodes     = NumericVector(nnodes);
    for(int i = 0; i<nnodes; i++){
      nodes[i] = i;
      }
    NumericVector edges;
    std::vector<int> neigh;
    int steps;
    int pos;
    int found;
    int j;
    int run;
    int sel;
    for(run = 0; run < nrun; run++){
      steps = 0;              
      pos   = from;
      found = 0;
      while(found == 0 & steps < stepmax){
        steps++;
        if(rand() % 100 < 5){
          pos = nodes[rand() % nnodes];
        } else {
          edges = adj(pos,_);
          neigh.clear();
          for(j = 0; j < nnodes; j++){
            if(edges(j) == 1){
              neigh.push_back(j);
              }
            }
          sel = rand() % neigh.size();
          pos = neigh[sel];
          }
        if(pos == to){
          found = 1;
          }
        }
      steplist(run) = steps;
      }
    return steplist;}


// [[Rcpp::export]]
bool notInVector(int item, std::vector<int> set){
  std::vector<int>::iterator it;
  it = find (set.begin(), set.end(), item);
  if (it != set.end()) return false;
  return true;
  }

// [[Rcpp::export]]
NumericVector staSearch2(GenericVector x){
  int from                = x[0];               // from node -> row
  int to                  = x[1];               // to node -> column
  int nrun                = x[2];               // how many runs
  int stepmax             = x[3];               // step limit
  NumericMatrix adj       = x[4];               // adjacanecy matrix
  NumericVector nogoNV    = x[5];
  NumericVector steplist  = NumericVector(nrun);// data container
  std::cout << from << "\n";
  std::vector<int> nogo;
  int nnodes              = adj.nrow();           // number of nodes
  NumericVector edges;
  std::vector<int> neigh;
  int steps, pos, found, j, run, sel;
  for(j = 0; j < nogoNV.size(); j++) nogo.push_back(int(nogoNV[j]));
  for(run = 0; run < nrun; run++){
    steps = 0;              
    pos   = from;
    found = 0;
    while(found == 0 & steps < stepmax){
      steps++;
      edges = adj(pos,_);
      neigh.clear();
      for(j = 0; j < nnodes; j++){
        if(edges(j) == 1){
          neigh.push_back(j);
          }
        }
      sel = rand() % neigh.size();
      pos = neigh[sel];
      if(pos == to && notInVector(pos,nogo)){
        found = 1;
      }
    }
    steplist(run) = steps;
  }
  return steplist;}

