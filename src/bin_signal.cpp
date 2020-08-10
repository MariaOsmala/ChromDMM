#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix bin_signal(IntegerMatrix data, int bin_width ) {
  int window_size = data.ncol();
  int numrow = data.nrow();
  int numbin = window_size/bin_width;

  if (window_size % bin_width != 0) {
    Rprintf("bin_signal: Window width is not a multiple of bin width!\n");
  }

  IntegerMatrix res(numrow, numbin);
  Function rseq("seq");
  Function paste0("paste0");
  Function setcolnames("colnames<-");

  IntegerVector sequence = rseq(0, window_size-1, bin_width);


  for (int row = 0; row < numrow; row++) {
    IntegerVector vec = data(row, _);

    for (int bin = 0; bin < sequence.length(); bin++) {
      res(row, bin) = sum(vec[seq(sequence[bin], sequence[bin]+bin_width-1)]);
    }
  }

  if (data.nrow() > 1)
    res = setcolnames(res, paste0("bin", seq_len(res.ncol())));

  return res;
}

/* If given the data as binned (contains the unbinned version), extracts the
 * data in the shifted window for sample(s) indicated by indices
 * This function is required when shifting the profiles
 * If we have learned that the profile should be flipped, it is flipped
 * List binned, list of M of matrix of size N x S
 * indices are the sample numbers between 1:N
 */

// [[Rcpp::export]]
List extract_binned_signal(List binned, IntegerVector indices) {

  List result = clone(binned);

  SEXP inner_windows = binned.attr("windows"); //the shifted window
  int bin_width = binned.attr("bin.width"); //40
  List nonbinned = binned.attr("nonbinned"); //list of length M, N x 2000 matrixes
  LogicalVector flips = binned.attr("flips");

  Environment IRanges = Environment::namespace_env("BiocGenerics");
  Function width = IRanges["width"];
  Function start = IRanges["start"];
  Function end = IRanges["end"];
  Function paste0("paste0");

  List windows = as<List>(width(inner_windows));
  List starts = as<List>(start(inner_windows));
  List ends = as<List>(end(inner_windows));

  int window = as<int>(windows[0]);

  if (window % bin_width != 0) {
    Rprintf("extract_signal: Window width is not a multiple of bin width!\n");
  }

  for(int l = 0; l < binned.length(); l++) { //loop of chromatin features

    IntegerMatrix nonbinned_datatype = as<IntegerMatrix>(nonbinned[l]); //N x 2000
    IntegerMatrix binned_datatype = as<IntegerMatrix>(result[l]);       //N x S

    for (int i = 0; i < indices.length(); i++) { //only one
      int update_index = indices[i]-1; //convert to 0-based index, required by c++
      IntegerMatrix row = nonbinned_datatype(seq(update_index, update_index),
                                             seq(as<int>(starts[update_index])-1,
                                                 as<int>(ends[update_index])-1));
      row = bin_signal(row, bin_width); //bin the feature vector
      if(flips[update_index]) { //flip
        IntegerVector revtmp = rev(row(0, _));
        row(0, _) = revtmp;
      }
      //shifted, flipped and binned
      binned_datatype(update_index, _) = row;
    }

    List dimnames = binned_datatype.attr("dimnames");

    if (dimnames.size() > 1) {
      binned_datatype.attr("dimnames") = List::create(dimnames[0],
                                          paste0("bin", seq_len(binned_datatype.ncol())));
    }
    result[l] = binned_datatype;
  }

  return result;

}