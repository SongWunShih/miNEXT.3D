#' mixture iNterpolation and EXTrapolation of Hill number
#'
#' \code{miNEXT3D}: Mixture interpolation and extrapolation of Hill number with order q
#'
#' @param data a \code{matrix}, \code{data.frame} (species by assemblages). The first column is the original assemblage, the second column is the transform assemblage. \cr
#' @param diversity selection of diversity type: \code{'TD'} = Taxonomic diversity, \code{'PD'} = Phylogenetic diversity, and \code{'FD'} = Functional diversity.
#' @param knots an integer specifying the number of \code{knots} (say K, default is 11).
#' each knot represents a particular sample size for which mixture diversity estimate will be calculated.
#' @param size an integer vector of sample sizes (number of individuals of original assemblage) for which mixture diversity estimates will be computed.
#' If NULL, then diversity estimates will be computed for those sample sizes determined by the specified/default \code{knots}.
#' @param q a numerical vector specifying the diversity orders. Default is c(0, 1, 2).
#' @param PDtree (required only when \code{diversity = "PD"}), a phylogenetic tree in Newick format for all observed species in the pooled assemblage.
#' @param PDreftime (required only when \code{diversity = "PD"}), a vector of numerical values specifying reference times for PD. Default is \code{NULL} (i.e., the age of the root of PDtree).
#' @param FDdistM (required only when \code{diversity = "FD"}), a species pairwise distance matrix for all species in the pooled assemblage.
#' @param FDtau (required only when \code{diversity = "FD"} and \code{FDtype = "tau_values"}), a numerical vector between 0 and 1 specifying tau values (threshold levels).
#' @param FDtype (required only when \code{diversity = "FD"}), select FD type: \code{FDtype = "tau_values"} for FD under specified threshold values, or \code{FDtype = "AUC"} (area under the curve of tau-profile) for an overall FD which integrates all threshold values between zero and one. Default is \code{"AUC"}
#' @param nboot a positive integer specifying the number of bootstrap replications when assessing sampling uncertainty and constructing confidence intervals. Enter 0 to skip the bootstrap procedures. Default is 0.


#' @import reshape2
#' @import ape
#' @import ggplot2
#' @import dplyr
#' @import iNEXT.3D
#' @import magrittr
#' @import future.apply
#' @import phytools
#' @import iNEXT.beta3D
#' @useDynLib miNEXT.3D, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


#' @return a list of five objects: \code{$out} for showing mixture diversity estimates for different original assemblage probabilitu along with related statistics;
#' \code{$q0_ana} for showing composition information in a mixed sample;
#' \code{$ori.prop} for showing different original assemblage probabilitu;
#' \code{$m.v} for showing different number of individuals of original assemblage 'm1.v' and transform assemblage 'm2.v';
#' and \code{$line_type} for showing the type of mixture sample. \cr
#' 'Rarefaction' means mixture of two rarefaction curves. \cr
#' 'Extrapolation' means mixture of one rarefaction curve and one extrapolation curve.

#' @examples
#' # diversity = 'TD'
#' data(dunes)
#' data = dunes$data
#' data1 = data[,c(2,1)]
#' output1 <- miNEXT3D(data1, diversity = 'TD', nboot = 20)
#' output1
#'
#' # diversity = 'PD'
#' data(dunes)
#' data <- dunes$data
#' data1 = data[,c(2,1)]
#' tree <- dunes$tree
#' output2 <- miNEXT3D(data1, diversity = 'PD', PDtree = tree, nboot = 20)
#' output2
#'
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(dunes)
#' data <- dunes$data
#' data1 = data[,c(2,1)]
#' distM = as.matrix(dunes$dist)
#' output3 <- miNEXT3D(data1, diversity = 'FD', FDdistM = distM, FDtau = 0.3, FDtype = 'tau_values', nboot = 10)
#' output3
#'
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(dunes)
#' data <- dunes$data
#' data1 = data[,c(2,1)]
#' distM = as.matrix(dunes$dist)
#' output4 <- miNEXT3D(data1, diversity = 'FD', FDdistM = distM, FDtype = 'AUC', nboot = 0)
#' output4
#'
#' @export


miNEXT3D = function(data, diversity, knots = 11, size = NULL, q = c(0,1,2),
                     PDtree, PDreftime = NULL, FDdistM, FDtau = NULL, FDtype = 'AUC', nboot = 0){

  est = NULL

  if(diversity == 'TD'){
    est = RTD.est(data, knots = knots, size = size, q = q, nboot = nboot)
  }else if(diversity == 'PD'){
    est = RPD.est(data, PDtree = PDtree, PDreftime = PDreftime,
                  knots = knots, size = size, q = q, nboot = nboot)
  }else if(diversity == 'FD'){
    if(FDtype == 'AUC'){
      est = RFD.est(data, FDdistM = FDdistM, FDtau = FDtau,
                    knots = knots, size = size, q = q, nboot = nboot)
    }else if(FDtype == 'tau_values'){
      est = RFD.singletau.est(data, FDdistM = FDdistM, tau = FDtau,
                              knots = knots, size = size, q = q, nboot = nboot)
    }else{
      print("'FDtype' must be 'AUC' or 'tau_values")
    }
  }else{
    print("The diversity must be 'TD', 'PD' or 'FD'.")
  }

  return(est)
}

#' ggplot2 extension for an iNEXT object
#'
#' \code{ggmiNEXT3D}: the \code{\link[ggplot2]{ggplot}} extension for \code{\link{iNEXT3D}} Object to plot sample-size- and coverage-based rarefaction/extrapolation curves along with a bridging sample completeness curve
#' @param x an \code{iNEXT} object computed by \code{\link{iNEXT}}.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1});
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).
#' @param facet.var create a separate plot for each value of a specified variable:
#'  no separation \cr (\code{facet.var="None"});
#'  a separate plot for each diversity order (\code{facet.var="Order.q"});
#'  a separate plot for each assemblage (\code{facet.var="Assemblage"});
#'  a separate plot for each combination of order x assemblage (\code{facet.var="Both"}).
#' @param color.var create curves in different colors for values of a specified variable:
#'  all curves are in the same color (\code{color.var="None"});
#'  use different colors for diversity orders (\code{color.var="Order.q"});
#'  use different colors for sites (\code{color.var="Assemblage"});
#'  use different colors for combinations of order x assemblage (\code{color.var="Both"}).
#' @param ... other arguments passed on to methods. Not currently used.
#' @return a ggplot2 object

















