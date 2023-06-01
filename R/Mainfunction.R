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



#' @return a list of five objects: \code{$Mixture} for showing mixture diversity estimates for different size of Original and Transform assemblagt along with related statistics;
#' \code{$Original} for showing Original diversity estimates for the size of Original assemblage ;
#' \code{$Transform} for showing Transform diversity estimates for the size of Transform assemblage ;
#' \code{$q0_ana} for showing composition information in a mixed sample;


#' @examples
#' # diversity = 'TD'
#' data(dunes)
#' data <- dunes$data
#' data1 <- data[,c(2,1)]
#' output1 <- miNEXT3D(data1, diversity = 'TD', nboot = 10)
#' output1
#'
#' # diversity = 'PD'
#' data(dunes)
#' data <- dunes$data
#' data1 <- data[,c(2,1)]
#' tree <- dunes$tree
#' output2 <- miNEXT3D(data1, diversity = 'PD', PDtree = tree, nboot = 10)
#' output2
#'
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(dunes)
#' data <- dunes$data
#' data1 <- data[,c(2,1)]
#' distM <- as.matrix(dunes$dist)
#' output3 <- miNEXT3D(data1, diversity = 'FD', FDdistM = distM, FDtau = 0.3, FDtype = 'tau_values', nboot = 10)
#' output3
#'
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(dunes)
#' data <- dunes$data
#' data1 <- data[,c(2,1)]
#' distM <- as.matrix(dunes$dist)
#' output4 <- miNEXT3D(data1, diversity = 'FD', FDdistM = distM, FDtype = 'AUC', nboot = 0)
#' output4
#'
#' @export

miNEXT3D = function(data, diversity, knots = 11, size = NULL, q = c(0,1,2),
                     PDtree, PDreftime = NULL, FDdistM, FDtau = NULL, FDtype = 'AUC', nboot = 0){

  est = NULL

  if(diversity == 'TD'){
    est = RTD.est(data, knots = knots, size = size, q = q, nboot = nboot)

    m_v = est$m.v
    m_v[m_v == 0] = 1 #避免iNEXT3D跑 m = 0的情況

    est.ori = iNEXT3D(apply(data, 1, as.vector)[1,], diversity = 'TD',
                      size = m_v[,1],nboot = nboot)
    est.trans = iNEXT3D(apply(data, 1, as.vector)[2,], diversity = 'TD',
                      size = m_v[,2],nboot = nboot)

  }else if(diversity == 'PD'){
    est = RPD.est(data, PDtree = PDtree, PDreftime = PDreftime,
                  knots = knots, size = size, q = q, nboot = nboot)
    m_v = est$m.v
    m_v[m_v == 0] = 1 #避免iNEXT3D跑 m = 0的情況
    est.ori = iNEXT3D(apply(data, 1, as.vector)[1,], diversity = 'PD',PDtree = PDtree,
                      size = m_v[,1],nboot = nboot)
    est.trans = iNEXT3D(apply(data, 1, as.vector)[2,], diversity = 'PD',PDtree = PDtree,
                        size = m_v[,2],nboot = nboot)
  }else if(diversity == 'FD'){
    if(FDtype == 'AUC'){
      est = RFD.est(data, FDdistM = FDdistM, FDtau = FDtau,
                    knots = knots, size = size, q = q, nboot = nboot)
      m_v = est$m.v
      m_v[m_v == 0] = 1 #避免iNEXT3D跑 m = 0的情況
      est.ori = iNEXT3D(apply(data, 1, as.vector)[1,], diversity = 'FD',FDdistM = FDdistM,
                        size = m_v[,1],nboot = nboot)
      est.trans = iNEXT3D(apply(data, 1, as.vector)[2,], diversity = 'FD',FDdistM = FDdistM,
                          size = m_v[,2],nboot = nboot)
    }else if(FDtype == 'tau_values'){
      est = RFD.singletau.est(data, FDdistM = FDdistM, tau = FDtau,
                              knots = knots, size = size, q = q, nboot = nboot)
      m_v = est$m.v
      m_v[m_v == 0] = 1 #避免iNEXT3D跑 m = 0的情況
      est.ori = iNEXT3D(apply(data, 1, as.vector)[1,], diversity = 'FD', FDdistM = FDdistM,
                        FDtype = FDtype, FDtau = FDtau,
                        size = m_v[,1],nboot = nboot)
      est.trans = iNEXT3D(apply(data, 1, as.vector)[2,], diversity = 'FD', FDdistM = FDdistM,
                          FDtype = FDtype, FDtau = FDtau,
                          size = m_v[,2],nboot = nboot)
    }else{
      print("'FDtype' must be 'AUC' or 'tau_values")
      break
    }
  }else{
    print("The diversity must be 'TD', 'PD' or 'FD'.")
    break
  }

  clean.real.fn = function(est,est.ori,est.trans,data,diversity,FDtype = 'AUC'){
    p.v = est$ori.prop
    est.final = est$out[order(est$out$Order.q),]
    rownames(est.final) = NULL
    est.final$m1 = rep(est$m.v[,1],length(unique(est.final$Order.q)))
    est.final$m2 = rep(est$m.v[,2],length(unique(est.final$Order.q)))
    if(sum(data[,1])>=sum(data[,2])){
      checkpoint = round(((sum(data[,1])-sum(data[,2]))/sum(data[,1]))*100,3)
      est.final$linetype = est.final$points = NA
      est.final$linetype[which(est.final$prop.v >= checkpoint)] = 'Rarefaction'
      est.final$linetype[which(est.final$prop.v < checkpoint)] = 'Extrapolation'
      est.final$points[which(est.final$prop.v > checkpoint)] = 'Rarefaction'
      est.final$points[which(est.final$prop.v < checkpoint)] = 'Extrapolation'
      est.final$points[which(est.final$prop.v == checkpoint)] = 'Observed'
      est.final$Method = "Mixture"
    }else{
      est.final$linetype = est.final$points = NA
      est.final$linetype = 'Rarefaction'
      est.final$points = 'Rarefaction'
      est.final$Method = "Mixture"
    }


    # if(sum(data[,1]) > sum(data[,2])){
    #   est.final = data.frame(est.final,
    #                          linetype = linetype,
    #                          points = points)
    #
    # }else{
    #   est.final = data.frame(est.final,linetype = "rarefaction")
    # }
    if(diversity == 'TD'){
      est.final = select(est.final,Order.q,m1,m2,prop.v,qmiNEXT_TD,LCL,UCL,points,linetype,Method)
      tmp = est.ori$iNextEst$size_based
      type = tmp$Method
      type[type == 'Observed'] = "Rarefaction"
      est.ori.df = data.frame(Order.q = tmp$Order.q,
                              m1 = tmp$m,
                              prop.v = round(tmp$m/sum(data[,1])*100,3),
                              qmiNEXT_TD = tmp$qD,
                              LCL = tmp$qD.LCL,
                              UCL = tmp$qD.UCL,
                              points = tmp$Method,
                              linetype = type,
                              Method = 'Original')

      tmp = est.trans$iNextEst$size_based
      type = tmp$Method
      type[type == 'Observed'] = "Rarefaction"
      est.trans.df = data.frame(Order.q = tmp$Order.q,
                                m2 = tmp$m,
                                prop.v = round(tmp$m/sum(data[,1])*100,3),
                                qmiNEXT_TD = tmp$qD,
                                LCL = tmp$qD.LCL,
                                UCL = tmp$qD.UCL,
                                points = tmp$Method,
                                linetype = type,
                                Method = 'Transform')
    }else if(diversity == 'PD'){
      est.final = select(est.final,Order.q,m1,m2,prop.v,qmiNEXT_PD,LCL,UCL,points,linetype,Method)
      tmp = est.ori$PDiNextEst$size_based
      type = tmp$Method
      type[type == 'Observed'] = "Rarefaction"
      est.ori.df = data.frame(Order.q = tmp$Order.q,
                              m1 = tmp$m,
                              prop.v = round(tmp$m/sum(data[,1])*100,3),
                              qmiNEXT_PD = tmp$qPD,
                              LCL = tmp$qPD.LCL,
                              UCL = tmp$qPD.UCL,
                              points = tmp$Method,
                              linetype = type,
                              Method = 'Original')

      tmp = est.trans$PDiNextEst$size_based
      type = tmp$Method
      type[type == 'Observed'] = "Rarefaction"
      est.trans.df = data.frame(Order.q = tmp$Order.q,
                                m2 = tmp$m,
                                prop.v = round(tmp$m/sum(data[,1])*100,3),
                                qmiNEXT_PD = tmp$qPD,
                                LCL = tmp$qPD.LCL,
                                UCL = tmp$qPD.UCL,
                                points = tmp$Method,
                                linetype = type,
                                Method = 'Transform')
    }else if(diversity == 'FD' & FDtype == 'AUC'){
      est.final = select(est.final,Order.q,m1,m2,prop.v,qmiNEXT_FD,LCL,UCL,points,linetype,Method)
      tmp = est.ori$AUCiNextEst$size_based
      type = tmp$Method
      type[type == 'Observed'] = "Rarefaction"
      est.ori.df = data.frame(Order.q = tmp$Order.q,
                              m1 = tmp$m,
                              prop.v = round(tmp$m/sum(data[,1])*100,3),
                              qmiNEXT_FD = tmp$qAUC,
                              LCL = tmp$qAUC.LCL,
                              UCL = tmp$qAUC.UCL,
                              points = tmp$Method,
                              linetype = type,
                              Method = 'Original')

      tmp = est.trans$AUCiNextEst$size_based
      type = tmp$Method
      type[type == 'Observed'] = "Rarefaction"
      est.trans.df = data.frame(Order.q = tmp$Order.q,
                                m2 = tmp$m,
                                prop.v = round(tmp$m/sum(data[,1])*100,3),
                                qmiNEXT_FD = tmp$qAUC,
                                LCL = tmp$qAUC.LCL,
                                UCL = tmp$qAUC.UCL,
                                points = tmp$Method,
                                linetype = type,
                                Method = 'Transform')
    }else if(diversity == 'FD' & FDtype == 'tau_values'){
      est.final = select(est.final,Order.q,m1,m2,prop.v,qmiNEXT_FD_singletau,LCL,UCL,points,linetype,Method)
      tmp = est.ori$FDiNextEst$size_based
      type = tmp$Method
      type[type == 'Observed'] = "Rarefaction"
      est.ori.df = data.frame(Order.q = tmp$Order.q,
                              m1 = tmp$m,
                              prop.v = round(tmp$m/sum(data[,1])*100,3),
                              qmiNEXT_FD_singletau = tmp$qFD,
                              LCL = tmp$qFD.LCL,
                              UCL = tmp$qFD.UCL,
                              points = tmp$Method,
                              linetype = type,
                              Method = 'Original')

      tmp = est.trans$FDiNextEst$size_based
      type = tmp$Method
      type[type == 'Observed'] = "Rarefaction"
      est.trans.df = data.frame(Order.q = tmp$Order.q,
                                m2 = tmp$m,
                                prop.v = round(tmp$m/sum(data[,1])*100,3),
                                qmiNEXT_FD_singletau = tmp$qFD,
                                LCL = tmp$qFD.LCL,
                                UCL = tmp$qFD.UCL,
                                points = tmp$Method,
                                linetype = type,
                                Method = 'Transform')
    }
    est.trans.df = dplyr::filter(est.trans.df,prop.v <= 100)


    return(list(Mixture = est.final, Orignal = est.ori.df, Transform = est.trans.df))

  }
  final = clean.real.fn(est,est.ori,est.trans,data,diversity,FDtype = FDtype)
  return(final)
}



#' ggplot2 extension for an miNEXT3D object
#'
#' \code{ggmiNEXT3D}: the \code{\link[ggplot2]{ggplot}} extension for \code{\link{miNEXT3D}} Object to plot sized-based rarefaction/extrapolation curves
#' @param out an \code{miNEXT3D} object computed by \code{\link{miNEXT3D}}.
#' @return a ggplot2 object
#'
#' @examples
#' # diversity = 'TD'
#' data(dunes)
#' data <- dunes$data
#' data1 <- data[,c(2,1)]
#' output1 <- miNEXT3D(data1, diversity = 'TD', nboot = 10)
#' plot1 <- ggmiNEXT3D(output1)
#' plot1
#'
#' # diversity = 'PD'
#' data(dunes)
#' data <- dunes$data
#' data1 <- data[,c(2,1)]
#' tree <- dunes$tree
#' output2 <- miNEXT3D(data1, diversity = 'PD', PDtree = tree, nboot = 10)
#' plot2 <- ggmiNEXT3D(output2)
#' plot2
#'
#' # diversity = 'FD' & FDtype = 'tau_values'
#' data(dunes)
#' data <- dunes$data
#' data1 <- data[,c(2,1)]
#' distM <- as.matrix(dunes$dist)
#' output3 <- miNEXT3D(data1, diversity = 'FD', FDdistM = distM, FDtau = 0.3, FDtype = 'tau_values', nboot = 10)
#' plot3 <- ggmiNEXT3D(output3)
#' plot3
#'
#' # diversity = 'FD' & FDtype = 'AUC'
#' data(dunes)
#' data <- dunes$data
#' data1 <- data[,c(2,1)]
#' distM <- as.matrix(dunes$dist)
#' output4 <- miNEXT3D(data1, diversity = 'FD', FDdistM = distM, FDtype = 'AUC', nboot = 0)
#' plot4 <- ggmiNEXT3D(output4)
#' plot4
#'
#' @export
#'
#'
#'
ggmiNEXT3D = function(out){
  tmp1 = subset(out$Mixture,select = -c(m1,m2))
  tmp2 = subset(out$Orignal,select = -c(m1))
  tmp3 = subset(out$Transform,select = -c(m2))
  final = rbind(tmp1,tmp2,tmp3)
  final$Method = factor(final$Method,levels = c('Mixture','Original','Transform'))
  ylabtmp = colnames(final)[3]
  colnames(final)[3] = "miNEXT"
  data.sub = final[which(final$points == "Observed"),]

  library(ggplot2)
  g = ggplot(final,aes(x = prop.v, y = miNEXT)) +
    geom_line(aes(col = Method,lty = linetype),size = 1.5)+
    geom_ribbon(aes(x = prop.v, ymin = LCL, ymax = UCL,fill = Method), alpha = 0.3)+
    geom_point(aes(col = Method),size = 3,data = data.sub)+
    facet_wrap(~Order.q,scales = "free_y",labeller = labeller(Order.q = c(`0` = "q = 0", `1` = "q = 1", `2` = "q = 2"))) +
    scale_linetype_manual(values=c("Rarefaction" = "solid","Extrapolation" = "dashed"))+
    scale_color_manual(values=c("Mixture"= "#C0392B","Original" = "black", "Transform" = "#3498DB"))+
    scale_fill_manual(values=c("Mixture"= "#C0392B","Original" = "black", "Transform" = "#3498DB"))+
    scale_x_continuous(name = paste0('Number of individuals',"\n","(Proportion % in original habitat)"),
                       breaks = round(seq(0,100,length.out = 5)),
                       labels = paste0(round(seq(0,max(out$Mixture$m1),length.out = 5)), "\n(", round(seq(0,100,length.out = 5), 1), "%)"))+
    theme(legend.position = "bottom",
          legend.box = "vertical",
          legend.title = element_blank(),
          legend.key.width = unit(1.2, "cm"),
          legend.margin = margin(0, 0, 0, 0),
          legend.box.margin = margin(-10, -10, -5, -10),
          text = element_text(size = 12),
          plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
          strip.text = element_text(size = 14, face = "bold"))+
    ylab(ylabtmp)

  return(g)
}












