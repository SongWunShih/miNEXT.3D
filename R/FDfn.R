DetAbu <- function(x, zero=FALSE){
  x <- unlist(x)
  n <- sum(x)
  f1 <- sum(x==1)
  f2 <- sum(x==2)
  f3 <- sum(x==3)
  if(f2==0){
    f1 <- max(f1 - 1, 0)
    f2 <- 1
  }
  A1 <- f1 / n * ((n-1)*f1 / ((n-1)*f1 + 2*max(f2,1)))
  A2 <- f2 / choose(n, 2) * ((n-2)*f2 / ((n-2)*f2 + 3*max(f3,1)))^2
  if(zero==FALSE) x <- x[x>0]
  q.solve <- function(q){
    e <- A1 / sum(x/n*exp(-q*x))
    out <- sum((x/n * (1 - e * exp(-q*x)))^2) - sum(choose(x,2)/choose(n,2)) + A2
    abs(out)
  }
  q <- tryCatch(optimize(q.solve, c(0,1))$min, error = function(e) {1})
  e <- A1 / sum(x/n*exp(-q*x))
  o <- x/n * (1 - e * exp(-q*x))
  o
}

checkFD = function (data, distM) {
  distM = as.matrix(distM)

  if (is.null(rownames(data)) | is.null(rownames(distM))) {
    warning("The species names are not provided in data or distance matrix.",
            call. = FALSE)
    rownames(data) <- rownames(distM) <- colnames(distM) <- paste0("Species",
                                                                   1:nrow(data))
  }
  else {
    if (sum(rownames(data) %in% rownames(distM)) != nrow(data))
      stop("Data and distance matrix contain unmatched species",
           call. = FALSE)
  }
  distM = distM[rownames(distM) %in% rownames(data), colnames(distM) %in%
                  rownames(data)]
  if (nrow(data) != nrow(distM))
    stop("The number of species in data should equal to that in distance matrix",
         call. = FALSE)
  order_sp <- match(rownames(data), rownames(distM))
  distM <- distM[order_sp, order_sp]
  distM <- distM[rowSums(data) > 0, rowSums(data) > 0]
  data <- data[rowSums(data) > 0, , drop = FALSE]
  return(list(checkdistM = distM, checkdata =  data))
}

Bootstrap_distance_matrix <- function (data, distance_matrix, f0.hat, datatype){
  if (datatype == "incidence_freq") {
    n = data[1]
    X = data[-1]
    u = sum(data)
  }
  else if (datatype == "abundance") {
    n = sum(data)
    X = data
  }
  distance = as.matrix(distance_matrix)
  dij = distance
  F.1 <- sum(dij[, X == 1])
  F.2 <- sum(dij[, X == 2])
  F11 <- sum(dij[X == 1, X == 1])
  F22 <- sum(dij[X == 2, X == 2])
  if (datatype == "abundance") {
    F.0hat <- ifelse(F.2 > 0, ((n - 1)/n) * (F.1^2/(2 *
                                                      F.2)), ((n - 1)/n) * (F.1 * (F.1 - 0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n - 2) * (n - 3) * (F11^2)/(4 *
                                                              n * (n - 1) * F22)), ((n - 2) * (n - 3) * (F11 *
                                                                                                           (F11 - 0.01))/(4 * n * (n - 1))))
  }
  else if (datatype == "incidence_freq") {
    F.0hat <- ifelse(F.2 > 0, ((n - 1)/n) * (F.1^2/(2 *
                                                      F.2)), ((n - 1)/n) * (F.1 * (F.1 - 0.01)/(2)))
    F00hat <- ifelse(F22 > 0, ((n - 1)^2 * (F11^2)/(4 *
                                                      n * n * F22)), ((n - 1) * (n - 1) * (F11 * (F11 -
                                                                                                    0.01))/(4 * n * n)))
  }
  if (f0.hat == 0) {
    d = dij
  }
  else if (f0.hat == 1) {
    d.0bar <- matrix(rep(F.0hat/length(X)/f0.hat, length(X) *
                           f0.hat), length(X), f0.hat)
    d00 = matrix(0, f0.hat, f0.hat)
    d <- cbind(dij, d.0bar)
    aa <- cbind(t(d.0bar), d00)
    d <- rbind(d, aa)
    diag(d) = 0
  }
  else {
    d.0bar <- matrix(rep(F.0hat/length(X)/f0.hat, length(X) *
                           f0.hat), length(X), f0.hat)
    fo.num = (f0.hat * (f0.hat - 1))/2
    d00 = matrix(0, f0.hat, f0.hat)
    d00[upper.tri(d00)] = (F00hat/2)/fo.num
    d00 <- pmax(d00, t(d00))
    d <- cbind(dij, d.0bar)
    aa <- cbind(t(d.0bar), d00)
    d <- rbind(d, aa)
    diag(d) = 0
  }
  return(d)
}

FDq2<-function(x1,x2,vi.all,m1,m2,n1,n2){
  x1 = x1[vi.all>0]
  x2 = x2[vi.all>0]
  vi.all = vi.all[vi.all>0]

  D2 = 0
  if(m1 !=0 || m2!=0){
    part1 = 1/(m1+m2)^2*sum(m1/n1*vi.all*x1 + m2/n2*vi.all*x2)
    part2 = (m1*(m1-1))/(m1+m2)^2*sum(vi.all*x1*(x1-1))/(n1*(n1-1))
    part3 = (m2*(m2-1))/(m1+m2)^2*sum(vi.all*x2*(x2-1))/(n2*(n2-1))
    part4 = 2*m1*m2/(m1+m2)^2*sum(vi.all*x1*x2)/(n1*n2)
    D2 = 1/(part1 + part2 + part3 + part4)

  }
  else if(m1==0 && m2==0)
    D2<-0

  return(D2)
}

# FD ===============================================================================
FD.h.theo_fn = function(m1,m2,ai1,ai2,vi,S){
  gtheo.tmp = matrix(0,m1+1,m2+1)
  for (i in 1:S) {
    for (k1 in 0:m1) {
      for (k2 in 0:m2) {
        if((k1+k2) == 0){
          gtheo.tmp[k1+1,k2+1] = gtheo.tmp[k1+1,k2+1]
        }else{
          gtheo.tmp[k1+1,k2+1] = gtheo.tmp[k1+1,k2+1] + vi[i] * exp(lchoose(m1,k1)+lchoose(m2,k2)) * (ai1[i]^k1) * ((1-ai1[i])^(m1-k1)) * (ai2[i]^k2) * ((1-ai2[i])^(m2-k2))
        }
      }
    }
  }

  return(gtheo.tmp)
}
## auc
mix.aivi.popu.FD = function(data, m1, m2, FDdistM, tau){

  FDdistM = as.matrix(FDdistM)

  data.list = list()
  for (i in 1:ncol(data)) {
    data.list[[i]] = data[,i]
    names(data.list[[i]]) = rownames(data)
  }
  data.sum = (m1/(m1+m2))*data[,1] + (m2/(m1+m2))*data[,2]
  names(data.sum) = rownames(data)

  ai_vi.popu.fn = function (data, dij, tau1, filt_zero = FALSE) {
    # if (filt_zero) {
    #   dij <- dij[data > 0, data > 0]
    #   data <- data[data > 0]
    # }
    out <- lapply(tau1, function(tau_) {
      dij_ <- dij
      if (tau_ == 0) {
        dij_[dij_ > 0] <- 1
        a <- as.vector((1 - dij_/1) %*% data)
      }
      else {
        dij_[which(dij_ > tau_, arr.ind = T)] <- tau_
        a <- as.vector((1 - dij_/tau_) %*% data)
      }
      if (filt_zero) {
        data <- data[a != 0]
        a <- a[a != 0]
      }

      v <- data/a
      data.frame(ai = a, vi = v)
    })
    out_a <- matrix(sapply(out, function(x) x[, 1]), ncol = length(tau1))
    out_v <- matrix(sapply(out, function(x) x[, 2]), ncol = length(tau1))
    colnames(out_a) <- colnames(out_v) <- paste0("tau_", round(tau1,
                                                               3))
    output = list(ai = out_a, vi = out_v)

    return(output)
  }
  aivi_1 = ai_vi.popu.fn(data = data.list[[1]],dij = FDdistM, tau1 = tau, filt_zero = F)
  aivi_2 = ai_vi.popu.fn(data = data.list[[2]],dij = FDdistM, tau1 = tau, filt_zero = F)
  aivi_sum = ai_vi.popu.fn(data = data.sum,dij = FDdistM, tau1 = tau,filt_zero = F)
  aivi_sum$vi[is.nan(aivi_sum$vi)] = 0


  return(list(ai1_mat = aivi_1$ai, ai2_mat = aivi_2$ai, vi_all_mat = aivi_sum$vi))
}
RFD.obs <- function(data, FDdistM, FDtau = NULL, knots = 11, size = NULL, q = c(0,1,2), n1, n2) {



  if (is.null(FDtau)) {
    tau <- seq(0, 1, length.out = 30)
  }else{
    tau <- FDtau
  }

  if(n1>n2){
    a = round(seq(n2,n1, length.out = knots - round(knots/2) + 1))[-1]
    m2.v = c(round(seq(0, n2, length.out = round(knots/2))),a)
    m1.v = n1 - m2.v
    if(!is.null(size)){
      m1.v <- size
      m2.v <- n1- m1.v
    }
    m.v = cbind(m1.v, m2.v)
    inc.idx = which(m2.v<=n2)
    obs.idx = which(m2.v==n2)
    ext.idx = which(m2.v>n2)
    line_type = rep(NA,length(m1.v))
    line_type[inc.idx] = 'Rarefaction'
    line_type[obs.idx] = 'Observe'
    line_type[ext.idx] = 'Extrapolation'
    prop.v = round(((m1.v/n1) * 100),3)

  }else{
    m2.v = round(seq(0, n1, length.out = knots))
    m1.v = n1 - m2.v
    if(!is.null(size)){
      m1.v <- size
      m2.v <- n1- m1.v
    }
    m.v = cbind(m1.v, m2.v)
    prop.v = round(((m1.v/n1) * 100),3)
    line_type = rep("Rarefaction",knots)
  }



  # 不同的m1, m2有不同的ai_all, vi_all。 ai1, ai2則相同
  aivi.popu.data = lapply(1:nrow(m.v), function (k) mix.aivi.popu.FD(data, m1.v[k], m2.v[k], FDdistM = FDdistM,
                                                                     tau = tau))

  output.list = lapply(1:length(tau), function(tt) {
    # Li.v = apply(Lis, 1, function(x) x[1]*x[2])
    S <- nrow(data)
    qlength <- length(q)
    out = list()
    ghat_pt2 = list()


    for (zz in 1:length(m1.v)) {
      ai1.v = aivi.popu.data[[zz]]$ai1_mat[,tt]
      ai2.v = aivi.popu.data[[zz]]$ai2_mat[,tt]
      vi.all.v = aivi.popu.data[[zz]]$vi_all_mat[,tt]
      # Vbar = sum(ai.all.v * vi.all.v)/(n1 + n2)
      Vbar = 1
      out[[zz]] <- matrix(0, qlength, 1)
      m1 = m1.v[zz]
      m2 = m2.v[zz]

      ghat_pt2[[zz]] <- FD.h.theo_fn(m1, m2, ai1.v, ai2.v ,vi.all.v, S)




      for (j in 1:qlength) {
        for (k1 in 0:m1) {
          for (k2 in 0:m2) {
            if (q[j] == 0) {
              out[[zz]][j,1] <- out[[zz]][j,1] + ghat_pt2[[zz]][k1+1,k2+1]
            }else if (q[j] == 1) {
              if(k1 == 0 & k2 == 0){
                aa = 0
              }else{
                aa = ((k1+k2)/((m1+m2)*Vbar)) * log((k1+k2)/((m1+m2)*Vbar))
              }
              out[[zz]][j,1] <- out[[zz]][j,1] - aa  * ghat_pt2[[zz]][k1+1,k2+1]
            }else if (q[j] == 2) {
              out[[zz]][j,1] <- out[[zz]][j,1] + ((k1+k2)/((m1+m2)*Vbar))^2 * ghat_pt2[[zz]][k1+1,k2+1]
            }else {
              out[[zz]][j,1] <- out[[zz]][j,1] + ((k1+k2)/((m1+m2)*Vbar))^q[j] * ghat_pt2[[zz]][k1+1,k2+1]
            }
          }

        }

        if (q[j] == 0) {
          out[[zz]][j,1] <- out[[zz]][j,1]
        } else if (q[j] == 1) {
          out[[zz]][j,1] <- exp(out[[zz]][j,1])
        } else if (q[j] == 2) {
          out[[zz]][j,1] <- 1 / out[[zz]][j,1]
        } else {
          out[[zz]][j,1] <- out[[zz]][j,1]^(1/(1-q[j]))
        }
      }
    }

    output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
    rownames(output) = q
    colnames(output) = prop.v
    output = reshape2::melt(output)
    output$tau = tau[tt]
    colnames(output) = c("Order.q","prop.v","qFD_mix","tau")
    output
  })
  output.bind = output.list[[1]]
  for (j in 2:length(output.list)) {
    output.bind = rbind(output.bind,output.list[[j]])
  }
  auc.table = matrix(NA,length(q),length(prop.v))
  p.v = prop.v
  for (i in 1:length(q)) {
    for (j in 1:length(prop.v)) {
      tmp1 = dplyr::filter(output.bind,Order.q == q[i] & prop.v == p.v[j]) # filter的名字不能重複，而且一定要加dplyr
      auc_1 = sum(tmp1$qFD_mix[seq_along(tmp1$qFD_mix[-1])] * diff(tau))
      auc_2 = sum(tmp1$qFD_mix[-1] * diff(tau))
      auc.table[i,j] = (auc_1 + auc_2)/2
    }
  }
  rownames(auc.table) = q
  colnames(auc.table) = p.v
  auc.table = reshape2::melt(auc.table)
  colnames(auc.table) = c("Order.q","prop.v","qmiNEXT_FD")
  return(list(out = auc.table,ori.prop = prop.v,m.v = m.v,line_type = line_type))
}
mix.aivi.sample.FD = function(data, m1, m2, FDdistM, tau){

  FDdistM = as.matrix(FDdistM)

  data.list = list()
  for (i in 1:ncol(data)) {
    data.list[[i]] = data[,i]
    names(data.list[[i]]) = rownames(data)
  }
  data.sum = (m1/(m1+m2))*data[,1] + (m2/(m1+m2))*data[,2]
  names(data.sum) = rownames(data)

  ai_vi.sample.fn = function (data, dij, tau1, filt_zero = FALSE, integer = TRUE) {
    # if (filt_zero) {
    #   dij <- dij[data > 0, data > 0]
    #   data <- data[data > 0]
    # }
    out <- lapply(tau1, function(tau_) {
      dij_ <- dij
      if (tau_ == 0) {
        dij_[dij_ > 0] <- 1
        a <- as.vector((1 - dij_/1) %*% data)
      }
      else {
        dij_[which(dij_ > tau_, arr.ind = T)] <- tau_
        a <- as.vector((1 - dij_/tau_) %*% data)
      }
      if (filt_zero) {
        data <- data[a != 0]
        a <- a[a != 0]
      }

      if (integer){
        a[a < 1 & a > 0] <- 1
        a <- round(a)
      }

      v <- data/a
      data.frame(ai = a, vi = v)
    })
    out_a <- matrix(sapply(out, function(x) x[, 1]), ncol = length(tau1))
    out_v <- matrix(sapply(out, function(x) x[, 2]), ncol = length(tau1))
    colnames(out_a) <- colnames(out_v) <- paste0("tau_", round(tau1,
                                                               3))
    output = list(ai = out_a, vi = out_v)

    return(output)
  }
  aivi_1 = ai_vi.sample.fn(data = data.list[[1]],dij = FDdistM, tau1 = tau, filt_zero = F)
  aivi_2 = ai_vi.sample.fn(data = data.list[[2]],dij = FDdistM, tau1 = tau, filt_zero = F)
  aivi_sum = ai_vi.sample.fn(data = data.sum,dij = FDdistM, tau1 = tau, filt_zero = F, integer = F)
  aivi_sum$vi[is.nan(aivi_sum$vi)] = 0
  # ai1_vec = list()
  # ai2_vec = list()
  # vi_all_vec = list()
  # for (i in 1:length(aivi_1)) {
  #   ai1_vec[[i]] = aivi_1[[i]]$ai
  #   ai2_vec[[i]] = aivi_2[[i]]$ai
  #   vi_all_vec[[i]] = aivi_sum[[i]]$vi
  # }
  # lapply(aivi_sum, function(k) sum(k[,1]*k[,2]))

  return(list(ai1_mat = aivi_1$ai, ai2_mat = aivi_2$ai, vi_all_mat = aivi_sum$vi))
}
RFD.est <- function(data, FDdistM, FDtau = NULL, knots = 11, size = NULL, q = c(0,1,2), conf = 0.95, nboot = 0) {

  # check = checkFD(data,FDdistM)
  #
  # FDdistM = check$checkdistM
  # data = check$checkdata
  n1 = sum(data[,1])
  n2 = sum(data[,2])

  if (is.null(FDtau)) {
    tau <- seq(0, 1, length.out = 30)
  }else{
    tau <- FDtau
  }

  if(n1>n2){
    a = round(seq(n2,n1, length.out = knots - round(knots/2) + 1))[-1]
    m2.v = c(round(seq(0, n2, length.out = round(knots/2))),a)
    m1.v = n1 - m2.v
    if(!is.null(size)){
      m1.v <- size
      m2.v <- n1- m1.v
    }
    m.v = cbind(m1.v, m2.v)
    prop.v = round(((m1.v/n1) * 100),3)
    inc.idx = which(m2.v<=n2)
    obs.idx = which(m2.v==n2)
    ext.idx = which(m2.v>n2)
    line_type = rep(NA,length(m1.v))
    line_type[inc.idx] = 'Rarefaction'
    line_type[obs.idx] = 'Observe'
    line_type[ext.idx] = 'Extrapolation'

    m.inc.v = cbind(m1.v[inc.idx],m2.v[inc.idx])
    prop.inc = prop.v[inc.idx]
    m.ext.v = cbind(m1.v[-inc.idx],m2.v[-inc.idx])
    prop.ext = prop.v[-inc.idx]
    p1 = DetAbu(data[,1], zero = T)
    p2 = DetAbu(data[,2], zero = T)
    p_bind = cbind(p1,p2)
    rownames(p_bind) = rownames(data)
  }else{
    m2.v = round(seq(0, n1, length.out = knots))
    m1.v = n1 - m2.v
    if(!is.null(size)){
      m1.v <- size
      m2.v <- n1- m1.v
    }
    m.v = cbind(m1.v, m2.v)
    prop.v = round(((m1.v/n1) * 100),3)
    line_type = rep("Rarefaction",length(m1.v))
  }

  # 不同的m1, m2有不同的ai_all, vi_all。 ai1, ai2則相同
  aivi.sample.data = lapply(1:nrow(m.v), function (k) mix.aivi.sample.FD(data, m1.v[k], m2.v[k],FDdistM = FDdistM, tau = tau))



  if(n1<=n2){

    output.list = lapply(1:length(tau), function(tt) {
      # Li.v = apply(Lis, 1, function(x) x[1]*x[2])
      # S <- nrow(data)
      qlength <- length(q)
      out = list()
      un1 = un2 = sh12 = c()

      for (zz in 1:length(m1.v)) {
        ai1.v = aivi.sample.data[[zz]]$ai1_mat[,tt]
        ai2.v = aivi.sample.data[[zz]]$ai2_mat[,tt]
        # ai.all.v = aivi.sample.data[[zz]]$ai_all_mat[,tt]
        vi.all.v = aivi.sample.data[[zz]]$vi_all_mat[,tt]

        out[[zz]] <- matrix(0, qlength, 1)
        m1 = m1.v[zz]
        m2 = m2.v[zz]

        for (j in 1:qlength) {
          if (q[j] == 0) {
            out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])
          } else if (q[j] == 1) {
            out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]))
          } else if (q[j] == 2) {
            out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
          } else {
            out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
          }
        }

        ## q0 ana =========================================
        datanew = cbind(ai1.v,ai2.v)
        rownames(datanew) = names(vi.all.v) = rownames(data)
        datash = datanew[(data[,1] > 0 & data[,2]>0), , drop=F]
        dataun1 = datanew[(data[,1] > 0 & data[,2] == 0), , drop=F]
        dataun2 = datanew[(data[,1] == 0 & data[,2] > 0), , drop=F]
        ## 都用 sh_abun_FD因為在單群落下ai1, ai2都不一定非0
        un1[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun1)],xi1 = dataun1[,1],xi2 = dataun1[,2], n1 = n1,
                                 m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
        un2[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun2)],xi1 = dataun2[,1],xi2 = dataun2[,2], n1 = n1,
                                 m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
        sh12[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(datash)],xi1 = datash[,1],xi2 = datash[,2], n1 = n1,
                                  m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
      }

      output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
      rownames(output) = q
      colnames(output) = prop.v
      output = reshape2::melt(output)
      output$tau = tau[tt]
      colnames(output) = c("Order.q","prop.v","qFD_mix","tau")

      # q0_ana ===============
      q0_ana = data.frame(prop.v = prop.v ,
                          q0_un1 = un1, q0_un2 = un2,
                          q0_sh = sh12)
      q0_ana = reshape2::melt(q0_ana,id = c('prop.v'))
      colnames(q0_ana)[2:3] = c('Type','Value')
      q0_ana$tau = tau[tt]

      list(output = output,q0_ana_ = q0_ana)
    })
    ## auc ================================================
    output.bind = output.list[[1]]$output
    for (j in 2:length(output.list)) {
      output.bind = rbind(output.bind,output.list[[j]]$output)
    }
    auc.table = matrix(NA,length(q),length(prop.v))
    p.v = prop.v
    for (i in 1:length(q)) {
      for (j in 1:length(prop.v)) {
        tmp1 = dplyr::filter(output.bind,Order.q == q[i] & prop.v == p.v[j]) # filter的名字不能重複，而且一定要加dplyr
        auc_1 = sum(tmp1$qFD_mix[seq_along(tmp1$qFD_mix[-1])] * diff(tau))
        auc_2 = sum(tmp1$qFD_mix[-1] * diff(tau))
        auc.table[i,j] = (auc_1 + auc_2)/2
      }
    }
    rownames(auc.table) = q
    colnames(auc.table) = p.v
    auc.table = reshape2::melt(auc.table)
    colnames(auc.table) = c("Order.q","prop.v","qmiNEXT_FD")

    ## q0_ana auc ===============
    q0_ana.bind = output.list[[1]]$q0_ana_
    for (j in 2:length(output.list)) {
      q0_ana.bind = rbind(q0_ana.bind,output.list[[j]]$q0_ana_)
    }
    com_type = unique(q0_ana.bind$Type)
    q0_ana_auc.table = matrix(NA,length(com_type),length(prop.v))
    p.v = prop.v
    for (i in 1:length(com_type)) {
      for (j in 1:length(prop.v)) {
        tmp1 = dplyr::filter(q0_ana.bind,Type == com_type[i],prop.v == p.v[j]) # filter的名字不能重複，而且一定要加dplyr
        auc_1 = sum(tmp1$Value[seq_along(tmp1$Value[-1])] * diff(tau))
        auc_2 = sum(tmp1$Value[-1] * diff(tau))
        q0_ana_auc.table[i,j] = (auc_1 + auc_2)/2
      }
    }
    rownames(q0_ana_auc.table) = com_type
    colnames(q0_ana_auc.table) = p.v
    q0_ana_auc.table = reshape2::melt(q0_ana_auc.table)
    colnames(q0_ana_auc.table) = c("Type","prop.v","Value")
    q0_ana_auc.table = dplyr::arrange(q0_ana_auc.table,Type)
    q0_ana_auc.table = q0_ana_auc.table[,c("prop.v","Type","Value")]

    if(nboot>0){
      data_gamma = rowSums(data)
      ## 開平行 ======================================================
      # future::plan(multisession,workers = parallel::detectCores()-1)
      # =======================================================

      se = future.apply::future_lapply(1:nboot, function(boottime){

        p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, "abundance")
        ## 新增
        if (nrow(p_bt) > nrow(data)){
          unseen_p = p_bt[-(1:nrow(data)), ] %>% matrix(ncol = ncol(data))
          unseen_name = sapply(1:nrow(unseen_p), function(i) paste0("unseen_",
                                                                    i))
          rownames(p_bt) = c(rownames(data), unseen_name)
          f0_hat = nrow(p_bt) - nrow(data)
          distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), FDdistM, f0_hat, "abundance")
          rownames(distance_matrix_bt) = colnames(distance_matrix_bt) = rownames(p_bt)
          data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1,size = sum(data[, k]), prob = p_bt[, k]))
          rownames(data_bt) = rownames(p_bt)
        }else {
          distance_matrix_bt = FDdistM
          data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[, k]), prob = p_bt[,k]))
          rownames(data_bt) = rownames(data)
        }

        aivi.sample.data = lapply(1:nrow(m.v), function (k) mix.aivi.sample.FD(data_bt, m1.v[k], m2.v[k], FDdistM = distance_matrix_bt,
                                                                               tau = tau))
        output.list = lapply(1:length(tau), function(tt) {
          # S <- nrow(data_bt)
          qlength <- length(q)
          out = list()
          un1 = un2 = sh12 = c()

          for (zz in 1:length(m1.v)) {
            ai1.v = aivi.sample.data[[zz]]$ai1_mat[,tt]
            ai2.v = aivi.sample.data[[zz]]$ai2_mat[,tt]
            # ai.all.v = aivi.sample.data[[zz]]$ai_all_mat[,tt]
            vi.all.v = aivi.sample.data[[zz]]$vi_all_mat[,tt]

            out[[zz]] <- matrix(0, qlength, 1)
            m1 = m1.v[zz]
            m2 = m2.v[zz]

            for (j in 1:qlength) {
              if (q[j] == 0) {
                out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])
              } else if (q[j] == 1) {
                out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]))
              } else if (q[j] == 2) {
                out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
              } else {
                out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
              }
            }
            ## q0 ana =========================================
            datanew = cbind(ai1.v,ai2.v)
            rownames(datanew) = names(vi.all.v) = rownames(data_bt)

            datash = datanew[(data_bt[,1] > 0 & data_bt[,2]>0), , drop=F]
            dataun1 = datanew[(data_bt[,1] > 0 & data_bt[,2] == 0), , drop=F]
            dataun2 = datanew[(data_bt[,1] == 0 & data_bt[,2] > 0), , drop=F]

            ## 都用 sh_abun_FD因為在單群落下ai1, ai2都不一定非0
            un1[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun1)],xi1 = dataun1[,1],xi2 = dataun1[,2], n1 = n1,
                                     m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
            un2[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun2)],xi1 = dataun2[,1],xi2 = dataun2[,2], n1 = n1,
                                     m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
            sh12[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(datash)],xi1 = datash[,1],xi2 = datash[,2], n1 = n1,
                                      m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
          }

          output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
          rownames(output) = q
          colnames(output) = prop.v
          output = reshape2::melt(output)
          output$tau = tau[tt]
          colnames(output) = c("Order.q","prop.v","qFD_mix","tau")

          # q0_ana ===============
          q0_ana = data.frame(prop.v = prop.v ,
                              q0_un1 = un1, q0_un2 = un2,
                              q0_sh = sh12)
          q0_ana = reshape2::melt(q0_ana,id = c('prop.v'))
          colnames(q0_ana)[2:3] = c('Type','Value')
          q0_ana$tau = tau[tt]

          list(output = output,q0_ana_ = q0_ana)
        })

        ## auc =======================================================
        output.bind = output.list[[1]]$output
        for (j in 2:length(output.list)) {
          output.bind = rbind(output.bind,output.list[[j]]$output)
        }
        auc.table = matrix(NA,length(q),length(prop.v))
        p.v = prop.v
        for (i in 1:length(q)) {
          for (j in 1:length(prop.v)) {
            tmp1 = dplyr::filter(output.bind,Order.q == q[i] & prop.v == p.v[j]) # filter的名字不能重複，而且一定要加dplyr
            auc_1 = sum(tmp1$qFD_mix[seq_along(tmp1$qFD_mix[-1])] * diff(tau))
            auc_2 = sum(tmp1$qFD_mix[-1] * diff(tau))
            auc.table[i,j] = (auc_1 + auc_2)/2
          }
        }
        rownames(auc.table) = q
        colnames(auc.table) = p.v
        auc.table = reshape2::melt(auc.table)
        colnames(auc.table) = c("Order.q","prop.v","qmiNEXT_FD")
        res = auc.table$qmiNEXT_FD

        ## q0_ana auc ====================================
        q0_ana.bind = output.list[[1]]$q0_ana_
        for (j in 2:length(output.list)) {
          q0_ana.bind = rbind(q0_ana.bind,output.list[[j]]$q0_ana_)
        }
        com_type = unique(q0_ana.bind$Type)
        q0_ana_auc.table = matrix(NA,length(com_type),length(prop.v))

        for (i in 1:length(com_type)) {
          for (j in 1:length(prop.v)) {
            tmp1 = dplyr::filter(q0_ana.bind,Type == com_type[i],prop.v == p.v[j]) # filter的名字不能重複，而且一定要加dplyr
            auc_1 = sum(tmp1$Value[seq_along(tmp1$Value[-1])] * diff(tau))
            auc_2 = sum(tmp1$Value[-1] * diff(tau))
            q0_ana_auc.table[i,j] = (auc_1 + auc_2)/2
          }
        }
        rownames(q0_ana_auc.table) = com_type
        colnames(q0_ana_auc.table) = p.v
        q0_ana_auc.table = reshape2::melt(q0_ana_auc.table)
        colnames(q0_ana_auc.table) = c("Type","prop.v","Value")
        q0_ana_auc.table = dplyr::arrange(q0_ana_auc.table,Type)

        res_q0ana = q0_ana_auc.table$Value
        return(list(allana = res,q0ana = res_q0ana))

      },future.seed = NULL)

      ## 開平行 ======================================================
      # future::plan(NULL)
      # =======================================================
      se_all = do.call(cbind,lapply(se, function(k) k$allana))
      se_q0ana = do.call(cbind,lapply(se, function(k) k$q0ana))
      se_all = apply(se_all, 1, sd)
      se_q0ana = apply(se_q0ana, 1, sd)

      qtile <- qnorm(1 - (1 - conf)/2)
      auc.table$LCL = auc.table$qmiNEXT_FD - qtile*se_all
      auc.table$UCL = auc.table$qmiNEXT_FD + qtile*se_all

      q0_ana_auc.table$LCL = q0_ana_auc.table$Value - qtile*se_q0ana
      q0_ana_auc.table$UCL = q0_ana_auc.table$Value + qtile*se_q0ana

      q0_ana_auc.table$LCL[q0_ana_auc.table$LCL<0] = 0


    }else{
      auc.table$LCL = NA
      auc.table$UCL = NA
      q0_ana_auc.table$LCL = NA
      q0_ana_auc.table$UCL = NA

    }
    return(list(out = auc.table,q0_ana = q0_ana_auc.table,
                ori.prop = prop.v,m.v = m.v,line_type = line_type))

  }else{
    ## Extrapolation 改成全部m都跑，不然後面for loop的idx對不上
    if(nrow(m.ext.v) != 0){
      aivi.ext = lapply(1:nrow(m.v), function (k) mix.aivi.popu.FD(p_bind, m.v[k,1], m.v[k,2], FDdistM = FDdistM,
                                                                   tau = tau))
    }


    output.list = lapply(1:length(tau), function(tt) {

      # S <- nrow(data)
      qlength <- length(q)
      out = list()
      un1 = un2 = sh12 = c()

      for (zz in inc.idx) {
        ai1.v = aivi.sample.data[[zz]]$ai1_mat[,tt]
        ai2.v = aivi.sample.data[[zz]]$ai2_mat[,tt]
        # ai.all.v = aivi.sample.data[[zz]]$ai_all_mat[,tt]
        vi.all.v = aivi.sample.data[[zz]]$vi_all_mat[,tt]

        out[[zz]] <- matrix(0, qlength, 1)
        m1 = m1.v[zz]
        m2 = m2.v[zz]

        for (j in 1:qlength) {
          if (q[j] == 0) {
            out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])
          } else if (q[j] == 1) {
            out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]))
          } else if (q[j] == 2) {
            out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
          } else {
            out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
          }
        }
        ## q0 ana =========================================
        datanew = cbind(ai1.v,ai2.v)
        rownames(datanew) = names(vi.all.v) = rownames(data)
        datash = datanew[(data[,1] > 0 & data[,2]>0), , drop=F]
        dataun1 = datanew[(data[,1] > 0 & data[,2] == 0), , drop=F]
        dataun2 = datanew[(data[,1] == 0 & data[,2] > 0), , drop=F]

        ## 都用 sh_abun_FD因為在單群落下ai1, ai2都不一定非0
        un1[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun1)],xi1 = dataun1[,1],xi2 = dataun1[,2], n1 = n1,
                                 m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
        un2[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun2)],xi1 = dataun2[,1],xi2 = dataun2[,2], n1 = n1,
                                 m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
        sh12[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(datash)],xi1 = datash[,1],xi2 = datash[,2], n1 = n1,
                                  m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))

      }
      # 新增外插 ========================================

      for (zz in ext.idx) {
        ai1.v = aivi.sample.data[[zz]]$ai1_mat[,tt]
        ai2.v = aivi.sample.data[[zz]]$ai2_mat[,tt]
        # ai.all.v = aivi.sample.data[[zz]]$ai_all_mat[,tt]
        vi.all.v = aivi.sample.data[[zz]]$vi_all_mat[,tt]

        out[[zz]] <- matrix(0, qlength, 1)
        m1 = m1.v[zz]
        m2_all = m2.v[zz]
        m2_s = m2.v[zz] - n2

        # 下面的m2改n2
        ## new =============================================

        ai1.ext.v = aivi.ext[[zz]]$ai1_mat[,tt]
        ai2.ext.v = aivi.ext[[zz]]$ai2_mat[,tt]
        vi.all.ext.v = aivi.ext[[zz]]$vi_all_mat[,tt]

        h0.FD = h0_hat_cpp_FD(ai1.ext.v,ai2.ext.v,vi.all.ext.v,m1,m2_s,n1,n2)
        ## 思考h1_hat_FD的xi1, xi2要放啥
        h1.FD = h1_hat_FD(ai1.ext.v, ai2.ext.v, vi.all.ext.v, ai1.v, ai2.v,
                          m1, m2_all, n1, n2)

        ## ===================================================
        for (j in 1:qlength) {
          if (q[j] == 0) {
            out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]) + h0.FD
          } else if (q[j] == 1) {
            out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]) + h1.FD)
          } else if (q[j] == 2) {
            out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2_all, n1, n2)
          } else {
            out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
          }
        }
      }
      # ==============================================================
      output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
      rownames(output) = q
      colnames(output) = prop.v
      output = reshape2::melt(output)
      output$tau = tau[tt]
      colnames(output) = c("Order.q","prop.v","qFD_mix","tau")

      # q0_ana ===============
      if(nrow(m.inc.v) != 0){
        q0_ana = data.frame(prop.v = prop.inc ,
                            q0_un1 = un1, q0_un2 = un2,
                            q0_sh = sh12)
        q0_ana = reshape2::melt(q0_ana,id = c('prop.v'))
        colnames(q0_ana)[2:3] = c('Type','Value')
        q0_ana$tau = tau[tt]
      }else{
        q0_ana = NULL
      }

      list(output = output,q0_ana_ = q0_ana)
    })

    ## auc ================================================
    output.bind = output.list[[1]]$output
    for (j in 2:length(output.list)) {
      output.bind = rbind(output.bind,output.list[[j]]$output)
    }
    auc.table = matrix(NA,length(q),length(prop.v))
    p.v = prop.v
    for (i in 1:length(q)) {
      for (j in 1:length(prop.v)) {
        tmp1 = dplyr::filter(output.bind,Order.q == q[i] & prop.v == p.v[j]) # filter的名字不能重複，而且一定要加dplyr
        auc_1 = sum(tmp1$qFD_mix[seq_along(tmp1$qFD_mix[-1])] * diff(tau))
        auc_2 = sum(tmp1$qFD_mix[-1] * diff(tau))
        auc.table[i,j] = (auc_1 + auc_2)/2
      }
    }
    rownames(auc.table) = q
    colnames(auc.table) = p.v
    auc.table = reshape2::melt(auc.table)
    colnames(auc.table) = c("Order.q","prop.v","qmiNEXT_FD")

    ## q0_ana auc ===============
    if(nrow(m.inc.v) != 0){
      q0_ana.bind = output.list[[1]]$q0_ana_
      for (j in 2:length(output.list)) {
        q0_ana.bind = rbind(q0_ana.bind,output.list[[j]]$q0_ana_)
      }
      com_type = unique(q0_ana.bind$Type)
      q0_ana_auc.table = matrix(NA,length(com_type),length(prop.inc))
      p.v = prop.inc
      for (i in 1:length(com_type)) {
        for (j in 1:length(prop.inc)) {
          tmp1 = dplyr::filter(q0_ana.bind,Type == com_type[i],prop.v == p.v[j]) # filter的名字不能重複，而且一定要加dplyr
          auc_1 = sum(tmp1$Value[seq_along(tmp1$Value[-1])] * diff(tau))
          auc_2 = sum(tmp1$Value[-1] * diff(tau))
          q0_ana_auc.table[i,j] = (auc_1 + auc_2)/2
        }
      }
      rownames(q0_ana_auc.table) = com_type
      colnames(q0_ana_auc.table) = prop.inc
      q0_ana_auc.table = reshape2::melt(q0_ana_auc.table)
      colnames(q0_ana_auc.table) = c("Type","prop.v","Value")
      q0_ana_auc.table = dplyr::arrange(q0_ana_auc.table,Type)
      q0_ana_auc.table = q0_ana_auc.table[,c("prop.v","Type","Value")]
    }else{
      q0_ana_auc.table = NULL
    }

    if(nboot>0){
      data_gamma = rowSums(data)
      ## 開平行 ======================================================
      # future::plan(multisession,workers = parallel::detectCores()-1)
      # =======================================================

      se = future.apply::future_lapply(1:nboot, function(boottime){
        p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, "abundance")
        ## 新增
        if (nrow(p_bt) > nrow(data)){
          unseen_p = p_bt[-(1:nrow(data)), ] %>% matrix(ncol = ncol(data))
          unseen_name = sapply(1:nrow(unseen_p), function(i) paste0("unseen_",
                                                                    i))
          rownames(p_bt) = c(rownames(data), unseen_name)
          f0_hat = nrow(p_bt) - nrow(data)
          distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), FDdistM, f0_hat, "abundance")
          rownames(distance_matrix_bt) = colnames(distance_matrix_bt) = rownames(p_bt)
          data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1,size = sum(data[, k]), prob = p_bt[, k]))
          rownames(data_bt) = rownames(p_bt)
        }else {
          distance_matrix_bt = FDdistM
          data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[, k]), prob = p_bt[,k]))
          rownames(data_bt) = rownames(data)
        }

        aivi.sample.data = lapply(1:nrow(m.v), function (k) mix.aivi.sample.FD(data_bt, m1.v[k], m2.v[k], FDdistM = distance_matrix_bt,
                                                                               tau = tau))

        # new ===================================================================
        p_bind_bt = cbind(DetAbu(data_bt[,1], zero = T),DetAbu(data_bt[,2], zero = T))
        rownames(p_bind_bt) = rownames(data_bt)

        ## Extrapolation 改成全部m都跑，不然後面for迴圈idx對不上
        if(nrow(m.ext.v) != 0){
          aivi.ext = lapply(1:nrow(m.v), function (k) mix.aivi.popu.FD(p_bind_bt, m.v[k,1], m.v[k,2], FDdistM = distance_matrix_bt,
                                                                       tau = tau))
        }

        # ===================================================================================
        output.list = lapply(1:length(tau), function(tt) {
          # Li.v = apply(Lis, 1, function(x) x[1]*x[2])
          S <- nrow(data_bt)
          qlength <- length(q)
          out = list()
          un1 = un2 = sh12 = c()

          for (zz in inc.idx) {
            ai1.v = aivi.sample.data[[zz]]$ai1_mat[,tt]
            ai2.v = aivi.sample.data[[zz]]$ai2_mat[,tt]
            # ai.all.v = aivi.sample.data[[zz]]$ai_all_mat[,tt]
            vi.all.v = aivi.sample.data[[zz]]$vi_all_mat[,tt]

            out[[zz]] <- matrix(0, qlength, 1)
            m1 = m1.v[zz]
            m2 = m2.v[zz]

            for (j in 1:qlength) {
              if (q[j] == 0) {
                out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])
              } else if (q[j] == 1) {
                out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]))
              } else if (q[j] == 2) {
                out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
              } else {
                out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
              }
            }
            ## q0 ana =========================================
            datanew = cbind(ai1.v,ai2.v)
            rownames(datanew) = names(vi.all.v) = rownames(data_bt)
            datash = datanew[(data_bt[,1] > 0 & data_bt[,2]>0), , drop=F]
            dataun1 = datanew[(data_bt[,1] > 0 & data_bt[,2] == 0), , drop=F]
            dataun2 = datanew[(data_bt[,1] == 0 & data_bt[,2] > 0), , drop=F]

            ## 都用 sh_abun_FD因為在單群落下ai1, ai2都不一定非0
            un1[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun1)],xi1 = dataun1[,1],xi2 = dataun1[,2], n1 = n1,
                                     m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
            un2[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun2)],xi1 = dataun2[,1],xi2 = dataun2[,2], n1 = n1,
                                     m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
            sh12[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(datash)],xi1 = datash[,1],xi2 = datash[,2], n1 = n1,
                                      m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))

          }
          # 新增外插 ========================================

          for (zz in ext.idx) {
            ai1.v = aivi.sample.data[[zz]]$ai1_mat[,tt]
            ai2.v = aivi.sample.data[[zz]]$ai2_mat[,tt]
            # ai.all.v = aivi.sample.data[[zz]]$ai_all_mat[,tt]
            vi.all.v = aivi.sample.data[[zz]]$vi_all_mat[,tt]

            out[[zz]] <- matrix(0, qlength, 1)
            m1 = m1.v[zz]
            m2_all = m2.v[zz]
            m2_s = m2.v[zz] - n2

            # 下面的m2改n2
            ## new =============================================

            ai1.ext.v = aivi.ext[[zz]]$ai1_mat[,tt]
            ai2.ext.v = aivi.ext[[zz]]$ai2_mat[,tt]
            vi.all.ext.v = aivi.ext[[zz]]$vi_all_mat[,tt]

            h0.FD = h0_hat_cpp_FD(ai1.ext.v,ai2.ext.v,vi.all.ext.v,m1,m2_s,n1,n2)
            ## 思考h1_hat_FD的xi1, xi2要放啥
            h1.FD = h1_hat_FD(ai1.ext.v, ai2.ext.v, vi.all.ext.v, ai1.v, ai2.v,
                              m1, m2_all, n1, n2)

            ## ===================================================
            for (j in 1:qlength) {
              if (q[j] == 0) {
                out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]) + h0.FD
              } else if (q[j] == 1) {
                out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]) + h1.FD)
              } else if (q[j] == 2) {
                out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2_all, n1, n2)
              } else {
                out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
              }
            }
          }
          # ==============================================================
          output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
          rownames(output) = q
          colnames(output) = prop.v
          output = reshape2::melt(output)
          output$tau = tau[tt]
          colnames(output) = c("Order.q","prop.v","qFD_mix","tau")
          # q0_ana ===============
          if(nrow(m.inc.v) != 0){
            q0_ana = data.frame(prop.v = prop.inc ,
                                q0_un1 = un1, q0_un2 = un2,
                                q0_sh = sh12)
            q0_ana = reshape2::melt(q0_ana,id = c('prop.v'))
            colnames(q0_ana)[2:3] = c('Type','Value')
            q0_ana$tau = tau[tt]
          }else{
            q0_ana = NULL
          }
          list(output = output,q0_ana_ = q0_ana)
        })

        ## auc ================================================
        output.bind = output.list[[1]]$output
        for (j in 2:length(output.list)) {
          output.bind = rbind(output.bind,output.list[[j]]$output)
        }
        auc.table = matrix(NA,length(q),length(prop.v))
        p.v = prop.v
        for (i in 1:length(q)) {
          for (j in 1:length(prop.v)) {
            tmp1 = dplyr::filter(output.bind,Order.q == q[i] & prop.v == p.v[j]) # filter的名字不能重複，而且一定要加dplyr
            auc_1 = sum(tmp1$qFD_mix[seq_along(tmp1$qFD_mix[-1])] * diff(tau))
            auc_2 = sum(tmp1$qFD_mix[-1] * diff(tau))
            auc.table[i,j] = (auc_1 + auc_2)/2
          }
        }
        rownames(auc.table) = q
        colnames(auc.table) = p.v
        auc.table = reshape2::melt(auc.table)
        colnames(auc.table) = c("Order.q","prop.v","qmiNEXT_FD")
        res = auc.table$qmiNEXT_FD

        ## q0_ana auc ===============
        if(nrow(m.inc.v) != 0){
          q0_ana.bind = output.list[[1]]$q0_ana_
          for (j in 2:length(output.list)) {
            q0_ana.bind = rbind(q0_ana.bind,output.list[[j]]$q0_ana_)
          }
          com_type = unique(q0_ana.bind$Type)
          q0_ana_auc.table = matrix(NA,length(com_type),length(prop.inc))

          for (i in 1:length(com_type)) {
            for (j in 1:length(prop.inc)) {
              tmp1 = dplyr::filter(q0_ana.bind,Type == com_type[i],prop.v == prop.inc[j]) # filter的名字不能重複，而且一定要加dplyr
              auc_1 = sum(tmp1$Value[seq_along(tmp1$Value[-1])] * diff(tau))
              auc_2 = sum(tmp1$Value[-1] * diff(tau))
              q0_ana_auc.table[i,j] = (auc_1 + auc_2)/2
            }
          }
          rownames(q0_ana_auc.table) = com_type
          colnames(q0_ana_auc.table) = prop.inc
          q0_ana_auc.table = reshape2::melt(q0_ana_auc.table)
          colnames(q0_ana_auc.table) = c("Type","prop.v","Value")
          q0_ana_auc.table = dplyr::arrange(q0_ana_auc.table,Type)
          res_q0ana = q0_ana_auc.table$Value
        }else{
          res_q0ana = NULL
        }
        ## ===============================
        return(list(allana = res,q0ana = res_q0ana))
      },future.seed = NULL)

      ## 關平行 ======================================================
      # future::plan(NULL)
      # =======================================================
      se_all = do.call(cbind,lapply(se, function(k) k$allana))
      se_all = apply(se_all, 1, sd)

      qtile <- qnorm(1 - (1 - conf)/2)
      auc.table$LCL = auc.table$qmiNEXT_FD - qtile*se_all
      auc.table$UCL = auc.table$qmiNEXT_FD + qtile*se_all

      if(nrow(m.inc.v) != 0){
        se_q0ana = do.call(cbind,lapply(se, function(k) k$q0ana))
        se_q0ana = apply(se_q0ana, 1, sd)
        q0_ana_auc.table$LCL = q0_ana_auc.table$Value - qtile*se_q0ana
        q0_ana_auc.table$UCL = q0_ana_auc.table$Value + qtile*se_q0ana

        q0_ana_auc.table$LCL[q0_ana_auc.table$LCL<0] = 0
      }

    }else{
      auc.table$LCL = NA
      auc.table$UCL = NA
      if(nrow(m.inc.v) != 0){
        q0_ana_auc.table$LCL = NA
        q0_ana_auc.table$UCL = NA
      }
    }
    return(list(out = auc.table,q0_ana = q0_ana_auc.table,
                ori.prop = prop.v,m.v = m.v,line_type = line_type))
  }


}

## single tau weight by sample
mix.aivi.popu.FD.single = function(data, m1, m2, FDdistM, tau){
  FDdistM = as.matrix(FDdistM)
  data.list = list()
  for (i in 1:ncol(data)) {
    data.list[[i]] = data[,i]
    names(data.list[[i]]) = rownames(data)
  }
  data.sum = (m1/(m1+m2))*data[,1] + (m2/(m1+m2))*data[,2]
  names(data.sum) = rownames(data)

  ai_vi.popu.fn = function (data, dij_, tau_, filt_zero = FALSE) {
    # if (filt_zero) {
    #   dij <- dij[data > 0, data > 0]
    #   data <- data[data > 0]
    # }

    if (tau_ == 0) {
      dij_[dij_ > 0] <- 1
      a <- as.vector((1 - dij_/1) %*% data)
    }
    else {
      dij_[which(dij_ > tau_, arr.ind = T)] <- tau_
      a <- as.vector((1 - dij_/tau_) %*% data)
    }
    if (filt_zero) {
      data <- data[a != 0]
      a <- a[a != 0]
    }

    v <- data/a


    output = data.frame(ai = a, vi = v)

    return(output)
  }
  aivi_1 = ai_vi.popu.fn(data = data.list[[1]],dij_ = FDdistM, tau_ = tau, filt_zero = F)
  aivi_2 = ai_vi.popu.fn(data = data.list[[2]],dij_ = FDdistM, tau_ = tau, filt_zero = F)
  aivi_sum = ai_vi.popu.fn(data = data.sum,dij_ = FDdistM, tau_ = tau,filt_zero = F)
  aivi_sum$vi[is.nan(aivi_sum$vi)] = 0


  return(list(ai1_mat = aivi_1$ai, ai2_mat = aivi_2$ai,
              ai_all_mat = aivi_sum$ai, vi_all_mat = aivi_sum$vi))
}
mix.aivi.sample.FD.single = function(data, m1, m2, FDdistM, tau){
  FDdistM = as.matrix(FDdistM)
  data.list = list()
  for (i in 1:ncol(data)) {
    data.list[[i]] = data[,i]
    names(data.list[[i]]) = rownames(data)
  }
  data.sum = (m1/(m1+m2))*data[,1] + (m2/(m1+m2))*data[,2]
  names(data.sum) = rownames(data)
  ai_vi.sample.fn = function (data1, dij_, tau_, filt_zero = FALSE, integer = TRUE) {
    # if (filt_zero) {
    #   dij <- dij[data1 > 0, data1 > 0]
    #   data1 <- data1[data1 > 0]
    # }

    if (tau_ == 0) {
      dij_[dij_ > 0] <- 1
      a <- as.vector((1 - dij_/1) %*% data1)
    }else {
      dij_[which(dij_ > tau_, arr.ind = T)] <- tau_
      a <- as.vector((1 - dij_/tau_) %*% data1)
    }
    if (filt_zero) {
      data1 <- data1[a != 0]
      a <- a[a != 0]
    }

    if (integer){
      a[a < 1 & a > 0] <- 1
      a <- round(a)
    }

    v <- data1/a


    output = data.frame(ai = a, vi = v)

    return(output)
  }


  aivi_1 = ai_vi.sample.fn(data1 = data.list[[1]],dij_ = FDdistM, tau_ = tau, filt_zero = F)
  aivi_2 = ai_vi.sample.fn(data1 = data.list[[2]],dij_ = FDdistM, tau_ = tau, filt_zero = F)
  aivi_sum = ai_vi.sample.fn(data1 = data.sum,dij_ = FDdistM, tau_ = tau, filt_zero = F,integer = F)
  aivi_sum$vi[is.nan(aivi_sum$vi)] = 0
  # ai1_vec = list()
  # ai2_vec = list()
  # vi_all_vec = list()
  # for (i in 1:length(aivi_1)) {
  #   ai1_vec[[i]] = aivi_1[[i]]$ai
  #   ai2_vec[[i]] = aivi_2[[i]]$ai
  #   vi_all_vec[[i]] = aivi_sum[[i]]$vi
  # }
  # lapply(aivi_sum, function(k) sum(k[,1]*k[,2]))

  return(list(ai1_mat = aivi_1$ai, ai2_mat = aivi_2$ai,
              vi1_mat = aivi_1$vi, vi2_mat = aivi_2$vi,
              ai_all_mat = aivi_sum$ai, vi_all_mat = aivi_sum$vi))
}
RFD.singletau.obs <- function(data, FDdistM, tau, knots = 11, size = NULL, q = c(0,1,2), n1, n2) {

  if(n1>n2){
    a = round(seq(n2,n1, length.out = knots - round(knots/2) + 1))[-1]
    m2.v = c(round(seq(0, n2, length.out = round(knots/2))),a)
    m1.v = n1 - m2.v
    if(!is.null(size)){
      m1.v <- size
      m2.v <- n1- m1.v
    }
    m.v = cbind(m1.v, m2.v)
    inc.idx = which(m2.v<=n2)
    obs.idx = which(m2.v==n2)
    ext.idx = which(m2.v>n2)
    line_type = rep(NA,length(m1.v))
    line_type[inc.idx] = 'Rarefaction'
    line_type[obs.idx] = 'Observe'
    line_type[ext.idx] = 'Extrapolation'
    prop.v = round(((m1.v/n1) * 100),3)

  }else{
    m2.v = round(seq(0, n1, length.out = knots))
    m1.v = n1 - m2.v
    if(!is.null(size)){
      m1.v <- size
      m2.v <- n1- m1.v
    }
    m.v = cbind(m1.v, m2.v)
    prop.v = round(((m1.v/n1) * 100),3)
    line_type = rep("Rarefaction",knots)
  }



  # 不同的m1, m2有不同的ai_all, vi_all。 ai1, ai2則相同
  aivi.popu.data = lapply(1:nrow(m.v), function (k) mix.aivi.popu.FD.single(data, m1.v[k], m2.v[k], FDdistM = FDdistM,
                                                                            tau = tau))

  # Li.v = apply(Lis, 1, function(x) x[1]*x[2])
  S <- nrow(data)
  qlength <- length(q)
  out = list()
  ghat_pt2 = list()

  for (zz in 1:length(m1.v)) {
    ai1.v = aivi.popu.data[[zz]]$ai1_mat
    ai2.v = aivi.popu.data[[zz]]$ai2_mat
    # ai.all.v = aivi.popu.data[[zz]]$ai_all_mat
    vi.all.v = aivi.popu.data[[zz]]$vi_all_mat
    # Vbar = sum(ai.all.v * vi.all.v)/(n1 + n2)
    Vbar = 1
    out[[zz]] <- matrix(0, qlength, 1)
    m1 = m1.v[zz]
    m2 = m2.v[zz]

    ghat_pt2[[zz]] <- FD.h.theo_fn(m1, m2, ai1.v, ai2.v ,vi.all.v, S)




    for (j in 1:qlength) {
      for (k1 in 0:m1) {
        for (k2 in 0:m2) {
          if (q[j] == 0) {
            out[[zz]][j,1] <- out[[zz]][j,1] + ghat_pt2[[zz]][k1+1,k2+1]
          }else if (q[j] == 1) {
            if(k1 == 0 & k2 == 0){
              aa = 0
            }else{
              aa = ((k1+k2)/((m1+m2)*Vbar)) * log((k1+k2)/((m1+m2)*Vbar))
            }
            out[[zz]][j,1] <- out[[zz]][j,1] - aa  * ghat_pt2[[zz]][k1+1,k2+1]
          }else if (q[j] == 2) {
            out[[zz]][j,1] <- out[[zz]][j,1] + ((k1+k2)/((m1+m2)*Vbar))^2 * ghat_pt2[[zz]][k1+1,k2+1]
          }else {
            out[[zz]][j,1] <- out[[zz]][j,1] + ((k1+k2)/((m1+m2)*Vbar))^q[j] * ghat_pt2[[zz]][k1+1,k2+1]
          }
        }

      }

      if (q[j] == 0) {
        out[[zz]][j,1] <- out[[zz]][j,1]
      } else if (q[j] == 1) {
        out[[zz]][j,1] <- exp(out[[zz]][j,1])
      } else if (q[j] == 2) {
        out[[zz]][j,1] <- 1 / out[[zz]][j,1]
      } else {
        out[[zz]][j,1] <- out[[zz]][j,1]^(1/(1-q[j]))
      }
    }
  }

  output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
  rownames(output) = q
  colnames(output) = prop.v
  output = reshape2::melt(output)

  colnames(output) = c("Order.q","prop.v","qmiNEXT_FD_singletau")
  return(list(out = output,ori.prop = prop.v,m.v = m.v,line_type = line_type))
}
RFD.singletau.est <- function(data, FDdistM, tau = NULL, knots = 11, size = NULL, q = c(0,1,2), conf = 0.95, nboot = 0) {

  if (length(tau) > 1){
    stop("Threshold must be one number between 0 and 1. Check the length of tau.",
         call. = FALSE)
  }
  if (is.null(tau)) {
    tmp = rowSums(data)
    tmp <- matrix(tmp/sum(tmp), ncol = 1)
    tau <- sum((tmp %*% t(tmp)) * FDdistM)
  }else if (tau < 0 | tau > 1){
    stop("Threshold must be one number between 0 and 1. Use NULL to set it to dmean.",
         call. = FALSE)
  }else{
    tau <- tau
  }


  n1 = sum(data[,1])
  n2 = sum(data[,2])

  if(n1>n2){
    a = round(seq(n2,n1, length.out = knots - round(knots/2) + 1))[-1]
    m2.v = c(round(seq(0, n2, length.out = round(knots/2))),a)
    m1.v = n1 - m2.v
    if(!is.null(size)){
      m1.v <- size
      m2.v <- n1- m1.v
    }
    m.v = cbind(m1.v, m2.v)
    prop.v = round(((m1.v/n1) * 100),3)
    inc.idx = which(m2.v<=n2)
    obs.idx = which(m2.v==n2)
    ext.idx = which(m2.v>n2)
    line_type = rep(NA,length(m1.v))
    line_type[inc.idx] = 'Rarefaction'
    line_type[obs.idx] = 'Observe'
    line_type[ext.idx] = 'Extrapolation'

    m.inc.v = cbind(m1.v[inc.idx],m2.v[inc.idx])
    prop.inc = prop.v[inc.idx]
    m.ext.v = cbind(m1.v[-inc.idx],m2.v[-inc.idx])
    prop.ext = prop.v[-inc.idx]
    p1 = DetAbu(data[,1], zero = T)
    p2 = DetAbu(data[,2], zero = T)
    p_bind = cbind(p1,p2)
    rownames(p_bind) = rownames(data)
  }else{
    m2.v = round(seq(0, n1, length.out = knots))
    m1.v = n1 - m2.v
    if(!is.null(size)){
      m1.v <- size
      m2.v <- n1- m1.v
    }
    m.v = cbind(m1.v, m2.v)
    prop.v = round(((m1.v/n1) * 100),3)
    line_type = rep("Rarefaction",length(m1.v))
  }

  # 不同的m1, m2有不同的ai_all, vi_all。 ai1, ai2則相同
  aivi.sample.data = lapply(1:nrow(m.v), function (k) mix.aivi.sample.FD.single(data, m1.v[k], m2.v[k],FDdistM = FDdistM, tau = tau))



  if(n1<=n2){

    # Li.v = apply(Lis, 1, function(x) x[1]*x[2])
    S <- nrow(data)
    qlength <- length(q)
    out = list()
    un1 = un2 = sh12 = c()

    for (zz in 1:length(m1.v)) {
      ai1.v = aivi.sample.data[[zz]]$ai1_mat
      ai2.v = aivi.sample.data[[zz]]$ai2_mat
      # ai.all.v = aivi.sample.data[[zz]]$ai_all_mat
      vi.all.v = aivi.sample.data[[zz]]$vi_all_mat

      out[[zz]] <- matrix(0, qlength, 1)
      m1 = m1.v[zz]
      m2 = m2.v[zz]

      for (j in 1:qlength) {
        if (q[j] == 0) {
          out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])
        } else if (q[j] == 1) {
          out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]))
        } else if (q[j] == 2) {
          out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
        } else {
          out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
        }
      }

      ## q0 ana =========================================
      datanew = cbind(ai1.v,ai2.v)
      rownames(datanew) = names(vi.all.v) = rownames(data)
      datash = datanew[(data[,1] > 0 & data[,2]>0), , drop=F]
      dataun1 = datanew[(data[,1] > 0 & data[,2] == 0), , drop=F]
      dataun2 = datanew[(data[,1] == 0 & data[,2] > 0), , drop=F]

      ## 都用 sh_abun_FD因為在單群落下ai1, ai2都不一定非0
      un1[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun1)],xi1 = dataun1[,1],xi2 = dataun1[,2], n1 = n1,
                               m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
      un2[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun2)],xi1 = dataun2[,1],xi2 = dataun2[,2], n1 = n1,
                               m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
      sh12[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(datash)],xi1 = datash[,1],xi2 = datash[,2], n1 = n1,
                                m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))

    }

    output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
    rownames(output) = q
    colnames(output) = prop.v
    output = reshape2::melt(output)
    colnames(output) = c("Order.q","prop.v","qmiNEXT_FD_singletau")


    q0_ana = data.frame(prop.v = prop.v ,
                        q0_un1 = un1, q0_un2 = un2,
                        q0_sh = sh12)
    q0_ana = reshape2::melt(q0_ana,id = c('prop.v'))
    colnames(q0_ana)[2:3] = c('Type','Value')

    if(nboot>0){
      data_gamma = rowSums(data)
      se = future.apply::future_lapply(1:nboot, function(boottime){
        p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, "abundance")
        ## 新增
        if (nrow(p_bt) > nrow(data)){
          unseen_p = p_bt[-(1:nrow(data)), ] %>% matrix(ncol = ncol(data))
          unseen_name = sapply(1:nrow(unseen_p), function(i) paste0("unseen_",
                                                                  i))
          rownames(p_bt) = c(rownames(data), unseen_name)
          f0_hat = nrow(p_bt) - nrow(data)
          distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), FDdistM, f0_hat, "abundance")
          rownames(distance_matrix_bt) = colnames(distance_matrix_bt) = rownames(p_bt)
          data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1,size = sum(data[, k]), prob = p_bt[, k]))
          rownames(data_bt) = rownames(p_bt)
        }else {
          distance_matrix_bt = FDdistM
          data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[, k]), prob = p_bt[,k]))
          rownames(data_bt) = rownames(data)
        }

        aivi.sample.data = lapply(1:nrow(m.v), function (k) mix.aivi.sample.FD.single(data_bt, m1.v[k], m2.v[k], FDdistM = distance_matrix_bt,
                                                                                      tau = tau))

        S <- nrow(data_bt)
        qlength <- length(q)
        out = list()
        un1 = un2 = sh12 = c()

        for (zz in 1:length(m1.v)) {
          ai1.v = aivi.sample.data[[zz]]$ai1_mat
          ai2.v = aivi.sample.data[[zz]]$ai2_mat
          # ai.all.v = aivi.sample.data[[zz]]$ai_all_mat
          vi.all.v = aivi.sample.data[[zz]]$vi_all_mat

          out[[zz]] <- matrix(0, qlength, 1)
          m1 = m1.v[zz]
          m2 = m2.v[zz]

          for (j in 1:qlength) {
            if (q[j] == 0) {
              out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])
            } else if (q[j] == 1) {
              out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]))
            } else if (q[j] == 2) {
              out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
            } else {
              out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
            }
          }

          ## q0 ana =========================================
          datanew = cbind(ai1.v,ai2.v)
          rownames(datanew) = names(vi.all.v) = rownames(data_bt)

          datash = datanew[(data_bt[,1] > 0 & data_bt[,2]>0), , drop=F]
          dataun1 = datanew[(data_bt[,1] > 0 & data_bt[,2] == 0), , drop=F]
          dataun2 = datanew[(data_bt[,1] == 0 & data_bt[,2] > 0), , drop=F]

          ## 都用 sh_abun_FD因為在單群落下ai1, ai2都不一定非0
          un1[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun1)],xi1 = dataun1[,1],xi2 = dataun1[,2], n1 = n1,
                                   m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
          un2[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun2)],xi1 = dataun2[,1],xi2 = dataun2[,2], n1 = n1,
                                   m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
          sh12[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(datash)],xi1 = datash[,1],xi2 = datash[,2], n1 = n1,
                                    m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))


        }

        output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
        rownames(output) = q
        colnames(output) = prop.v
        output = reshape2::melt(output)

        colnames(output) = c("Order.q","prop.v","qmiNEXT_FD_singletau")

        res = output$qmiNEXT_FD_singletau

        ## q0_ana ==========================================
        q0_ana = data.frame(prop.v = prop.v ,
                            q0_un1 = un1, q0_un2 = un2,
                            q0_sh = sh12)
        q0_ana = reshape2::melt(q0_ana,id = c('prop.v'))
        colnames(q0_ana)[2:3] = c('Type','Value')
        res_q0ana = q0_ana$Value
        return(list(allana = res,q0ana = res_q0ana))

      },future.seed = NULL)
      se_all = do.call(cbind,lapply(se, function(k) k$allana))
      se_q0ana = do.call(cbind,lapply(se, function(k) k$q0ana))
      se_all = apply(se_all, 1, sd)
      se_q0ana = apply(se_q0ana, 1, sd)

      qtile <- qnorm(1 - (1 - conf)/2)
      output$LCL = output$qmiNEXT_FD_singletau - qtile*se_all
      output$UCL = output$qmiNEXT_FD_singletau + qtile*se_all

      q0_ana$LCL = q0_ana$Value - qtile*se_q0ana
      q0_ana$UCL = q0_ana$Value + qtile*se_q0ana

      q0_ana$LCL[q0_ana$LCL<0] = 0

    }else{
      output$LCL = NA
      output$UCL = NA
      q0_ana$LCL = NA
      q0_ana$UCL = NA
    }
    return(list(out = output,q0_ana = q0_ana,ori.prop = prop.v,m.v = m.v,line_type = line_type))

  }else{

    if(nrow(m.ext.v) != 0){
      aivi.ext = lapply(1:nrow(m.v), function (k) mix.aivi.popu.FD.single(p_bind, m.v[k,1], m.v[k,2], FDdistM = FDdistM,
                                                                          tau = tau))
    }

    S <- nrow(data)
    qlength <- length(q)
    out = list()
    un1 = un2 = sh12 = c()

    for (zz in inc.idx) {
      ai1.v = aivi.sample.data[[zz]]$ai1_mat
      ai2.v = aivi.sample.data[[zz]]$ai2_mat
      # ai.all.v = aivi.sample.data[[zz]]$ai_all_mat
      vi.all.v = aivi.sample.data[[zz]]$vi_all_mat

      out[[zz]] <- matrix(0, qlength, 1)
      m1 = m1.v[zz]
      m2 = m2.v[zz]
      for (j in 1:qlength) {
        if (q[j] == 0) {
          out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])
        } else if (q[j] == 1) {
          out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]))
        } else if (q[j] == 2) {
          out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
        } else {
          out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
        }
      }

      ## q0 ana =========================================
      datanew = cbind(ai1.v,ai2.v)
      rownames(datanew) = names(vi.all.v) = rownames(data)
      datash = datanew[(data[,1] > 0 & data[,2]>0), , drop=F]
      dataun1 = datanew[(data[,1] > 0 & data[,2] == 0), , drop=F]
      dataun2 = datanew[(data[,1] == 0 & data[,2] > 0), , drop=F]

      ## 都用 sh_abun_FD因為在單群落下ai1, ai2都不一定非0
      un1[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun1)],xi1 = dataun1[,1],xi2 = dataun1[,2], n1 = n1,
                               m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
      un2[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun2)],xi1 = dataun2[,1],xi2 = dataun2[,2], n1 = n1,
                               m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
      sh12[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(datash)],xi1 = datash[,1],xi2 = datash[,2], n1 = n1,
                                m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))


    }
    # 新增外插 ========================================
    for (zz in ext.idx) {
      ai1.v = aivi.sample.data[[zz]]$ai1_mat
      ai2.v = aivi.sample.data[[zz]]$ai2_mat
      # ai.all.v = aivi.sample.data[[zz]]$ai_all_mat
      vi.all.v = aivi.sample.data[[zz]]$vi_all_mat

      out[[zz]] <- matrix(0, qlength, 1)
      m1 = m1.v[zz]
      m2_all = m2.v[zz]
      m2_s = m2.v[zz] - n2

      # 下面的m2改n2
      ## new =============================================

      ai1.ext.v = aivi.ext[[zz]]$ai1_mat
      ai2.ext.v = aivi.ext[[zz]]$ai2_mat
      vi.all.ext.v = aivi.ext[[zz]]$vi_all_mat

      h0.FD = h0_hat_cpp_FD(ai1.ext.v,ai2.ext.v,vi.all.ext.v,m1,m2_s,n1,n2)
      ## 思考h1_hat_FD的xi1, xi2要放啥
      h1.FD = h1_hat_FD(ai1.ext.v, ai2.ext.v, vi.all.ext.v, ai1.v, ai2.v,
                        m1, m2_all, n1, n2)


      ## ===================================================
      for (j in 1:qlength) {
        if (q[j] == 0) {
          out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]) + h0.FD
        } else if (q[j] == 1) {
          out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]) + h1.FD)
        } else if (q[j] == 2) {
          # 記得放m2_all
          out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2_all, n1, n2)
        } else {
          out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
        }
      }
    }
    # ==============================================================
    output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
    rownames(output) = q
    colnames(output) = prop.v
    output = reshape2::melt(output)
    colnames(output) = c("Order.q","prop.v","qmiNEXT_FD_singletau")

    if(nrow(m.inc.v) != 0){
      q0_ana = data.frame(prop.v = prop.inc ,
                          q0_un1 = un1, q0_un2 = un2,
                          q0_sh = sh12)
      q0_ana = reshape2::melt(q0_ana,id = c('prop.v'))
      colnames(q0_ana)[2:3] = c('Type','Value')
    }else{
      q0_ana = NULL
    }


    if(nboot>0){
      data_gamma = rowSums(data)
      se = future.apply::future_lapply(1:nboot, function(boottime){
        p_bt = bootstrap_population_multiple_assemblage(data, data_gamma, "abundance")
        ## 新增
        if (nrow(p_bt) > nrow(data)){
          unseen_p = p_bt[-(1:nrow(data)), ] %>% matrix(ncol = ncol(data))
          unseen_name = sapply(1:nrow(unseen_p), function(i) paste0("unseen_",
                                                                    i))
          rownames(p_bt) = c(rownames(data), unseen_name)
          f0_hat = nrow(p_bt) - nrow(data)
          distance_matrix_bt = Bootstrap_distance_matrix(rowSums(data), FDdistM, f0_hat, "abundance")
          rownames(distance_matrix_bt) = colnames(distance_matrix_bt) = rownames(p_bt)
          data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1,size = sum(data[, k]), prob = p_bt[, k]))
          rownames(data_bt) = rownames(p_bt)
        }else {
          distance_matrix_bt = FDdistM
          data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[, k]), prob = p_bt[,k]))
          rownames(data_bt) = rownames(data)
        }


        aivi.sample.data = lapply(1:nrow(m.v), function (k) mix.aivi.sample.FD.single(data_bt, m1.v[k], m2.v[k], FDdistM = distance_matrix_bt,
                                                                                      tau = tau))

        # new ===================================================================
        p_bind_bt = cbind(DetAbu(data_bt[,1], zero = T),DetAbu(data_bt[,2], zero = T))
        rownames(p_bind_bt) = rownames(data_bt)

        if(nrow(m.ext.v) != 0){
          aivi.ext = lapply(1:nrow(m.v), function (k) mix.aivi.popu.FD.single(p_bind_bt, m.v[k,1], m.v[k,2], FDdistM = distance_matrix_bt,
                                                                              tau = tau))
        }

        # ===================================================================================

        # Li.v = apply(Lis, 1, function(x) x[1]*x[2])
        S <- nrow(data_bt)
        qlength <- length(q)
        out = list()
        un1 = un2 = sh12 = c()

        for (zz in inc.idx) {
          ai1.v = aivi.sample.data[[zz]]$ai1_mat
          ai2.v = aivi.sample.data[[zz]]$ai2_mat
          # ai.all.v = aivi.sample.data[[zz]]$ai_all_mat
          vi.all.v = aivi.sample.data[[zz]]$vi_all_mat

          out[[zz]] <- matrix(0, qlength, 1)
          m1 = m1.v[zz]
          m2 = m2.v[zz]

          for (j in 1:qlength) {
            if (q[j] == 0) {
              out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])
            } else if (q[j] == 1) {
              out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]))
            } else if (q[j] == 2) {
              out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
            } else {
              out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
            }
          }

          ## q0 ana =========================================
          datanew = cbind(ai1.v,ai2.v)
          rownames(datanew) = names(vi.all.v) = rownames(data_bt)
          datash = datanew[(data_bt[,1] > 0 & data_bt[,2]>0), , drop=F]
          dataun1 = datanew[(data_bt[,1] > 0 & data_bt[,2] == 0), , drop=F]
          dataun2 = datanew[(data_bt[,1] == 0 & data_bt[,2] > 0), , drop=F]

          ## 都用 sh_abun_FD因為在單群落下ai1, ai2都不一定非0
          un1[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun1)],xi1 = dataun1[,1],xi2 = dataun1[,2], n1 = n1,
                                   m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
          un2[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(dataun2)],xi1 = dataun2[,1],xi2 = dataun2[,2], n1 = n1,
                                   m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
          sh12[zz] = sum(sh_abun_FD(vi = vi.all.v[rownames(datash)],xi1 = datash[,1],xi2 = datash[,2], n1 = n1,
                                    m1 =  m1.v[zz], n2 = n2, m2 = m2.v[zz]))
        }
        # 新增外插 ========================================
        for (zz in ext.idx) {
          ai1.v = aivi.sample.data[[zz]]$ai1_mat
          ai2.v = aivi.sample.data[[zz]]$ai2_mat
          # ai.all.v = aivi.sample.data[[zz]]$ai_all_mat
          vi.all.v = aivi.sample.data[[zz]]$vi_all_mat

          out[[zz]] <- matrix(0, qlength, 1)
          m1 = m1.v[zz]
          m2_all = m2.v[zz]
          m2_s = m2.v[zz] - n2

          # 下面的m2改n2
          ## new =============================================
          ai1.ext.v = aivi.ext[[zz]]$ai1_mat
          ai2.ext.v = aivi.ext[[zz]]$ai2_mat
          # ai.all.v = aivi.ext[[zz]]$ai_all_mat
          vi.all.ext.v = aivi.ext[[zz]]$vi_all_mat

          h0.FD = h0_hat_cpp_FD(ai1.ext.v,ai2.ext.v,vi.all.ext.v,m1,m2_s,n1,n2)
          ## 思考h1_hat_FD的xi1, xi2要放啥
          h1.FD = h1_hat_FD(ai1.ext.v, ai2.ext.v, vi.all.ext.v, ai1.v, ai2.v,
                            m1, m2_all, n1, n2)
          ## ===================================================
          for (j in 1:qlength) {
            if (q[j] == 0) {
              out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]) + h0.FD
            } else if (q[j] == 1) {
              out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j]) + h1.FD)
            } else if (q[j] == 2) {
              out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2_all, n1, n2)
            } else {
              out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
            }
          }
        }
        # ==============================================================
        output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
        rownames(output) = q
        colnames(output) = prop.v
        output = reshape2::melt(output)

        colnames(output) = c("Order.q","prop.v","qmiNEXT_FD_singletau")


        res = output$qmiNEXT_FD_singletau

        ## q0_ana ==========================================
        if(nrow(m.inc.v) != 0){
          q0_ana = data.frame(prop.v = prop.inc ,
                              q0_un1 = un1, q0_un2 = un2,
                              q0_sh = sh12)
          q0_ana = reshape2::melt(q0_ana,id = c('prop.v'))
          colnames(q0_ana)[2:3] = c('Type','Value')
          res_q0ana = q0_ana$Value
        }else{
          res_q0ana = NULL
        }
        return(list(allana = res,q0ana = res_q0ana))
      },future.seed = NULL)
      se_all = do.call(cbind,lapply(se, function(k) k$allana))
      se_all = apply(se_all, 1, sd)

      qtile <- qnorm(1 - (1 - conf)/2)
      output$LCL = output$qmiNEXT_FD_singletau - qtile*se_all
      output$UCL = output$qmiNEXT_FD_singletau + qtile*se_all

      if(nrow(m.inc.v) != 0){
        se_q0ana = do.call(cbind,lapply(se, function(k) k$q0ana))
        se_q0ana = apply(se_q0ana, 1, sd)
        q0_ana$LCL = q0_ana$Value - qtile*se_q0ana
        q0_ana$UCL = q0_ana$Value + qtile*se_q0ana

        q0_ana$LCL[q0_ana$LCL<0] = 0
      }

    }else{
      output$LCL = NA
      output$UCL = NA
      if(nrow(m.inc.v) != 0){
        q0_ana$LCL = NA
        q0_ana$UCL = NA
      }
    }
    return(list(out = output,q0_ana = q0_ana,ori.prop = prop.v,m.v = m.v,line_type = line_type))
  }


}
