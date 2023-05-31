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

FDq2<-function(x1,x2,vi.all,m1,m2,n1,n2){
  
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
RFD.est <- function(data, FDdistM, FDtau = NULL, knots = 11, size = NULL, q = c(0,1,2), nboot = 0) {
  
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
      S <- nrow(data)
      qlength <- length(q)
      out = list()
      
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
            out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)
          } else if (q[j] == 1) {
            out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S))
          } else if (q[j] == 2) {
            out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
          } else {
            out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)^(1/(1-q[j]))
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
    
    
    if(nboot>0){
      data_gamma = rowSums(data)
      ## 開平行 ======================================================
      # future::plan(multisession,workers = parallel::detectCores()-1)
      # =======================================================
      
      se = future.apply::future_lapply(1:nboot, function(boottime){
        p_bt = iNEXT.beta3D:::bootstrap_population_multiple_assemblage(data, 
                                                                       data_gamma, "abundance")
        f0_hat = nrow(p_bt) - nrow(data)
        distance_matrix_bt = iNEXT.beta3D:::Bootstrap_distance_matrix(rowSums(data), 
                                                                      FDdistM, f0_hat, "abundance")
        data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, 
                                                             size = sum(data[, k]), prob = p_bt[, k]))
        
        aivi.sample.data = lapply(1:nrow(m.v), function (k) mix.aivi.sample.FD(data_bt, m1.v[k], m2.v[k], FDdistM = distance_matrix_bt,
                                                                               tau = tau))
        output.list = lapply(1:length(tau), function(tt) {
          # Li.v = apply(Lis, 1, function(x) x[1]*x[2])
          S <- nrow(data_bt)
          qlength <- length(q)
          out = list()
          
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
                out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)
              } else if (q[j] == 1) {
                out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S))
              } else if (q[j] == 2) {
                out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
              } else {
                out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)^(1/(1-q[j]))
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
        res = auc.table$qmiNEXT_FD
        res
      },future.seed = NULL) %>% cbind()
      
      ## 開平行 ======================================================
      # future::plan(NULL)
      # =======================================================
      se = abind::abind(se,along = 2)
      se = apply(se, 1, sd)
      
      auc.table$LCL = auc.table$qmiNEXT_FD - 1.96*se
      auc.table$UCL = auc.table$qmiNEXT_FD + 1.96*se
      
    }else{
      auc.table$LCL = NA
      auc.table$UCL = NA
    }
    return(list(out = auc.table,ori.prop = prop.v,m.v = m.v,line_type = line_type))
    
  }else{
    ## Extrapolation 改成全部m都跑，不然後面for loop的idx對不上
    if(nrow(m.ext.v) != 0){
      aivi.ext = lapply(1:nrow(m.v), function (k) mix.aivi.popu.FD(p_bind, m.v[k,1], m.v[k,2], FDdistM = FDdistM,
                                                                   tau = tau))
    }
    
    
    output.list = lapply(1:length(tau), function(tt) {
      # Li.v = apply(Lis, 1, function(x) x[1]*x[2])
      S <- nrow(data)
      qlength <- length(q)
      out = list()
      
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
            out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)
          } else if (q[j] == 1) {
            out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S))
          } else if (q[j] == 2) {
            out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
          } else {
            out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)^(1/(1-q[j]))
          }
        }
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
        h1.FD = h1_hat_FD(ai1.ext.v, ai2.ext.v, vi.all.ext.v, p_bind[,1], p_bind[,2],
                          m1, m2_all, n1, n2)
        
        ## ===================================================
        for (j in 1:qlength) {
          if (q[j] == 0) {
            out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S) + h0.FD
          } else if (q[j] == 1) {
            out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S) + h1.FD)
          } else if (q[j] == 2) {
            out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2_all, n1, n2)
          } else {
            out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)^(1/(1-q[j]))
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
        tmp1 = dplyr::filter(output.bind,Order.q == q[i] & prop.v == p.v[j])
        auc_1 = sum(tmp1$qFD_mix[seq_along(tmp1$qFD_mix[-1])] * diff(tau))
        auc_2 = sum(tmp1$qFD_mix[-1] * diff(tau))
        auc.table[i,j] = (auc_1 + auc_2)/2
      }
    }
    rownames(auc.table) = q
    colnames(auc.table) = p.v
    auc.table = reshape2::melt(auc.table)
    colnames(auc.table) = c("Order.q","prop.v","qmiNEXT_FD")
    
    
    if(nboot>0){
      data_gamma = rowSums(data)
      ## 開平行 ======================================================
      # future::plan(multisession,workers = parallel::detectCores()-1)
      # =======================================================
      
      se = future.apply::future_lapply(1:nboot, function(boottime){
        p_bt = iNEXT.beta3D:::bootstrap_population_multiple_assemblage(data, 
                                                                       data_gamma, "abundance")
        f0_hat = nrow(p_bt) - nrow(data)
        distance_matrix_bt = iNEXT.beta3D:::Bootstrap_distance_matrix(rowSums(data), 
                                                                      FDdistM, f0_hat, "abundance")
        
        data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, 
                                                             size = sum(data[, k]), prob = p_bt[, k]))
        
        
        
        
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
                out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)
              } else if (q[j] == 1) {
                out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S))
              } else if (q[j] == 2) {
                out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
              } else {
                out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)^(1/(1-q[j]))
              }
            }
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
            h1.FD = h1_hat_FD(ai1.ext.v, ai2.ext.v, vi.all.ext.v, p_bind_bt[,1], p_bind_bt[,2],
                              m1, m2_all, n1, n2)
            
            ## ===================================================
            for (j in 1:qlength) {
              if (q[j] == 0) {
                out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S) + h0.FD
              } else if (q[j] == 1) {
                out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S) + h1.FD)
              } else if (q[j] == 2) {
                out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2_all, n1, n2)
              } else {
                out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)^(1/(1-q[j]))
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
            tmp1 = dplyr::filter(output.bind,Order.q == q[i] & prop.v == p.v[j])
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
        res
      },future.seed = NULL) %>% cbind()
      
      ## 關平行 ======================================================
      # future::plan(NULL)
      # =======================================================
      se = abind::abind(se,along = 2)
      se = apply(se, 1, sd)
      
      auc.table$LCL = auc.table$qmiNEXT_FD - 1.96*se
      auc.table$UCL = auc.table$qmiNEXT_FD + 1.96*se
      
    }else{
      auc.table$LCL = NA
      auc.table$UCL = NA
    }
    return(list(out = auc.table,ori.prop = prop.v,m.v = m.v,line_type = line_type))
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
RFD.singletau.est <- function(data, FDdistM, tau, knots = 11, size = NULL, q = c(0,1,2), nboot = 0) {
  
  
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
          out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)
        } else if (q[j] == 1) {
          out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S))
        } else if (q[j] == 2) {
          out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
        } else {
          out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)^(1/(1-q[j]))
        }
      }
    }
    
    output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
    rownames(output) = q
    colnames(output) = prop.v
    output = reshape2::melt(output)
    colnames(output) = c("Order.q","prop.v","qmiNEXT_FD_singletau")
    
    if(nboot>0){
      data_gamma = rowSums(data)
      se = future.apply::future_lapply(1:nboot, function(boottime){
        p_bt = iNEXT.beta3D:::bootstrap_population_multiple_assemblage(data, 
                                                                       data_gamma, "abundance")
        f0_hat = nrow(p_bt) - nrow(data)
        distance_matrix_bt = iNEXT.beta3D:::Bootstrap_distance_matrix(rowSums(data), 
                                                                      FDdistM, f0_hat, "abundance")
        data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1,size = sum(data[, k]), prob = p_bt[, k]))
        
        aivi.sample.data = lapply(1:nrow(m.v), function (k) mix.aivi.sample.FD.single(data_bt, m1.v[k], m2.v[k], FDdistM = distance_matrix_bt,
                                                                                      tau = tau))
        
        # Li.v = apply(Lis, 1, function(x) x[1]*x[2])
        S <- nrow(data_bt)
        qlength <- length(q)
        out = list()
        
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
              out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)
            } else if (q[j] == 1) {
              out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S))
            } else if (q[j] == 2) {
              out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
            } else {
              out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)^(1/(1-q[j]))
            }
          }
        }
        
        output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
        rownames(output) = q
        colnames(output) = prop.v
        output = reshape2::melt(output)
        
        colnames(output) = c("Order.q","prop.v","qmiNEXT_FD_singletau")
        
        res = output$qmiNEXT_FD_singletau
        res
      },future.seed = NULL) %>% cbind()
      se = abind::abind(se,along = 2)
      se = apply(se, 1, sd)
      
      output$LCL = output$qmiNEXT_FD_singletau - 1.96*se
      output$UCL = output$qmiNEXT_FD_singletau + 1.96*se
      
    }else{
      output$LCL = NA
      output$UCL = NA
    }
    return(list(out = output,ori.prop = prop.v,m.v = m.v,line_type = line_type))
    
  }else{
    
    if(nrow(m.ext.v) != 0){
      aivi.ext = lapply(1:nrow(m.v), function (k) mix.aivi.popu.FD.single(p_bind, m.v[k,1], m.v[k,2], FDdistM = FDdistM,
                                                                          tau = tau))
    }
    
    
    # Li.v = apply(Lis, 1, function(x) x[1]*x[2])
    S <- nrow(data)
    qlength <- length(q)
    out = list()
    
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
          out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)
        } else if (q[j] == 1) {
          out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S))
        } else if (q[j] == 2) {
          out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
        } else {
          out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)^(1/(1-q[j]))
        }
      }
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
      h1.FD = h1_hat_FD(ai1.ext.v, ai2.ext.v, vi.all.ext.v, p_bind[,1], p_bind[,2],
                        m1, m2_all, n1, n2)
      
      
      ## ===================================================
      for (j in 1:qlength) {
        if (q[j] == 0) {
          out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S) + h0.FD
        } else if (q[j] == 1) {
          out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S) + h1.FD)
        } else if (q[j] == 2) {
          # 記得放m2_all
          out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2_all, n1, n2)
        } else {
          out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)^(1/(1-q[j]))
        }
      }
    }
    # ==============================================================
    output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
    rownames(output) = q
    colnames(output) = prop.v
    output = reshape2::melt(output)
    colnames(output) = c("Order.q","prop.v","qmiNEXT_FD_singletau")
    
    # if(nrow(m.inc.v) != 0){
    #   aivi.inc = lapply(1:nrow(m.inc.v),
    #                     function (k) mix.aivi.sample.FD.single(data, m.inc.v[k,1], m.inc.v[k,2],
    #                                                            FDdistM = FDdistM, tau = tau))
    #   un1 = un2 = sh12 = c()
    #   for (i in 1:length(aivi.inc)) {
    #     datanew = cbind(aivi.inc[[i]]$ai1_mat,aivi.inc[[i]]$ai2_mat)
    #     datash = datanew[(datanew[,1]>0 & datanew[,2]>0), , drop=F]
    #     dataun1 = datanew[(datanew[,1]>0 & datanew[,2]==0), , drop=F]
    #     dataun2 = datanew[(datanew[,1]==0 & datanew[,2]>0), , drop=F]
    #     un1[i] = sum(un_abun_FD(vi = aivi.inc[[i]]$vi1_mat,xi = aivi.inc[[i]]$ai1_mat, n = n1, m = m.inc.v[i,1]))
    #     un2[i] = sum(un_abun_FD(vi = aivi.inc[[i]]$vi2_mat,xi = aivi.inc[[i]]$ai2_mat, n = n2, m = m.inc.v[i,2]))
    #     sh12[i] =  sum(sh_abun_FD(vi = aivi.inc[[i]]$vi_all_mat,xi1 = aivi.inc[[i]]$ai1_mat, xi2 = aivi.inc[[i]]$ai2_mat,
    #                               n1 = n1, m1 =  m.inc.v[i,1],
    #                               n2 = n2, m2 =  m.inc.v[i,2]))
    #     
    #   }
    #   
    #   
    #   q0_ana = data.frame(prop.v = prop.inc ,
    #                       q0_un1 = un1, q0_un2 = un2,
    #                       q0_sh = sh12)
    #   q0_ana = reshape2::melt(q0_ana,id = c('prop.v'))
    #   colnames(q0_ana)[2:3] = c('Type','Value')  
    # }else{
    #   q0_ana = NULL
    # }
    
    if(nboot>0){
      data_gamma = rowSums(data)
      se = future.apply::future_lapply(1:nboot, function(boottime){
        p_bt = iNEXT.beta3D:::bootstrap_population_multiple_assemblage(data, 
                                                                       data_gamma, "abundance")
        
        f0_hat = nrow(p_bt) - nrow(data)
        distance_matrix_bt = iNEXT.beta3D:::Bootstrap_distance_matrix(rowSums(data), 
                                                                      FDdistM, f0_hat, "abundance")
        
        data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[, k]), prob = p_bt[, k]))
        
        
        
        
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
              out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)
            } else if (q[j] == 1) {
              out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S))
            } else if (q[j] == 2) {
              out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2, n1, n2)
            } else {
              out[[zz]][j,1] <- FD_h_hat_fn(m1, m2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)^(1/(1-q[j]))
            }
          }
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
          h1.FD = h1_hat_FD(ai1.ext.v, ai2.ext.v, vi.all.ext.v, p_bind_bt[,1], p_bind_bt[,2],
                            m1, m2_all, n1, n2)
          ## ===================================================
          for (j in 1:qlength) {
            if (q[j] == 0) {
              out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S) + h0.FD
            } else if (q[j] == 1) {
              out[[zz]][j,1] <- exp(FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S) + h1.FD)
            } else if (q[j] == 2) {
              out[[zz]][j,1] <- FDq2(ai1.v, ai2.v, vi.all.v, m1, m2_all, n1, n2)
            } else {
              out[[zz]][j,1] <- FD_h_hat_fn(m1, n2, n1, n2, vi.all.v, ai1.v, ai2.v, q[j], S)^(1/(1-q[j]))
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
        res
      },future.seed = NULL) %>% cbind()
      se = abind::abind(se,along = 2)
      se = apply(se, 1, sd)
      
      output$LCL = output$qmiNEXT_FD_singletau - 1.96*se
      output$UCL = output$qmiNEXT_FD_singletau + 1.96*se
      
    }else{
      output$LCL = NA
      output$UCL = NA
    }
    return(list(out = output,ori.prop = prop.v,m.v = m.v,line_type = line_type))
  }
  
  
}