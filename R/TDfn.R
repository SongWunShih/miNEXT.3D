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

bootstrap_population_multiple_assemblage <- function (data, data_gamma, datatype){
  if (datatype == "abundance") {
    S_obs = sum(data_gamma > 0)
    n = sum(data_gamma)
    f1 = sum(data_gamma == 1)
    f2 = sum(data_gamma == 2)
    f0_hat = ifelse(f2 == 0, (n - 1)/n * f1 * (f1 - 1)/2,
                    (n - 1)/n * f1^2/2/f2) %>% ceiling()
    output = apply(data, 2, function(x) {
      p_i_hat = iNEXT.3D:::EstiBootComm.Ind(Spec = x)
      if (length(p_i_hat) != length(x)) {
        p_i_hat_unobs = p_i_hat[(length(x) + 1):length(p_i_hat)]
        p_i_hat_obs = p_i_hat[1:length(x)]
        p_i_hat = c(p_i_hat_obs, rep(0, f0_hat))
        candidate = which(p_i_hat == 0)
        chosen = sample(x = candidate, size = min(length(p_i_hat_unobs),
                                                  length(candidate)), replace = F)
        p_i_hat[chosen] = (1 - sum(p_i_hat))/length(chosen)
        p_i_hat
      }
      else {
        p_i_hat = c(p_i_hat, rep(0, f0_hat))
        p_i_hat
      }
    })
  }
  if (datatype == "incidence") {
    S_obs = sum(data_gamma > 0)
    t = data_gamma[1]
    Q1 = sum(data_gamma == 1)
    Q2 = sum(data_gamma == 2)
    Q0_hat = if (Q2 == 0) {
      ((t - 1)/t) * (Q1 * (Q1 - 1)/2)
    }
    else {
      ((t - 1)/t) * ((Q1^2)/(2 * Q2))
    } %>% ceiling
    output = apply(data, 2, function(x) {
      pi_i_hat = iNEXT.3D:::EstiBootComm.Sam(Spec = x)
      if (length(pi_i_hat) != (length(x) - 1)) {
        pi_i_hat_unobs = pi_i_hat[length(x):length(pi_i_hat)]
        pi_i_hat_obs = pi_i_hat[1:(length(x) - 1)]
        pi_i_hat = c(pi_i_hat_obs, rep(0, Q0_hat))
        candidate = which(pi_i_hat == 0)
        chosen = sample(x = candidate, size = min(length(pi_i_hat_unobs),
                                                  length(candidate)), replace = F)
        pi_i_hat[chosen] = pi_i_hat_unobs
        pi_i_hat
      }
      else {
        pi_i_hat = c(pi_i_hat, rep(0, Q0_hat))
        pi_i_hat
      }
    })
  }
  return(output)
}

TDq2<-function(x1,x2,m1,m2,n1,n2){

  D2 = 0
  if(m1 !=0 || m2!=0){
    D2 = 1/(m1+m2)+sum(x1*(x1-1)/(n1*(n1-1)))*m1*(m1-1)/(m1+m2)^2+sum(x2*(x2-1)/(n2*(n2-1)))*m2*(m2-1)/(m1+m2)^2+2*sum(x1*x2/(n1*n2))*m1*m2/(m1+m2)^2
    D2<-1/D2
  }
  else if(m1==0 && m2==0)
    D2<-0
  return(D2)
}

## TD ==================================================================
TD.f.theo_fn = function(m1,m2,ai1,ai2,S){
  # tmp = matrix(0, length(q), 1)
  gtheo.tmp = matrix(0,m1+1,m2+1)
  for (i in 1:S) {
    for (k1 in 0:m1) {
      for (k2 in 0:m2) {
        if((k1+k2) == 0){
          gtheo.tmp[k1+1,k2+1] = gtheo.tmp[k1+1,k2+1]
        }else{
          gtheo.tmp[k1+1,k2+1] = gtheo.tmp[k1+1,k2+1] + exp(lchoose(m1,k1)+lchoose(m2,k2)) * (ai1[i]^k1) * ((1-ai1[i])^(m1-k1)) * (ai2[i]^k2) * ((1-ai2[i])^(m2-k2))
        }
      }
    }
  }

  return(gtheo.tmp)
}
RTD.obs <- function(data, knots = 11, size = NULL, q = c(0,1,2), n1, n2) {

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

  S = nrow(data) ## 注意是B不是S
  qlength <- length(q)
  out = list()
  ghat_pt2 = list()
  ai1.v = data[,1]
  ai2.v = data[,2]

  for (zz in 1:length(m1.v)) {


    out[[zz]] <- matrix(0, qlength, 1)
    m1 = m1.v[zz]
    m2 = m2.v[zz]

    ghat_pt2[[zz]] <- TD.f.theo_fn(m1, m2, ai1.v, ai2.v , S)




    for (j in 1:qlength) {
      for (k1 in 0:m1) {
        for (k2 in 0:m2) {
          if (q[j] == 0) {
            out[[zz]][j,1] <- out[[zz]][j,1] + ghat_pt2[[zz]][k1+1,k2+1]
          }else if (q[j] == 1) {
            if(k1 == 0 & k2 == 0){
              aa = 0
            }else{
              aa = ((k1+k2)/((m1+m2))) * log((k1+k2)/((m1+m2)))
            }
            out[[zz]][j,1] <- out[[zz]][j,1] - aa  * ghat_pt2[[zz]][k1+1,k2+1]
          }else if (q[j] == 2) {
            out[[zz]][j,1] <- out[[zz]][j,1] + ((k1+k2)/((m1+m2)))^2 * ghat_pt2[[zz]][k1+1,k2+1]
          }else {
            out[[zz]][j,1] <- out[[zz]][j,1] + ((k1+k2)/((m1+m2)))^q[j] * ghat_pt2[[zz]][k1+1,k2+1]
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
  colnames(output) = c("Order.q","prop.v","qmiNEXT_TD")
  output

  return(list(out = output,ori.prop = prop.v,m.v = m.v,line_type = line_type))
}
RTD.est <- function(data, knots = 11, size = NULL, q = c(0,1,2), conf = 0.95, nboot = 0) {

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

  if(n1<=n2){
    # S <- nrow(data)
    qlength <- length(q)
    out = list()
    ai1.v = data[,1]
    ai2.v = data[,2]


    for (zz in 1:length(m1.v)) {


      out[[zz]] <- matrix(0, qlength, 1)
      m1 = m1.v[zz]
      m2 = m2.v[zz]

      for (j in 1:qlength) {
        if (q[j] == 0) {
          out[[zz]][j,1] <- TD_f_hat_fn(m1, m2, n1, n2, ai1.v, ai2.v, q[j])
        } else if (q[j] == 1) {
          out[[zz]][j,1] <- exp(TD_f_hat_fn(m1, m2, n1, n2, ai1.v, ai2.v, q[j]))
        } else if (q[j] == 2) {
          out[[zz]][j,1] <- TDq2(ai1.v, ai2.v, m1, m2, n1, n2)
        } else {
          out[[zz]][j,1] <- TD_f_hat_fn(m1, m2, n1, n2, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
        }
      }

    }

    output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
    rownames(output) = q
    colnames(output) = prop.v
    output = reshape2::melt(output)
    colnames(output) = c("Order.q","prop.v","qmiNEXT_TD")
    output
    ## q0 ana =========================================
    datash = data[(data[,1]>0 & data[,2]>0), , drop=F]
    dataun1 = data[(data[,1]>0 & data[,2]==0), , drop=F]
    dataun2 = data[(data[,1]==0 & data[,2]>0), , drop=F]

    un1 = sapply(m1.v,function(i){
      sum(un_abun(xi = dataun1[,1], n = n1, m = i))
    })
    un2 = sapply(m2.v,function(i){
      sum(un_abun(xi = dataun2[,2], n = n2, m = i))
    })
    sh12 = sapply(1:length(m1.v),function(i){
      sum(sh_abun(xi1 = datash[,1],xi2 = datash[,2], n1 = n1, m1 =  m1.v[i],
                  n2 = n2, m2 = m2.v[i]))
    })
    q0_ana = data.frame(prop.v = prop.v ,
                        q0_un1 = un1, q0_un2 = un2,
                        q0_sh = sh12)
    q0_ana = reshape2::melt(q0_ana,id = c('prop.v'))
    colnames(q0_ana)[2:3] = c('Type','Value')
    # ====================================================================
    if(nboot > 0){
      data_gamma = rowSums(data)
      se = future.apply::future_lapply(1:nboot, function(boottime) {
        p_bt = bootstrap_population_multiple_assemblage(data,data_gamma, "abundance")
        data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1,
                                                             size = sum(data[, k]), prob = p_bt[, k]))

        # S <- nrow(data_bt)
        qlength <- length(q)
        out = list()

        ai1.v = data_bt[,1]
        ai2.v = data_bt[,2]
        for (zz in 1:length(m1.v)) {


          out[[zz]] <- matrix(0, qlength, 1)
          m1 = m1.v[zz]
          m2 = m2.v[zz]

          for (j in 1:qlength) {
            if (q[j] == 0) {
              out[[zz]][j,1] <- TD_f_hat_fn(m1, m2, n1, n2, ai1.v, ai2.v, q[j])
            } else if (q[j] == 1) {
              out[[zz]][j,1] <- exp(TD_f_hat_fn(m1, m2, n1, n2, ai1.v, ai2.v, q[j]))
            } else if (q[j] == 2) {
              out[[zz]][j,1] <- TDq2(ai1.v, ai2.v, m1, m2, n1, n2)
            } else {
              out[[zz]][j,1] <- TD_f_hat_fn(m1, m2, n1, n2, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
            }
          }

        }



        output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
        rownames(output) = q
        colnames(output) = prop.v
        output = reshape2::melt(output)
        colnames(output) = c("Order.q","prop.v","qmiNEXT_TD")
        res = output$qmiNEXT_TD

        ## q0 ana =========================================
        datash = data_bt[(data_bt[,1]>0 & data_bt[,2]>0), , drop=F]
        dataun1 = data_bt[(data_bt[,1]>0 & data_bt[,2]==0), , drop=F]
        dataun2 = data_bt[(data_bt[,1]==0 & data_bt[,2]>0), , drop=F]

        un1 = sapply(m1.v,function(i){
          sum(un_abun(xi = dataun1[,1], n = n1, m = i))
        })
        un2 = sapply(m2.v,function(i){
          sum(un_abun(xi = dataun2[,2], n = n2, m = i))
        })
        sh12 = sapply(1:length(m1.v),function(i){
          sum(sh_abun(xi1 = datash[,1],xi2 = datash[,2], n1 = n1, m1 =  m1.v[i],
                      n2 = n2, m2 = m2.v[i]))
        })
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
      output$LCL = output$qmiNEXT_TD - qtile*se_all
      output$UCL = output$qmiNEXT_TD + qtile*se_all

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

    # S <- nrow(data)
    qlength <- length(q)
    out = list()
    ai1.v = data[,1]
    ai2.v = data[,2]


    for (zz in inc.idx) {


      out[[zz]] <- matrix(0, qlength, 1)
      m1 = m1.v[zz]
      m2 = m2.v[zz]

      for (j in 1:qlength) {
        if (q[j] == 0) {
          out[[zz]][j,1] <- TD_f_hat_fn(m1, m2, n1, n2, ai1.v, ai2.v, q[j])
        } else if (q[j] == 1) {
          out[[zz]][j,1] <- exp(TD_f_hat_fn(m1, m2, n1, n2, ai1.v, ai2.v, q[j]))
        } else if (q[j] == 2) {
          out[[zz]][j,1] <- TDq2(ai1.v, ai2.v, m1, m2, n1, n2)
        } else {
          out[[zz]][j,1] <- TD_f_hat_fn(m1, m2, n1, n2, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
        }
      }
    }

    for (zz in ext.idx) {

      out[[zz]] <- matrix(0, qlength, 1)
      m1 = m1.v[zz]
      m2_all = m2.v[zz]
      m2_s = m2.v[zz] - n2

      ## new =============================================
      h0.TD = h0_hat_cpp(p_bind[,1],p_bind[,2],m1,m2_s,n1,n2)
      h1.TD = h1_hat(p_bind[,1], p_bind[,2], ai1.v, ai2.v,
                     m1, m2_all, n1, n2)



      for (j in 1:qlength) {
        if (q[j] == 0) {
          out[[zz]][j,1] <- TD_f_hat_fn(m1, n2, n1, n2, ai1.v, ai2.v, q[j]) + h0.TD
        } else if (q[j] == 1) {
          out[[zz]][j,1] <- exp(TD_f_hat_fn(m1, n2, n1, n2, ai1.v, ai2.v, q[j]) + h1.TD)
        } else if (q[j] == 2) {
          out[[zz]][j,1] <- TDq2(ai1.v, ai2.v, m1, m2_all, n1, n2)
        } else {
          out[[zz]][j,1] <- TD_f_hat_fn(m1, n2, n1, n2, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
        }
      }
    }

    output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
    rownames(output) = q
    colnames(output) = prop.v
    output = reshape2::melt(output)
    colnames(output) = c("Order.q","prop.v","qmiNEXT_TD")

    ## q0 ana =========================================
    if(nrow(m.inc.v) != 0){
      datash = data[(data[,1]>0 & data[,2]>0), , drop=F]
      dataun1 = data[(data[,1]>0 & data[,2]==0), , drop=F]
      dataun2 = data[(data[,1]==0 & data[,2]>0), , drop=F]

      un1 = sapply(m.inc.v[,1],function(i){
        sum(un_abun(xi = dataun1[,1], n = n1, m = i))
      })
      un2 = sapply(m.inc.v[,2],function(i){
        sum(un_abun(xi = dataun2[,2], n = n2, m = i))
      })
      sh12 = sapply(1:length(m.inc.v[,1]),function(i){
        sum(sh_abun(xi1 = datash[,1],xi2 = datash[,2], n1 = n1, m1 =  m.inc.v[i,1],
                    n2 = n2, m2 = m.inc.v[i,2]))
      })
      q0_ana = data.frame(prop.v = prop.inc ,
                          q0_un1 = un1, q0_un2 = un2,
                          q0_sh = sh12)
      q0_ana = reshape2::melt(q0_ana,id = c('prop.v'))
      colnames(q0_ana)[2:3] = c('Type','Value')
    }else{
      q0_ana = NULL
    }


    if(nboot > 0){
      data_gamma = rowSums(data)
      se = future.apply::future_lapply(1:nboot, function(boottime) {
        p_bt = bootstrap_population_multiple_assemblage(data,data_gamma, "abundance")
        data_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1,
                                                             size = sum(data[, k]), prob = p_bt[, k]))

        # new ==================================
        p_bind_bt = cbind(DetAbu(data_bt[,1], zero = T),DetAbu(data_bt[,2], zero = T))
        # ==========================================
        # S <- nrow(data_bt)
        qlength <- length(q)
        out = list()

        ai1.v = data_bt[,1]
        ai2.v = data_bt[,2]


        for (zz in inc.idx) {


          out[[zz]] <- matrix(0, qlength, 1)
          m1 = m1.v[zz]
          m2 = m2.v[zz]

          for (j in 1:qlength) {
            if (q[j] == 0) {
              out[[zz]][j,1] <- TD_f_hat_fn(m1, m2, n1, n2, ai1.v, ai2.v, q[j])
            } else if (q[j] == 1) {
              out[[zz]][j,1] <- exp(TD_f_hat_fn(m1, m2, n1, n2, ai1.v, ai2.v, q[j]))
            } else if (q[j] == 2) {
              out[[zz]][j,1] <- TDq2(ai1.v, ai2.v, m1, m2, n1, n2)
            } else {
              out[[zz]][j,1] <- TD_f_hat_fn(m1, m2, n1, n2, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
            }
          }
        }

        for (zz in ext.idx) {

          out[[zz]] <- matrix(0, qlength, 1)
          m1 = m1.v[zz]
          m2_all = m2.v[zz]
          m2_s = m2.v[zz] - n2

          ## new =============================================
          h0.TD = h0_hat_cpp(p_bind_bt[,1],p_bind_bt[,2],m1,m2_s,n1,n2)
          h1.TD = h1_hat(p_bind_bt[,1], p_bind_bt[,2], ai1.v, ai2.v,
                         m1, m2_all, n1, n2)



          for (j in 1:qlength) {
            if (q[j] == 0) {
              out[[zz]][j,1] <- TD_f_hat_fn(m1, n2, n1, n2, ai1.v, ai2.v, q[j]) + h0.TD
            } else if (q[j] == 1) {
              out[[zz]][j,1] <- exp(TD_f_hat_fn(m1, n2, n1, n2, ai1.v, ai2.v, q[j]) + h1.TD)
            } else if (q[j] == 2) {
              out[[zz]][j,1] <- TDq2(ai1.v, ai2.v, m1, m2_all, n1, n2)
            } else {
              out[[zz]][j,1] <- TD_f_hat_fn(m1, n2, n1, n2, ai1.v, ai2.v, q[j])^(1/(1-q[j]))
            }
          }
        }

        output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
        rownames(output) = q
        colnames(output) = prop.v
        output = reshape2::melt(output)
        colnames(output) = c("Order.q","prop.v","qmiNEXT_TD")
        res = output$qmiNEXT_TD

        ## q0 ana =========================================
        if(nrow(m.inc.v) != 0){
          datash = data_bt[(data_bt[,1]>0 & data_bt[,2]>0), , drop=F]
          dataun1 = data_bt[(data_bt[,1]>0 & data_bt[,2]==0), , drop=F]
          dataun2 = data_bt[(data_bt[,1]==0 & data_bt[,2]>0), , drop=F]

          un1 = sapply(m.inc.v[,1],function(i){
            sum(un_abun(xi = dataun1[,1], n = n1, m = i))
          })
          un2 = sapply(m.inc.v[,2],function(i){
            sum(un_abun(xi = dataun2[,2], n = n2, m = i))
          })
          sh12 = sapply(1:length(m.inc.v[,1]),function(i){
            sum(sh_abun(xi1 = datash[,1],xi2 = datash[,2], n1 = n1, m1 =  m.inc.v[i,1],
                        n2 = n2, m2 = m.inc.v[i,2]))
          })
          q0_ana = data.frame(prop.v = prop.inc,
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
      output$LCL = output$qmiNEXT_TD - qtile*se_all
      output$UCL = output$qmiNEXT_TD + qtile*se_all

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
