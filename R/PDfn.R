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

PDq2<-function(x1,x2,Li.all,Tbar,m1,m2,n1,n2){
  x1 = x1[Li.all>0]
  x2 = x2[Li.all>0]
  Li.all = Li.all[Li.all>0]
  D2 = 0
  if(m1 !=0 || m2!=0){
    part1 = 1/((m1+m2)*Tbar)^2*sum(m1/n1*Li.all*x1 + m2/n2*Li.all*x2)
    part2 = (m1*(m1-1))/((m1+m2)*Tbar)^2*sum(Li.all*x1*(x1-1))/(n1*(n1-1))
    part3 = (m2*(m2-1))/((m1+m2)*Tbar)^2*sum(Li.all*x2*(x2-1))/(n2*(n2-1))
    part4 = 2*m1*m2/((m1+m2)*Tbar)^2*sum(Li.all*x1*x2)/(n1*n2)
    D2 = 1/(part1 + part2 + part3 + part4)
  }
  else if(m1==0 && m2==0)
    D2<-0
  return(D2)
}

# PD ===========================================================================
PD.g.theo_fn = function(m1,m2,ai1,ai2,Li,S){
  # tmp = matrix(0, length(q), 1)
  gtheo.tmp = matrix(0,m1+1,m2+1)
  for (i in 1:S) {
    for (k1 in 0:m1) {
      for (k2 in 0:m2) {
        if((k1+k2) == 0){
          gtheo.tmp[k1+1,k2+1] = gtheo.tmp[k1+1,k2+1]
        }else{
          gtheo.tmp[k1+1,k2+1] = gtheo.tmp[k1+1,k2+1] + Li[i] * exp(lchoose(m1,k1)+lchoose(m2,k2)) * (ai1[i]^k1) * ((1-ai1[i])^(m1-k1)) * (ai2[i]^k2) * ((1-ai2[i])^(m2-k2))
        }
      }
    }
  }

  return(gtheo.tmp)
}
mix.Liai.popu.PD = function (data, oritree, datatype, PDreftime, nT) {
  rm.idx = which(rowSums(data) == 0)
  if(sum(rm.idx) == 0){
    data1 = data
  }else{
    data1 = data[-rm.idx,]
  }

  # data.list = list()
  # for (i in 1:ncol(data1)) {
  #   data.list[[i]] = data1[,i]
  #   names(data.list[[i]]) = rownames(data1)
  # }
  # data.sum = rowSums(data1)
  # names(data.sum) = rownames(data1)

  checktree = iNEXT.3D:::check.tree(data, datatype, oritree, PDreftime,
                                    nT)
  PDreftime1 = checktree[[1]]
  mytree = checktree[[2]]
  mydata = checktree[[3]]
  names(mydata[[1]]) = names(mydata[[2]]) = rownames(data1)

  aL1 = iNEXT.3D:::phyBranchAL_Abu(phylo = mytree, data = mydata[[1]],
                                   datatype, refT = PDreftime1,remove0 = F)

  aL2 = iNEXT.3D:::phyBranchAL_Abu(phylo = mytree, data = mydata[[2]],
                                   datatype, refT = PDreftime1,remove0 = F)
  # 母體row sum後要除2
  aL = iNEXT.3D:::phyBranchAL_Abu(phylo = mytree, data = rowSums(data1)/2,
                                  datatype, refT = PDreftime1)

  # 幫Li加名字，確保有對上ai
  names(aL$treeNabu$branch.length) = names(aL$treeNabu$branch.abun)
  return(list(ai1.v = aL1$treeNabu$branch.abun, ai2.v = aL2$treeNabu$branch.abun,
              ai_all.v = aL$treeNabu$branch.abun, Li_all.v = aL$treeNabu$branch.length,
              Tbar = aL$treeH))
}
mix.Liai.sample.PD = function (data, oritree, datatype, PDreftime, nT) {
  rm.idx = which(rowSums(data) == 0)
  if(sum(rm.idx) == 0){
    data1 = data
  }else{
    data1 = data[-rm.idx,]
  }

  # data.list = list()
  # for (i in 1:ncol(data1)) {
  #   data.list[[i]] = data1[,i]
  #   names(data.list[[i]]) = rownames(data1)
  # }
  # data.sum = rowSums(data1)
  # names(data.sum) = rownames(data1)

  checktree = iNEXT.3D:::check.tree(data, datatype, oritree, PDreftime,
                                    nT)
  PDreftime1 = checktree[[1]]
  mytree = checktree[[2]]
  mydata = checktree[[3]]
  names(mydata[[1]]) = names(mydata[[2]]) = rownames(data1)

  aL1 = iNEXT.3D:::phyBranchAL_Abu(phylo = mytree, data = mydata[[1]],
                                   datatype, refT = PDreftime1,remove0 = F)

  aL2 = iNEXT.3D:::phyBranchAL_Abu(phylo = mytree, data = mydata[[2]],
                                   datatype, refT = PDreftime1,remove0 = F)

  aL = iNEXT.3D:::phyBranchAL_Abu(phylo = mytree, data = rowSums(data1),
                                  datatype, refT = PDreftime1)

  # 幫Li加名字，確保有對上ai
  names(aL$treeNabu$branch.length) = names(aL$treeNabu$branch.abun)
  return(list(ai1.v = aL1$treeNabu$branch.abun, ai2.v = aL2$treeNabu$branch.abun,
              ai_all.v = aL$treeNabu$branch.abun, Li_all.v = aL$treeNabu$branch.length,
              Tbar = aL$treeH,refT = PDreftime1))
}
RPD.obs <- function(data, PDtree, PDreftime = NULL, knots = 11, size = NULL, q = c(0,1,2), n1, n2) {

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

  Livi.popu.data = mix.Liai.popu.PD(data = data, oritree = PDtree, datatype = "abundance",
                                    PDreftime = PDreftime, nT = "NULL")

  B <- length(Livi.popu.data$ai_all.v) ## 注意是B不是S
  qlength <- length(q)
  out = list()
  ghat_pt2 = list()
  ai1.v = Livi.popu.data$ai1.v
  ai2.v = Livi.popu.data$ai2.v
  ai.all.v = Livi.popu.data$ai_all.v
  Li.all.v = Livi.popu.data$Li_all.v
  Tbar = Livi.popu.data$Tbar

  for (zz in 1:length(m1.v)) {


    out[[zz]] <- matrix(0, qlength, 1)
    m1 = m1.v[zz]
    m2 = m2.v[zz]

    ghat_pt2[[zz]] <- PD.g.theo_fn(m1, m2, ai1.v, ai2.v ,Li.all.v, B)




    for (j in 1:qlength) {
      for (k1 in 0:m1) {
        for (k2 in 0:m2) {
          if (q[j] == 0) {
            out[[zz]][j,1] <- out[[zz]][j,1] + ghat_pt2[[zz]][k1+1,k2+1]
          }else if (q[j] == 1) {
            if(k1 == 0 & k2 == 0){
              aa = 0
            }else{
              aa = ((k1+k2)/((m1+m2)*Tbar)) * log((k1+k2)/((m1+m2)*Tbar))
            }
            out[[zz]][j,1] <- out[[zz]][j,1] - aa  * ghat_pt2[[zz]][k1+1,k2+1]
          }else if (q[j] == 2) {
            out[[zz]][j,1] <- out[[zz]][j,1] + ((k1+k2)/((m1+m2)*Tbar))^2 * ghat_pt2[[zz]][k1+1,k2+1]
          }else {
            out[[zz]][j,1] <- out[[zz]][j,1] + ((k1+k2)/((m1+m2)*Tbar))^q[j] * ghat_pt2[[zz]][k1+1,k2+1]
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

      out[[zz]][j,1] = out[[zz]][j,1]/Tbar

    }
  }

  output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
  rownames(output) = q
  colnames(output) = prop.v
  output = reshape2::melt(output)
  colnames(output) = c("Order.q","prop.v","qmiNEXT_PD")
  output
  return(list(out = output,ori.prop = prop.v,m.v = m.v,line_type = line_type))

}
RPD.est <- function(data, PDtree, PDreftime = NULL, knots = 11, size = NULL, q = c(0,1,2), conf = 0.95, nboot = 0) {

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

  # checktree = iNEXT.3D:::check.tree(data, datatype, PDtree, PDreftime,
  #                                   nT)

  Liai.sample.data = mix.Liai.sample.PD(data = data, oritree = PDtree, datatype = "abundance",
                                        PDreftime = PDreftime, nT = NULL)

  reft = Liai.sample.data$refT

  if(n1<=n2){
    # B <- length(Liai.sample.data$Li_all.v) ## 注意是B不是S
    qlength <- length(q)
    out = list()
    ai1.v = Liai.sample.data$ai1.v
    ai2.v = Liai.sample.data$ai2.v
    # ai.all.v = Liai.sample.data$ai_all.v
    Li.all.v = Liai.sample.data$Li_all.v
    Tbar = Liai.sample.data$Tbar

    for (zz in 1:length(m1.v)) {
      out[[zz]] <- matrix(0, qlength, 1)
      m1 = m1.v[zz]
      m2 = m2.v[zz]

      for (j in 1:qlength) {
        if (q[j] == 0) {
          out[[zz]][j,1] <- PD_g_hat_fn(m1, m2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar)
        } else if (q[j] == 1) {
          out[[zz]][j,1] <- exp(PD_g_hat_fn(m1, m2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar))
        } else if (q[j] == 2) {
          out[[zz]][j,1] <- PDq2(ai1.v, ai2.v, Li.all.v, Tbar, m1, m2, n1, n2)
        } else {
          out[[zz]][j,1] <- PD_g_hat_fn(m1, m2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar)^(1/(1-q[j]))
        }

        out[[zz]][j,1] = out[[zz]][j,1]/Tbar
      }
    }

    output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
    rownames(output) = q
    colnames(output) = prop.v
    output = reshape2::melt(output)
    colnames(output) = c("Order.q","prop.v","qmiNEXT_PD")
    output

    ## q0 ana =========================================
    datanew = cbind(ai1.v,ai2.v)
    datash = datanew[(datanew[,1]>0 & datanew[,2]>0), , drop=F]
    dataun1 = datanew[(datanew[,1]>0 & datanew[,2]==0), , drop=F]
    dataun2 = datanew[(datanew[,1]==0 & datanew[,2]>0), , drop=F]

    un1 = sapply(m1.v,function(i){
      sum(un_abun_PD(Li = Li.all.v[rownames(dataun1)],xi = dataun1[,1], n = n1, m = i))
    })
    un2 = sapply(m2.v,function(i){
      sum(un_abun_PD(Li = Li.all.v[rownames(dataun2)],xi = dataun2[,2], n = n2, m = i))
    })
    sh12 = sapply(1:length(m1.v),function(i){
      sum(sh_abun_PD(Li = Li.all.v[rownames(datash)],xi1 = datash[,1],xi2 = datash[,2], n1 = n1, m1 =  m1.v[i],
                     n2 = n2, m2 = m2.v[i]))
    })

    un1 = un1/Tbar
    un2 = un2/Tbar
    sh12 = sh12/Tbar

    q0_ana = data.frame(prop.v = prop.v ,
                        q0_un1 = un1, q0_un2 = un2,
                        q0_sh = sh12)
    q0_ana = reshape2::melt(q0_ana,id = c('prop.v'))
    colnames(q0_ana)[2:3] = c('Type','Value')
    # ============================================================================
    if(nboot > 0){
      data_gamma = rowSums(data)
      se = future.apply::future_lapply(1:nboot, function(boottime) {
        tree_bt = PDtree
        p_bt = bootstrap_population_multiple_assemblage(data,data_gamma, "abundance")
        unseen_p = p_bt[-(1:nrow(data)), ] %>% matrix(ncol = ncol(data))
        #  =========================================================================
        if (nrow(p_bt) > nrow(data) & sum(unseen_p) > 0) {
          unseen = unseen_p[which(rowSums(unseen_p) >
                                    0), ]
          unseen = matrix(unseen, ncol = ncol(unseen_p),
                          byrow = T)
          p_bt = rbind(p_bt[(1:nrow(data)), ], unseen)
          unseen_name = sapply(1:nrow(unseen), function(i) paste0("unseen_",
                                                                  i))
          rownames(p_bt) = c(rownames(data), unseen_name)
          x_bt = sapply(1:ncol(data),
                        function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[, k]))

          rownames(x_bt) = rownames(p_bt)
          if (sum(x_bt[-(1:nrow(data)), ]) > 0) {
            g0_hat = apply(data, 2, function(x) {
              n = sum(x)
              f1 = sum(x == 1)
              f2 = sum(x == 2)
              # remove0 = F 看要不要加 (有沒有加一樣)
              aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree,remove0 = F,
                                              data = x, rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,
                                                   1]
              aL = aL$treeNabu %>% select(branch.abun,
                                          branch.length)
              g1 = aL$branch.length[aL$branch.abun ==
                                      1] %>% sum
              g2 = aL$branch.length[aL$branch.abun ==
                                      2] %>% sum
              g0_hat = ifelse(g2 > ((g1 * f2)/(2 *
                                                 f1)), ((n - 1)/n) * (g1^2/(2 * g2)),
                              ((n - 1)/n) * (g1 * (f1 - 1)/(2 *
                                                              (f2 + 1))))
              g0_hat[is.na(g0_hat)] = 0
              g0_hat
            })
            # te表示新的bootstrap樣本中有看到但原本的樣本沒看到的物種
            te = (x_bt[1:nrow(data), ] * (data == 0)) > 0
            # used_length表示bootstrap中有看到但原本沒看到的branch length和
            used_length = sapply(1:ncol(data), function(i) {
              if (sum(te[, i]) == 0)
                return(0)
              else {
                # remove0 = F 看要不要加
                iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree,
                                           data = x_bt[1:nrow(data), i],
                                           rootExtend = T, refT = reft)$treeNabu %>%
                  subset(label %in% names(which(te[,
                                                   i] == TRUE))) %>% select(branch.length) %>%
                  sum
              }
            })
            g0_hat = g0_hat - used_length
            g0_hat[g0_hat < 0] = 0
            unseen_sample = x_bt[-(1:nrow(data)),]
            if (is.vector(unseen_sample))
              unseen_sample = matrix(unseen_sample,
                                     ncol = ncol(x_bt), byrow = T)
            #L0_hat為g0_hat平分給所有沒看到的unseen sample
            L0_hat = sapply(1:length(g0_hat), function(i) if (sum(unseen_sample[,i] > 0) > 0)
              (g0_hat[i]/nrow(unseen))
              else 0)
            L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample),
                                     ncol(unseen_sample), byrow = T) *
                                unseen_sample))/rowSums(unseen_sample)
            L0_hat[which(rowSums(unseen_sample) == 0)] = 0
            for (i in 1:length(L0_hat)) {
              tip = list(edge = matrix(c(2, 1),
                                       1, 2), tip.label = unseen_name[i],
                         edge.length = L0_hat[i], Nnode = 1)
              class(tip) = "phylo"
              tree_bt = tree_bt + tip
            }
          }
          else {
            x_bt = x_bt[1:nrow(data), ]
            p_bt = p_bt[1:nrow(data), ]
          }
        }else {
          p_bt = p_bt[1:nrow(data), ]
          x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[, k]), prob = p_bt[,k]))
          rownames(x_bt) = rownames(data)
        }
        # ===============================================================================
        # p_bt
        # x_bt
        # tree_bt

        Liai.sample.data = mix.Liai.sample.PD(data = x_bt, oritree = tree_bt, datatype = "abundance",
                                              PDreftime = PDreftime, nT = NULL)

        B <- length(Liai.sample.data$Li_all.v) ## 注意是B不是S
        qlength <- length(q)
        out = list()
        ai1.v = Liai.sample.data$ai1.v
        ai2.v = Liai.sample.data$ai2.v
        ai.all.v = Liai.sample.data$ai_all.v
        Li.all.v = Liai.sample.data$Li_all.v
        # Tbar要改，因為bootstrap的tree為non-ultrametric
        Tbar = sum(Li.all.v*ai.all.v)/(n1 + n2)

        for (zz in 1:length(m1.v)) {
          out[[zz]] <- matrix(0, qlength, 1)
          m1 = m1.v[zz]
          m2 = m2.v[zz]

          for (j in 1:qlength) {
            if (q[j] == 0) {
              out[[zz]][j,1] <- PD_g_hat_fn(m1, m2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar)
            } else if (q[j] == 1) {
              out[[zz]][j,1] <- exp(PD_g_hat_fn(m1, m2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar))
            } else if (q[j] == 2) {
              out[[zz]][j,1] <- PDq2(ai1.v, ai2.v, Li.all.v, Tbar, m1, m2, n1, n2)
            } else {
              out[[zz]][j,1] <- PD_g_hat_fn(m1, m2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar)^(1/(1-q[j]))
            }

            out[[zz]][j,1] = out[[zz]][j,1]/Tbar
          }
        }


        output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
        rownames(output) = q
        colnames(output) = prop.v
        output = reshape2::melt(output)
        colnames(output) = c("Order.q","prop.v","qmiNEXT_PD")
        res = output$qmiNEXT_PD

        ## q0 ana =========================================
        datanew = cbind(ai1.v,ai2.v)
        datash = datanew[(datanew[,1]>0 & datanew[,2]>0), , drop=F]
        dataun1 = datanew[(datanew[,1]>0 & datanew[,2]==0), , drop=F]
        dataun2 = datanew[(datanew[,1]==0 & datanew[,2]>0), , drop=F]

        un1 = sapply(m1.v,function(i){
          sum(un_abun_PD(Li = Li.all.v[rownames(dataun1)],xi = dataun1[,1], n = n1, m = i))
        })
        un2 = sapply(m2.v,function(i){
          sum(un_abun_PD(Li = Li.all.v[rownames(dataun2)],xi = dataun2[,2], n = n2, m = i))
        })
        sh12 = sapply(1:length(m1.v),function(i){
          sum(sh_abun_PD(Li = Li.all.v[rownames(datash)],xi1 = datash[,1],xi2 = datash[,2], n1 = n1, m1 =  m1.v[i],
                         n2 = n2, m2 = m2.v[i]))
        })

        un1 = un1/Tbar
        un2 = un2/Tbar
        sh12 = sh12/Tbar

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
      output$LCL = output$qmiNEXT_PD - qtile*se_all
      output$UCL = output$qmiNEXT_PD + qtile*se_all

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

    ## new ================================
    if(nrow(m.ext.v) != 0){
      Liai.ext = mix.Liai.popu.PD(data = p_bind, oritree = PDtree, datatype = "abundance",
                                  PDreftime = PDreftime, nT = NULL)

      ai1.ext.v = Liai.ext$ai1.v
      ai2.ext.v = Liai.ext$ai2.v
      Li.all.ext.v = Liai.ext$Li_all.v
    }

    # =================================================================
    # B <- length(Liai.sample.data$Li_all.v) ## 注意是B不是S
    qlength <- length(q)
    out = list()
    ai1.v = Liai.sample.data$ai1.v
    ai2.v = Liai.sample.data$ai2.v
    # ai.all.v = Liai.sample.data$ai_all.v
    Li.all.v = Liai.sample.data$Li_all.v
    Tbar = Liai.sample.data$Tbar

    for (zz in inc.idx) {


      out[[zz]] <- matrix(0, qlength, 1)
      m1 = m1.v[zz]
      m2 = m2.v[zz]

      for (j in 1:qlength) {
        if (q[j] == 0) {
          out[[zz]][j,1] <- PD_g_hat_fn(m1, m2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar)
        } else if (q[j] == 1) {
          out[[zz]][j,1] <- exp(PD_g_hat_fn(m1, m2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar))
        } else if (q[j] == 2) {
          out[[zz]][j,1] <- PDq2(ai1.v, ai2.v, Li.all.v, Tbar, m1, m2, n1, n2)
        } else {
          out[[zz]][j,1] <- PD_g_hat_fn(m1, m2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar)^(1/(1-q[j]))
        }

        out[[zz]][j,1] = out[[zz]][j,1]/Tbar
      }
    }

    for (zz in ext.idx) {
      out[[zz]] <- matrix(0, qlength, 1)
      m1 = m1.v[zz]
      m2_all = m2.v[zz]
      m2_s = m2.v[zz] - n2
      # 下面的m2改n2
      ## new =============================================
      h0.PD = h0_hat_cpp_PD(ai1.ext.v,ai2.ext.v,Li.all.ext.v,m1,m2_s,n1,n2)
      ## 思考h1_hat_FD的xi1, xi2要放啥
      h1.PD = h1_hat_PD(ai1.ext.v, ai2.ext.v, Li.all.ext.v, Tbar, ai1.ext.v, ai2.ext.v,
                        m1, m2_all, n1, n2)
      ## ===================================================
      for (j in 1:qlength) {
        if (q[j] == 0) {
          out[[zz]][j,1] <- PD_g_hat_fn(m1, n2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar) + h0.PD
        } else if (q[j] == 1) {
          out[[zz]][j,1] <- exp(PD_g_hat_fn(m1, n2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar) + h1.PD)
        } else if (q[j] == 2) {
          out[[zz]][j,1] <- PDq2(ai1.v, ai2.v, Li.all.v, Tbar, m1, m2_all, n1, n2)
        } else {
          out[[zz]][j,1] <- PD_g_hat_fn(m1, n2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar)^(1/(1-q[j]))
        }

        out[[zz]][j,1] = out[[zz]][j,1]/Tbar
      }
    }

    output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
    rownames(output) = q
    colnames(output) = prop.v
    output = reshape2::melt(output)
    colnames(output) = c("Order.q","prop.v","qmiNEXT_PD")

    ## q0 ana =========================================
    if(nrow(m.inc.v) != 0){
      datanew = cbind(ai1.v,ai2.v)

      datash = datanew[(datanew[,1]>0 & datanew[,2]>0), , drop=F]
      dataun1 = datanew[(datanew[,1]>0 & datanew[,2]==0), , drop=F]
      dataun2 = datanew[(datanew[,1]==0 & datanew[,2]>0), , drop=F]


      un1 = sapply(m.inc.v[,1],function(i){
        sum(un_abun_PD(Li = Li.all.v[rownames(dataun1)],xi = dataun1[,1], n = n1, m = i))
      })
      un2 = sapply(m.inc.v[,2],function(i){
        sum(un_abun_PD(Li = Li.all.v[rownames(dataun2)],xi = dataun2[,2], n = n2, m = i))
      })
      sh12 = sapply(1:length(m.inc.v[,1]),function(i){
        sum(sh_abun_PD(Li = Li.all.v[rownames(datash)],xi1 = datash[,1],xi2 = datash[,2], n1 = n1, m1 =  m.inc.v[i,1],
                       n2 = n2, m2 = m.inc.v[i,2]))
      })

      un1 = un1/Tbar
      un2 = un2/Tbar
      sh12 = sh12/Tbar

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
        tree_bt = PDtree
        p_bt = bootstrap_population_multiple_assemblage(data,data_gamma, "abundance")
        unseen_p = p_bt[-(1:nrow(data)), ] %>% matrix(ncol = ncol(data))
        #  =========================================================================
        if (nrow(p_bt) > nrow(data) & sum(unseen_p) > 0) {
          unseen = unseen_p[which(rowSums(unseen_p) >
                                    0), ]
          unseen = matrix(unseen, ncol = ncol(unseen_p),
                          byrow = T)
          p_bt = rbind(p_bt[(1:nrow(data)), ], unseen)
          unseen_name = sapply(1:nrow(unseen), function(i) paste0("unseen_",i))
          rownames(p_bt) = c(rownames(data), unseen_name)
          x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1, size = sum(data[,k]), prob = p_bt[, k]))

          rownames(x_bt) = rownames(p_bt)
          if (sum(x_bt[-(1:nrow(data)), ]) > 0) {
            g0_hat = apply(data, 2, function(x) {
              n = sum(x)
              f1 = sum(x == 1)
              f2 = sum(x == 2)
              # remove0 = F 看要不要加
              aL = iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree, remove0 = F,
                                              data = x, rootExtend = T, refT = reft)
              aL$treeNabu$branch.length = aL$BLbyT[,
                                                   1]
              aL = aL$treeNabu %>% select(branch.abun,
                                          branch.length)
              g1 = aL$branch.length[aL$branch.abun ==
                                      1] %>% sum
              g2 = aL$branch.length[aL$branch.abun ==
                                      2] %>% sum
              g0_hat = ifelse(g2 > ((g1 * f2)/(2 *
                                                 f1)), ((n - 1)/n) * (g1^2/(2 * g2)),
                              ((n - 1)/n) * (g1 * (f1 - 1)/(2 * (f2 + 1))))
              g0_hat[is.na(g0_hat)] = 0
              g0_hat
            })
            # te表示新的bootstrap樣本中有看到但原本的樣本沒看到的物種
            te = (x_bt[1:nrow(data), ] * (data == 0)) > 0
            # used_length表示bootstrap中有看到但原本沒看到的branch length和
            used_length = sapply(1:ncol(data), function(i) {
              if (sum(te[, i]) == 0)
                return(0)
              else {
                # remove0 = F 看要不要加
                iNEXT.3D:::phyBranchAL_Abu(phylo = PDtree,
                                           data = x_bt[1:nrow(data), i],
                                           rootExtend = T, refT = reft)$treeNabu %>%
                  subset(label %in% names(which(te[, i] == TRUE))) %>% select(branch.length) %>%
                  sum
              }
            })
            g0_hat = g0_hat - used_length
            g0_hat[g0_hat < 0] = 0
            unseen_sample = x_bt[-(1:nrow(data)),]
            if (is.vector(unseen_sample))
              unseen_sample = matrix(unseen_sample,
                                     ncol = ncol(x_bt), byrow = T)
            #L0_hat為g0_hat平分給所有沒看到的unseen sample
            L0_hat = sapply(1:length(g0_hat), function(i) if (sum(unseen_sample[,i] > 0) > 0)
              (g0_hat[i]/nrow(unseen))
              else 0)
            L0_hat = rowSums((matrix(L0_hat, nrow(unseen_sample),
                                     ncol(unseen_sample), byrow = T) *
                                unseen_sample))/rowSums(unseen_sample)
            L0_hat[which(rowSums(unseen_sample) ==
                           0)] = 0
            for (i in 1:length(L0_hat)) {
              tip = list(edge = matrix(c(2, 1), 1, 2), tip.label = unseen_name[i],
                         edge.length = L0_hat[i], Nnode = 1)
              class(tip) = "phylo"
              tree_bt = tree_bt + tip
            }
          }
          else {
            x_bt = x_bt[1:nrow(data), ]
            p_bt = p_bt[1:nrow(data), ]
          }
        }else {
          p_bt = p_bt[1:nrow(data), ]
          x_bt = sapply(1:ncol(data), function(k) rmultinom(n = 1,
                                                            size = sum(data[, k]), prob = p_bt[, k]))
          rownames(x_bt) = rownames(data)
        }
        # ===============================================================================
        # p_bt
        # x_bt
        # tree_bt

        Liai.sample.data = mix.Liai.sample.PD(data = x_bt, oritree = tree_bt, datatype = "abundance",
                                              PDreftime = PDreftime, nT = NULL)

        ## new ===================================================================
        p_bind_bt = cbind(DetAbu(x_bt[,1], zero = T),DetAbu(x_bt[,2], zero = T))

        if(nrow(m.ext.v) != 0){
          Liai.ext = mix.Liai.popu.PD(data = p_bind_bt, oritree = tree_bt, datatype = "abundance",
                                      PDreftime = PDreftime, nT = NULL)

          ai1.ext.v = Liai.ext$ai1.v
          ai2.ext.v = Liai.ext$ai2.v

          Li.all.ext.v = Liai.ext$Li_all.v
        }

        # =================================================================
        # B <- length(Liai.sample.data$Li_all.v) ## 注意是B不是S
        qlength <- length(q)
        out = list()
        ai1.v = Liai.sample.data$ai1.v
        ai2.v = Liai.sample.data$ai2.v
        ai.all.v = Liai.sample.data$ai_all.v
        Li.all.v = Liai.sample.data$Li_all.v
        # Tbar要改，因為bootstrap的tree為non-ultrametric
        Tbar = sum(Li.all.v*ai.all.v)/(n1 + n2)

        for (zz in inc.idx) {


          out[[zz]] <- matrix(0, qlength, 1)
          m1 = m1.v[zz]
          m2 = m2.v[zz]

          for (j in 1:qlength) {
            if (q[j] == 0) {
              out[[zz]][j,1] <- PD_g_hat_fn(m1, m2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar)
            } else if (q[j] == 1) {
              out[[zz]][j,1] <- exp(PD_g_hat_fn(m1, m2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar))
            } else if (q[j] == 2) {
              out[[zz]][j,1] <- PDq2(ai1.v, ai2.v, Li.all.v, Tbar, m1, m2, n1, n2)
            } else {
              out[[zz]][j,1] <- PD_g_hat_fn(m1, m2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar)^(1/(1-q[j]))
            }

            out[[zz]][j,1] = out[[zz]][j,1]/Tbar
          }
        }

        for (zz in ext.idx) {
          out[[zz]] <- matrix(0, qlength, 1)
          m1 = m1.v[zz]
          m2_all = m2.v[zz]
          m2_s = m2.v[zz] - n2
          # 下面的m2改n2
          ## new =============================================
          h0.PD = h0_hat_cpp_PD(ai1.ext.v,ai2.ext.v,Li.all.ext.v,m1,m2_s,n1,n2)
          ## 思考h1_hat_FD的xi1, xi2要放啥
          h1.PD = h1_hat_PD(ai1.ext.v, ai2.ext.v, Li.all.ext.v, Tbar, ai1.ext.v, ai2.ext.v,
                            m1, m2_all, n1, n2)
          ## ===================================================
          for (j in 1:qlength) {
            if (q[j] == 0) {
              out[[zz]][j,1] <- PD_g_hat_fn(m1, n2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar) + h0.PD
            } else if (q[j] == 1) {
              out[[zz]][j,1] <- exp(PD_g_hat_fn(m1, n2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar) + h1.PD)
            } else if (q[j] == 2) {
              out[[zz]][j,1] <- PDq2(ai1.v, ai2.v, Li.all.v, Tbar, m1, m2_all, n1, n2)
            } else {
              out[[zz]][j,1] <- PD_g_hat_fn(m1, n2, n1, n2, Li.all.v, ai1.v, ai2.v, q[j], Tbar)^(1/(1-q[j]))
            }

            out[[zz]][j,1] = out[[zz]][j,1]/Tbar
          }
        }

        output = matrix(sapply(out, function(x) x[, 1]), ncol = length(m2.v))
        rownames(output) = q
        colnames(output) = prop.v
        output = reshape2::melt(output)
        colnames(output) = c("Order.q","prop.v","qmiNEXT_PD")

        res = output$qmiNEXT_PD
        ## q0 ana =========================================
        if(nrow(m.inc.v) != 0){
          datanew = cbind(ai1.v,ai2.v)

          datash = datanew[(datanew[,1]>0 & datanew[,2]>0), , drop=F]
          dataun1 = datanew[(datanew[,1]>0 & datanew[,2]==0), , drop=F]
          dataun2 = datanew[(datanew[,1]==0 & datanew[,2]>0), , drop=F]


          un1 = sapply(m.inc.v[,1],function(i){
            sum(un_abun_PD(Li = Li.all.v[rownames(dataun1)],xi = dataun1[,1], n = n1, m = i))
          })
          un2 = sapply(m.inc.v[,2],function(i){
            sum(un_abun_PD(Li = Li.all.v[rownames(dataun2)],xi = dataun2[,2], n = n2, m = i))
          })
          sh12 = sapply(1:length(m.inc.v[,1]),function(i){
            sum(sh_abun_PD(Li = Li.all.v[rownames(datash)],xi1 = datash[,1],xi2 = datash[,2], n1 = n1, m1 =  m.inc.v[i,1],
                           n2 = n2, m2 = m.inc.v[i,2]))
          })

          un1 = un1/Tbar
          un2 = un2/Tbar
          sh12 = sh12/Tbar

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
      output$LCL = output$qmiNEXT_PD - qtile*se_all
      output$UCL = output$qmiNEXT_PD + qtile*se_all

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
