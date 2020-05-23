#############################################################################################
######################Algorithm 2: ChordlessCycles(G)########################################
#############################################################################################
`ChordlessCycles`<-function(amat){#amat: adjacency matrix of an undirected graph
  counter<-0
  blocked<-c()
  Adj<-list()
  label<-c()
  ccycles<-list()
  ############################################################################################
  ###########Algorithm 3: DegreeLabeling(G)###################################################
  ############################################################################################
  `DegreeLabeling`<-function(amat){
    degree<-c()
    color<-c()
    valency<-c()
    # Adj<-list()
    # label<-c()
    n<-nrow(amat)
    for (i in 1:n) {
      neighbor<-c()
      degree[i]<-valency[i]<-0
      blocked[i]<<-0
      #color[i]=0 means color[i]=white#
      color[i]<-0
      for (j in 1:nrow(amat)) {
        if(amat[i,j]==1){
          degree[i]<-degree[i]+1
          valency[i]<-degree[i]
          neighbor<-c(neighbor,j)
        }
      }
      Adj[[i]]<<-neighbor
    }
    min_degree<-v<-0
    for (i in 1:n) {
      min_degree<-n
      for (j in 1:n) {
        if((color[j]==0) && (degree[j]<min_degree)){
          v<-j
          min_degree<-degree[j]
        }
      }
      label[v]<<-i
      #color[i]=1 means color[i]=black#
      color[v]<-1
      for (k in 1:n) {
        if((amat[v,k]==1) && (color[k]==0)){
          degree[k]<-degree[k]-1
        }
      }
    }
    return(list(label=label,valency=valency,Adj=Adj))
  }
  ##########################################################################################
  ################Algorithm 4: Triplets(G)##################################################
  ##########################################################################################
  `Triplets`<-function(amat){
    triplet<-list()
    n<-nrow(amat)
    valency<-DegreeLabeling(amat)$valency
    Adj<-DegreeLabeling(amat)$Adj
    label<-DegreeLabeling(amat)$label
    t<-1
    c<-1
    for (i in 1:n) {
      #Generate all triplets on form <x, u, y>.
      if(valency[i]>=2){
        sub2Adj<-combn(Adj[[i]],2)
        for (j in 1:ncol(sub2Adj)) {
          if(label[i]<label[sub2Adj[1,j]] && label[sub2Adj[1,j]]<label[sub2Adj[2,j]]){
            if(amat[sub2Adj[1,j],sub2Adj[2,j]]==1){
              ccycles[[c]]<<-c(sub2Adj[1,j],i,sub2Adj[2,j])
              c<-c+1
            }else{
              triplet[[t]]<-c(sub2Adj[1,j],i,sub2Adj[2,j])
              t<-t+1
            }
          }
          if(label[i]<label[sub2Adj[2,j]] && label[sub2Adj[2,j]]<label[sub2Adj[1,j]]){
            if(amat[sub2Adj[1,j],sub2Adj[2,j]]==1){
              ccycles[[c]]<<-c(sub2Adj[2,j],i,sub2Adj[1,j])
              c<-c+1
            }else{
              triplet[[t]]<-c(sub2Adj[2,j],i,sub2Adj[1,j])
              t<-t+1
            }
          }
        }
      } 
    }
    counter<<-length(ccycles)+1
    return(list(triplet=triplet,ccycles=ccycles))
  }
  
  ############################################################################################
  ################Algorithm 6: BlockNeighbors(v, blocked)#####################################
  ############################################################################################
  `BlockNeighbors`<-function(v){
    for (i in 1:length(Adj[[v]])) {
      blocked[Adj[[v]][i]]<<-blocked[Adj[[v]][i]]+1
    }
  }
  
  ############################################################################################
  ########################Algorithm 7: UnblockNeighbors(v, blocked)###########################
  ############################################################################################
  `UnblockNeighbors`<-function(v){
    for (i in 1:length(Adj[[v]])){
      if(blocked[Adj[[v]][i]]>0){
        blocked[Adj[[v]][i]]<<-blocked[Adj[[v]][i]]-1
      }
    }
  }
  
  ############################################################################################
  #######################Algorithm 5: CC_Visit(p, C, key, blocked)############################
  ############################################################################################
  `CC_Visit`<-function(cpath,ccycles){
    key<-label[cpath[2]]
    BlockNeighbors(cpath[length(cpath)])
    for (i in Adj[[cpath[length(cpath)]]]) {
      if(label[i]>key && blocked[i]==1){
        cpath1<-c(cpath,i)
        if(cpath[1] %in% Adj[[i]]){
          ccycles[[counter]]<<-cpath1
          counter<<-counter+1
        }else{
          CC_Visit(cpath1,ccycles)
        }
      }
    }
    UnblockNeighbors(cpath[length(cpath)])
  }
  ###################################################Algorithm 2: continue...#################
  G<-DegreeLabeling(amat)
  triplet<-Triplets(amat)$triplet
  t<-1
  while (length(triplet)>0) {
    cpath<-triplet[[t]]
    BlockNeighbors(cpath[2])
    ccycles[counter]<-CC_Visit(cpath,ccycles)
    UnblockNeighbors(cpath[2])
    triplet[[t]]<-NULL
  }
  return(ccycles)
}
################################################################################################
################################################################################################
################################################################################################
`studeny_rules`<-function(amat){
  ###A function that computes the parents of given node v
  parents<-function(v,amat){
    pv<-c()
    for (i in 1:ncol(amat)) {
      if(amat[i,v]==1 && amat[v,i]==0){
        pv<-c(pv,i)
      }
    }
    return(pv)
  }
  ###A function that computes the children of given node v
  children<-function(v,amat){
    pv<-c()
    for (i in 1:ncol(amat)) {
      if(amat[i,v]==0 && amat[v,i]==1){
        pv<-c(pv,i)
      }
    }
    return(pv)
  }
  ###A function that computes the neighbors of given node v
  neighbours<-function(v,amat){
    pv<-c()
    for (i in 1:ncol(amat)) {
      if(amat[i,v]==1 && amat[v,i]==1){
        pv<-c(pv,i)
      }
    }
    return(pv)
  }
  ##########################################################################################################
  ###A function that determines the type of chordless cycles: undirected, directed, or partially directed###
  ##########################################################################################################
  ccycles_type<-function(vect,amat){
    cu<-cd<-co<-0
    l<-length(vect)
    for (i in 1:(l-1)) {
      if(amat[vect[i],vect[i+1]]==1 &&  amat[vect[i+1],vect[i]]==1){
        cu<-cu+1
      }
      if(amat[vect[i],vect[i+1]]==1 &&  amat[vect[i+1],vect[i]]==0){
        cd<-cd+1
      }
      if(amat[vect[i],vect[i+1]]==0 &&  amat[vect[i+1],vect[i]]==1){
        co<-co+1
      }
    }
    if(amat[vect[l],vect[1]]==1 &&  amat[vect[1],vect[l]]==1){
      cu<-cu+1
    }
    if(amat[vect[l],vect[1]]==1 &&  amat[vect[1],vect[l]]==0){
      cd<-cd+1
    }
    if(amat[vect[l],vect[1]]==0 &&  amat[vect[1],vect[l]]==1){
      co<-co+1
    }
    ##message("cu= ",cu," cd= ",cd," co= ", co)
    if(cu==l){
      #undirected cycle
      return(0)
    }else if(cd==1 && co==1){
      m<-amat[vect,vect]
      b<-which(m-t(m)>0,arr.ind = TRUE)
      if(b[1,1]==b[2,1]){
        #not a partially directed cycle with one colider -><-
        return(3)
      }else{
        #not a partially directed cycle not include a colider
        return(4)
      }
    }else if((cd==0 && co==2) || (co==0 && cd==2)){
      m<-amat[vect,vect]
      b<-which(m-t(m)>0,arr.ind = TRUE)
      if(b[1,2]==b[2,1] || b[1,1]==b[2,2]){
        #partially directed cycle with 2 sequential directed edge ->->
        return(2) 
      }else{
        #partially directed cycle with 2 non-sequential directed edge
        return(4)
      }
    }else if((cd==0 && co==1) || (co==0 && cd==1)){
      #partially directed cycle with exactly 1 directed edge
      return(1)
    }else{
      #others
      return(4)
    }
  }
  ###################comute the skeleton of the pattern############
  `skelet` <- function(amat)
  {
    0 + (amat + t(amat) > 0)
  }
  #################################################################
  ###################necessity rule for ccycles of length > 3######
  #################################################################
  necessity_rule<-function(lst,vec,amat){
    amat<-amat
    counter<-length(lst)
    while (counter>0) {
      l<-lst[[1]]
      #print(l)
      m<-amat[l,l]
      b<-which(m-t(m)>0,arr.ind = TRUE)
      if(vec[1]==1){
        cs<-parents(l[b[1,2]],amat)
        s<-cs[which(cs %in% l)]
        #print(s)
        cr<-neighbours(s,amat)
        r<-cr[which(cr %in% l)]
        #print(r)
        amat[r,s]<-0
      }else{
        if(b[1,1]==b[2,2]){
          cs<-parents(l[b[1,1]],amat)
          s<-cs[which(cs %in% l)]
          #print(s)
          cr<-neighbours(s,amat)
          r<-cr[which(cr %in% l)]
          #print(r)
          amat[r,s]<-0
        }else{
          cs<-parents(l[b[1,2]],amat)
          s<-cs[which(cs %in% l)]
          #print(s)
          cr<-neighbours(s,amat)
          r<-cr[which(cr %in% l)]
          #print(r)
          amat[r,s]<-0
        }
      }
      lst<-lst[-1]
      vec<-vec[-1]
      if(length(lst)>0){
        for (i in 1:length(lst)) {
          vec[i]<-ccycles_type(lst[[i]],amat)
        }
        ctemp<-intersect(which(vec!=1),which(vec!=2))
        if(length(ctemp)>0){
          vec<-vec[-ctemp]
          lst<-lst[-ctemp]
        } 
      }
      #message("*********")
      #print(amat)
      #print(vec)
      #print(lst)
      #message("*********")
      counter<-length(lst)
    }
    return(amat)
  }
  #################################################################
  ###################Find Double-cycle rule candidate##############
  #################################################################
  find_double_cycle_candidate<-function(r,a,t,lst,amat){
    for (i in 1:length(lst)) {
      if((r %in% lst[[i]]) && (t %in% lst[[i]])){
        #adjacents of vertex r in amat
        radjacents<-union(union(parents(r,amat),children(r,amat)),neighbours(r,amat))
        cs<-setdiff(radjacents, t)
        s<-cs[which(cs %in% lst[[i]])]
        schild<-c()
        schildren<-children(s,amat)
        if(length(schildren)>0){
          schild<-schildren[which(schildren %in% lst[[i]])] 
        }
        if((length(s)==1) && (s!=a) && (length(schild)==1)){
          return(1)
        }
      }
    }
    return(0)
  }
  
  #################################################################
  #finding chordless cycles of the learnt pattern
  umat<-skelet(amat)
  ccycles<-ChordlessCycles(umat)
  #################################################################
  possible<-TRUE
  while(possible){
    possible<-FALSE
    vectype<-c()
    #for ccycles of size >3 & of type 1 or 2
    ctemp<-list()
    vtemp<-c()
    vect<-c()
    #for ccycles of type 0
    uv<-c()
    uvec<-c()
    
    feasible<-TRUE
    while(feasible){
      feasible<-FALSE
      lcc<-length(ccycles)
      if(lcc>0){
        for (i in 1:lcc) {
          vectype[i]<-ccycles_type(ccycles[[i]],amat)
          if((vectype[i]==1 || vectype[i]==2) && (length(ccycles[[i]])>3)){
            vtemp<-c(vtemp,i)
            vect<-c(vect,vectype[i])
          }
          if(vectype[i]==0){
            uvec<-c(uvec,i)
          }
        }
      }
      
      ###if length(vtemp)>0, apply necessity rule
      if(length(vtemp)>0){
        ctemp<-ccycles[vtemp]
        ###print(ctemp)
        utemp<-ccycles[uvec]
        ###print(utemp)
        amat<-necessity_rule(ctemp,vect,amat)
        # vect<-NULL
        # ctemp<-NULL
        #update undirected ccycles if needed
        if(length(uvec)>0){
          for (i in 1:length(uvec)) {
            uv[i]<-ccycles_type(utemp[[i]],amat)
          }
          temp0<-union(which(uv==1),which(uv==2))
          ##print(uv)
          if(length(temp0)>0){
            # for (j in 1:length(temp0)) {
            #   vect[j]<-uv[temp0[j]]
            #   ctemp[[j]]<-utemp[[temp0[j]]]
            # }
            feasible<-TRUE
          }
        }
      }
    }
    
    
    ###prepare for Double-cycle rule
    ltype0<-list()
    vtype0<-c()
    ltype1<-list()
    vtype1<-c()
    ltype2<-list()
    vtype2<-c()
    ltype3<-list()
    vtype3<-c()
    if(lcc>0){
      for (i in 1:lcc) {
        vectype[i]<-ccycles_type(ccycles[[i]],amat)
        if(vectype[i]==1){
          vtype1<-c(vtype1,i)
        }
        if(vectype[i]==2){
          vtype2<-c(vtype2,i)
        }
        if(vectype[i]==3){
          vtype3<-c(vtype3,i)
        }
        if(vectype[i]==0){
          vtype0<-c(vtype0,i)
        }
      }
    }
    
    ltype0<-ccycles[vtype0]
    ltype1<-ccycles[vtype1]
    ltype2<-ccycles[vtype2]
    ltype3<-ccycles[vtype3]
    ltype123<-union(union(ltype1,ltype2),ltype3)
    ###if length(ltype1)>0, try Double-cycle rule
    counter<-length(ltype1)
    while (counter>0) {
      l<-ltype1[[1]]
      ##message("l=",l)
      m<-amat[l,l]
      ##print(m)
      b<-which(m-t(m)>0,arr.ind = TRUE)
      ##print(b)
      ca<-parents(l[b[1,2]],amat)
      ##message("ca=",ca)
      a<-ca[which(ca %in% l)]
      ##message("a=",a)
      cr<-neighbours(a,amat)
      r<-cr[which(cr %in% l)]
      ##message("r=",r)
      ct<-setdiff(neighbours(r,amat), a)
      t<-ct[which(ct %in% l)]
      #if find double_cycle candidate, apply Double-cycle rule
      if(find_double_cycle_candidate(r,a,t,ltype123,amat)==1){
        amat[t,r]<-0
        ###update cycles' type if needed###
        ltype1<-ltype1[-1]
        vtype1<-vtype1[-1]
        if(length(ltype1)>0){
          for (i in 1:(counter-1)) {
            vtype1[i]<-ccycles_type(ltype1[[i]],amat)
          }
          temp<-which(vtype1!=1)
          if(length(temp)>0){
            vtype1<-vtype1[-temp]
            ltype1<-ltype1[-temp]
          } 
        }
        counter<-length(ltype1)
      }else{
        ltype1<-ltype1[-1]
        vtype1<-vtype1[-1]
        counter<-length(ltype1)
      }
    }
    ###number of possible ccycles of type 0 converted to type 1 or 2
    temp0<-0
    if(length(ltype0)>0){
      for (i in 1:length(ltype0)) {
        typ<-ccycles_type(ltype0[[i]],amat)
        if(typ==1 || typ==2){
          temp0<-1
          break;
        }
      }
    }
    ### if possible==TRUE, repeat the procedure
    if(temp0>0){
      possible<-TRUE
    }
  }
  return(list(matrix=amat,num_ccycles=length(ccycles)))
}
#################################################################################################################
#################################################################################################################
#################################################################################################################
`learn.lwf.norm`<-function(data,p.value,method ="stable",LCG=FALSE,...){
  #### First, load the package pcalg and the data set. ####
  V <- colnames(data)
  covariance<-cov(data)
  #print(covariance)
  n<-nrow(data)
  suffStat<-list(C=cor(data),n=n)
  mm <- switch(method,
               stable = 1,
               original = 2,
               stable.fast = 3,
               0)
  if (mm == 0) stop("Invalid method!")
  if(mm==1){
    skel<-skeleton(suffStat,
                   indepTest = gaussCItest, ## (partial correlations)
                   labels =V,alpha=p.value, verbose = FALSE,...)
  }
  if(mm==2){
    skel<-skeleton(suffStat,labels =V,method="original",
                   indepTest = gaussCItest, ## (partial correlations)
                   alpha=p.value, verbose = FALSE,...)
  }
  if(mm==3){
    skel<-skeleton(suffStat,labels =V,method="stable.fast",
                   indepTest = gaussCItest, ## (partial correlations)
                   alpha=p.value, verbose = FALSE,...)
  }
  wmat<-as(skel@graph,"matrix")####gives adjacency matrix of pcalg
  #wmat<-wmat[nrow(wmat):1,ncol(wmat):1]
  # rownames(wmat)<-V
  # colnames(wmat)<-V
  vset <- rownames(wmat)
  #print(wmat)
  #draw(wmat)
  zmat<-wmat
  p <- nrow(wmat)
  if(p>2){
    ##################################################### 
    ##############Rule0(R0): Complex Recovery
    #####################################################  
    sep.pairs <-skel@sepset
    #n.sep <- length(sep.pairs)
    #if (n.sep == 0) return(wmat)
    for (i in 1:(p-1)) {#message("i: ",i)
      for (k in (i+1):p) {#message("k: ",k)
        for (turn in 1:2) {
          u <- if(turn == 1) i else k
          v <- if(turn == 1) k else i
          #message("u: ",u," ",class(u)," v: ",v)
          sep <- sep.pairs[[i]][[k]]
          if(!is.null(sep)){
            #message("sep: ",sep," ",class(sep))
            nb.u <- which(zmat[u,] == 1)
            #message("nb.u: ",nb.u," ",class(nb.u))
            nb.u.size <- length(nb.u)
            if (nb.u.size > 0) {
              for (j in 1:nb.u.size) {
                w <- nb.u[j]
                #message("w: ",w," ",class(w))
                newsep <- unique(append(V[sep], V[w]))
                #message("newsep",newsep," ", class(newsep))
                idx <- c(V[u], V[v], newsep)
                #message("idx",idx," ",class(idx))
                #print(idx)
                res <- norm.ci.test(covariance[idx, idx], n, V[u], V[v], newsep)
                if (res$p.value < p.value &&
                    (-1 - res$deviance) < wmat[w,u]) {
                  wmat[w, u] <-  -1 - res$deviance
                }
              }
            }
          }
        }
      }
    }
    idx <- which(wmat - t(wmat) < 0)
    wmat[idx] <- 0
    wmat <- 0 + (wmat != 0)
    # rownames(wmat)<-V
    # colnames(wmat)<-V
    wmat<-pattern(wmat)
    
    #################################################################
    ######## Converting the learnt pattern to a LCG if needed########
    #################################################################
    
    if(LCG){
      #if(!is.chaingraph(wmat)){
      
      #########################################################################################
      ###########Then apply the rules: Transitivity, NECESSITY, and Double-cycle rule##########
      #######Here we use LEMMA 4.5 of 4.2. Convergence of the Algorithm (Studeny 1997)#########
      #########################################################################################
      wmat<-studeny_rules(wmat)$matrix
      num_ccycles<-studeny_rules(wmat)$num_ccycles
      #message("num_ccycles=",num_ccycles)
      #} 
    }
  }
  return(list(matrix=wmat,num_ccycles=num_ccycles))
  #return(wmat)
}
#######################################################################################################
#######################################################################################################
#######################################################################################################
`learn.lwf.multinom`<-function(data,suffStat,p.value,Test,LCG=FALSE,...){
  #### First, load the package pcalg and the data set. ####
  V <- colnames(data)
  freq.table<-as.freq.tb(data.matrix(data, rownames.force = NA))
  #print(covariance)
  #n<-nrow(data)
  # suffStat<-list(C=cor(data),n=n)
  # mm <- switch(method,
  #              stable = 1,
  #              original = 2,
  #              stable.fast = 3,
  #              0)
  # if (mm == 0) stop("Invalid method!")
  # if(mm==1){
  #   skel<-skeleton(suffStat,
  #                  indepTest = gaussCItest, ## (partial correlations)
  #                  labels =V,alpha=p.value, verbose = FALSE,...)
  # }
  # if(mm==2){
  #   skel<-skeleton(suffStat,labels =V,method="original",
  #                  indepTest = gaussCItest, ## (partial correlations)
  #                  alpha=p.value, verbose = FALSE,...)
  # }
  # if(mm==3){
  #   skel<-skeleton(suffStat,labels =V,method="stable.fast",
  #                  indepTest = gaussCItest, ## (partial correlations)
  #                  alpha=p.value, verbose = FALSE,...)
  # }
  skel<-pcalg::skeleton(suffStat,
                        indepTest=Test, ## G square Test for (Conditional) Independence of Binary Variables
                        alpha=p.value,labels=V, verbose = FALSE)
  wmat<-as(skel@graph,"matrix")####gives adjacency matrix of pcalg
  #wmat<-wmat[nrow(wmat):1,ncol(wmat):1]
  # rownames(wmat)<-V
  # colnames(wmat)<-V
  #vset <- rownames(wmat)
  #print(wmat)
  #draw(wmat)
  zmat<- wmat
  p <- nrow(wmat)
  if(p>2){
    ##################################################### 
    ##############Rule0(R0): Complex Recovery
    #####################################################  
    sep.pairs <-skel@sepset
    #n.sep <- length(sep.pairs)
    #if (n.sep == 0) return(wmat)
    for (i in 1:(p-1)) {#message("i: ",i)
      for (k in (i+1):p) {#message("k: ",k)
        for (turn in 1:2) {
          u <- if(turn == 1) i else k
          v <- if(turn == 1) k else i
          #message("u: ",u," ",class(u)," v: ",v)
          sep <- sep.pairs[[i]][[k]]
          if(!is.null(sep)){
            #message("sep: ",sep," ",class(sep))
            nb.u <- which(zmat[u,] == 1)
            #message("nb.u: ",nb.u," ",class(nb.u))
            nb.u.size <- length(nb.u)
            if (nb.u.size > 0) {
              for (j in 1:nb.u.size) {
                w <- nb.u[j]
                #message("w: ",w," ",class(w))
                newsep <- unique(append(V[sep], V[w]))
                #message("newsep",newsep," ", class(newsep))
                idx <- c(V[u], V[v], newsep)
                #message("idx",idx," ",class(idx))
                #print(idx)
                res <- multinom.ci.test(freq.table, V[u], V[v], newsep)
                if (res$p.value < p.value &&
                    (-1 - res$deviance) < wmat[w,u]) {
                  wmat[w, u] <-  -1 - res$deviance
                }
              }
            }
          }
        }
      }
    }
    idx <- which(wmat - t(wmat) < 0)
    wmat[idx] <- 0
    wmat <- 0 + (wmat != 0)
    # rownames(wmat)<-V
    # colnames(wmat)<-V
    wmat<-pattern(wmat)
    
    
    #################################################################
    ######## Converting the learnt pattern to a LCG if needed########
    #################################################################
    
    if(LCG){
      #if(!is.chaingraph(wmat)){
      
      #########################################################################################
      ###########Then apply the rules: NECESSITY, and Double-cycle rule##########
      #######Here we use LEMMA 4.5 of 4.2. Convergence of the Algorithm (Studeny 1997)#########
      #########################################################################################
      wmat<-studeny_rules(wmat)$matrix
      num_ccycles<-studeny_rules(wmat)$num_ccycles
      #message("num_ccycles=",num_ccycles)
      #} 
    }
  }
  return(list(matrix=wmat,num_ccycles=num_ccycles))
  #return(wmat)
}
#############################################################################################################
#############################################################################################################
#############################################################################################################
`comp.cgs`<-function(truecg,lcg){
  vset <- rownames(truecg)
  ########Skeleton comparison#############
  `skelet` <- function(amat)
  {
    0 + (amat + t(amat) > 0)
  }
  
  skel<-skelet(lcg)
  trueskel<-skelet(truecg)
  a <-  trueskel- skel[vset, vset]
  e.missing <- length(which(a == 1))/2
  e.extra <- length(which(a == -1))/2
  e.total = length(which(trueskel == 1))/2
  t.total<-length(which(skel == 1))/2
  #######Pattern comparison##############
  truepat<-pattern(truecg)
  pat<-pattern(lcg)
  truearr <- which(truepat - t(truepat) == 1)
  pat <- pat[vset, vset]
  arr <- which(pat - t(pat) == 1)
  a.missing <- length(truearr)-length(which(match(truearr, arr)>0))
  a.extra <- length(arr)-length(which(match(arr, truearr)>0))
  a.total = length(truearr)
  #########Structural Hamming distance#################
  shd<-0
  target<-lcg[vset,vset]
  for (i in 1:(length(vset)-1)) {
    for (j in (i+1):length(vset)) {
      if((target[i,j]!=truecg[i,j]) || (target[j,i]!=truecg[j,i])){
        shd<-shd+1
      }
    }
  }
  ###########Error measures###############
  tp<-t.total-e.extra
  N=choose(length(vset),2)-e.total
  tn<-N-e.extra
  tpr<-tp/e.total
  fpr<-e.extra/N
  acc<-(tp+tn)/(e.total+N)
  #######################################
  return(list(TPR=tpr,
              FPR=fpr,
              ACC=acc,
              a.total=a.total,
              a.missing=a.missing,
              a.extra=a.extra,
              SHD=shd))
}
###################################################################################################
###################################################################################################
###################################################################################################