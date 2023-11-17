###############################################################################L

#             Multi-state model specification for movement between refuges in 
#               western Tennessee and the LMAV by mallards



# This code is the model specification for the multi-state model to estimate 
# movement between refuges for mallards in western Tennessee. The project is 
# part of a mallard project at Tennessee Tech University in collaboration with 
# Tennessee Wildlife Resources Agency. 

# Dependencies: none



# Allison C Keever and Nick Masto
# Tennessee Tech University
# akeever@tntech.edu & nmasto1214@gmail.com
# GitHub: akeever2 & nmasto
# Created: 5/5/21
# Last updated: 8/30/23  # Change number of sanctuary nodes



###############################################################################L



# Constant Full model ----------------------------------------------------------
sink("fullconstant_MallardMSM.txt")
cat("
    model {
    
    #-----------------------------------------------------------------------
    # States and Observation of states (S; O):
    # 1 - alive at BB; seen at BB
    # 2 - alive at LINWR; seen at LINWR
    # 3 - alive at BS; seen at BS
    # 4 - alive at HI; seen at HI
    # 5 - alive at M; seen at M
    # 6 - alive at WL; seen at WL
    # 7 - alive at CNWR; seen at CNWR
    # 8 - alive at BLNWR; seen at BLNWR
    # 9 - alive at LL; seen at LL
    # 10 - alive at RLNWR_N; seen at RLNWR_N
    # 11 - alive at HB; seen at HB
    # 12 - alive at RLNWR_S; seen at RLNWR_S
    # 13 - alive at P; seen at P
    # 14 - dead or lost forever; not seen
    
    #-----------------------------------------------------------------------
    # Parameters
    # phi: survival probability for all sites. Because GPS transmitters can
    #       fail, this is more accurately interpreted as transmitter survival
    #       probability, which may or may not be related to actual the animals'
    #       survival. We assume phi is constant across sites. 
    #
    # p: recapture probability for all sites. Because we are using GPS transmitters 
    #     with potentially different fix rates and only using fixes on refuge, 
    #     this (p) is essentially the probability a fix occured on refuge between 
    #     10:00 and 19:00. We will also assume a constant probability
    #
    # psi.BB[1:14]:       probability of transition from BB to any other state 
    #                     that was observed in data
    # psi.LINWR[1:14]:    probability of transition from LINWR to any other state 
    #                     that was observed in data
    # psi.BS[1:14]:       probability of transition from BS to any other state 
    #                     that was observed in data
    # psi.HI[1:14]:       probability of transition from HI to any other state 
    #                     that was observed in data
    # psi.M[1:14]:        probability of transition from M to any other state 
    #                     that was observed in data
    # psi.WL[1:14]:       probability of transition from WL to any other state 
    #                     that was observed in data
    # psi.CNWR[1:14]:     probability of transition from CNWR to any other state 
    #                     that was observed in data
    # psi.BLNWR[1:14]:    probability of transition from BLNWR to any other state 
    #                     that was observed in data
    # psi.LL[1:14]:       probability of transition from LL to any other state 
    #                     that was observed in data
    # psi.RLNWR_N[1:14]:  probability of transition from RLNWR_N to any other state 
    #                     that was observed in data
    # psi.HB[1:14]:       probability of transition from HB to any other state 
    #                     that was observed in data
    # psi.RLNWR_S[1:14]:  probability of transition from RLNWR_S to any other state 
    #                     that was observed in data
    # psi.P[1:14]:        probability of transition from P to any other state 
    #                     that was observed in data
    #-----------------------------------------------------------------------
    
    
    
    ### Priors and constraints -----------
    
    # Prior for recapture probability / fix probability and survival
    # Uniform distribution from 0-1 assumes anywhere in parameter space
    # Is there a way to come up with a more informative prior? Probably but good for now.
    
    p ~ dunif(0, 1)
    phi ~ dunif(0, 1)
    
    
    # Priors for transition probabilities - multinomial logit. We will put normal
    # priors on the logit scale for an intercept and any covariates for the probabilty 
    # of transitioning to any other refuge. Then we will constrain the transitions 
    # so they sum to < 1, and calculate the probability of staying on refuge.
    
    # Priors for intercepts on logit scale for transitioning from refuge to any 
    # other refuge. It will be constant for each refuge and only vary based on 
    # covariates. If we were to add intercept for each refuge it would be too
    # system-specific so we are assuming that refuge transitions depend solely 
    # on island biogeo. and individual characteristics.
    
    b0 ~ dnorm(0, 0.04)
    
    # Priors for effect sizes
    b.dist ~ dnorm(0, 0.04)
    b.size1 ~ dnorm(0, 0.04)
    b.size2 ~ dnorm(0, 0.04)
    b.sex ~ dnorm(0, 0.04)
    b.age ~ dnorm(0, 0.04)
    
    
    # Transition probabilities
    for(i in 1:ninds){     # loop over no. of individuals (rows)
      for(s in 2:13){      # And possible transitions 2-13 possibilities
        # GLM on the logit scale for transition probabilities to any other refuge
        lpsi.BB[i,s] <- b0 + b.dist * Distance[1,s-1] + b.size1 * Size1[1] + 
                        b.size2 * Size2[s-1,1] + b.sex * Sex[i] + b.age * Age[i]
        lpsi.LINWR[i,s] <- b0 + b.dist * Distance[2,s-1] + b.size1 * Size1[2] + 
                        b.size2 * Size2[s-1,2] + b.sex * Sex[i] + b.age * Age[i]
        lpsi.BS[i,s] <- b0 + b.dist * Distance[3,s-1] + b.size1 * Size1[3] + 
                        b.size2 * Size2[s-1,3] + b.sex * Sex[i] + b.age * Age[i]
        lpsi.HI[i,s] <- b0 + b.dist * Distance[4,s-1] + b.size1 * Size1[4] + 
                        b.size2 * Size2[s-1,4] + b.sex * Sex[i] + b.age * Age[i]
        lpsi.M[i,s] <- b0 + b.dist * Distance[5,s-1] + b.size1 * Size1[5] + 
                        b.size2 * Size2[s-1,5] + b.sex * Sex[i] + b.age * Age[i]
        lpsi.WL[i,s] <- b0 + b.dist * Distance[6,s-1] + b.size1 * Size1[6] + 
                        b.size2 * Size2[s-1,6] + b.sex * Sex[i] + b.age * Age[i]
        lpsi.CNWR[i,s] <- b0 + b.dist * Distance[7,s-1] + b.size1 * Size1[7] + 
                        b.size2 * Size2[s-1,7] + b.sex * Sex[i] + b.age * Age[i]
        lpsi.BLNWR[i,s] <- b0 + b.dist * Distance[8,s-1] + b.size1 * Size1[8] + 
                        b.size2 * Size2[s-1,8] + b.sex * Sex[i] + b.age * Age[i]
        lpsi.LL[i,s] <- b0 + b.dist * Distance[9,s-1] + b.size1 * Size1[9] + 
                        b.size2 * Size2[s-1,9] + b.sex * Sex[i] + b.age * Age[i]
        lpsi.RLNWR_N[i,s] <- b0 + b.dist * Distance[10,s-1] + b.size1 * Size1[10] + 
                        b.size2 * Size2[s-1,10] + b.sex * Sex[i] + b.age * Age[i]
        lpsi.HB[i,s] <- b0 + b.dist * Distance[11,s-1] + b.size1 * Size1[11] + 
                        b.size2 * Size2[s-1,11] + b.sex * Sex[i] + b.age * Age[i]
        lpsi.RLNWR_S[i,s] <- b0 + b.dist * Distance[12,s-1] + b.size1 * Size1[12] + 
                        b.size2 * Size2[s-1,12] + b.sex * Sex[i] + b.age * Age[i]
        lpsi.P[i,s] <- b0 + b.dist * Distance[13,s-1] + b.size1 * Size1[13] + 
                        b.size2 * Size2[s-1,13] + b.sex * Sex[i] + b.age * Age[i]
        
        
        # Transform to probability scale and constrain probabilites < 1
        
        psi.BB[i,s] <- exp(lpsi.BB[i,s]) / (1 + exp(lpsi.BB[i,2]) + exp(lpsi.BB[i,3]) + 
          exp(lpsi.BB[i,4]) + exp(lpsi.BB[i,5]) + exp(lpsi.BB[i,6]) + exp(lpsi.BB[i,7]) + 
          exp(lpsi.BB[i,8]) + exp(lpsi.BB[i,9]) + exp(lpsi.BB[i,10]) + exp(lpsi.BB[i,11]) +
          exp(lpsi.BB[i,12]) + exp(lpsi.BB[i,13]))
        
        psi.LINWR[i,s] <- exp(lpsi.LINWR[i,s]) / (1 + exp(lpsi.LINWR[i,2]) + 
          exp(lpsi.LINWR[i,3]) + exp(lpsi.LINWR[i,4]) + exp(lpsi.LINWR[i,5]) + 
          exp(lpsi.LINWR[i,6]) + exp(lpsi.LINWR[i,7]) + exp(lpsi.LINWR[i,8]) + 
          exp(lpsi.LINWR[i,9]) + exp(lpsi.LINWR[i,10]) + exp(lpsi.LINWR[i,11]) +
          exp(lpsi.LINWR[i,12]) + exp(lpsi.LINWR[i,13]))
          
        psi.BS[i,s] <- exp(lpsi.BS[i,s]) / (1 + exp(lpsi.BS[i,2]) + exp(lpsi.BS[i,3]) + 
          exp(lpsi.BS[i,4]) + exp(lpsi.BS[i,5]) + exp(lpsi.BS[i,6]) + exp(lpsi.BS[i,7]) + 
          exp(lpsi.BS[i,8]) + exp(lpsi.BS[i,9]) + exp(lpsi.BS[i,10]) + exp(lpsi.BS[i,11]) +
          exp(lpsi.BS[i,12]) + exp(lpsi.BS[i,13]))
          
        psi.HI[i,s] <- exp(lpsi.HI[i,s]) / (1 + exp(lpsi.HI[i,2]) + exp(lpsi.HI[i,3]) + 
          exp(lpsi.HI[i,4]) + exp(lpsi.HI[i,5]) + exp(lpsi.HI[i,6]) + exp(lpsi.HI[i,7]) + 
          exp(lpsi.HI[i,8]) + exp(lpsi.HI[i,9]) + exp(lpsi.HI[i,10]) + exp(lpsi.HI[i,11]) +
          exp(lpsi.HI[i,12]) + exp(lpsi.HI[i,13]))
          
        psi.M[i,s] <- exp(lpsi.M[i,s]) / (1 + exp(lpsi.M[i,2]) + exp(lpsi.M[i,3]) + 
          exp(lpsi.M[i,4]) + exp(lpsi.M[i,5]) + exp(lpsi.M[i,6]) + exp(lpsi.M[i,7]) + 
          exp(lpsi.M[i,8]) + exp(lpsi.M[i,9]) + exp(lpsi.M[i,10]) + exp(lpsi.M[i,11]) +
          exp(lpsi.M[i,12]) + exp(lpsi.M[i,13]))
          
        psi.WL[i,s] <- exp(lpsi.WL[i,s]) / (1 + exp(lpsi.WL[i,2]) + exp(lpsi.WL[i,3]) + 
          exp(lpsi.WL[i,4]) + exp(lpsi.WL[i,5]) + exp(lpsi.WL[i,6]) + exp(lpsi.WL[i,7]) + 
          exp(lpsi.WL[i,8]) + exp(lpsi.WL[i,9]) + exp(lpsi.WL[i,10]) + exp(lpsi.WL[i,11]) +
          exp(lpsi.WL[i,12]) + exp(lpsi.WL[i,13]))
    
        psi.CNWR[i,s] <- exp(lpsi.CNWR[i,s]) / (1 + exp(lpsi.CNWR[i,2]) + exp(lpsi.CNWR[i,3]) + 
          exp(lpsi.CNWR[i,4]) + exp(lpsi.CNWR[i,5]) + exp(lpsi.CNWR[i,6]) + exp(lpsi.CNWR[i,7]) + 
          exp(lpsi.CNWR[i,8]) + exp(lpsi.CNWR[i,9]) + exp(lpsi.CNWR[i,10]) + exp(lpsi.CNWR[i,11]) +
          exp(lpsi.CNWR[i,12]) + exp(lpsi.CNWR[i,13]))
    
        psi.BLNWR[i,s] <- exp(lpsi.BLNWR[i,s]) / (1 + exp(lpsi.BLNWR[i,2]) + 
          exp(lpsi.BLNWR[i,3]) + exp(lpsi.BLNWR[i,4]) + exp(lpsi.BLNWR[i,5]) + 
          exp(lpsi.BLNWR[i,6]) + exp(lpsi.BLNWR[i,7]) + exp(lpsi.BLNWR[i,8]) + 
          exp(lpsi.BLNWR[i,9]) + exp(lpsi.BLNWR[i,10]) + exp(lpsi.BLNWR[i,11]) +
          exp(lpsi.BLNWR[i,12]) + exp(lpsi.BLNWR[i,13]))
          
        psi.LL[i,s] <- exp(lpsi.LL[i,s]) / (1 + exp(lpsi.LL[i,2]) + exp(lpsi.LL[i,3]) + 
          exp(lpsi.LL[i,4]) + exp(lpsi.LL[i,5]) + exp(lpsi.LL[i,6]) + exp(lpsi.LL[i,7]) + 
          exp(lpsi.LL[i,8]) + exp(lpsi.LL[i,9]) + exp(lpsi.LL[i,10]) + exp(lpsi.LL[i,11]) +
          exp(lpsi.LL[i,12]) + exp(lpsi.LL[i,13]))
          
        psi.RLNWR_N[i,s] <- exp(lpsi.RLNWR_N[i,s]) / (1 + exp(lpsi.RLNWR_N[i,2]) + exp(lpsi.RLNWR_N[i,3]) + 
          exp(lpsi.RLNWR_N[i,4]) + exp(lpsi.RLNWR_N[i,5]) + exp(lpsi.RLNWR_N[i,6]) + exp(lpsi.RLNWR_N[i,7]) + 
          exp(lpsi.RLNWR_N[i,8]) + exp(lpsi.RLNWR_N[i,9]) + exp(lpsi.RLNWR_N[i,10]) + exp(lpsi.RLNWR_N[i,11]) +
          exp(lpsi.RLNWR_N[i,12]) + exp(lpsi.RLNWR_N[i,13]))
          
        psi.HB[i,s] <- exp(lpsi.HB[i,s]) / (1 + exp(lpsi.HB[i,2]) + exp(lpsi.HB[i,3]) + 
          exp(lpsi.HB[i,4]) + exp(lpsi.HB[i,5]) + exp(lpsi.HB[i,6]) + exp(lpsi.HB[i,7]) + 
          exp(lpsi.HB[i,8]) + exp(lpsi.HB[i,9]) + exp(lpsi.HB[i,10]) + exp(lpsi.HB[i,11]) +
          exp(lpsi.HB[i,12]) + exp(lpsi.HB[i,13]))
    
        psi.RLNWR_S[i,s] <- exp(lpsi.RLNWR_S[i,s]) / (1 + exp(lpsi.RLNWR_S[i,2]) + exp(lpsi.RLNWR_S[i,3]) + 
          exp(lpsi.RLNWR_S[i,4]) + exp(lpsi.RLNWR_S[i,5]) + exp(lpsi.RLNWR_S[i,6]) + exp(lpsi.RLNWR_S[i,7]) + 
          exp(lpsi.RLNWR_S[i,8]) + exp(lpsi.RLNWR_S[i,9]) + exp(lpsi.RLNWR_S[i,10]) + exp(lpsi.RLNWR_S[i,11]) +
          exp(lpsi.RLNWR_S[i,12]) + exp(lpsi.RLNWR_S[i,13]))
    
        psi.P[i,s] <- exp(lpsi.P[i,s]) / (1 + exp(lpsi.P[i,2]) + exp(lpsi.P[i,3]) + 
          exp(lpsi.P[i,4]) + exp(lpsi.P[i,5]) + exp(lpsi.P[i,6]) + exp(lpsi.P[i,7]) + 
          exp(lpsi.P[i,8]) + exp(lpsi.P[i,9]) + exp(lpsi.P[i,10]) + exp(lpsi.P[i,11]) +
          exp(lpsi.P[i,12]) + exp(lpsi.P[i,13]))
          
      } # end s
          
      # Calculate probability of staying
      # Simply 1 - sum of transition probs
    
      psi.BB[i,1] <- 1 - sum(psi.BB[i,2:13])
      psi.LINWR[i,1] <- 1 - sum(psi.LINWR[i,2:13])
      psi.BS[i,1] <- 1 - sum(psi.BS[i,2:13])
      psi.HI[i,1] <- 1 - sum(psi.HI[i,2:13])
      psi.M[i,1] <- 1 - sum(psi.M[i,2:13])
      psi.WL[i,1] <- 1 - sum(psi.WL[i,2:13])
      psi.CNWR[i,1] <- 1 - sum(psi.CNWR[i,2:13])
      psi.BLNWR[i,1] <- 1 - sum(psi.BLNWR[i,2:13])
      psi.LL[i,1] <- 1 - sum(psi.LL[i,2:13])
      psi.RLNWR_N[i,1] <- 1 - sum(psi.RLNWR_N[i,2:13])
      psi.HB[i,1] <- 1 - sum(psi.HB[i,2:13])
      psi.RLNWR_S[i,1] <- 1 - sum(psi.RLNWR_S[i,2:13])
      psi.P[i,1] <- 1 - sum(psi.P[i,2:13])
      
    } # end i
    

    ### Define state-transition and observer matrices --------
    for(i in 1:ninds){
      for(t in f[i]:(noccs - 1)){
      
        # Movement (ps) from BB (1) to any other refuge and observation (po)
        ps[1,i,t,1] <- psi.BB[i,1] * phi
        ps[1,i,t,2] <- psi.BB[i,2] * phi
        ps[1,i,t,3] <- psi.BB[i,3] * phi
        ps[1,i,t,4] <- psi.BB[i,4] * phi
        ps[1,i,t,5] <- psi.BB[i,5] * phi
        ps[1,i,t,6] <- psi.BB[i,6] * phi
        ps[1,i,t,7] <- psi.BB[i,7] * phi
        ps[1,i,t,8] <- psi.BB[i,8] * phi
        ps[1,i,t,9] <- psi.BB[i,9] * phi
        ps[1,i,t,10] <- psi.BB[i,10] * phi
        ps[1,i,t,11] <- psi.BB[i,11] * phi
        ps[1,i,t,12] <- psi.BB[i,12] * phi
        ps[1,i,t,13] <- psi.BB[i,13] * phi
        ps[1,i,t,14] <- 1 - phi
        
        po[1,i,t,1] <- p
        po[1,i,t,2] <- 0
        po[1,i,t,3] <- 0
        po[1,i,t,4] <- 0
        po[1,i,t,5] <- 0
        po[1,i,t,6] <- 0
        po[1,i,t,7] <- 0
        po[1,i,t,8] <- 0
        po[1,i,t,9] <- 0
        po[1,i,t,10] <- 0
        po[1,i,t,11] <- 0
        po[1,i,t,12] <- 0
        po[1,i,t,13] <- 0
        po[1,i,t,14] <- 1 - p
        
        # Movement from LINWR (2) to any other refuge and observation (po)
        ps[2,i,t,1] <- psi.LINWR[i,2] * phi
        ps[2,i,t,2] <- psi.LINWR[i,1] * phi
        ps[2,i,t,3] <- psi.LINWR[i,3] * phi
        ps[2,i,t,4] <- psi.LINWR[i,4] * phi
        ps[2,i,t,5] <- psi.LINWR[i,5] * phi
        ps[2,i,t,6] <- psi.LINWR[i,6] * phi
        ps[2,i,t,7] <- psi.LINWR[i,7] * phi
        ps[2,i,t,8] <- psi.LINWR[i,8] * phi
        ps[2,i,t,9] <- psi.LINWR[i,9] * phi
        ps[2,i,t,10] <- psi.LINWR[i,10] * phi
        ps[2,i,t,11] <- psi.LINWR[i,11] * phi
        ps[2,i,t,12] <- psi.LINWR[i,12] * phi
        ps[2,i,t,13] <- psi.LINWR[i,13] * phi
        ps[2,i,t,14] <- 1 - phi
        
        po[2,i,t,1] <- 0
        po[2,i,t,2] <- p
        po[2,i,t,3] <- 0
        po[2,i,t,4] <- 0
        po[2,i,t,5] <- 0
        po[2,i,t,6] <- 0
        po[2,i,t,7] <- 0
        po[2,i,t,8] <- 0
        po[2,i,t,9] <- 0
        po[2,i,t,10] <- 0
        po[2,i,t,11] <- 0
        po[2,i,t,12] <- 0
        po[2,i,t,13] <- 0
        po[2,i,t,14] <- 1 - p
        
        # Movement from BS (3) to any other refuge and observation (po)
        ps[3,i,t,1] <- psi.BS[i,2] * phi
        ps[3,i,t,2] <- psi.BS[i,3] * phi
        ps[3,i,t,3] <- psi.BS[i,1] * phi
        ps[3,i,t,4] <- psi.BS[i,4] * phi
        ps[3,i,t,5] <- psi.BS[i,5] * phi
        ps[3,i,t,6] <- psi.BS[i,6] * phi
        ps[3,i,t,7] <- psi.BS[i,7] * phi
        ps[3,i,t,8] <- psi.BS[i,8] * phi
        ps[3,i,t,9] <- psi.BS[i,9] * phi
        ps[3,i,t,10] <- psi.BS[i,10] * phi
        ps[3,i,t,11] <- psi.BS[i,11] * phi
        ps[3,i,t,12] <- psi.BS[i,12] * phi
        ps[3,i,t,13] <- psi.BS[i,13] * phi
        ps[3,i,t,14] <- 1 - phi
        
        po[3,i,t,1] <- 0
        po[3,i,t,2] <- 0
        po[3,i,t,3] <- p
        po[3,i,t,4] <- 0
        po[3,i,t,5] <- 0
        po[3,i,t,6] <- 0
        po[3,i,t,7] <- 0
        po[3,i,t,8] <- 0
        po[3,i,t,9] <- 0
        po[3,i,t,10] <- 0
        po[3,i,t,11] <- 0
        po[3,i,t,12] <- 0
        po[3,i,t,13] <- 0
        po[3,i,t,14] <- 1 - p
        
        # Movement from HI (4) to any other refuge and observation (po)
        ps[4,i,t,1] <- psi.HI[i,2] * phi
        ps[4,i,t,2] <- psi.HI[i,3] * phi
        ps[4,i,t,3] <- psi.HI[i,4] * phi
        ps[4,i,t,4] <- psi.HI[i,1] * phi
        ps[4,i,t,5] <- psi.HI[i,5] * phi
        ps[4,i,t,6] <- psi.HI[i,6] * phi
        ps[4,i,t,7] <- psi.HI[i,7] * phi
        ps[4,i,t,8] <- psi.HI[i,8] * phi
        ps[4,i,t,9] <- psi.HI[i,9] * phi
        ps[4,i,t,10] <- psi.HI[i,10] * phi
        ps[4,i,t,11] <- psi.HI[i,11] * phi
        ps[4,i,t,12] <- psi.HI[i,12] * phi
        ps[4,i,t,13] <- psi.HI[i,13] * phi
        ps[4,i,t,14] <- 1 - phi
        
        po[4,i,t,1] <- 0
        po[4,i,t,2] <- 0
        po[4,i,t,3] <- 0
        po[4,i,t,4] <- p
        po[4,i,t,5] <- 0
        po[4,i,t,6] <- 0
        po[4,i,t,7] <- 0
        po[4,i,t,8] <- 0
        po[4,i,t,9] <- 0
        po[4,i,t,10] <- 0
        po[4,i,t,11] <- 0
        po[4,i,t,12] <- 0
        po[4,i,t,13] <- 0
        po[4,i,t,14] <- 1 - p
        
        # Movement from M (5) to any other refuge and observation (po)
        ps[5,i,t,1] <- psi.M[i,2] * phi
        ps[5,i,t,2] <- psi.M[i,3] * phi
        ps[5,i,t,3] <- psi.M[i,4] * phi
        ps[5,i,t,4] <- psi.M[i,5] * phi
        ps[5,i,t,5] <- psi.M[i,1] * phi
        ps[5,i,t,6] <- psi.M[i,6] * phi
        ps[5,i,t,7] <- psi.M[i,7] * phi
        ps[5,i,t,8] <- psi.M[i,8] * phi
        ps[5,i,t,9] <- psi.M[i,9] * phi
        ps[5,i,t,10] <- psi.M[i,10] * phi
        ps[5,i,t,11] <- psi.M[i,11] * phi
        ps[5,i,t,12] <- psi.M[i,12] * phi
        ps[5,i,t,13] <- psi.M[i,13] * phi
        ps[5,i,t,14] <- 1 - phi
        
        po[5,i,t,1] <- 0
        po[5,i,t,2] <- 0
        po[5,i,t,3] <- 0
        po[5,i,t,4] <- 0
        po[5,i,t,5] <- p
        po[5,i,t,6] <- 0
        po[5,i,t,7] <- 0
        po[5,i,t,8] <- 0
        po[5,i,t,9] <- 0
        po[5,i,t,10] <- 0
        po[5,i,t,11] <- 0
        po[5,i,t,12] <- 0
        po[5,i,t,13] <- 0
        po[5,i,t,14] <- 1 - p
        
        # Movement from WL (6) to any other refuge and observation (po)
        ps[6,i,t,1] <- psi.WL[i,2] * phi
        ps[6,i,t,2] <- psi.WL[i,3] * phi
        ps[6,i,t,3] <- psi.WL[i,4] * phi
        ps[6,i,t,4] <- psi.WL[i,5] * phi
        ps[6,i,t,5] <- psi.WL[i,6] * phi
        ps[6,i,t,6] <- psi.WL[i,1] * phi
        ps[6,i,t,7] <- psi.WL[i,7] * phi
        ps[6,i,t,8] <- psi.WL[i,8] * phi
        ps[6,i,t,9] <- psi.WL[i,9] * phi
        ps[6,i,t,10] <- psi.WL[i,10] * phi
        ps[6,i,t,11] <- psi.WL[i,11] * phi
        ps[6,i,t,12] <- psi.WL[i,12] * phi
        ps[6,i,t,13] <- psi.WL[i,13] * phi
        ps[6,i,t,14] <- 1 - phi
        
        po[6,i,t,1] <- 0
        po[6,i,t,2] <- 0
        po[6,i,t,3] <- 0
        po[6,i,t,4] <- 0
        po[6,i,t,5] <- 0
        po[6,i,t,6] <- p
        po[6,i,t,7] <- 0
        po[6,i,t,8] <- 0
        po[6,i,t,9] <- 0
        po[6,i,t,10] <- 0
        po[6,i,t,11] <- 0
        po[6,i,t,12] <- 0
        po[6,i,t,13] <- 0
        po[6,i,t,14] <- 1 - p
        
        # Movement from CNWR (7) to any other refuge and observation (po)
        ps[7,i,t,1] <- psi.CNWR[i,2] * phi
        ps[7,i,t,2] <- psi.CNWR[i,3] * phi
        ps[7,i,t,3] <- psi.CNWR[i,4] * phi
        ps[7,i,t,4] <- psi.CNWR[i,5] * phi
        ps[7,i,t,5] <- psi.CNWR[i,6] * phi
        ps[7,i,t,6] <- psi.CNWR[i,7] * phi
        ps[7,i,t,7] <- psi.CNWR[i,1] * phi
        ps[7,i,t,8] <- psi.CNWR[i,8] * phi
        ps[7,i,t,9] <- psi.CNWR[i,9] * phi
        ps[7,i,t,10] <- psi.CNWR[i,10] * phi
        ps[7,i,t,11] <- psi.CNWR[i,11] * phi
        ps[7,i,t,12] <- psi.CNWR[i,12] * phi
        ps[7,i,t,13] <- psi.CNWR[i,13] * phi
        ps[7,i,t,14] <- 1 - phi
        
        po[7,i,t,1] <- 0
        po[7,i,t,2] <- 0
        po[7,i,t,3] <- 0
        po[7,i,t,4] <- 0
        po[7,i,t,5] <- 0
        po[7,i,t,6] <- 0
        po[7,i,t,7] <- p
        po[7,i,t,8] <- 0
        po[7,i,t,9] <- 0
        po[7,i,t,10] <- 0
        po[7,i,t,11] <- 0
        po[7,i,t,12] <- 0
        po[7,i,t,13] <- 0
        po[7,i,t,14] <- 1 - p
        
        # Movement from BLNWR (8) to any other refuge and observation (po)
        ps[8,i,t,1] <- psi.BLNWR[i,2] * phi
        ps[8,i,t,2] <- psi.BLNWR[i,3] * phi
        ps[8,i,t,3] <- psi.BLNWR[i,4] * phi
        ps[8,i,t,4] <- psi.BLNWR[i,5] * phi
        ps[8,i,t,5] <- psi.BLNWR[i,6] * phi
        ps[8,i,t,6] <- psi.BLNWR[i,7] * phi
        ps[8,i,t,7] <- psi.BLNWR[i,8] * phi
        ps[8,i,t,8] <- psi.BLNWR[i,1] * phi
        ps[8,i,t,9] <- psi.BLNWR[i,9] * phi
        ps[8,i,t,10] <- psi.BLNWR[i,10] * phi
        ps[8,i,t,11] <- psi.BLNWR[i,11] * phi
        ps[8,i,t,12] <- psi.BLNWR[i,12] * phi
        ps[8,i,t,13] <- psi.BLNWR[i,13] * phi
        ps[8,i,t,14] <- 1 - phi
        
        po[8,i,t,1] <- 0
        po[8,i,t,2] <- 0
        po[8,i,t,3] <- 0
        po[8,i,t,4] <- 0
        po[8,i,t,5] <- 0
        po[8,i,t,6] <- 0
        po[8,i,t,7] <- 0
        po[8,i,t,8] <- p
        po[8,i,t,9] <- 0
        po[8,i,t,10] <- 0
        po[8,i,t,11] <- 0
        po[8,i,t,12] <- 0
        po[8,i,t,13] <- 0
        po[8,i,t,14] <- 1 - p
        
        # Movement from LL (9) to any other refuge and observation (po)
        ps[9,i,t,1] <- psi.LL[i,2] * phi
        ps[9,i,t,2] <- psi.LL[i,3] * phi
        ps[9,i,t,3] <- psi.LL[i,4] * phi
        ps[9,i,t,4] <- psi.LL[i,5] * phi
        ps[9,i,t,5] <- psi.LL[i,6] * phi
        ps[9,i,t,6] <- psi.LL[i,7] * phi
        ps[9,i,t,7] <- psi.LL[i,8] * phi
        ps[9,i,t,8] <- psi.LL[i,9] * phi
        ps[9,i,t,9] <- psi.LL[i,1] * phi
        ps[9,i,t,10] <- psi.LL[i,10] * phi
        ps[9,i,t,11] <- psi.LL[i,11] * phi
        ps[9,i,t,12] <- psi.LL[i,12] * phi
        ps[9,i,t,13] <- psi.LL[i,13] * phi
        ps[9,i,t,14] <- 1 - phi
        
        po[9,i,t,1] <- 0
        po[9,i,t,2] <- 0
        po[9,i,t,3] <- 0
        po[9,i,t,4] <- 0
        po[9,i,t,5] <- 0
        po[9,i,t,6] <- 0
        po[9,i,t,7] <- 0
        po[9,i,t,8] <- 0
        po[9,i,t,9] <- p
        po[9,i,t,10] <- 0
        po[9,i,t,11] <- 0
        po[9,i,t,12] <- 0
        po[9,i,t,13] <- 0
        po[9,i,t,14] <- 1 - p
        
        # Movement from RLNWR_N (10) to any other refuge and observation (po)
        ps[10,i,t,1] <- psi.RLNWR_N[i,2] * phi
        ps[10,i,t,2] <- psi.RLNWR_N[i,3] * phi
        ps[10,i,t,3] <- psi.RLNWR_N[i,4] * phi
        ps[10,i,t,4] <- psi.RLNWR_N[i,5] * phi
        ps[10,i,t,5] <- psi.RLNWR_N[i,6] * phi
        ps[10,i,t,6] <- psi.RLNWR_N[i,7] * phi
        ps[10,i,t,7] <- psi.RLNWR_N[i,8] * phi
        ps[10,i,t,8] <- psi.RLNWR_N[i,9] * phi
        ps[10,i,t,9] <- psi.RLNWR_N[i,10] * phi
        ps[10,i,t,10] <- psi.RLNWR_N[i,1] * phi
        ps[10,i,t,11] <- psi.RLNWR_N[i,11] * phi
        ps[10,i,t,12] <- psi.RLNWR_N[i,12] * phi
        ps[10,i,t,13] <- psi.RLNWR_N[i,13] * phi
        ps[10,i,t,14] <- 1 - phi
        
        po[10,i,t,1] <- 0
        po[10,i,t,2] <- 0
        po[10,i,t,3] <- 0
        po[10,i,t,4] <- 0
        po[10,i,t,5] <- 0
        po[10,i,t,6] <- 0
        po[10,i,t,7] <- 0
        po[10,i,t,8] <- 0
        po[10,i,t,9] <- 0
        po[10,i,t,10] <- p
        po[10,i,t,11] <- 0
        po[10,i,t,12] <- 0
        po[10,i,t,13] <- 0
        po[10,i,t,14] <- 1 - p
        
        # Movement from HB (11) to any other refuge and observation (po)
        ps[11,i,t,1] <- psi.HB[i,2] * phi
        ps[11,i,t,2] <- psi.HB[i,3] * phi
        ps[11,i,t,3] <- psi.HB[i,4] * phi
        ps[11,i,t,4] <- psi.HB[i,5] * phi
        ps[11,i,t,5] <- psi.HB[i,6] * phi
        ps[11,i,t,6] <- psi.HB[i,7] * phi
        ps[11,i,t,7] <- psi.HB[i,8] * phi
        ps[11,i,t,8] <- psi.HB[i,9] * phi
        ps[11,i,t,9] <- psi.HB[i,10] * phi
        ps[11,i,t,10] <- psi.HB[i,11] * phi
        ps[11,i,t,11] <- psi.HB[i,1] * phi
        ps[11,i,t,12] <- psi.HB[i,12] * phi
        ps[11,i,t,13] <- psi.HB[i,13] * phi
        ps[11,i,t,14] <- 1 - phi
        
        po[11,i,t,1] <- 0
        po[11,i,t,2] <- 0
        po[11,i,t,3] <- 0
        po[11,i,t,4] <- 0
        po[11,i,t,5] <- 0
        po[11,i,t,6] <- 0
        po[11,i,t,7] <- 0
        po[11,i,t,8] <- 0
        po[11,i,t,9] <- 0
        po[11,i,t,10] <- 0
        po[11,i,t,11] <- p
        po[11,i,t,12] <- 0
        po[11,i,t,13] <- 0
        po[11,i,t,14] <- 1 - p
    
        # Movement from RLNWR_S (12) to any other refuge and observation (po)
        ps[12,i,t,1] <- psi.RLNWR_S[i,2] * phi
        ps[12,i,t,2] <- psi.RLNWR_S[i,3] * phi
        ps[12,i,t,3] <- psi.RLNWR_S[i,4] * phi
        ps[12,i,t,4] <- psi.RLNWR_S[i,5] * phi
        ps[12,i,t,5] <- psi.RLNWR_S[i,6] * phi
        ps[12,i,t,6] <- psi.RLNWR_S[i,7] * phi
        ps[12,i,t,7] <- psi.RLNWR_S[i,8] * phi
        ps[12,i,t,8] <- psi.RLNWR_S[i,9] * phi
        ps[12,i,t,9] <- psi.RLNWR_S[i,10] * phi
        ps[12,i,t,10] <- psi.RLNWR_S[i,11] * phi
        ps[12,i,t,11] <- psi.RLNWR_S[i,12] * phi
        ps[12,i,t,12] <- psi.RLNWR_S[i,1] * phi
        ps[12,i,t,13] <- psi.RLNWR_S[i,13] * phi
        ps[12,i,t,14] <- 1 - phi
        
        po[12,i,t,1] <- 0
        po[12,i,t,2] <- 0
        po[12,i,t,3] <- 0
        po[12,i,t,4] <- 0
        po[12,i,t,5] <- 0
        po[12,i,t,6] <- 0
        po[12,i,t,7] <- 0
        po[12,i,t,8] <- 0
        po[12,i,t,9] <- 0
        po[12,i,t,10] <- 0
        po[12,i,t,11] <- 0
        po[12,i,t,12] <- p
        po[12,i,t,13] <- 0
        po[12,i,t,14] <- 1 - p
    
        # Movement from P (13) to any other refuge and observation (po)
        ps[13,i,t,1] <- psi.P[i,2] * phi
        ps[13,i,t,2] <- psi.P[i,3] * phi
        ps[13,i,t,3] <- psi.P[i,4] * phi
        ps[13,i,t,4] <- psi.P[i,5] * phi
        ps[13,i,t,5] <- psi.P[i,6] * phi
        ps[13,i,t,6] <- psi.P[i,7] * phi
        ps[13,i,t,7] <- psi.P[i,8] * phi
        ps[13,i,t,8] <- psi.P[i,9] * phi
        ps[13,i,t,9] <- psi.P[i,10] * phi
        ps[13,i,t,10] <- psi.P[i,11] * phi
        ps[13,i,t,11] <- psi.P[i,12] * phi
        ps[13,i,t,12] <- psi.P[i,13] * phi
        ps[13,i,t,13] <- psi.P[i,1] * phi
        ps[13,i,t,14] <- 1 - phi
        
        po[13,i,t,1] <- 0
        po[13,i,t,2] <- 0
        po[13,i,t,3] <- 0
        po[13,i,t,4] <- 0
        po[13,i,t,5] <- 0
        po[13,i,t,6] <- 0
        po[13,i,t,7] <- 0
        po[13,i,t,8] <- 0
        po[13,i,t,9] <- 0
        po[13,i,t,10] <- 0
        po[13,i,t,11] <- 0
        po[13,i,t,12] <- 0
        po[13,i,t,13] <- p
        po[13,i,t,14] <- 1 - p
    
    
        # Movement from dead/lost (14) to any other refuge and observation (po)
        ps[14,i,t,1] <- 0
        ps[14,i,t,2] <- 0
        ps[14,i,t,3] <- 0
        ps[14,i,t,4] <- 0
        ps[14,i,t,5] <- 0
        ps[14,i,t,6] <- 0
        ps[14,i,t,7] <- 0
        ps[14,i,t,8] <- 0
        ps[14,i,t,9] <- 0
        ps[14,i,t,10] <- 0
        ps[14,i,t,11] <- 0
        ps[14,i,t,12] <- 0
        ps[14,i,t,13] <- 0
        ps[14,i,t,14] <- 1
        
        po[14,i,t,1] <- 0
        po[14,i,t,2] <- 0
        po[14,i,t,3] <- 0
        po[14,i,t,4] <- 0
        po[14,i,t,5] <- 0
        po[14,i,t,6] <- 0
        po[14,i,t,7] <- 0
        po[14,i,t,8] <- 0
        po[14,i,t,9] <- 0
        po[14,i,t,10] <- 0
        po[14,i,t,11] <- 0
        po[14,i,t,12] <- 0
        po[14,i,t,13] <- 0
        po[14,i,t,14] <- 1
        
      } # end t
    } # end i

    
    
    ### Likelihood -----------
    for(i in 1:ninds){
      # Define latent state at first capture
      z[i,f[i]] <- init_state[i]
      
      for(t in (f[i]+1):noccs){
        
        # Ecological / state process: draw S(t) given S(t-1)
        z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
        
        # Observation process: draw O(t) given S(t)
        y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } # end t
    } # end i
    
    
    }", fill = TRUE)
sink()