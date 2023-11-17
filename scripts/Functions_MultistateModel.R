###############################################################################L

#   Functions for multi-state model to estimate movement between refuges in 
#   western Tennessee and the LMAV by mallards



# This code  specifies a variety of custom functions necessary to run the analysis
# for the multi-state mallard movement transition model to estimate movement 
# probabilities among waterfowl refuges for mallards in western Tennessee. 
# The project is part of a mallard project at Tennessee Tech University 
# in collaboration with Tennessee Wildlife Resources Agency. 

# Dependencies: none

# Allison C Keever
# Tennessee Tech University
# akeever@tntech.edu
# GitHub: akeever2
# Created: 5/8/21
# Last updated: 5/12/21

###############################################################################L


# Function to create known latent states
known_state <- function(ms, notseen){
  state <- ms
  state[state == notseen] <- NA
  for(i in 1:nrow(ms)){
    m <- min(which(!is.na(state[i,])))
    state[i,m] <- NA
  }
  return(as.matrix(state))
}

# Function to create good starting values for unknown z
init_z <- function(ms, f){           
  states <- max(ms, na.rm = TRUE)
  bs <- ms
  
  for(i in 1:nrow(ms)){
    if(f[i] != ncol(ms)){       
      for(t in (f[i]):ncol(ms)){  # In essense, this function just tells
        if(ms[i,t] == states){    # model "instead of searching param space",
          ms[i,t] <- ms[i,t-1]    # use last refuge location if GPS fix
        }                         # is not there; provides good and likely start 
      }                           # value based on previous occasion thereby 
                                  # speeding convergence.
      ms[i,1:f[i]] <- NA_integer_
      bs[i,1:f[i]] <- NA_integer_
      
    } else {
      ms[i,] <- NA_integer_
    }
  }
  
  for(i in 1:nrow(ms)){
    for(t in 1:ncol(ms)){
      ms[i,t] <- ifelse(bs[i,t] == states, ms[i,t], NA_integer_)
    }
  }
  
  return(as.matrix(ms))
}

# Get initial state for each individual
first_state <- function(ms, f){
  init <- NULL
  for(i in 1:nrow(ms)){
    x <- as.numeric(ms[i, f[i]]) # self-explanatory
    init <- c(init, x)
  }
  return(init)
}


# Function to calculate probabilities for our refuges based on results
calc_refuge_probs <- function(from.ref = NULL, to.ref = NULL, dist_datum = dist_datum, 
                              Refs = Refs){ 
  D <- 1 - diag(13) # diag number of rows of Refs
  #D <- 1 - diag(nrows(Refs)) # diag number of rows of Refs
  D[D == 1] <- t(dist_datum) # transpose distance df
  
  refuge_datum <- data.frame(From = Refs$Name, Size1 = Refs$Area_sqkm) %>%
    expand_grid(., 
                data.frame(To = Refs$Name, Size2 = Refs$Area_sqkm)) %>%
    mutate(dist = as.vector(D)) %>% 
    expand_grid(., 
                expand_grid(age = c(0, 1), 
                            sex = c(0, 1)))
  
  
  if(!is.null(from.ref) & is.null(to.ref)){ # b0 b/c we decided constant probs very general application
    x <- out %>% spread_draws(b.dist, b.size1, b.size2, b.sex, b.age, b0) %>% # 
      expand_grid(., 
                  refuge_datum %>%
                    filter(From == from.ref, To != from.ref)) %>%
      mutate(#b0 = case_when(
      #   From == "BB" ~ b0.BB, 
      #   From == "LINWR" ~ b0.LINWR,
      #   From == "BS" ~ b0.BS, 
      #   From == "HI" ~ b0.HI, 
      #   From == "M" ~ b0.M, 
      #   From == "WL" ~ b0.WL,
      #   From == "CNWR" ~ b0.CNWR,
      #   From == "BLNWR" ~ b0.BLNWR,
      #   From == "LL" ~ b0.LL,
      #   From == "RLNWR_N" ~ b0.RLNWR_N,    # Fixed refs
      #   From == "HB" ~ b0.HB,
      #   From == "RLNWR_S" ~ b0.RLNWR_S,    # 2
      #   From == "P" ~ b0.P),               # 3
        lp = exp(b0 + b.dist * dist + b.size1 * Size1 + b.size2 * Size2 + 
                   b.sex * sex + b.age * age)) %>%
      group_by(age, sex, .draw) %>%
      mutate(psi = lp / (1 + sum(lp))) %>%
      bind_rows(expand_grid(From = from.ref, To = from.ref, age = c(0, 1), sex = c(0, 1), 
                            lp = NA, psi = NA, 
                            data.frame(.draw = 1:max(.$.draw)))) %>%
      group_by(age, sex, .draw) %>% 
      mutate(psi = case_when(
        To == from.ref ~ 1 - sum(psi, na.rm = TRUE),
        TRUE ~ psi)) %>% 
      group_by(From, To, age, sex) %>% 
      select(psi) %>%
      mean_qi(.width = 0.9)
  }
  
  if(!is.null(to.ref) & is.null(from.ref)){
    rs <- c("BB","LINWR","BS","HI","M","WL","CNWR","BLNWR","LL","RLNWR_N","HB", "RLNWR_S", "P") # added
    
    for(r in rs){
      y <- out %>% spread_draws(b.dist, b.size1, b.size2, b.sex, b.age, b0) %>%
        expand_grid(., 
                    refuge_datum %>%
                      filter(From == r, To != r)) %>%
        mutate(#b0 = case_when(
          # From == "BB" ~ b0.BB, 
          # From == "LINWR" ~ b0.LINWR,
          # From == "BS" ~ b0.BS, 
          # From == "HI" ~ b0.HI, 
          # From == "M" ~ b0.M, 
          # From == "WL" ~ b0.WL,
          # From == "CNWR" ~ b0.CNWR,
          # From == "BLNWR" ~ b0.BLNWR,
          # From == "LL" ~ b0.LL,
          # From == "RLNWR_N" ~ b0.RLNWR_N,    # Fixed refs
          # From == "HB" ~ b0.HB,
          # From == "RLNWR_S" ~ b0.RLNWR_S,    # 2
          # From == "P" ~ b0.P),               # 3
          lp = exp(b0 + b.dist * dist + b.size1 * Size1 + b.size2 * Size2 + 
                     b.sex * sex + b.age * age)) %>%
        group_by(age, sex, .draw) %>%
        mutate(psi = lp / (1 + sum(lp))) %>% # log p to psi
        bind_rows(expand_grid(From = r, To = r, age = c(0, 1), sex = c(0, 1), 
                              lp = NA, psi = NA,  
                              data.frame(.draw = 1:max(.$.draw)))) %>%
        mutate(psi = case_when(
          To == r ~ 1 - sum(psi, na.rm = TRUE),
          TRUE ~ psi)) %>% 
        filter(To == to.ref)
      
      if(r == "BB"){
        x <- y %>% 
          group_by(From, To, age, sex) %>% 
          select(psi) %>%
          mean_qi(.width = 0.9)
      } else {
        x <- x %>% 
          bind_rows(y %>% 
                      group_by(From, To, age, sex) %>% 
                      select(psi) %>%
                      mean_qi(.width = 0.9))
      }
    }
  }
  
  if(!is.null(to.ref) & !is.null(from.ref)){
    x <-  out %>% spread_draws(b.dist, b.size1, b.size2, b.sex, b.age, b0) %>%
      expand_grid(., 
                  refuge_datum %>%
                    filter(From == from.ref, To != from.ref)) %>%
      mutate(#b0 = case_when(
        # From == "BB" ~ b0.BB, 
        # From == "LINWR" ~ b0.LINWR,
        # From == "BS" ~ b0.BS, 
        # From == "HI" ~ b0.HI, 
        # From == "M" ~ b0.M, 
        # From == "WL" ~ b0.WL,
        # From == "CNWR" ~ b0.CNWR,
        # From == "BLNWR" ~ b0.BLNWR,
        # From == "LL" ~ b0.LL,
        # From == "RLNWR_N" ~ b0.RLNWR_N,    # Fixed refs
        # From == "HB" ~ b0.HB,
        # From == "RLNWR_S" ~ b0.RLNWR_S,    # 2
        # From == "P" ~ b0.P),               # 3
        lp = exp(b0 + b.dist * dist + b.size1 * Size1 + b.size2 * Size2 + 
                   b.sex * sex + b.age * age)) %>%
      group_by(age, sex, .draw) %>%
      mutate(psi = lp / (1 + sum(lp))) %>%
      bind_rows(expand_grid(From = from.ref, To = from.ref, age = c(0, 1), sex = c(0, 1), 
                            lp = NA, psi = NA, 
                            data.frame(.draw = 1:max(.$.draw)))) %>%
      group_by(age, sex, .draw) %>% 
      mutate(psi = case_when(
        To == from.ref ~ 1 - sum(psi, na.rm = TRUE),
        TRUE ~ psi)) %>% 
      group_by(From, To, age, sex) %>%
      filter(To == to.ref) %>% 
      select(psi, .draw) %>%
      mean_qi(.width = 0.9)
  }
  
  return(x)
}

# Calculate probability of direction (i.e., direction of the effect)
# This is really cool - I remember Ally showing me this a while back

calc_prob <- function(model = base, param = "b0.l", relationship = "positive") {
  out <- model
  sims <- eval(parse(text = paste("out$sims.list$", param, sep = "")))
  partial <- ifelse(relationship == "positive", length(sims[which(sims > 0)]), 
                    length(sims[which(sims < 0)]))
  prob <- partial / length(sims)
  
  return(prob)
}