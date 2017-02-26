#' Multisite Causal Mediation Analysis
#'
#' Performs causal mediation analysis in multisite trials. It is used to estimate both the population average and between-site variance of direct and indirect effects. 
#'
#' @param data The data set for analysis.
#' @param y The name of the outcome variable (string).
#' @param treatment The name of the treatment variable (string).
#' @param mediator The name of the mediator variable (string).
#' @param X A vector of variable names (string) of pretreatment confounders, which will be included in the propensity score model. For now, the multilevel propensity score model only allows for one random intercept.
#' @param site The variable name for the site ID (string).
#' @return A list contains the estimates of the between-site variance of direct effect, that of indirect effect, and the correlation between the direct and indirect effects across sites ($Random_effects), and the population average direct and indirect effect estimates along with their hypothesis testing results ($Fixed_effects).
#' @author Xu Qin and Guanglei Hong
#' @references Qin, X., & Hong, G (in press). A weighting method for assessing between-site heterogeneity in causal mediation mechanism. Journal of Educational and Behavioral Statistics.
#' @export
#' @importFrom stats as.formula binomial coef fitted lm model.matrix pnorm
#' @importFrom lme4 VarCorr fixef glmer ranef
#' @importFrom Matrix bdiag
#' @importFrom statmod gauss.quad.prob
#' @examples 
#' data(sim)
#'
#' msmediate(data = sim, y = "y", treatment = "tr", mediator = "me", X = c("x1", "x2", "x3"), site = "site")
#'                     
msmediate = function(data, y, treatment, mediator, X, site) {
  data = as.data.frame(data)
  colnames(data)[which(colnames(data) == treatment)] = "tr"
  colnames(data)[which(colnames(data) == mediator)] = "me"
  colnames(data)[which(colnames(data) == site)] = "site"
  colnames(data)[which(colnames(data) == y)] = "y"
  ranX = 1
  data = data[order(data$site), ]
  data1 = data[data$tr == 1, ]
  data0 = data[data$tr == 0, ]
  
  l_tr = suppressWarnings({glmer(as.formula(paste("me", "~", paste(X, collapse="+"), "+(", ranX, "|site)")), data = data1, family = binomial, nAGQ = 10)})
  data$p1[data$tr == 1] = fitted(l_tr)
  l_ctrl = suppressWarnings({glmer(as.formula(paste("me", "~", paste(X, collapse="+"), "+(", ranX, "|site)")), data = data0, family = binomial, nAGQ = 10)})
  data$p0[data$tr == 0] = fitted(l_ctrl)
  
  predict.lmer = function(l, data, X, ranX) {
    nj = NULL
    for (j in as.numeric(rownames(ranef(l)$site))) {
      nj = c(nj, sum(data$site == j))
    }
    ranef = NULL
    for (k in 1:dim(ranef(l)$site)[2]) {
      ranef = cbind(ranef, rep(ranef(l)$site[, k], nj))
    }
    pred_logit = X %*% fixef(l) + apply(ranX * ranef, 1, sum)
    pred = exp(pred_logit)/(exp(pred_logit) + 1)
  }
  data$p0[data$tr == 1] = predict.lmer(l_ctrl, data1, model.matrix(l_tr), rep(1, nrow(data1)))
  data$p1[data$tr == 0] = predict.lmer(l_tr, data0, model.matrix(l_ctrl), rep(1, nrow(data0)))
  
  data$rmpw[data$tr == 1 & data$me == 1] = (data$p0/data$p1)[data$tr == 1 & data$me == 1]
  data$rmpw[data$tr == 1 & data$me == 0] = ((1 - data$p0)/(1 - data$p1))[data$tr == 1 & data$me == 0]
  data$rmpw[data$tr == 0] = 1

  data = data[order(data$site),]
  data0 = data[data$tr == 0, ]
  data1 = data[data$tr == 1, ]
  
  f = function(x, alpha, s, F, B){
    return(1 / (1 + exp(-(x %*% alpha + s * F * B)))) # When the estimated variance of both two or one of the two random effects is 1
  }
  
  nnodes = 10
  temp = gauss.quad.prob(nnodes, "normal") # Here, random effects are assumed normal
  nodes = temp$nodes
  weights = temp$weights
  # When there is one random effect, i.e. when the estimated variance of one random effect is 0
  B_1 = t(as.matrix(nodes))
  A_1 = weights
  
  x = model.matrix(as.formula(paste("me", "~", paste(X, collapse="+"))), data = data) # Here we have the same covariates in the models for both the treatment and control group
  h1_fun = function(data0, data){
    if(data0$tr[1] == 0){
      l0 = l_ctrl
    } else{
      l0 = l_tr
    }        
    alpha0 = fixef(l0)
    if(round(VarCorr(l0)$site[1], 5) == 0){
      F0 = 0
      v0 = 0
      s0 = matrix(0, nrow(x), 1)
      J0 = 0
      A0 = 1
      B0 = as.matrix(0)
    } else {
      F0 = sqrt(VarCorr(l0)$site[1])
      v0 = ranef(l0)$site
      s0 = matrix(1, nrow(x), 1)
      J0 = 1
      A0 = A_1
      B0 = B_1
    }
    p0 = NULL
    for(q in 1:length(A0)){
      p0 = cbind(p0, f(x = x, alpha = alpha0, s = s0, F = F0, B = B0[, q]))
    }
    p0 = as.matrix(p0)
    
    h_alpha0 = NULL
    h_F0 = NULL
    for(j in unique(data$site)){
      l0j = 0
      h_alpha0j = 0
      h_F0j = 0
      for(q in 1:length(A0)){
        p0q = p0[, q]
        if(data0$tr[1] == 0){
          f0jq = 1
          for(i in which(data$site == j)){
            f0jq = f0jq * (p0q[i]^data$me[i] * (1 - p0q[i])^(1 - data$me[i]))^(1 - data$tr[i])
          }
          l0j = l0j + f0jq * A0[q]
          h_alpha0j = h_alpha0j + as.numeric(f0jq * A0[q] * (1 - data$tr[data$site == j]) * (data$me[data$site == j] - p0q[data$site == j])) * x[data$site == j,]      
          h_F0j =  h_F0j + as.numeric(f0jq * A0[q] * (1 - data$tr[data$site == j]) * (data$me[data$site == j] - p0q[data$site == j])) * kronecker(t(B0[, q]), s0[data$site == j,]) %*% t(J0)    
        } else {
          f0jq = 1
          for(i in which(data$site == j)){
            f0jq = f0jq * (p0q[i]^data$me[i] * (1 - p0q[i])^(1 - data$me[i]))^(data$tr[i])
          }
          l0j = l0j + f0jq * A0[q]
          h_alpha0j = h_alpha0j + as.numeric(f0jq * A0[q] * data$tr[data$site == j] * (data$me[data$site == j] - p0q[data$site == j])) * x[data$site == j,]              
          h_F0j =  h_F0j + as.numeric(f0jq * A0[q] * data$tr[data$site == j] * (data$me[data$site == j] - p0q[data$site == j])) * kronecker(t(B0[, q]), s0[data$site == j,]) %*% t(J0)    
        }        
      }
      h_alpha0 = rbind(h_alpha0, h_alpha0j/l0j)
      h_F0 = rbind(h_F0, h_F0j/l0j)
    }
    
    if(round(VarCorr(l0)$site[1], 5) == 0){
      h1_0 = h_alpha0
    } else {
      h1_0 = cbind(h_alpha0, h_F0)
    }
    
    if(round(VarCorr(l0)$site[1], 5) == 0){
      theta0 = matrix(0, 1, length(unique(data0$site)))
    } else {
      theta0 = solve(F0) %*% t(v0)
    }
    
    return(list(h1 = h1_0, theta = theta0, s = s0, J = J0))
  }
  
  h1_list0 = h1_fun(data0, data)
  h1_list1 = h1_fun(data1, data)
  h1_0 = h1_list0$h1 
  h1_1 = h1_list1$h1
  h1 = cbind(h1_0, h1_1)
  p = ncol(h1)
  theta0 = h1_list0$theta
  theta1 = h1_list1$theta
  s0 = h1_list0$s
  s1 = h1_list1$s
  J0 = h1_list0$J
  J1 = h1_list1$J
  
  N = nrow(data)
  J = length(unique(data$site))
  h2 = matrix(0, N, J * 3)
  mu = NULL
  nj = NULL
  for (j in unique(data$site)) {
    dataj = data[data$site == j, ]
    nj = c(nj, nrow(dataj))
    mu_0 = mean(dataj$y[dataj$tr == 0])
    mu_counter = sum(dataj$y * dataj$rmpw * dataj$tr)/sum(dataj$rmpw * dataj$tr)
    mu_1 = mean(dataj$y[dataj$tr == 1])
    h_0 = (data$y - mu_0) * (1 - data$tr) * (data$site == j)
    h_counter = (data$y - mu_counter) * data$rmpw * data$tr * (data$site == j)
    h_1 = (data$y - mu_1) * data$tr * (data$site == j)
    h2[, (3 * (which(unique(data$site) == j) - 1) + 1):(3 * which(unique(data$site) == j))] = cbind(h_0, 
                                                                                                    h_counter, h_1)
    mu = c(mu, mu_0, mu_counter, mu_1)
  }
  
  h = cbind(h1, h2)  # Stack the moment functions from both steps
  H = 1/N * t(as.matrix(h)) %*% as.matrix(h)  # See Equation S13
  ## The derivation of R (See Equation S14)
  R11 = as.matrix(bdiag(-1/N * t(h1_0) %*% as.matrix(h1_0), -1/N * t(h1_1) %*% as.matrix(h1_1)))
  
  R22_diag = NULL
  for (j in unique(data$site)) {
    R22_diag = c(R22_diag, mean(-(1 - data$tr) * (data$site == j)), mean(-data$rmpw * data$tr * (data$site == j)), mean(-data$tr * (data$site == j)))
  }
  R22 = diag(R22_diag)
  
  p0 = data$p0
  p1 = data$p1
  R21 = NULL
  for (j in 1:J) {
    dh_counter_dalpha0 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (data$me/p1 - (1 - data$me)/(1 - p1)) *   
                                 p0 * (1 - p0) * x, 2, mean)
    dh_counter_dF0 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (data$me/p1 - (1 - data$me)/(1 - p1)) * p0 
                           * (1 - p0) * kronecker(t(theta0[, j]), s0) %*% t(J0), 2, mean)
    dh_counter_dalpha1 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (-data$me * p0/p1 * (1 - p1) + (1 - data$me) 
                                                                  * (1 - p0)/(1 - p1) * p1) * x, 2, mean)
    dh_counter_dF1 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (-data$me * p0/p1 * (1 - p1) + (1 - data$me) 
                                                              * (1 - p0)/(1 - p1) * p1) * kronecker(t(theta1[, j]), s1) %*% t(J1), 2, mean)
    if (ncol(h1_0) - ncol(x) == 0 & ncol(h1_1) - ncol(x) == 0) {
      R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dalpha1), rep(0, p)))
    } else if (ncol(h1_0) - ncol(x) == 0) {
      R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dalpha1, dh_counter_dF1), 
                             rep(0, p)))
    } else if (ncol(h1_1) - ncol(x) == 0) {
      R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dF0, dh_counter_dalpha1), 
                             rep(0, p)))
    } else {
      R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dF0, dh_counter_dalpha1, 
                                          dh_counter_dF1), rep(0, p)))
    }
  }
  
  R12 = matrix(0, p, 3 * J)
  
  R = matrix(0, p + 3 * J, p + 3 * J)
  R[1:p, 1:p] = R11
  R[1:p, (p + 1):(p + 3 * J)] = R12
  R[(p + 1):(p + 3 * J), 1:p] = R21
  R[(p + 1):(p + 3 * J), (p + 1):(p + 3 * J)] = R22
  
  V_all = 1/N * solve(R) %*% H %*% t(solve(R))  # See Appendix A
  V_mu = V_all[(ncol(V_all) - 3 * J + 1):ncol(V_all), (ncol(V_all) - 3 * J + 1):ncol(V_all)]  
  
  Phi = matrix(c(-1, 0, 1, -1, 0, 1), 2, 3)
  for (j in 1:(J - 1)) {
    Phi = as.matrix(bdiag(Phi, matrix(c(-1, 0, 1, -1, 0, 1), 2, 3)))
  }
  beta = Phi %*% mu # Estimated site-specific direct and indirect effects
  V_beta = Phi %*% V_mu %*% t(Phi)
  dej = beta[seq(1, 2 * J, 2)]
  iej = beta[seq(2, 2 * J, 2)]
  
  Psi = cbind(rep(c(1, 0), J), rep(c(0, 1), J))
  gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% beta
  de = gamma[1]
  ie = gamma[2]
  G = 0
  V = 0
  for(j in 1:J){
    G = G + (beta[(2 * j - 1): (2 * j)] - gamma) %*% t(beta[(2 * j - 1): (2 * j)] - gamma)    
    Vj = V_beta[(2 * j - 1):(2 * j), (2 * j - 1):(2 * j)]
    V = V + Vj
  }
  tau = 1/(J - 1) * (G - V + 1/J * t(Psi) %*% V_beta %*% Psi)
  
  V_gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% (V_beta + kronecker(diag(rep(1,J)), tau)) %*% Psi %*%   
    solve(t(Psi) %*% Psi)
  V_de = V_gamma[1, 1]
  V_ie = V_gamma[2, 2]   
  
  se_de = sqrt(V_de)
  se_ie = sqrt(V_ie)
  
  z_de = de/se_de
  z_ie = ie/se_ie
  
  p_de = (1 - pnorm(abs(z_de))) * 2
  p_ie = (1 - pnorm(abs(z_ie))) * 2
  
  est = cbind(round(c(de, ie), 5), round(c(se_de, se_ie), 5), round(c(z_de, z_ie), 3), round(c(p_de, p_ie), 3))
  
  sig = NULL
  sig[est[, 4] <= 0.001] = "**"
  sig[est[, 4] > 0.001 & est[, 4] <= 0.01] = "*"
  sig[est[, 4] > 0.01 & est[, 4] <= 0.05] = "."
  sig[est[, 4] > 0.05] = ""
  
  est = cbind(est, sig)
  est[est[, 4] < 0.001, 4] = "<0.001"
  est = as.data.frame(est)
  colnames(est) = c("Estimate", "Std.Error", "t value", "Pr(>|t|)", "")
  rownames(est) = c("Natrual Direct Effect", "Natrual Indirect Effect")
  
  cor_de_ie = tau[1, 2]/sqrt(tau[1, 1] * tau[2, 2])
  tau[tau < 0] = 0
  cor_de_ie[cor_de_ie > 1] = 1
  cor_de_ie[cor_de_ie < -1] = -1
  
  var = cbind(round(c(tau[1, 1], tau[2, 2]), 3), round(c(sqrt(tau[1, 1]), sqrt(tau[2, 2])), 3), c("", round(cor_de_ie, 3)))
  var = as.data.frame(var)
  colnames(var) = c("Variance", "Std.Dev.", "Corr")
  rownames(var) = c("Natrual Direct Effect", "Natrual Indirect Effect")
  
  result = list(Random_effects = var, Fixed_effects = est)
  
  return(result)
}

#' Include Variance Testing for Multisite Causal Mediation Analysis
#'
#' Performs hypothesis testing for the between-site variance of direct effect and that of indirect effect, besides providing the same output as given by the function msmediate().
#'
#' @param data The data set for analysis.
#' @param y The name of the outcome variable (string).
#' @param treatment The name of the treatment variable (string).
#' @param mediator The name of the mediator variable (string).
#' @param X A vector of variable names (string) of pretreatment confounders, which will be included in the propensity score model. For now, the multilevel propensity score model only allows for one random intercept.
#' @param site The variable name for the site ID (string).
#' @param npermute The number of permutations for the permutation test. The default value is 200. It may take a long time, depending on the sample size and the length of X.
#' @return A list contains the hypothesis testing results of the between-site variance of the causal effects, besides the same output as given by the function msmediate().
#' @author Xu Qin and Guanglei Hong
#' @references Qin, X., & Hong, G (in press). A weighting method for assessing between-site heterogeneity in causal mediation mechanism. Journal of Educational and Behavioral Statistics.
#' @export
#' @importFrom stats as.formula binomial coef fitted lm model.matrix pnorm
#' @importFrom lme4 VarCorr fixef glmer ranef
#' @importFrom Matrix bdiag
#' @importFrom statmod gauss.quad.prob
#' @examples 
#' data(sim)
#'
#' vartest.msmediate(data = sim, y = "y", treatment = "tr", mediator = "me", X = c("x1", "x2", "x3"), site = "site", npermute = 2)
#'
vartest.msmediate = function(data, y, treatment, mediator, X, site, npermute = 200) {
  data = as.data.frame(data)
  colnames(data)[which(colnames(data) == treatment)] = "tr"
  colnames(data)[which(colnames(data) == mediator)] = "me"
  colnames(data)[which(colnames(data) == site)] = "site"
  colnames(data)[which(colnames(data) == y)] = "y"
  est = function(data, y, treatment, mediator, X, site){
    ranX = 1
    data = data[order(data$site), ]
    data1 = data[data$tr == 1, ]
    data0 = data[data$tr == 0, ]
    
    l_tr = suppressWarnings({glmer(as.formula(paste("me", "~", paste(X, collapse="+"), "+(", ranX, "|site)")), data = data1, family = binomial, nAGQ = 10)})
    data$p1[data$tr == 1] = fitted(l_tr)
    l_ctrl = suppressWarnings({glmer(as.formula(paste("me", "~", paste(X, collapse="+"), "+(", ranX, "|site)")), data = data0, family = binomial, nAGQ = 10)})
    data$p0[data$tr == 0] = fitted(l_ctrl)
    
    predict.lmer = function(l, data, X, ranX) {
      nj = NULL
      for (j in as.numeric(rownames(ranef(l)$site))) {
        nj = c(nj, sum(data$site == j))
      }
      ranef = NULL
      for (k in 1:dim(ranef(l)$site)[2]) {
        ranef = cbind(ranef, rep(ranef(l)$site[, k], nj))
      }
      pred_logit = X %*% fixef(l) + apply(ranX * ranef, 1, sum)
      pred = exp(pred_logit)/(exp(pred_logit) + 1)
    }
    data$p0[data$tr == 1] = predict.lmer(l_ctrl, data1, model.matrix(l_tr), rep(1, nrow(data1)))
    data$p1[data$tr == 0] = predict.lmer(l_tr, data0, model.matrix(l_ctrl), rep(1, nrow(data0)))
    
    data$rmpw[data$tr == 1 & data$me == 1] = (data$p0/data$p1)[data$tr == 1 & data$me == 1]
    data$rmpw[data$tr == 1 & data$me == 0] = ((1 - data$p0)/(1 - data$p1))[data$tr == 1 & data$me == 0]
    data$rmpw[data$tr == 0] = 1
    
    data = data[order(data$site),]
    data0 = data[data$tr == 0, ]
    data1 = data[data$tr == 1, ]
    
    f = function(x, alpha, s, F, B){
      return(1 / (1 + exp(-(x %*% alpha + s * F * B)))) # When the estimated variance of both two or one of the two random effects is 1
    }
    
    nnodes = 10
    temp = gauss.quad.prob(nnodes, "normal") # Here, random effects are assumed normal
    nodes = temp$nodes
    weights = temp$weights
    # When there is one random effect, i.e. when the estimated variance of one random effect is 0
    B_1 = t(as.matrix(nodes))
    A_1 = weights
    
    x = model.matrix(as.formula(paste("me", "~", paste(X, collapse="+"))), data = data) # Here we have the same covariates in the models for both the treatment and control group
    h1_fun = function(data0, data){
      if(data0$tr[1] == 0){
        l0 = l_ctrl
      } else{
        l0 = l_tr
      }        
      alpha0 = fixef(l0)
      if(round(VarCorr(l0)$site[1], 5) == 0){
        F0 = 0
        v0 = 0
        s0 = matrix(0, nrow(x), 1)
        J0 = 0
        A0 = 1
        B0 = as.matrix(0)
      } else {
        F0 = sqrt(VarCorr(l0)$site[1])
        v0 = ranef(l0)$site
        s0 = matrix(1, nrow(x), 1)
        J0 = 1
        A0 = A_1
        B0 = B_1
      }
      p0 = NULL
      for(q in 1:length(A0)){
        p0 = cbind(p0, f(x = x, alpha = alpha0, s = s0, F = F0, B = B0[, q]))
      }
      p0 = as.matrix(p0)
      
      h_alpha0 = NULL
      h_F0 = NULL
      for(j in unique(data$site)){
        l0j = 0
        h_alpha0j = 0
        h_F0j = 0
        for(q in 1:length(A0)){
          p0q = p0[, q]
          if(data0$tr[1] == 0){
            f0jq = 1
            for(i in which(data$site == j)){
              f0jq = f0jq * (p0q[i]^data$me[i] * (1 - p0q[i])^(1 - data$me[i]))^(1 - data$tr[i])
            }
            l0j = l0j + f0jq * A0[q]
            h_alpha0j = h_alpha0j + as.numeric(f0jq * A0[q] * (1 - data$tr[data$site == j]) * (data$me[data$site == j] - p0q[data$site == j])) * x[data$site == j,]      
            h_F0j =  h_F0j + as.numeric(f0jq * A0[q] * (1 - data$tr[data$site == j]) * (data$me[data$site == j] - p0q[data$site == j])) * kronecker(t(B0[, q]), s0[data$site == j,]) %*% t(J0)    
          } else {
            f0jq = 1
            for(i in which(data$site == j)){
              f0jq = f0jq * (p0q[i]^data$me[i] * (1 - p0q[i])^(1 - data$me[i]))^(data$tr[i])
            }
            l0j = l0j + f0jq * A0[q]
            h_alpha0j = h_alpha0j + as.numeric(f0jq * A0[q] * data$tr[data$site == j] * (data$me[data$site == j] - p0q[data$site == j])) * x[data$site == j,]              
            h_F0j =  h_F0j + as.numeric(f0jq * A0[q] * data$tr[data$site == j] * (data$me[data$site == j] - p0q[data$site == j])) * kronecker(t(B0[, q]), s0[data$site == j,]) %*% t(J0)    
          }        
        }
        h_alpha0 = rbind(h_alpha0, h_alpha0j/l0j)
        h_F0 = rbind(h_F0, h_F0j/l0j)
      }
      
      if(round(VarCorr(l0)$site[1], 5) == 0){
        h1_0 = h_alpha0
      } else {
        h1_0 = cbind(h_alpha0, h_F0)
      }
      
      if(round(VarCorr(l0)$site[1], 5) == 0){
        theta0 = matrix(0, 1, length(unique(data0$site)))
      } else {
        theta0 = solve(F0) %*% t(v0)
      }
      
      return(list(h1 = h1_0, theta = theta0, s = s0, J = J0))
    }
    
    h1_list0 = h1_fun(data0, data)
    h1_list1 = h1_fun(data1, data)
    h1_0 = h1_list0$h1 
    h1_1 = h1_list1$h1
    h1 = cbind(h1_0, h1_1)
    p = ncol(h1)
    theta0 = h1_list0$theta
    theta1 = h1_list1$theta
    s0 = h1_list0$s
    s1 = h1_list1$s
    J0 = h1_list0$J
    J1 = h1_list1$J
    
    N = nrow(data)
    J = length(unique(data$site))
    h2 = matrix(0, N, J * 3)
    mu = NULL
    nj = NULL
    for (j in unique(data$site)) {
      dataj = data[data$site == j, ]
      nj = c(nj, nrow(dataj))
      mu_0 = mean(dataj$y[dataj$tr == 0])
      mu_counter = sum(dataj$y * dataj$rmpw * dataj$tr)/sum(dataj$rmpw * dataj$tr)
      mu_1 = mean(dataj$y[dataj$tr == 1])
      h_0 = (data$y - mu_0) * (1 - data$tr) * (data$site == j)
      h_counter = (data$y - mu_counter) * data$rmpw * data$tr * (data$site == j)
      h_1 = (data$y - mu_1) * data$tr * (data$site == j)
      h2[, (3 * (which(unique(data$site) == j) - 1) + 1):(3 * which(unique(data$site) == j))] = cbind(h_0, 
                                                                                                      h_counter, h_1)
      mu = c(mu, mu_0, mu_counter, mu_1)
    }
    
    h = cbind(h1, h2)  # Stack the moment functions from both steps
    H = 1/N * t(as.matrix(h)) %*% as.matrix(h)  # See Equation S13
    ## The derivation of R (See Equation S14)
    R11 = as.matrix(bdiag(-1/N * t(h1_0) %*% as.matrix(h1_0), -1/N * t(h1_1) %*% as.matrix(h1_1)))
    
    R22_diag = NULL
    for (j in unique(data$site)) {
      R22_diag = c(R22_diag, mean(-(1 - data$tr) * (data$site == j)), mean(-data$rmpw * data$tr * (data$site == j)), mean(-data$tr * (data$site == j)))
    }
    R22 = diag(R22_diag)
    
    p0 = data$p0
    p1 = data$p1
    R21 = NULL
    for (j in 1:J) {
      dh_counter_dalpha0 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (data$me/p1 - (1 - data$me)/(1 - p1)) *   
                                   p0 * (1 - p0) * x, 2, mean)
      dh_counter_dF0 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (data$me/p1 - (1 - data$me)/(1 - p1)) * p0 
                             * (1 - p0) * kronecker(t(theta0[, j]), s0) %*% t(J0), 2, mean)
      dh_counter_dalpha1 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (-data$me * p0/p1 * (1 - p1) + (1 - data$me) 
                                                                    * (1 - p0)/(1 - p1) * p1) * x, 2, mean)
      dh_counter_dF1 = apply(h2[, 3 * (j - 1) + 2]/data$rmpw * (-data$me * p0/p1 * (1 - p1) + (1 - data$me) 
                                                                * (1 - p0)/(1 - p1) * p1) * kronecker(t(theta1[, j]), s1) %*% t(J1), 2, mean)
      if (ncol(h1_0) - ncol(x) == 0 & ncol(h1_1) - ncol(x) == 0) {
        R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dalpha1), rep(0, p)))
      } else if (ncol(h1_0) - ncol(x) == 0) {
        R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dalpha1, dh_counter_dF1), 
                               rep(0, p)))
      } else if (ncol(h1_1) - ncol(x) == 0) {
        R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dF0, dh_counter_dalpha1), 
                               rep(0, p)))
      } else {
        R21 = rbind(R21, rbind(rep(0, p), c(dh_counter_dalpha0, dh_counter_dF0, dh_counter_dalpha1, 
                                            dh_counter_dF1), rep(0, p)))
      }
    }
    
    R12 = matrix(0, p, 3 * J)
    
    R = matrix(0, p + 3 * J, p + 3 * J)
    R[1:p, 1:p] = R11
    R[1:p, (p + 1):(p + 3 * J)] = R12
    R[(p + 1):(p + 3 * J), 1:p] = R21
    R[(p + 1):(p + 3 * J), (p + 1):(p + 3 * J)] = R22
    
    V_all = 1/N * solve(R) %*% H %*% t(solve(R))  # See Appendix A
    V_mu = V_all[(ncol(V_all) - 3 * J + 1):ncol(V_all), (ncol(V_all) - 3 * J + 1):ncol(V_all)]  
    
    Phi = matrix(c(-1, 0, 1, -1, 0, 1), 2, 3)
    for (j in 1:(J - 1)) {
      Phi = as.matrix(bdiag(Phi, matrix(c(-1, 0, 1, -1, 0, 1), 2, 3)))
    }
    beta = Phi %*% mu # Estimated site-specific direct and indirect effects
    V_beta = Phi %*% V_mu %*% t(Phi)
    dej = beta[seq(1, 2 * J, 2)]
    iej = beta[seq(2, 2 * J, 2)]
    
    Psi = cbind(rep(c(1, 0), J), rep(c(0, 1), J))
    gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% beta
    de = gamma[1]
    ie = gamma[2]
    
    G = 0
    V = 0
    for(j in 1:J){
      G = G + (beta[(2 * j - 1): (2 * j)] - gamma) %*% t(beta[(2 * j - 1): (2 * j)] - gamma)    
      Vj = V_beta[(2 * j - 1):(2 * j), (2 * j - 1):(2 * j)]
      V = V + Vj
    }
    tau = 1/(J - 1) * (G - V + 1/J * t(Psi) %*% V_beta %*% Psi)
    
    V_gamma = solve(t(Psi) %*% Psi) %*% t(Psi) %*% (V_beta + kronecker(diag(rep(1,J)), tau)) %*% Psi %*%   
      solve(t(Psi) %*% Psi)
    V_de = V_gamma[1, 1]
    V_ie = V_gamma[2, 2]   
    
    se_de = sqrt(V_de)
    se_ie = sqrt(V_ie)
    
    z_de = de/se_de
    z_ie = ie/se_ie
    
    p_de = (1 - pnorm(abs(z_de))) * 2
    p_ie = (1 - pnorm(abs(z_ie))) * 2
    
    est = cbind(round(c(de, ie), 5), round(c(se_de, se_ie), 5), round(c(z_de, z_ie), 3), round(c(p_de, p_ie), 3))
    
    sig = NULL
    sig[est[, 4] <= 0.001] = "**"
    sig[est[, 4] > 0.001 & est[, 4] <= 0.01] = "*"
    sig[est[, 4] > 0.01 & est[, 4] <= 0.05] = "."
    sig[est[, 4] > 0.05] = ""
    
    est = cbind(est, sig)
    est[est[, 4] < 0.001, 4] = "<0.001"
    est = as.data.frame(est)
    colnames(est) = c("Estimate", "Std.Error", "t value", "Pr(>|t|)", "")
    rownames(est) = c("Natrual Direct Effect", "Natrual Indirect Effect")
    
    cor_de_ie = tau[1, 2]/sqrt(tau[1, 1] * tau[2, 2])
    tau[tau < 0] = 0
    cor_de_ie[cor_de_ie > 1] = 1
    cor_de_ie[cor_de_ie < -1] = -1
    
    var = cbind(round(c(tau[1, 1], tau[2, 2]), 3), round(c(sqrt(tau[1, 1]), sqrt(tau[2, 2])), 3), c("", round(cor_de_ie, 3)))
    var = as.data.frame(var)
    colnames(var) = c("Variance", "Std.Dev.", "Corr")
    rownames(var) = c("Natrual Direct Effect", "Natrual Indirect Effect")
    
    chisq_de = 0
    chisq_ie = 0
    for (j in 1:J) {
      Vj = V_beta[(2 * (j - 1) + 1):(2 * j), (2 * (j - 1) + 1):(2 * j)]
      chisq_de = chisq_de + (dej[j] - de)^2/Vj[1, 1]
      chisq_ie = chisq_ie + (iej[j] - ie)^2/Vj[2, 2]
    }
    
    result = list(Random_effects = var, Fixed_effects = est, chisq = c(chisq_de = chisq_de, chisq_ie = chisq_ie))
    
    return(result)
  }
  
  result = suppressWarnings({est(data, y, treatment, mediator, X, site)})
  
  siteID = data$site
  chisq_de_pm = NULL
  chisq_ie_pm = NULL
  for (i in 1:npermute) {
    data$site = sample(siteID, length(siteID), replace = F)
    est_list_pm = suppressWarnings({est(data, y, treatment, mediator, X, site)$chisq})
    chisq_de_pm = c(chisq_de_pm, est_list_pm["chisq_de"])
    chisq_ie_pm = c(chisq_ie_pm, est_list_pm["chisq_ie"])
  }
  
  pvalue_var_de = sum(chisq_de_pm >= result$chisq["chisq_de"])/npermute
  pvalue_var_ie = sum(chisq_ie_pm >= result$chisq["chisq_ie"])/npermute
  
  pvalue = round(c(pvalue_var_de, pvalue_var_ie), 5)
  sig = NULL
  sig[pvalue <= 0.001] = "**"
  sig[pvalue > 0.001 & pvalue <= 0.01] = "*"
  sig[pvalue > 0.01 & pvalue <= 0.05] = "."
  sig[pvalue > 0.05] = ""
  result$Random_effects = cbind(result$Random_effects, pvalue, sig)
  colnames(result$Random_effects)[4] = "p-value"
  colnames(result$Random_effects)[5] = ""
  
  return(list(Random_effects = result$Random_effects, Fixed_effects = result$Fixed_effects))
}