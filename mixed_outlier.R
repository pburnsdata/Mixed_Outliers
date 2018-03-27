
#################################################################################
#
# Combine Studentized residuals and Cook's D for Mixed Model Outlier Detection
#   
#   Get rid of some of the false positives from residuals resulting from obs. hijacking 
#   influence of real outliers without the ridiculous wait time of removing every obs.
#
#################################################################################


influence_modified <- function (model, cutoff=3) {
  
  
  #####################################
  #  Parameters
  #
  #  model: Linear Mixed Model (lmerMod class)
  #  cutoff: cutoff for studentized residuals
  #
  #####################################
  
  
  ############ Studentized Residuals to narrow down
  
  
  res <- residuals(model)
  H <- hatvalues(model)
  sigma <- summary(model)$sigm
  sres <- map_dbl(1:length(res), ~ res[[.]]/(sigma*sqrt(1-H[[.]]))) 
  test = abs(sres) > cutoff
  
  
  ############ Cook's D to slice off some false positives
  #Started from influence.ME
  
  ifelse(as.character(model@call)[3] == "data.update", data.adapted <- model.frame(model), 
         data.adapted <- get(as.character(model@call)[3]))
  
  data.adapted$outlier<-as.factor(test)  #Add column with outlier status from residuals for later subsetting
  
  original.no.estex <- which(substr(names(fixef(model)), 1, 
                                    6) != "estex.")
  n.pred <- length(fixef(model)[original.no.estex])
  if ("(weights)" %in% names(data.adapted)) {
    names(data.adapted)[names(data.adapted) == "(weights)"] <- as.character(model@call$weights)
  }
  if ("(offset)" %in% names(data.adapted)) {
    names(data.adapted)[names(data.adapted) == "(offset)"] <- as.character(model@call$offset)
  }
  if (sum(grepl("offset", names(data.adapted))) > 0) {
    names(data.adapted)[grep("offset", names(data.adapted))] <- gsub("offset\\\\(|\\\\)", 
                                                                     "", names(data.adapted)[grep("offset", names(data.adapted))])
  }
  
  data.adapted <- data.adapted[which(data.adapted$outlier==TRUE), , drop=FALSE] #subset so we're only checking those with big residuals
  
  n.obs <- nrow(data.adapted)
  or.fixed <- matrix(ncol = n.pred, nrow = 1, data = fixef(model)[original.no.estex])
  dimnames(or.fixed) <- list(NULL, names(fixef(model))[original.no.estex])
  or.se <- matrix(ncol = n.pred, nrow = 1, data = se.fixef(model)[original.no.estex])
  dimnames(or.se) <- list(NULL, names(fixef(model))[original.no.estex])
  or.vcov <- as.matrix(vcov(model)[original.no.estex, original.no.estex])
  dimnames(or.vcov) <- list(names(fixef(model)[original.no.estex]), 
                            names(fixef(model)[original.no.estex]))
  or.test <- coef(summary(model))[original.no.estex, 3]

      alt.fixed <- matrix(ncol = n.pred, nrow = n.obs, 
                          data = NA)
      dimnames(alt.fixed) <- list(1:n.obs, names(fixef(model))[original.no.estex])
      alt.se <- matrix(ncol = n.pred, nrow = n.obs, data = NA)
      dimnames(alt.se) <- list(1:n.obs, names(fixef(model))[original.no.estex])
      alt.vcov <- list()
      alt.test <- matrix(ncol = n.pred, nrow = n.obs, data = NA)
      dimnames(alt.test) <- list(1:n.obs, names(fixef(model))[original.no.estex])
    for (i in 1:n.obs) {
        model.updated <- exclude.influence(model, obs = i)
        altered.no.estex <- which(substr(names(fixef(model.updated)), 
                                         1, 6) != "estex.")
        alt.fixed[i, ] <- as.matrix(fixef(model.updated)[altered.no.estex])
        alt.se[i, ] <- as.matrix(se.fixef(model.updated)[altered.no.estex])
        alt.vcov[[i]] <- as.matrix(vcov(model.updated)[altered.no.estex, 
                                                       altered.no.estex])
        alt.test[i, ] <- as.matrix(coef(summary(model.updated))[, 
                                                                3][altered.no.estex])
      }
  
  estex <- list(or.fixed = or.fixed, or.se = or.se, or.vcov = or.vcov, 
                or.test = or.test, alt.fixed = alt.fixed, alt.se = alt.se, 
                alt.vcov = alt.vcov, alt.test = alt.test)
  class(estex) <- "estex"
  
  
  n.groups <- dim(estex$alt.fixed)[1]
  n.parameters <- dim(estex$alt.fixed)[2]
  
  sel <- 1:n.parameters
  a <- as.matrix(estex$or.fixed[, sel])
  b <- as.matrix(estex$alt.fixed[, sel])
  e <- NA
  if (n.groups == 1) {
    d <- as.matrix(estex$alt.vcov[[1]][sel, sel])
    e <- t(a - b) %*% solve(d) %*% (a - b)/length(a)
  }
  if (n.groups > 1) {
    for (i in 1:n.groups) {
      d <- as.matrix(estex$alt.vcov[[i]][sel, sel])
      e[i] <- t(a - b[i, ]) %*% solve(d) %*% (a - b[i, 
                                                    ])/length(a)
    }
    e <- as.matrix(e)
    rownames(e) <- rownames(b)
  }

  return(e)
  
}

