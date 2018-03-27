






############ Studentized Residuals to narrow down quickly

res <- residuals(model)
H <- hatvalues(model)
sigma <- summary(model)$sigm
sres <- map_dbl(1:length(res), ~ res[[.]]/(sigma*sqrt(1-H[[.]]))) 


############ Cook's D to account for influence 

function (model) 
{
  
  ifelse(as.character(model@call)[3] == "data.update", data.adapted <- model.frame(model), 
         data.adapted <- get(as.character(model@call)[3]))
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
  n.obs <- nrow(data.adapted)
  or.fixed <- matrix(ncol = n.pred, nrow = 1, data = fixef(model)[original.no.estex])
  dimnames(or.fixed) <- list(NULL, names(fixef(model))[original.no.estex])
  or.se <- matrix(ncol = n.pred, nrow = 1, data = se.fixef(model)[original.no.estex])
  dimnames(or.se) <- list(NULL, names(fixef(model))[original.no.estex])
  or.vcov <- as.matrix(vcov(model)[original.no.estex, original.no.estex])
  dimnames(or.vcov) <- list(names(fixef(model)[original.no.estex]), 
                            names(fixef(model)[original.no.estex]))
  or.test <- coef(summary(model))[original.no.estex, 3]
  if (obs) {
    if (is.null(select)) {
      alt.fixed <- matrix(ncol = n.pred, nrow = n.obs, 
                          data = NA)
      dimnames(alt.fixed) <- list(1:n.obs, names(fixef(model))[original.no.estex])
      alt.se <- matrix(ncol = n.pred, nrow = n.obs, data = NA)
      dimnames(alt.se) <- list(1:n.obs, names(fixef(model))[original.no.estex])
      alt.vcov <- list()
      alt.test <- matrix(ncol = n.pred, nrow = n.obs, data = NA)
      dimnames(alt.test) <- list(1:n.obs, names(fixef(model))[original.no.estex])
      for (i in 1:n.obs) {
        if (count == TRUE) {
          print(n.obs + 1 - i)
        }
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
    }
    if (!is.null(select)) {
      model.updated <- exclude.influence(model, obs = select)
      altered.no.estex <- which(substr(names(fixef(model.updated)), 
                                       1, 6) != "estex.")
      alt.fixed <- matrix(ncol = n.pred, nrow = 1, data = fixef(model.updated)[altered.no.estex])
      dimnames(alt.fixed) <- list("Altered model", names(fixef(model.updated))[altered.no.estex])
      alt.se <- matrix(ncol = n.pred, nrow = 1, data = se.fixef(model.updated)[altered.no.estex])
      dimnames(alt.se) <- list("Altered model", names(fixef(model.updated))[altered.no.estex])
      alt.vcov <- list()
      alt.vcov[[1]] <- as.matrix(vcov(model.updated)[altered.no.estex, 
                                                     altered.no.estex])
      dimnames(alt.vcov[[1]]) <- list(names(fixef(model.updated)[altered.no.estex]), 
                                      names(fixef(model.updated)[altered.no.estex]))
      alt.test <- matrix(ncol = n.pred, nrow = 1, data = coef(summary(model.updated))[, 
                                                                                      3][altered.no.estex])
      dimnames(alt.test) <- list("Altered model", names(fixef(model.updated))[altered.no.estex])
    }
  }
  estex <- list(or.fixed = or.fixed, or.se = or.se, or.vcov = or.vcov, 
                or.test = or.test, alt.fixed = alt.fixed, alt.se = alt.se, 
                alt.vcov = alt.vcov, alt.test = alt.test)
  class(estex) <- "estex"
  return(estex)
}

