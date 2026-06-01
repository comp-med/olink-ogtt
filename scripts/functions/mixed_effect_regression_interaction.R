###########################################
###########################################
#### function to compute linear regression
#### with categories and additional inter-
#### action and adjustment for age and sex

mixed.regression.v2 <- function(
  data,
  feat,
  formel,
  name,
  log.y = F,
  scale = F
) {
  ## 'data' -- data frame containing all needed values
  ## 'feat' -- which feature to use
  ## 'form' -- string of the left sid of the formula
  ## 'name' -- name for the residual plots
  ## 'log.y'-- logical whether to perform log-transformation prior regression
  ## 'scale'-- indicator whether to scale the data to ease comparison

  ## load packages for mixed models
  require("lmerTest")

  ## log-transform the data if needed
  if (log.y) {
    ii <- grep(paste(feat, collapse = "|"), names(data))
    data[, ii] <- log2(data[, ii])
  }

  ## scale if required
  if (scale) {
    ii <- grep(paste(feat, collapse = "|"), names(data))
    data[, ii] <- scale(data[, ii])
  }

  ## how many feature
  n <- length(feat)

  ## start pdf for residual plots
  pdf(paste("graphics/residuals_", name, ".pdf", sep = ""))
  par(
    mfrow = c(3, 3),
    mar = c(2.5, 2.5, 1.5, .5),
    mgp = c(1.3, .4, 0),
    tck = -.02,
    cex.lab = .7,
    cex.axis = .8,
    cex.main = .8
  )

  ## rename
  cat("calculate model :\n")

  ## loop over all variables
  for (j in 1:n) {
    ## define formula --> include baseline values in the adjustment set
    formula <- paste(feat[j], formel, sep = " ~ ")
    cat(formula, "\n")

    ## run mixed linear model
    foo.m <- lmer(as.formula(formula), data = data)
    foo.s <- summary(foo.m)
    #CI      <- confint(foo.m)

    ## use anova to get the fixed effect values for the variables
    foo.a <- anova(foo.m)

    car::qqPlot(
      foo.s$residuals,
      dist = "norm",
      col = "black",
      ylab = "Residual Quantiles",
      main = "Normal Probability Plot",
      pch = 19
    )

    ## perpare storage
    if (j == 1) {
      ## how many terms (same in anova and lm as only binary categories)
      dd <- nrow(foo.s$coefficients)
      # Name f?r die Spalten
      col.nam <- rownames(foo.s$coefficients)
      col.nam[1] <- "Intercept"
      # Arrays zum speichern der Daten --> nicht besetzte Zeilen, um sp?ter Regressionen zu filtern
      pvalue <- array(data = NA, dim = c(n, dd))
      colnames(pvalue) <- paste("pvalue_", col.nam, sep = "")
      # Spaltennamen verteilen, um später einen Datenframe als Ausgabe zu erzeugen
      stat <- array(data = NA, dim = c(n, dd))
      colnames(stat) <- paste("stderr_", col.nam, sep = "")
      beta <- array(data = NA, dim = c(n, dd))
      colnames(beta) <- paste("beta_", col.nam, sep = "")
      ## store the anova p values
      anov <- array(data = NA, dim = c(n, nrow(foo.a)))
      colnames(anov) <- paste("pvalue.anova", rownames(foo.a), sep = "_")
    }

    ## store the results
    for (l in 1:dd) {
      pvalue[j, l] <- foo.s$"coefficients"[l, 5]
      beta[j, l] <- unlist(foo.s$"coefficients"[l, 1])
      stat[j, l] <- unlist(foo.s$"coefficients"[l, 2])
    }

    ## store lambda used
    anov[j, ] <- unlist(foo.a[, 6])
  }

  ## return data frame
  results <- data.frame(feat, beta, pvalue, stat, anov)

  ## add FDR
  if (nrow(results) > 15) {
    for (j in grep("pvalue", names(results), value = T)) {
      results[, gsub("pvalue", "fdr", j)] <- p.adjust(
        results[, j],
        method = "BH"
      )
    }
  }

  dev.off()

  return(results)
}
