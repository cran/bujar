bjboost.fit <- function(y, ypred) {
  #if(!existsFunction('survfit.km'))
      #if (requireNamespace("survival", quietly = TRUE)) {
      #   survfit.km <- survival::survfitKM
      #} else {
      #	      stop("survival package is needed\n")
      ## do something else not involving survival.
  # }
	#  survfit.km <- getFromNamespace('survfitKM','survival')
  
  if(ncol(y) != 2)
	stop("y is not a right-censored Surv object")
  status <- y[, 2]
  yy <- y[, 1]
  N <- length(yy)
  timeorig <- yy
  order.orig <- 1:N
  dummystrat <- factor(rep(1, N))
        ehat <- timeorig - ypred
	state <- status
	state[ehat == max(ehat)] <- 1
	S <- structure(cbind(ehat, state), class = "Surv", type = "right")
	KM.ehat <- survfitKM(dummystrat, S, conf.type = "none", se.fit = FALSE)
        #if (requireNamespace("survival", quietly = TRUE)) {
	#KM.ehat <- survival::survfitKM(dummystrat, S, conf.type = "none", se.fit = FALSE)
      #} else {
      # 	      stop("survival package is needed\n")
      #  }
	#KM.ehat <- survfit.km(dummystrat, S, conf.type = "none", se.fit = FALSE)
	n.risk <- KM.ehat$n.risk
	surv <- KM.ehat$surv
	repeats <- c(diff( - n.risk), n.risk[length(n.risk)])
	surv <- rep(surv, repeats)
	w <-  - diff(c(1, surv))
	m <- order(ehat,  - status)
	bla <- cumsum((w * ehat[m]))
	bla <- (bla[length(bla)] - bla)/(surv + state[m])	## Put bla back into original order
	bl <- bla
	bl[(1:N)[m]] <- bla
                yhat <- ypred + bl
	yy[state == 0] <- yhat[state == 0]
  y.imputed <- yy
  list(y=y,y.imputed=y.imputed)
}
