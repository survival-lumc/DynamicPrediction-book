semicure <- function(formula = formula(data), cureform,
    data = sys.parent(), link = "logit", model = c("glm", "gam"),
    emmax=100, eps=1e-4, debug = T, tail = c("zero", "none"), savedata = T, ...)
{
    call <- match.call()
    model <- match.arg(model)
    tail <- match.arg(tail)

    nodata <- missing(data)
    m <- if (nodata) model.frame(formula)
        else model.frame(formula, data = data)
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    else {
        datalen <- nrow(Y)
        large <- max(Y[Y[, 2]==0, 1])
        small <- min(Y[Y[, 2]==0, 1])
        delta <- (large-small)/100
        cens <- Y[, 2]
        uncureprob <- Y[, 1]
        uncureprob[Y[, 2]==0] <-
            (large + delta - Y[Y[, 2]==0, 1])/(large - small + 2*delta)
        uncureprob[Y[, 2]==1] <- 1
    }

    coxformula <- paste(paste(as.character(as.expression(formula)),
        collapse = " "), "+ offset(log(uncureprob))")
    logisform <- paste("uncureprob", paste(as.character(as.expression(cureform)),
        collapse = " ")) 

    attach(data)
    on.exit(if (is.element("data", search())) detach("data"))
    if (debug) cat("Program is running ")
    for (em in 1:emmax) {
        if (debug) cat(".")
        logisfit <- eval(parse(text = paste(model, "(", logisform,
            ", family = binomial(link='", link, "'), ...)", sep = "")))
#       assign("uncureprob", uncureprob, frame = 1, immediate = T)
        assign("uncureprob", uncureprob)
        coxfit <- eval(parse(text = paste("coxph(", coxformula,
            ", subset = uncureprob != 0, method = 'breslow')")))
        survprob <- sapply(survival(coxfit, newdata =
            data.frame(uncureprob = rep(1, datalen), data), tail = tail,
            type = "kaplan-meier"), function(x)x$surv)
        coef2 <- c(coef(coxfit), coef(logisfit))
        if(em > 1) {
            test <- max(abs(1.0-coef1/coef2))
            if (is.na(test) || test <= eps) break
            else coef1 <- coef2
        }
        else coef1 <- coef2
        puncure <- fitted(logisfit)
        uncureprob <- ifelse(cens==1, 1,
            puncure*survprob/(1-puncure+puncure*survprob))
    }
    if (debug) cat(" done.\n")
    structure(list(call = call, coxfit = coxfit, logisfit = logisfit,
        tail = tail, data = if (savedata && !nodata)
        data.frame(data, uncureprob = uncureprob) else NULL),
        class = "semicure")
}

coef.semicure <- function(x)
{
    coef <- c(coef(x$coxfit), coef(x$logisfit))

    if(any(nas <- is.na(coef))) {
        if(is.null(names(coef)))
            names(coef) <- paste("b", 1:length(coef), sep = "")
            cat("\nCoefficients: (", sum(nas), 
                " not defined because of singularities)\n", sep = "")
    }
    coef
}

print.semicure <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    
    cat("\nFailure time distribution model:\n")
    print(x$coxfit, ...)

    cat("\nCure probability model:\n")
    print(x$logisfit, ...)
    invisible(x)
}

summary.semicure <- function(x)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    
    if (!is.null(x$data)) {
        cat("\nSummary of estimated uncured probabilities:\n")
        print(summary(x$data$uncureprob))
    }

    cat("\nSummary of failure time distribution model:\n")
    print(summary(x$coxfit))

    cat("\nSummary of cure probability model:\n")
    print(summary(x$logisfit))

    invisible(x)
}

predict.semicure <- function(object, newdata, newtime = NULL, overall = F,
    cumhaz = F, ...)
{
    call <- match.call()
    if(!inherits(object, "semicure")) stop("Object must be results of semicure")

    if (is.null(object$data))
        stop("Object must be fitted with data arguments")
    attach(object$data)
    on.exit(if (is.element("object$data", search()))
            detach("object$data"))

    pred <- survival(object$coxfit, newdata = data.frame(uncureprob =
        rep(1, nrow(newdata)), newdata), newtime = newtime, ...)
    cure <- 1 - predict.glm(object$logisfit, newdata = newdata,
            type = "response")

    for (i in 1:length(pred)) {
        pred[[i]]$curerate <- cure[i]
        if (overall)
            pred[[i]]$surv <- pred[[i]]$surv*(1 - pred[[i]]$curerate) + pred[[i]]$curerate
        if (cumhaz) pred[[i]]$surv <- -log(pred[[i]]$surv)
    }
    
    structure(list(call = call, cumhaz = cumhaz, overall = overall,
        prediction = pred), class = "predict.semicure")
}

print.predict.semicure <- function(x, digits = 3, rdigits = NULL)
{
    if(is.null(digits))
        digits <- options()$digits
    else {
        old.digits <- options(digits = digits)
        on.exit(options(old.digits))
    }
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
    }
    
    cat("\nPredicted ")

    if (x$overall) cat("overall ")
    else cat("uncured ")

    if (x$cumhaz) cat("cumulative hazard ")
    else cat("survival probabilities ")

    cat("and cure rates\n")

    x <- do.call("cbind", lapply(x$prediction, function(x)
        cbind(time = c(NA, NA, NA, x$time),
            surv = c(x$curerate, x$median, NA, x$surv))))
    px <- apply(x, 2, function(y, digits, rdigits)
    {
        yna <- is.na(y)
        tmpy <- y
        if(is.null(rdigits)) {
            y[2] <- format(zapsmall(tmpy[2], digits))
            y[-2] <- format(zapsmall(tmpy[-2], digits))
        }
        else {
            y[2] <- format(round(zapsmall(tmpy[2], digits), digits = rdigits))
            y[-2] <- format(round(zapsmall(tmpy[-2], digits), digits = rdigits))
        }
        y[yna] <- ""
        y
    }
    , digits, rdigits)
    px[1, seq(1, ncol(px), 2)] <- "cure rate"
    px[2, seq(1, ncol(px), 2)] <- "median"
    px[3, ] <- dimnames(px)[[2]]
    px <- data.frame(px, row.names = rep("", nrow(px)), dup.row.names = T)
    names(px)[seq(1, ncol(px), 2)] <- as.character(1:(ncol(px)/2))
    names(px)[seq(2, ncol(px), 2)] <- ""

    print.data.frame(px, digits = digits)

    invisible(x)
}

plot.predict.semicure <- function(object, type = "S",
    which = 1:length(object$curerate), xlab, ylab, ...)
{
    if (missing(xlab)) xlab <- "Time"
    if (missing(ylab)) ylab <- "Survival Probability"

    if (length(object$curerate)==1)
        plot(object$time, object$surv, type = type, xlab = xlab, ylab = ylab,
            ...)
    else 
        matplot(object$time, object$surv[, which], type = type, xlab = xlab,
            ylab = ylab, ...)
    invisible()
}

lines.predict.semicure <- function(object, type = "S",
    which = 1:length(object$curerate), ...)
{
    if (length(object$curerate)==1)
        lines(object$time, object$surv, type = type, ...)
    else 
        matlines(object$time, object$surv[, which], type = type, ...)
    invisible()
}

survival <- function(object, newdata, newtime = NULL,
    type = c("kaplan-meier", "aalen"), tail = c("none", "zero", "NA"))
{
    if(!inherits(object, "coxph")) stop("Object must be results of coxph.")

    tail <- match.arg(tail)
    type <- match.arg(type)
    nullcox <- inherits(object, "coxph.null")
    if (nullcox) sfit <- survfit(object, type = type)
    else sfit <- survfit(object, newdata = newdata, type = type)
    Terms <- terms(object)
    m <- model.frame(Terms, data = newdata)
    n <- nrow(newdata)
    strata <- NULL
    if(!is.null(sfit$strata)) {
        strata <- m[[untangle.specials(Terms, "strata", 1)$vars]]
        strata <- charmatch(as.character(strata), names(sfit$strata))
    }
    if (is.null(newtime))
        newtime <- as.list(model.extract(m, "response")[, 1])
    else {
        newtime <- sort(newtime)
        newtime <- eval(parse(text = paste("list(", paste(rep("newtime", n),
            collapse = ", "), ")", sep = "")))
    }

    oldtime <- object$y[, 1]
    small <- min(diff(sort(unique(c(0, oldtime)))))/2
    surv <- list()
    for (i in 1:n) {
        fit <- if (is.null(strata)) sfit[i]
            else sfit[strata[i], i]

        ti <- newtime[[i]]
        maxtime <- max(fit$time)
        large <- ti > maxtime
        ti[!large] <- ti[!large]-small
        ti[large] <- maxtime

        prob <- summary(fit, times = ti)$surv
        if (tail=="zero") prob[large] <- 0
        if (tail=="NA") prob[large] <- NA

##      med <- if (n > 1)
##              print.survfit.computations(fit)[1, "median"]
##          else print.survfit.computations(fit)["median"]
##      names(med) <- NULL
med <- NULL

        surv[[i]] <- list(time = newtime[[i]], surv = prob, median = med)
    }
    surv
}
