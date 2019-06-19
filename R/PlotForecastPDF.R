#'Plot one or multiple ensemble forecast pdfs for the same event
#'
#'@author Llorenç Lledó \email{llledo@bsc.es}
#'@description This function plots the probability distribution function of several ensemble forecasts for the same event, either initialized at different moments or by different models. Probabilities for tercile categories are computed, plotted in colors and annotated. An asterisk marks the tercile with higher probabilities. Probabilities for extreme categories (above P90 and below P10) can also be included as hatched areas. Individual ensemble members can be plotted as jittered points. The observed value is optionally shown as a diamond.
#'
#'@param fcst a dataframe or array containing all the ensember members for each frecast. If \code{'fcst'} is an array, it should have two labelled dimensions, and one of them should be \code{'members'}. If \code{'fcsts'} is a data.frame, each column shoul be a separate forecast, with the rows beeing the different ensemble members.
#'@param tercile.limits an array with P33 and P66 values that define the tercile categories.
#'@param extreme.limits (optional) an array with P10 and P90 values that define the extreme categories. (Default: extreme categories are not shown).
#'@param obs (optional) a number with the observed value. (Default: observation is not shown).
#'@param plotfile (optional) a filename (pdf, png...) where the plot will be saved. (Default: the plot is not saved).
#'@param title a string with the plot title.
#'@param var.name a string with the variable name and units.
#'@param fcst.names (optional) an array of strings with the titles of each individual forecast.
#'@param add.ensmemb either to add the ensemble members \code{'above'} (default) or \code{'below'} the pdf, or not (\code{'no'}).
#'@param color.set a selection of predefined color sets: use \code{'ggplot'} (default) for blue/green/red, \code{'s2s4e'} for blue/grey/orange, or \code{'hydro'} for yellow/gray/blue (suitable for precipitation and inflows).
#'
#'@return a ggplot object containing the plot.
#'
#'@importFrom data.table data.table
#'@importFrom data.table CJ
#'@importFrom data.table setkey
#'@import ggplot2
#'@importFrom reshape2 melt
#'@importFrom plyr .
#'@importFrom plyr dlply
#'
#'@examples
#'fcsts <- data.frame(fcst1 = rnorm(10), fcst2 = rnorm(10, 0.5, 1.2), 
#'                    fcst3 = rnorm(10, -0.5, 0.9))
#'PlotForecastPDF(fcsts,c(-1,1))
#'\donttest{
#'fcsts2 <- array(rnorm(100), dim = c(members = 20, fcst = 5))
#'PlotForecastPDF(fcsts2, c(-0.66, 0.66), extreme.limits = c(-1.2, 1.2), 
#'                fcst.names = paste0('random fcst ', 1 : 5), obs = 0.7)
#'}
#'@export
PlotForecastPDF <- function(fcst, tercile.limits, extreme.limits = NULL, obs = NULL, 
    plotfile = NULL, title = "Set a title", var.name = "Varname (units)", fcst.names = NULL, 
    add.ensmemb = c("above", "below", "no"), color.set = c("ggplot", "s2s4e", "hydro")) {
    value <- init <- extremes <- x <- ymin <- ymax <- tercile <- y <- xend <- yend <- yjitter <- MLT <- NULL
    ggColorHue <- function(n) {
      hues <- seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    #------------------------
    # Define color sets
    #------------------------
    color.set <- match.arg(color.set)
    if (color.set == "s2s4e") {
        colorFill <- rev(c("#FF764D", "#b5b5b5", "#33BFD1"))
        colorHatch <- c("deepskyblue3", "indianred3")
        colorMember <- c("#ffff7f")
        colorObs <- "purple"
        colorLab <- c("blue", "red")
    } else if (color.set == "hydro") {
        colorFill <- rev(c("#41CBC9", "#b5b5b5", "#FFAB38"))
        colorHatch <- c("darkorange1", "deepskyblue3")
        colorMember <- c("#ffff7f")
        colorObs <- "purple"
        colorLab <- c("darkorange3", "blue")
    } else if (color.set == "ggplot") {
        colorFill <- rev(ggColorHue(3))
        colorHatch <- c("deepskyblue3", "indianred1")
        colorMember <- c("#ffff7f")
        colorObs <- "purple"
        colorLab <- c("blue", "red")
    } else {
        stop("Parameter 'color.set' should be one of ggplot/s2s4e/hydro")
    }
    #------------------------
    # Check input arguments
    #------------------------
    add.ensmemb <- match.arg(add.ensmemb)
    if (length(tercile.limits) != 2) {
        stop("Parameter 'tercile.limits' should be an array with two limits for delimiting tercile categories")
    }
    if (tercile.limits[1] >= tercile.limits[2]) {
        stop("The provided tercile limits are in the wrong order")
    }
    if (!is.null(extreme.limits)) {
        if (length(extreme.limits) != 2) {
            stop("Parameter 'extreme.limits' should be an array with two limits for delimiting extreme categories")
        }
        if (extreme.limits[1] >= tercile.limits[1] | extreme.limits[2] <= tercile.limits[2]) {
            stop("The provided extreme limits are not consistent with tercile limits")
        }
    }
    #------------------------
    # Check fcst type and convert to data.frame if needed
    #------------------------
    if (is.array(fcst)) {
        if (!"members" %in% names(dim(fcst)) | length(dim(fcst)) != 2) {
            stop("Parameter 'fcst' should be a two-dimensional array with labelled dimensions and one of them should be 'members'")
        }
        dim.members <- which(names(dim(fcst)) == "members")
        if (dim.members == 1) {
            fcst.df <- data.frame(fcst)
        } else {
            fcst.df <- data.frame(t(fcst))
        }
    } else if (is.data.frame(fcst)) {
        fcst.df <- fcst
    } else {
        stop("Parameter 'fcst' should be an array or a data.frame")
    }
    #------------------------
    # Set proper fcst names
    #------------------------
    if (!is.null(fcst.names)) {
        if (length(fcst.names) != dim(fcst.df)[2]) {
            stop("Parameter 'fcst.names' should be an array with as many names as distinct forecasts provided")
        }
        colnames(fcst.df) <- factor(fcst.names, levels = fcst.names)
    }
    #------------------------
    # Produce a first plot with the pdf for each init in a panel
    #------------------------
    melt.df <- melt(fcst.df, variable.name = "init", id.vars = NULL)
    plot <- ggplot(melt.df, aes(x = value)) + geom_density(alpha = 1, na.rm = T) + 
        coord_flip() + facet_wrap(~init, strip.position = "top", nrow = 1) + xlim(range(c(obs, 
        density(melt.df$value, na.rm = T)$x)))
    ggp <- ggplot_build(plot)
    #------------------------
    # Gather the coordinates of the plots together with init and corresponding
    # terciles
    #------------------------
    tmp.df <- ggp$data[[1]][, c("x", "ymin", "ymax", "PANEL")]
    if (!is.null(ggp$layout$layout)) {
        tmp.df$init <- ggp$layout$layout$init[as.numeric(tmp.df$PANEL)]
    } else if (!is.null(ggp$layout$panel_layout)) {
        tmp.df$init <- ggp$layout$panel_layout$init[as.numeric(tmp.df$PANEL)]
    } else {
        stop("Cannot find PANELS in ggp object")
    }
    tmp.df$tercile <- factor(ifelse(tmp.df$x < tercile.limits[1], "Below normal", 
        ifelse(tmp.df$x < tercile.limits[2], "Normal", "Above normal")), levels = c("Below normal", 
        "Normal", "Above normal"))
    #------------------------
    # Get the height and width of a panel
    #------------------------
    pan.width <- diff(range(tmp.df$x))
    pan.height <- max(tmp.df$ymax)
    magic.ratio <- 9 * pan.height/pan.width
    #------------------------
    # Compute hatch coordinates for extremes
    #------------------------
    if (!is.null(extreme.limits)) {
        tmp.df$extremes <- factor(ifelse(tmp.df$x < extreme.limits[1], "Below P10", 
            ifelse(tmp.df$x < extreme.limits[2], "Normal", "Above P90")), levels = c("Below P10", 
            "Normal", "Above P90"))
        hatch.ls <- dlply(tmp.df, .(init, extremes), function(x) {
            # close the polygon
            tmp.df2 <- data.frame(x = c(x$x, max(x$x), min(x$x)), y = c(x$ymax, 0, 
                0))
            # compute the hatches for this polygon
            hatches <- .polygon.fullhatch(tmp.df2$x, tmp.df2$y, angle = 60, density = 10, 
                width.units = pan.width, height.units = pan.height)
            # add bottom segment
            end1 <- data.frame(x = x$x[1], y = x$ymax[1], xend = x$x[1], yend = 0)
            # add top segment
            end2 <- data.frame(x = x$x[length(x$x)], y = x$ymax[length(x$x)], xend = x$x[length(x$x)], 
                yend = 0)
            return(rbind(hatches, end1, end2))
        })
        attr <- attr(hatch.ls, "split_labels")
        for (i in 1:length(hatch.ls)) {
            hatch.ls[[i]] <- cbind(hatch.ls[[i]], attr[i, ])
        }
        hatch.df <- do.call("rbind", hatch.ls)
        # Compute max y for each extreme category
        max.ls <- dlply(tmp.df, .(init, extremes), function(x) {
            data.frame(y = min(0.85 * pan.height, max(x$ymax)))
        })
        attr <- attr(max.ls, "split_labels")
        for (i in 1:length(max.ls)) {
            max.ls[[i]] <- cbind(max.ls[[i]], attr[i, ])
        }
        max.df <- do.call("rbind", max.ls)
    }
    #------------------------
    # Compute jitter space for ensemble members
    #------------------------
    if (add.ensmemb != "no") {
        jitter.df <- melt(data.frame(dlply(melt.df, .(init), function(x) {
            .jitter.ensmemb(sort(x$value, na.last = T), pan.width/100)
        }), check.names = F), value.name = "yjitter", variable.name = "init", id.vars = NULL)
        jitter.df$x <- melt(data.frame(dlply(melt.df, .(init), function(x) {
            sort(x$value, na.last = T)
        })), value.name = "x", id.vars = NULL)$x
    }
    #------------------------
    # Get y coordinates for observed x values, using a cool data.table feature: merge
    # to nearest value
    #------------------------
    if (!is.null(obs)) {
        tmp.dt <- data.table(tmp.df, key = c("init", "x"))
        obs.dt <- data.table(init = factor(colnames(fcst.df), levels = colnames(fcst.df)), 
            value = rep(obs, dim(fcst.df)[2]))
        setkey(obs.dt, init, value)
        obs.xy <- tmp.dt[obs.dt, roll = "nearest"]
    }
    #------------------------
    # Fill each pdf with different colors for the terciles
    #------------------------
    plot <- plot + geom_ribbon(data = tmp.df, aes(x = x, ymin = ymin, ymax = ymax, 
        fill = tercile), alpha = 0.7)
    #------------------------
    # Add hatches for extremes
    #------------------------
    if (!is.null(extreme.limits)) {
        if (nrow(hatch.df[hatch.df$extremes != "Normal", ]) == 0) {
            warning("The provided extreme categories are outside the plot bounds. The extremes will not be drawn.")
            extreme.limits <- NULL
        } else {
            plot <- plot + geom_segment(data = hatch.df[hatch.df$extremes != "Normal", 
                ], aes(x = x, y = y, xend = xend, yend = yend, color = extremes))
        }
    }
    #------------------------
    # Add obs line
    #------------------------
    if (!is.null(obs)) {
        plot <- plot + geom_vline(data = obs.dt, aes(xintercept = value), linetype = "dashed", 
            color = colorObs)
    }
    #------------------------
    # Add ensemble members
    #------------------------
    if (add.ensmemb == "below") {
        plot <- plot + # this adds a grey box for ensmembers
        geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -pan.height/10), 
            fill = "gray95", color = "black", width = 0.2) + # this adds the ensemble members
        geom_point(data = jitter.df, color = "black", fill = colorMember, alpha = 1, 
            aes(x = x, y = -pan.height/10 - magic.ratio * yjitter, shape = "Ensemble members"))
    } else if (add.ensmemb == "above") {
        plot <- plot + geom_point(data = jitter.df, color = "black", fill = colorMember, 
            alpha = 1, aes(x = x, y = 0.7 * magic.ratio * yjitter, shape = "Ensemble members"))
    }
    #------------------------
    # Add obs diamond
    #------------------------
    if (!is.null(obs)) {
        plot <- plot + # this adds the obs diamond
        geom_point(data = obs.xy, aes(x = x, y = ymax, size = "Observation"), shape = 23, 
            color = "black", fill = colorObs)
    }
    #------------------------
    # Compute probability for each tercile and identify MLT
    #------------------------
    tmp.dt <- data.table(tmp.df)
    pct <- tmp.dt[, .(pct = integrate(approxfun(x, ymax), lower = min(x), upper = max(x))$value), 
        by = .(init, tercile)]
    tot <- pct[, .(tot = sum(pct)), by = init]
    pct <- merge(pct, tot, by = "init")
    pct$pct <- round(100 * pct$pct/pct$tot, 0)
    pct$MLT <- pct[, .(MLT = pct == max(pct)), by = init]$MLT
    lab.pos <- c(tercile.limits[1], mean(tercile.limits), tercile.limits[2])
    #------------------------
    # Compute probability for extremes
    #------------------------
    if (!is.null(extreme.limits)) {
        pct2 <- tmp.dt[, .(pct = integrate(approxfun(x, ymax), lower = min(x), upper = max(x))$value), 
            by = .(init, extremes)]
        tot2 <- pct2[, .(tot = sum(pct)), by = init]
        pct2 <- merge(pct2, tot2, by = "init")
        pct2$pct <- round(100 * pct2$pct/pct2$tot, 0)
        lab.pos2 <- c(extreme.limits[1], NA, extreme.limits[2])
        pct2 <- merge(pct2, max.df, by = c("init", "extremes"))
        # include potentially missing groups
        pct2 <- pct2[CJ(levels(pct2$init), factor(c("Below P10", "Normal", "Above P90"), 
            levels = c("Below P10", "Normal", "Above P90"))), ]
    }
    #------------------------
    # Add probability labels for terciles
    #------------------------
    if (add.ensmemb == "above") {
        labpos <- -0.2 * pan.height
        vjust <- 0
    } else {
        labpos <- 0
        vjust <- -0.5
    }
    plot <- plot + geom_text(data = pct, aes(x = lab.pos[as.integer(tercile)], y = labpos, 
        label = paste0(pct, "%"), hjust = as.integer(tercile) * 1.5 - 2.5), vjust = vjust, 
        angle = -90, size = 3.2) + geom_text(data = pct[MLT == T, ], aes(x = lab.pos[as.integer(tercile)], 
        y = labpos, label = "*", hjust = as.integer(tercile) * 3.5 - 5), vjust = 0.1, 
        angle = -90, size = 7, color = "black")
    #------------------------
    # Add probability labels for extremes
    #------------------------
    if (!is.null(extreme.limits)) {
        plot <- plot + geom_text(data = pct2[extremes != "Normal", ], aes(x = lab.pos2[as.integer(extremes)], 
            y = 0.9 * y, label = paste0(pct, "%"), hjust = as.integer(extremes) * 
                1.5 - 2.5), vjust = -0.5, angle = -90, size = 3.2, color = rep(colorLab, 
            dim(fcst.df)[2]))
    }
    #------------------------
    # Finish all theme and legend details
    #------------------------
    plot <- plot + theme_minimal() + scale_fill_manual(name = "Probability of\nterciles", 
        breaks = c("Above normal", "Normal", "Below normal"), values = colorFill, 
        drop = F) + scale_color_manual(name = "Probability of\nextremes", values = colorHatch) + 
        scale_shape_manual(name = "Ensemble\nmembers", values = c(21)) + scale_size_manual(name = "Observation", 
        values = c(3)) + labs(x = var.name, y = "Probability density\n(total area=1)", 
        title = title) + theme(axis.text.x = element_blank(), panel.grid.minor.x = element_blank(), 
        legend.key.size = unit(0.3, "in"), panel.border = element_rect(fill = NA, 
            color = "black"), strip.background = element_rect(colour = "black", fill = "gray80"), 
        panel.spacing = unit(0.2, "in"), panel.grid.major.x = element_line(color = "grey93")) + 
        guides(fill = guide_legend(order = 1), color = guide_legend(order = 2, reverse = T), 
            shape = guide_legend(order = 3, label = F), size = guide_legend(order = 4, 
                label = F))
    #------------------------
    # Save to plotfile if needed, and return plot
    #------------------------
    if (!is.null(plotfile)) {
        ggsave(plotfile, plot)
    }
    return(plot)
}
.jitter.ensmemb <- function(x, thr = 0.1) {
    # Idea: start with first level. Loop all points, and if distance to last point in
    # the level is more than a threshold, include the point to the level.  Otherwise
    # keep the point for another round.  Do one round in each direction to avoid
    # uggly patterns.
    if (is.unsorted(x, na.rm = T)) {
        stop("Provide a sorted array!")
    }
    lev <- x * 0
    lev[is.na(lev)] <- -1
    level <- 1
    while (any(lev == 0)) {
        last <- -1/0
        for (i in 1:length(x)) {
            if (lev[i] != 0) {
                next
            }
            if (x[i] - last > thr) {
                lev[i] <- level
                last <- x[i]
            }
        }
        level <- level + 1
        last <- 1/0
        for (i in seq(length(x), 1, -1)) {
            if (lev[i] != 0) {
                next
            }
            if (last - x[i] > thr) {
                lev[i] <- level
                last <- x[i]
            }
        }
        level <- level + 1
    }
    lev[lev == -1] <- NA
    return(lev * thr * sqrt(3)/2)
}
.polygon.onehatch <- function(x, y, x0, y0, xd, yd, fillOddEven = F) {
    halfplane <- as.integer(xd * (y - y0) - yd * (x - x0) <= 0)
    cross <- halfplane[-1L] - halfplane[-length(halfplane)]
    does.cross <- cross != 0
    if (!any(does.cross)) {
        return()
    }
    x1 <- x[-length(x)][does.cross]
    y1 <- y[-length(y)][does.cross]
    x2 <- x[-1L][does.cross]
    y2 <- y[-1L][does.cross]
    t <- (((x1 - x0) * (y2 - y1) - (y1 - y0) * (x2 - x1))/(xd * (y2 - y1) - yd * 
        (x2 - x1)))
    o <- order(t)
    tsort <- t[o]
    crossings <- cumsum(cross[does.cross][o])
    if (fillOddEven) {
        crossings <- crossings%%2
    }
    drawline <- crossings != 0
    lx <- x0 + xd * tsort
    ly <- y0 + yd * tsort
    lx1 <- lx[-length(lx)][drawline]
    ly1 <- ly[-length(ly)][drawline]
    lx2 <- lx[-1L][drawline]
    ly2 <- ly[-1L][drawline]
    return(data.frame(x = lx1, y = ly1, xend = lx2, yend = ly2))
}
.polygon.fullhatch <- function(x, y, density, angle, width.units, height.units, inches = c(5, 
    1)) {
    x <- c(x, x[1L])
    y <- c(y, y[1L])
    angle <- angle%%180
    upi <- c(width.units, height.units)/inches
    if (upi[1L] < 0) {
        angle <- 180 - angle
    }
    if (upi[2L] < 0) {
        angle <- 180 - angle
    }
    upi <- abs(upi)
    xd <- cos(angle/180 * pi) * upi[1L]
    yd <- sin(angle/180 * pi) * upi[2L]
    hatch.ls <- list()
    i <- 1
    if (angle < 45 || angle > 135) {
        if (angle < 45) {
            first.x <- max(x)
            last.x <- min(x)
        } else {
            first.x <- min(x)
            last.x <- max(x)
        }
        y.shift <- upi[2L]/density/abs(cos(angle/180 * pi))
        x0 <- 0
        y0 <- floor((min(y) - first.x * yd/xd)/y.shift) * y.shift
        y.end <- max(y) - last.x * yd/xd
        while (y0 < y.end) {
            hatch.ls[[i]] <- .polygon.onehatch(x, y, x0, y0, xd, yd)
            i <- i + 1
            y0 <- y0 + y.shift
        }
    } else {
        if (angle < 90) {
            first.y <- max(y)
            last.y <- min(y)
        } else {
            first.y <- min(y)
            last.y <- max(y)
        }
        x.shift <- upi[1L]/density/abs(sin(angle/180 * pi))
        x0 <- floor((min(x) - first.y * xd/yd)/x.shift) * x.shift
        y0 <- 0
        x.end <- max(x) - last.y * xd/yd
        while (x0 < x.end) {
            hatch.ls[[i]] <- .polygon.onehatch(x, y, x0, y0, xd, yd)
            i <- i + 1
            x0 <- x0 + x.shift
        }
    }
    return(do.call("rbind", hatch.ls))
}