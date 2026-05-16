## TEATIME embedded MAGOS = server-compatible OLD algorithm + fast mag.single
##  - All non-mag.single helpers verbatim from server's MAGOS install
##  - mag.single replaced with vectorised split+vapply version

cut.off.multiple <-
function (mag.var) 
{
    time0 <- as.numeric(Sys.time())
    mag.var.list <- mag.var$var.mean.all
    depths.data <- mag.var$mag.out$prep.data$depths
    vafs.data <- mag.var$mag.out$prep.data$vafs
    for (i in 1:length(mag.var.list)) {
        temp <- mag.var.list[[i]]
        temp <- data.frame(temp, row.names = NULL)
        temp[is.na(temp)] <- 0
        temp$var <- round(temp$var, 6)
        temp$var <- temp$var + 9.9999999999999995e-08
        mag.var.list[[i]] <- temp
    }
    names(mag.var.list) <- names(mag.var$var.mean.all)
    x.el.list <- mag.var$mag.out$x.el
    flag <- FALSE
    m <- c()
    str <- c()
    str.f <- c()
    cutstr <- c()
    cutstr2 <- c()
    i <- max(mag.var.list[[1]]$step)
    allPoints <- x.el.list[[i]]
    ok.points <- allPoints - allPoints
    n <- length(allPoints)
    id.ok <- c()
    temp.mv <- c()
    temp.points <- c()
    temp.step <- c()
    temp.point.step <- c()
    while (flag == FALSE) {
        v <- i
        stepI.x.el <- x.el.list[[i]]
        stepI.x.el.step <- stepI.x.el[nrow(stepI.x.el), ]
        stepII.x.el <- x.el.list[[i - 1]]
        stepII.breaks.x.el <- dplyr::setdiff(as.data.frame(stepII.x.el), 
            as.data.frame(stepI.x.el))
        temp.points <- rbind(temp.points, c(i, ok.points))
        temp.step <- rbind(temp.step, c(i, stepI.x.el.step))
        temp.point.step <- rbind(temp.point.step, c(i, sum(ok.points & 
            stepI.x.el.step)))
        keep <- c(0, 0)
        for (j in 1:length(mag.var.list)) {
            var.mean <- mag.var.list[[j]]
            stepI <- var.mean[var.mean$step == i, ]
            stepII <- var.mean[var.mean$step == i - 1, ]
            if (sum(ok.points & stepI.x.el.step) == 0) {
                stepII.breaks <- dplyr::setdiff(stepII[, -4], 
                  stepI[, -4])
                if (nrow(stepII.breaks) == 1) {
                  stepII.breaks <- rbind(stepII.breaks, stepII.breaks)
                }
                stepII.num <- unique(stepII[, 4])
                stepII.breaks <- cbind(stepII.breaks, step = stepII.num)
                for (ii in 1:2) {
                  var.sim <- c()
                  size <- ifelse(stepII.breaks$NumberOfPoints[ii] < 
                    5, 5, stepII.breaks$NumberOfPoints[ii])
                  depth.vec <- depths.data[depths.data$ID %in% 
                    names(stepII.breaks.x.el[ii, as.logical(stepII.breaks.x.el[ii, 
                      ]), drop = F]), j]
                  vaf.vec <- vafs.data[vafs.data$ID %in% names(stepII.breaks.x.el[ii, 
                    as.logical(stepII.breaks.x.el[ii, ]), drop = F]), 
                    j]
                  if (mean(depth.vec) > 150) {
                    type <- "H"
                  }
                  if (mean(depth.vec) <= 150) {
                    type <- "M"
                  }
                  if (mean(depth.vec) <= 40) {
                    type <- "L"
                  }
                  var.sim <- mag.exp.var.v3(efrq.vec = vaf.vec, 
                    edep.vec = depth.vec, num = size, type = type)
                  var.sim <- var.sim$exp.vat
                  m <- mean(var.sim) + 9.9999999999999995e-07
                  vv <- sd(var.sim)
                  vmax <- max(var.sim)
                  if (length(vaf.vec) > 50 & max(vaf.vec) - min(vaf.vec) > 
                    0.050000000000000003) {
                    var.check <- var(vaf.vec[-which.max(vaf.vec)[1]])
                  }
                  else {
                    var.check <- stepII.breaks[ii, 2]
                  }
                  var.threshold <- m + 2 * vv
                  if (type == "H") {
                    var.threshold <- m + 3 * vv
                  }
                  var.threshold <- m + 3 * vv
                  if (var.threshold > 0.001) {
                    var.threshold <- round(var.threshold, 3)
                  }
                  if (stepII.breaks[ii, "okay"] > 0 | floor(var.check * 
                    10^4)/10^4 < m + 3 * vv) {
                    keep[ii] <- keep[ii] + 1
                  }
                }
            }
        }
        for (ii in 1:2) {
            if (keep[ii] == length(mag.var.list)) {
                cutstr2 <- rbind(cutstr2, stepII.breaks.x.el[ii, 
                  ])
                str <- rbind(str, stepII.breaks[ii, ])
                ok.points <- ok.points + stepII.breaks.x.el[ii, 
                  ]
            }
        }
        if (sum(ok.points > 1) > 0) {
            print("ERROR")
        }
        if (sum(ok.points) == n) {
            flag <- TRUE
        }
        else {
            i <- (v - 1)
        }
    }
    time1 <- as.numeric(Sys.time())
    # cat("took ", (time1 - time0), " seconds.\n")  # silenced
    colors <- colSums(cutstr2 * c(1:nrow(cutstr2))) + 1
    data <- mag.var$mag.out$vaf.sorted
    final.data <- cbind(data, colors)
    t = final.data
    t.1 = c()
    for (i in 1:(ncol(t) - 2)) {
        temp = c()
        for (j in unique(t$colors)) {
            m.1 = mean(t[t$colors == j, i])
            v.1 = var(t[t$colors == j, i])
            a <- ((1 - m.1)/v.1 - 1/m.1) * m.1^2
            b <- a * (1/m.1 - 1)
            temp = rbind(temp, c(j, a, b))
        }
        colnames(temp) = c("colors", paste0("alpha.", i), paste0("beta.", 
            i))
        t.1 = cbind(t.1, temp[, -1])
    }
    t.1 = cbind(unique(t$colors), t.1)
    colnames(t.1)[1] = "colors"
    i = 1
    prob.all = c()
    for (i in 1:nrow(t)) {
        vf = t[i, 1:(ncol(t) - 2)]
        prob.s = c()
        for (j in 1:length(vf)) {
            probs = c()
            prob = c()
            for (ii in unique(t$colors)) {
                prob = c(prob, dbeta(vf[, j], shape1 = t.1[t.1[, 
                  1] == ii, 1 + (2 * j - 1)], shape2 = t.1[t.1[, 
                  1] == ii, 1 + 2 * j]))
            }
            prob.s = rbind(prob.s, prob)
            prob.s = round(apply(prob.s, 2, min), 2)
        }
        prob.all = rbind(prob.all, prob.s)
    }
    colnames(prob.all) = paste0("colors.", unique(t$colors))
    rownames(prob.all) = NULL
    final.probs = cbind(t, prob.all)
    results <- list(str = str, cutclust = cutstr, cutclust2 = cutstr2, 
        colors = colors, probs = final.probs, final.data = final.data)
    return(results)
}
cut.off.single <-
function (mag.output) 
{
    time0 <- as.numeric(Sys.time())
    var.mean <- mag.output$var.mean
    head(var.mean)
    var.mean <- data.frame(var.mean, row.names = NULL)
    var.mean[is.na(var.mean)] <- 0
    var.mean$var <- round(var.mean$var, 6)
    var.mean$var <- var.mean$var + 9.9999999999999995e-08
    flag <- FALSE
    m <- c()
    str <- c()
    str.f <- c()
    cutstr <- c()
    cutstr2 <- c()
    i <- max(var.mean$step)
    if (i == 1) {
        flag = T
        str <- var.mean
    }
    allPoints <- mag.output$x.el[[i]]
    allPoints <- as.data.frame(as.matrix(allPoints))
    ok.points <- allPoints - allPoints
    n <- length(allPoints)
    id.ok <- c()
    temp.mv <- c()
    temp.points <- c()
    temp.step <- c()
    temp.point.step <- c()
    var.sim.zero <- mag.exp.var.v3(efrq.vec = mag.output$vaf.sorted$vaf.1, 
        edep.vec = mag.output$depth.sorted$depth.1, num = n)
    var.sim.zero <- var.sim.zero$exp.vat
    m <- mean(var.sim.zero, na.rm = T) + 9.9999999999999995e-07
    vv <- sd(var.sim.zero, na.rm = T)
    if (var(mag.output$vaf.sorted$vaf.1) < m + 3 * vv) {
        flag <- TRUE
        print("Only one clone detected")
        str <- var.mean[var.mean$step == i, ]
        colors <- rep(1, n)
    }
    {
        while (flag == FALSE) {
            v <- i
            stepI <- var.mean[var.mean$step == i, ]
            stepI.x.el <- mag.output$x.el[[i]]
            stepI.x.el <- as.data.frame(as.matrix(stepI.x.el))
            stepI.x.el.step <- stepI.x.el[nrow(stepI.x.el), ]
            stepII <- var.mean[var.mean$step == i - 1, ]
            if (i > 1) {
                stepII.x.el <- mag.output$x.el[[i - 1]]
                stepII.x.el <- as.data.frame(as.matrix(stepII.x.el))
            }
            if (i == 1) {
                stepII.x.el <- mag.output$x.el.preprocess
                stepII.x.el <- as.data.frame(as.matrix(stepII.x.el))
            }
            temp.points <- rbind(temp.points, c(i, ok.points))
            temp.step <- rbind(temp.step, c(i, stepI.x.el.step))
            temp.point.step <- rbind(temp.point.step, c(i, sum(ok.points & 
                stepI.x.el.step)))
            if (sum(ok.points & stepI.x.el.step) == 0) {
                stepII.breaks <- dplyr::setdiff(stepII[, -4], 
                  stepI[, -4])
                if (nrow(stepII.breaks) == 1) {
                  stepII.breaks <- rbind(stepII.breaks, stepII.breaks)
                }
                stepII.num <- stepII[, 4]
                stepII.breaks <- cbind(stepII.breaks, stepII[1:2, 
                  4])
                stepII.breaks.x.el <- dplyr::setdiff(stepII.x.el, 
                  stepI.x.el)
                for (ii in 1:2) {
                  var.sim <- c()
                  size <- ifelse(stepII.breaks$NumberOfPoints[ii] < 
                    5, 5, stepII.breaks$NumberOfPoints[ii])
                  vafs <- mag.output$vaf.sorted$vaf.1[as.matrix(stepII.breaks.x.el[ii, 
                    ])]
                  depths <- mag.output$depth.sorted$depth.1[as.matrix(stepII.breaks.x.el[ii, 
                    ])]
                  var.sim <- mag.exp.var.v3(efrq.vec = (vafs), 
                    edep.vec = (depths), num = size)
                  var.sim <- var.sim$exp.vat
                  if (length(vafs) > 20 & max(vafs) - min(vafs) > 
                    0.050000000000000003) {
                    var.check <- var(vafs[-which.max(vafs)[1]])
                  }
                  else {
                    var.check <- stepII.breaks[ii, 2]
                  }
                  m <- mean(var.sim, na.rm = T) + 9.9999999999999995e-07
                  vv <- sd(var.sim, na.rm = T)
                  vmax <- max(var.sim, na.rm = T)
                  if (stepII.breaks[ii, "okay"] > 0 | var.check < 
                    m + 3 * vv) {
                    str.f <- c(str.f, stepII.breaks[ii, 1])
                    str <- rbind(str, stepII.breaks[ii, ])
                    cutstr2 <- rbind(cutstr2, stepII.breaks.x.el[ii, 
                      ])
                    ok.points <- ok.points + stepII.breaks.x.el[ii, 
                      ]
                  }
                }
            }
            if (sum(ok.points > 1) > 0) {
                print("ERROR")
            }
            if (sum(ok.points) == n) {
                flag <- TRUE
            }
            else {
                i <- (v - 1)
            }
        }
        for (i in 1:nrow(str)) {
            if (str$step[i] != 0) {
                xel <- mag.output$x.el[[str$step[i]]]
            }
            if (str$step[i] == 0) {
                xel <- mag.output$x.el.preprocess
            }
            for (j in 1:nrow(xel)) {
                m1 <- mean(mag.output$freq.s[xel[j, ] == 1])
                if (abs(m1 - str[i, 1]) < 0.001) {
                  cutstr <- rbind(cutstr, xel[j, ])
                }
            }
        }
        time1 <- as.numeric(Sys.time())
        # cat("took ", (time1 - time0), " seconds.\n")  # silenced
        colors <- colSums(cutstr * c(1:nrow(cutstr))) + 1
    }
    final.data <- cbind(mag.output$vaf.sorted, colors)
    t = final.data
    t.1 = c()
    for (i in 1:(ncol(t) - 2)) {
        temp = c()
        for (j in unique(t$colors)) {
            m.1 = mean(t[t$colors == j, i])
            v.1 = var(t[t$colors == j, i])
            a <- ((1 - m.1)/v.1 - 1/m.1) * m.1^2
            b <- a * (1/m.1 - 1)
            temp = rbind(temp, c(j, a, b))
        }
        colnames(temp) = c("colors", paste0("alpha.", i), paste0("beta.", 
            i))
        t.1 = cbind(t.1, temp[, -1, drop = FALSE])
    }
    t.1 = cbind(unique(t$colors), t.1)
    colnames(t.1)[1] = "colors"
    i = 1
    prob.all = c()
    for (i in 1:nrow(t)) {
        vf = t[i, 1:(ncol(t) - 2), drop = F]
        prob.s = c()
        for (j in 1:length(vf)) {
            probs = c()
            prob = c()
            for (ii in unique(t$colors)) {
                prob = c(prob, dbeta(vf[, j], shape1 = t.1[t.1[, 
                  1] == ii, 1 + (2 * j - 1)], shape2 = t.1[t.1[, 
                  1] == ii, 1 + 2 * j]))
            }
            prob.s = rbind(prob.s, prob)
            prob.s = round(apply(prob.s, 2, min), 2)
        }
        prob.all = rbind(prob.all, prob.s)
    }
    colnames(prob.all) = paste0("colors.", unique(t$colors))
    rownames(prob.all) = NULL
    final.probs = cbind(t, prob.all)
    results <- list(str = str, cutclust = cutstr, cutclust2 = cutstr2, 
        colors = colors, final.data = final.data, probs = final.probs)
    return(results)
}
fit.elements.reduce.v3 <-
function (elements, vaf.data) 
{
    result <- c()
    for (i in 1:(dim(vaf.data)[2] - 1)) {
        variants <- vaf.data[elements, i]
        dist <- max(variants) - min(variants)
        result <- rbind(result, dist)
    }
    return(max(result))
}
fit.elements.v3 <-
function (elements, vaf.data, depth.data) 
{
    result <- c()
    for (i in 1:(dim(vaf.data)[2] - 1)) {
        variants <- vaf.data[elements, i]
        depths <- depth.data[elements, i]
        dist <- mag.dist.v3(variants, depths)
        result <- rbind(result, dist)
    }
    return(result[which.max(result[, "dist"]), c("dist", "s1")])
}
get.var.v3.1 <-
function (x, x.nb, s) 
{
    x.nb <- x.nb[as.logical(x)]
    a <- mean(x.nb)
    b <- var(x.nb)
    c <- sum(as.logical(x))
    d <- s
    return(c(a, b, c, d))
}
mag.cn.call <-
function (cut, cn.input) 
{
    cn.output = c()
    temp = c()
    temp2 = data.frame(cut$final.data %>% group_by(colors) %>% 
        summarise_all(funs(mean)))
    for (m in 1:nrow(cn.input[[1]])) {
        res.m = c()
        for (c in unique(temp2$colors)) {
            res.c = 0
            for (s in 1:length(cn.input)) {
                res = c()
                cn.sample = cn.input[[s]]
                cn.tmp = cn.sample$cn[m]
                cn.vaf = cn.sample$vaf[m]
                for (k in 1:(10 * cn.tmp)) {
                  vaf.sc = temp2[temp2$colors == c, s + 1]
                  res = rbind(res, c(abs(cn.vaf - 2 * k * vaf.sc/cn.tmp), 
                    k))
                }
                res.c[1] = res.c[1] + res[which.min(res[, 1]), 
                  1]
                res.c = cbind(res.c, res[which.min(res[, 1]), 
                  2])
                colnames(res.c)[ncol(res.c)] = paste0("k.", s)
            }
            res.m = rbind(res.m, c(res.c, color = c))
            colnames(res.m) = c(colnames(res.c), "color")
        }
        cn.output = rbind(cn.output, res.m[which.min(res.m[, 
            1]), ])
        colnames(cn.output)[c(1, ncol(cn.output))] = c("error", 
            "cluster")
    }
    if (nrow(cut$final.data) < nrow(cn.input[[1]])) {
        print("Warning: number of CNV events is greater than SNVs. CNV assignment may not be reliable.")
    }
    return(data.frame(cn.output))
}
mag.dist.v3 <-
function (vafs, depths) 
{
    vaf <- mean(vafs[vafs > 0.001])
    if (sum(vafs > 0.001) == 0) {
        vaf <- mean(vafs)
    }
    depth <- mean(depths)
    s1 <- depth * vaf
    s2 <- depth - s1
    R <- c(dbeta(vafs[vafs != 0.001], shape1 = s1, shape2 = s2, 
        log = T), dbeta(vafs[vafs == 0.001], shape1 = s1, shape2 = s2, 
        log = T))
    v <- sum((vafs - vaf)^2)/length(vafs)
    v <- sqrt(v)
    range <- max(vafs) - min(vafs)
    v <- ifelse(v < 0.00050000000000000001, 0.00050000000000000001, 
        v)
    range <- ifelse(range < 0.0050000000000000001, 0.0050000000000000001, 
        range)
    if (mean(R) > 0) {
        R <- R/(range * v)
    }
    if (mean(R) < 0) {
        R <- -R * range * v
    }
    dist <- -mean(R)
    results <- c(dist = dist, s1 = s1, s2 = s2)
    return(results)
}
mag.exp.var.v2.1 <-
function (efrq, edep.vec, num, n = 1000) 
{
    expt.var <- c()
    depth.keep <- c()
    vafs <- c()
    if (length(edep.vec) > 50) {
        edep.vec <- edep.vec[edep.vec >= quantile(edep.vec, 0.14999999999999999) & 
            edep.vec <= quantile(edep.vec, 0.84999999999999998)]
    }
    else {
        edep.vec <- c(mean(edep.vec), mean(edep.vec))
    }
    for (i in 1:n) {
        edep <- sample(edep.vec, 1, replace = T)
        s1 <- efrq * edep
        s2 <- (1 - efrq) * edep
        x <- rbeta(num, shape1 = s1, shape2 = s2)
        vafs <- c(vafs, s1/(s1 + s2))
        expt.var <- c(expt.var, var(x))
        depth.keep <- c(depth.keep, edep)
        i <- i + 1
    }
    return(list(exp.vat = (expt.var), depths = depth.keep, vafs = vafs))
}
mag.exp.var.v3 <-
function (efrq.vec, edep.vec, num, n = 1000, type = "M") 
{
    expt.var <- c()
    depth.keep <- c()
    vafs <- c()
    a <- quantile(edep.vec, 0.050000000000000003)
    b <- quantile(edep.vec, 0.94999999999999996)
    lowerf.1 <- quantile(efrq.vec, 0.14999999999999999)
    higherf.1 <- quantile(efrq.vec, 0.84999999999999998)
    if (length(edep.vec) > 5) {
        edep.vec <- edep.vec[edep.vec >= a & edep.vec <= b]
    }
    else {
        edep.vec <- c(mean(edep.vec), mean(edep.vec))
    }
    if (length(efrq.vec) > 5) {
        efrq.vec <- efrq.vec[efrq.vec >= lowerf.1 & efrq.vec <= 
            higherf.1]
    }
    else {
        efrq.vec <- c(mean(efrq.vec), mean(efrq.vec))
    }
    for (i in 1:n) {
        edep <- sample(edep.vec, 1, replace = T)
        efrq <- mean(efrq.vec)
        s1 <- efrq * edep
        s2 <- (1 - efrq) * edep
        x <- rbeta(num, shape1 = s1, shape2 = s2)
        vafs <- c(vafs, s1/(s1 + s2))
        expt.var <- c(expt.var, var(x))
        depth.keep <- c(depth.keep, edep)
        i <- i + 1
    }
    return(list(exp.vat = (expt.var), depths = depth.keep, vafs = vafs))
}
mag.multiple <-
function (prep.data, x.el.cut = data.frame()) 
{
    if (1 > dim(prep.data$vafs)[2] - 1) {
        print("stop the algorithm and check the sampleNum; run mag.prepdata() first")
    }
    time0 <- as.numeric(Sys.time())
    params <- c()
    x.values <- c()
    x.steps <- c()
    row.loop <- c()
    x.step <- c()
    x.elements.list <- list()
    okay.list = list()
    if (ncol(prep.data$vafs) == 2) {
        print("executing single sample clustering....")
        flush.console()
    }
    order.var <- apply(prep.data$vafs[, -which(names(prep.data$vafs) == 
        "ID")], 1, prod)
    order.ind <- order(order.var)
    sorted.order.var <- order.var[order.ind]
    vaf.data <- prep.data$vafs
    depth.data <- prep.data$depths
    vaf.sorted <- vaf.data[order.ind, ]
    depth.sorted <- depth.data[order.ind, ]
    IDs <- vaf.sorted[, "ID"]
    n <- dim(vaf.data)[1]
    if (sum(dim(x.el.cut)) == 0) {
        print("Mag.reduce has not been run.")
        x.elements <- data.frame(diag(n))
        x.elements <- sapply(x.elements, as.logical)
        colnames(x.elements) <- IDs
    }
    else {
        x.elements <- x.el.cut
    }
    x.elements.list[[1]] = x.elements
    okay.id = c(1:nrow(x.elements))
    okay.list[[1]] = okay.id
    l <- nrow(x.elements)
    indx.cut <- c(1:nrow(x.elements))
    combn.cut <- t(combn(indx.cut, 2))
    cat("initial fitting ... ")
    flush.console()
    time1 <- as.numeric(Sys.time())
    initial.1 <- apply(combn.cut, 1, function(x) {
        fit.elements.v3(c(x.elements[x[1], ] | x.elements[x[2], 
            ]), vaf.sorted, depth.sorted)
    })
    time2 <- as.numeric(Sys.time())
    # cat("took ", time2 - time1, " seconds.\n")  # silenced
    flush.console()
    tempMatLoop <- Matrix(1000000, nrow = l, ncol = l, sparse = T)
    tempMatLoop.2 <- lower.tri(diag(1000000, nrow = l, ncol = l))
    tempMatLoop[tempMatLoop.2] <- initial.1[1, ]
    mat.loop <- t(tempMatLoop)
    tempPar1Loop <- lower.tri(diag(1000000, nrow = l, ncol = l))
    tempPar1Loop[tempPar1Loop] <- initial.1[2, ]
    par.1.loop <- t(tempPar1Loop)
    time2 <- as.numeric(Sys.time())
    s <- 2
    while (nrow(x.elements) > 1) {
        cat(s, "...", nrow(x.elements), "...\n")
        flush.console()
        x.values[s] <- min(mat.loop)
        el1 <- which(mat.loop == min(mat.loop), arr.ind = T)[1, 
            1]
        el2 <- which(mat.loop == min(mat.loop), arr.ind = T)[1, 
            2]
        el.rm <- c(el1, el2)
        params <- c(params, par.1.loop[el1, el2])
        x.step <- x.elements[el.rm[1], ] | x.elements[el.rm[2], 
            ]
        x.steps <- rbind(x.steps, x.step)
        x.elements <- rbind(x.elements[-el.rm, ], x.step)
        x.elements.list[[s]] <- x.elements
        okay.id = c(okay.id[-el.rm], 0)
        okay.list[[s]] = okay.id
        mat.loop <- mat.loop[-el.rm, -el.rm, drop = F]
        par.1.loop <- par.1.loop[-el.rm, -el.rm, drop = F]
        mat.column.update <- c()
        par1.column.update <- c()
        x1 <- x.elements[1:(nrow(x.elements) - 1), , drop = F]
        x2 <- x.elements[nrow(x.elements), , drop = F]
        fit.between <- apply(x1, 1, function(x) {
            fit.elements.v3(x | x2, vaf.sorted, depth.sorted)
        })
        mat.column.update <- fit.between[1, ]
        par1.column.update <- fit.between[2, ]
        if (nrow(mat.loop) != 0) {
            mat.loop <- cbind(mat.loop, mat.column.update)
            mat.loop <- rbind(mat.loop, rep(1000000, dim(mat.loop)[2]))
            par.1.loop <- cbind(par.1.loop, par1.column.update)
            par.1.loop <- rbind(par.1.loop, rep(1000000, dim(par.1.loop)[2]))
        }
        s <- s + 1
    }
    time3 <- as.numeric(Sys.time())
    cat("loop took ", time3 - time2, " seconds.\n")
    cat("total took ", time3 - time0, " seconds.\n")
    result <- list(x.el = x.elements.list, vaf.sorted = vaf.sorted, 
        depth.sorted = depth.sorted, prep.data = prep.data, okay.list = okay.list)
    return(result)
}
mag.multiple.run <-
function (input.data, reduce = T, threshold = 0.029999999999999999, 
    fold = F) 
{
    prep = mag.prepdata(input.data)
    purity = c()
    if (reduce) {
        red = mag.reduce.graph(prep, threshold = threshold)
        mag = mag.multiple(prep, x.el.cut = red$x.el.red)
        mv = mag.var(mag)
        cut = cut.off.multiple(mv)
    }
    if (!reduce) {
        mag = mag.multiple(prep)
        mv = mag.var(mag)
        cut = cut.off.multiple(mv)
    }
    temp = cut$final.data
    sum1 = temp %>% group_by(colors) %>% summarise_all(mean)
    purity = apply(sum1, 2, max) * 2
    purity = purity[c(-1, -length(purity))]
    fs = which(sum1[, -c(1, ncol(sum1))] > 0.5, arr.ind = T)
    print(purity)
    if (fold == F & any(purity > 1)) {
        purity = "There are clusters with frequency higher than 0.5- consider folding."
    }
    if (fold == T & nrow(fs) > 0) {
        fold.prep = prep
        for (i in 1:nrow(fs)) {
            bad_col = as.numeric(sum1[fs[i, 1], 1])
            bad_sam = fs[i, 2]
            bad_ID = temp$ID[temp$colors == bad_col]
            head(temp)
            fold.prep$vafs[fold.prep$vafs$ID %in% bad_ID, bad_sam] = 1 - 
                fold.prep$vafs[fold.prep$vafs$ID %in% bad_ID, 
                  bad_sam]
        }
        if (reduce) {
            red = mag.reduce.graph(fold.prep, threshold = threshold)
            mag = mag.multiple(fold.prep, x.el.cut = red$x.el.red)
            mv = mag.var(mag)
            cut = cut.off.multiple(mv)
        }
        if (!reduce) {
            mag = mag.multiple(fold.prep)
            mv = mag.var(mag)
            cut = cut.off.multiple(mv)
        }
    }
    temp = cut$final.data
    sum1 = temp %>% group_by(colors) %>% summarise_all(mean)
    purity = apply(sum1, 2, max) * 2
    purity = purity[c(-1, -length(purity))]
    res = list(results = temp, purity = purity, reduce = reduce, 
        mag = mag, cut = cut)
    return(res)
}
mag.prepdata <-
function (arg.data) 
{
    num <- dim(arg.data)[2]/2
    # cat("Number of samples: ", num)  # silenced
    arg.data <- data.frame(arg.data)
    vaf.data <- c()
    depth.data <- c()
    for (i in 1:num) {
        depth <- arg.data[, i * 2] + arg.data[, i * 2 - 1]
        vaf <- arg.data[, i * 2]/depth
        depth.data <- cbind(depth.data, depth)
        colnames(depth.data)[i] <- paste("depth.", i, sep = "")
        vaf.data <- cbind(vaf.data, vaf)
        colnames(vaf.data)[i] <- paste("vaf.", i, sep = "")
    }
    vaf.data <- as.matrix(vaf.data)
    vaf.data <- round(vaf.data, 3)
    vaf.data <- ifelse(vaf.data < 0.001, 0.001, vaf.data)
    vaf.data <- ifelse(vaf.data > 1 - (0.001), 1 - (0.001), vaf.data)
    counts.data <- data.frame(arg.data, row.names = NULL)
    vaf.data <- data.frame(vaf.data, row.names = NULL)
    depth.data <- data.frame(depth.data, row.names = NULL)
    ID <- c(1:dim(counts.data)[1])
    counts.data$ID <- ID
    vaf.data$ID <- ID
    depth.data$ID <- ID
    results <- list(counts = counts.data, vafs = vaf.data, depths = depth.data)
    return(results)
}
mag.reduce.graph <-
function (prep.data, x.el = matrix(), n = -100, threshold = 0.029999999999999999) 
{
    time0 <- as.numeric(Sys.time())
    order.var <- apply(prep.data$vafs[, -which(names(prep.data$vafs) == 
        "ID")], 1, prod)
    order.ind <- order(order.var)
    sorted.order.var <- order.var[order.ind]
    vaf.data <- prep.data$vafs
    depth.data <- prep.data$depths
    vaf.sorted <- vaf.data[order.ind, ]
    depth.sorted <- depth.data[order.ind, ]
    IDs <- vaf.sorted[, "ID"]
    if (dim(x.el)[1] == 1 & dim(x.el)[2] == 1) {
        l <- dim(vaf.data)[1]
        x.el <- data.frame(diag(l))
        x.el <- sapply(x.el, as.logical)
        colnames(x.el) <- IDs
    }
    mean.frqs <- t(apply(x.el, 1, function(x) {
        mean(sorted.order.var[x])
    }))
    x.el.sorted <- x.el[order(mean.frqs), ]
    x.el.final <- c()
    size.threshold <- 10
    l <- nrow(x.el)
    full <- floor(l/size.threshold)
    rem <- l%%size.threshold
    row.fold <- 0
    if (rem < 5) {
        full <- full - 1
        rem <- rem + size.threshold
    }
    optim.value.bound <- threshold
    ind.full <- c(1:size.threshold)
    combn.el <- t(combn(ind.full, 2))
    if (full > 0) {
        for (i in 1:full) {
            fold <- c((1 + (i - 1) * size.threshold):(i * size.threshold))
            x.el.fold <- x.el.sorted[fold, ]
            vafs.fold <- vaf.sorted[apply(x.el.fold, 2, function(x) sum(x) == 
                1), ]
            initial.full <- apply(combn.el, 1, function(x) {
                fit.elements.reduce.v3(x.el.fold[x[1], ] | x.el.fold[x[2], 
                  ], vaf.sorted)
            })
            tempMatLoop <- matrix(100, nrow = size.threshold, 
                ncol = size.threshold)
            tempMatLoop.2 <- lower.tri(diag(100, nrow = size.threshold, 
                ncol = size.threshold))
            tempMatLoop[tempMatLoop.2] <- initial.full
            mat.loop <- t(tempMatLoop)
            mat.loop.g <- mat.loop
            mat.loop.g[mat.loop <= optim.value.bound] <- 1
            mat.loop.g[mat.loop > optim.value.bound] <- 0
            sum(mat.loop.g == 1)
            g1 <- graph_from_adjacency_matrix(mat.loop.g, mode = "undirected")
            complete_points <- max_cliques(g1, min = 2)
            if (length(complete_points) > 0) {
                x.el.graph <- c()
                el.rm.graph <- c()
                for (i in 1:length(complete_points)) {
                  el.rm <- as.numeric(complete_points[[i]])
                  el.rm <- el.rm[!el.rm %in% el.rm.graph]
                  if (length(el.rm) > 0) {
                    x.el.graph <- rbind(x.el.graph, as.logical(colSums(x.el.fold[el.rm, 
                      , drop = F])))
                    el.rm.graph <- c(el.rm.graph, el.rm)
                  }
                }
                x.el.fold <- rbind(x.el.fold[-el.rm.graph, , 
                  drop = F], x.el.graph)
            }
            x.el.final <- rbind(x.el.final, x.el.fold)
            {
            }
        }
    }
    rems <- c((l - rem + 1):l)
    x.el.rem <- x.el.sorted[rems, ]
    ind.rem <- c(1:rem)
    combn.rem <- t(combn(ind.rem, 2))
    initial.rems <- apply(combn.rem, 1, function(x) {
        fit.elements.reduce.v3(x.el.rem[x[1], ] | x.el.rem[x[2], 
            ], vaf.sorted)
    })
    tempMatLoop <- matrix(100, nrow = rem, ncol = rem)
    tempMatLoop.2 <- lower.tri(diag(100, nrow = rem, ncol = rem))
    tempMatLoop[tempMatLoop.2] <- initial.rems
    mat.loop <- t(tempMatLoop)
    mat.loop.g <- mat.loop
    mat.loop.g[mat.loop <= optim.value.bound] <- 1
    mat.loop.g[mat.loop > optim.value.bound] <- 0
    sum(mat.loop.g == 1)
    g1 <- graph_from_adjacency_matrix(mat.loop.g, mode = "undirected")
    complete_points <- max_cliques(g1, min = 2)
    x.el.graph <- c()
    el.rm.graph <- c()
    if (length(complete_points) > 0) {
        for (i in 1:length(complete_points)) {
            el.rm <- as.numeric(complete_points[[i]])
            el.rm <- el.rm[!el.rm %in% el.rm.graph]
            if (length(el.rm) > 0) {
                x.el.graph <- rbind(x.el.graph, as.logical(colSums(x.el.rem[el.rm, 
                  , drop = F])))
                el.rm.graph <- c(el.rm.graph, el.rm)
            }
        }
        x.el.rem <- rbind(x.el.rem[-el.rm.graph, , drop = F], 
            x.el.graph)
        colSums(x.el.rem)
    }
    x.el.final <- rbind(x.el.final, x.el.rem)
    {
    }
    results <- list(x.el.red = x.el.final, vaf.data = vaf.data, 
        vaf.sorted = vaf.sorted, depth.data = depth.data, depth.sorted = depth.sorted)
    nn <- nrow(results$x.el.red)
    print(nn)
    time1 <- as.numeric(Sys.time())
    print("TIME:")
    time <- time1 - time0
    print(time)
    if (nn < size.threshold | nn == n) {
        return(results)
    }
    else {
        mag.reduce.graph(prep.data = prep.data, x.el = x.el.final, 
            n = nn, threshold = threshold)
    }
}
## --- mag.single (fast) ------------------------------------------------------
## Drop-in replacement for the SERVER's older mag.single (the one that gives
## 3 clusters on exampledata). All other functions in the package -- including
## mag.exp.var.v3, cut.off.single, fit.elements.v3, get.var.v3.1 -- are kept
## byte-for-byte from the server's source, so clustering decisions stay the
## same. Only the per-step stats-recompute is vectorised via membership +
## split(), and the N x N identity allocation is skipped.
mag.single <- function(prep.data) {
  if (is.null(prep.data$vafs$ID)) return(print('Run mag.prepdata on the data'))
  time0 <- as.numeric(Sys.time())

  vaf.data <- prep.data$vafs; vaf.data <- round(vaf.data, 3)
  vaf.data[, -dim(vaf.data)[2]] <- ifelse(vaf.data[, -dim(vaf.data)[2]] < 1e-3,   1e-3,   vaf.data[, -dim(vaf.data)[2]])
  vaf.data[, -dim(vaf.data)[2]] <- ifelse(vaf.data[, -dim(vaf.data)[2]] > 1-1e-3, 1-1e-3, vaf.data[, -dim(vaf.data)[2]])
  depth.data <- prep.data$depths

  ord <- order(vaf.data$vaf.1)
  sort.x          <- vaf.data$vaf.1[ord]
  sort.vaf.data   <- vaf.data[ord, ]
  sort.depth.data <- depth.data[ord, ]
  u.vafs <- unique(sort.x); N <- length(sort.x); l <- length(u.vafs)

  x.elements.preprc <- matrix(FALSE, nrow = l, ncol = N)
  membership <- integer(N)
  for (k in seq_along(u.vafs)) {
    mask <- (sort.x == u.vafs[k])
    x.elements.preprc[k, ] <- mask
    membership[mask] <- k
  }
  colnames(x.elements.preprc) <- sort.vaf.data$ID
  depth_col <- sort.depth.data$depth.1

  mat.loop   <- matrix(1000000, nrow = l, ncol = l)
  par.1.loop <- matrix( 100000, nrow = l, ncol = l)
  s.pairs <- cbind(1:(l-1), 2:l)
  s.likls <- t(apply(s.pairs, 1, function(x)
    fit.elements.v3(x.elements.preprc[x[1], ] | x.elements.preprc[x[2], ],
                    sort.vaf.data, sort.depth.data)))
  for (i in 2:l) {
    mat.loop[i-1, i]   <- s.likls[i-1, 1]
    par.1.loop[i-1, i] <- s.likls[i-1, 2]
  }
  time1 <- as.numeric(Sys.time())

  current_nrow <- l; okay <- seq_len(l)
  step_vm <- function(s) {
    f <- factor(membership, levels = seq_len(current_nrow))
    sx <- split(sort.x, f); sd <- split(depth_col, f)
    cbind(vapply(sx, mean, numeric(1)),
          vapply(sx, var,  numeric(1)),
          lengths(sx),
          rep(s, current_nrow),
          round(vapply(sd, mean, numeric(1))),
          okay)
  }

  vm.chunks <- vector("list", l); vm.chunks[[1]] <- step_vm(0)
  x.elements <- x.elements.preprc
  m.preprc   <- Matrix::Matrix(x.elements, sparse = TRUE)
  x.elements.list <- vector("list", l - 1)
  params <- c(); x.values <- numeric(l - 1); x.steps <- c()
  s <- 1; n <- N

  while (nrow(x.elements) >= 2) {
    mn  <- min(mat.loop)
    pos <- which(mat.loop == mn, arr.ind = TRUE)[1, ]
    el1 <- pos[1]; el2 <- pos[2]; el.rm <- c(el1, el2)
    x.values[s] <- mn
    params <- rbind(params, par.1.loop[el1, el2])

    x.step <- x.elements[el1, ] | x.elements[el2, ]
    x.steps <- rbind(x.steps, x.step)
    x.elements <- rbind(x.elements[-el.rm, ], x.step)
    x.elements.list[[s]] <- Matrix::Matrix(x.elements, sparse = TRUE)

    mn_el <- min(el1, el2); mx_el <- max(el1, el2)
    was_merged <- (membership == el1) | (membership == el2)
    membership <- membership - ifelse(membership > mx_el, 2L,
                              ifelse(membership > mn_el, 1L, 0L))
    new_idx <- current_nrow - 1L
    membership[was_merged] <- new_idx
    current_nrow <- new_idx
    okay <- c(okay[-el.rm], 0)
    vm.chunks[[s + 1]] <- step_vm(s)

    mat.loop   <- mat.loop  [-el.rm, -el.rm]
    par.1.loop <- par.1.loop[-el.rm, -el.rm]
    max_pos <- max(which(x.step)); min_pos <- min(which(x.step))
    neiIND <- integer(0)
    if      (min_pos == 1 & max_pos != n) neiIND <- which(x.elements[, max_pos + 1])
    else if (max_pos == n & min_pos != 1) neiIND <- which(x.elements[, min_pos - 1])
    else if (max_pos != n & min_pos != 1)
      neiIND <- c(which(x.elements[, min_pos - 1]), which(x.elements[, max_pos + 1]))

    if (length(neiIND) > 0) {
      nei <- x.elements[neiIND, , drop = FALSE]
      x.last <- x.elements[nrow(x.elements), , drop = FALSE]
      temp <- apply(nei, 1, function(x)
        fit.elements.v3(x | x.last, vaf.data = sort.vaf.data, depth.data = sort.depth.data))
      l2 <- sqrt(length(mat.loop))
      mat.column.update  <- rep(100000, l2); mat.column.update [neiIND] <- temp[1, ]
      par1.column.update <- rep(100000, l2); par1.column.update[neiIND] <- temp[2, ]
      mat.loop   <- cbind(mat.loop,   mat.column.update);   mat.loop   <- rbind(mat.loop,   rep(100000, ncol(mat.loop)))
      par.1.loop <- cbind(par.1.loop, par1.column.update);  par.1.loop <- rbind(par.1.loop, rep(100000, ncol(par.1.loop)))
    }
    s <- s + 1
  }
  var.mean <- do.call(rbind, vm.chunks[!sapply(vm.chunks, is.null)])
  colnames(var.mean) <- c("mean", "var", "NumberOfPoints", "step", "depth", "okay")
  x.elements.list <- x.elements.list[!sapply(x.elements.list, is.null)]
  time2 <- as.numeric(Sys.time())
  cat(sprintf("[MAGOS-fast]  N=%d  l=%d  init=%.1fs  merge=%.1fs  total=%.1fs\n",
              N, l, time1 - time0, time2 - time1, time2 - time0))
  flush.console()

  list(x.el = x.elements.list, x.el.preprocess = m.preprc,
       params = params, var.mean = var.mean, prep.data = prep.data,
       vaf.sorted = sort.vaf.data, depth.sorted = sort.depth.data,
       freq.s = sort.x, time = list(time0, time1, time2))
}
mag.single.run <-
function (input.data, fold = F) 
{
    data.prep <- mag.prepdata(input.data)
    purity <- 1
    mag <- mag.single(data.prep)
    cut <- cut.off.single(mag)
    temp <- merge(cut$final.data, data.prep$depths, by = "ID")
    sum1 <- temp %>% group_by(colors) %>% summarize(meanVAF = mean(vaf.1))
    purity <- 2 * max(sum1$meanVAF)
    if (purity > 1 & fold == F) {
        purity = "There are clusters with frequency higher than 0.5- consider folding."
    }
    if (fold == T & sum(sum1$meanVAF > 0.5) > 0) {
        fold.colors <- sum1$colors[sum1$meanVAF > 0.5]
        data.folded <- temp
        data.folded$vaf.1[data.folded$colors %in% fold.colors] <- 1 - 
            data.folded$vaf.1[data.folded$colors %in% fold.colors]
        data.prep.fold <- data.prep
        data.prep.fold$vafs$vaf.1 <- data.folded$vaf.1
        mag <- mag.single(data.prep.fold)
        cut <- cut.off.single(mag)
        temp <- merge(cut$final.data, data.prep$depths, by = "ID")
        sum2 <- temp %>% group_by(colors) %>% summarize(meanVAF = mean(vaf.1))
        purity <- 2 * max(sum2$meanVAF)
    }
    temp <- temp[, c(2, 4, 3, 1)]
    res <- list(mag = mag, cut = cut, results = temp, fold = fold, 
        purity = purity)
    return(res)
}
mag.var <-
function (mag.out) 
{
    x.el.list <- mag.out$x.el
    data <- data.frame(mag.out$vaf.sorted)
    num <- dim(data)[2] - 1
    var.mean.list <- list()
    okay.list = mag.out$okay.list
    for (i in 1:num) {
        var.mean <- c()
        x.data <- data[, i]
        names(data)
        depth.sorted.each <- mag.out$depth.sorted[, i]
        for (s in 2:length(x.el.list)) {
            x.el <- x.el.list[[s]]
            okay.2 = okay.list[[s]]
            temp <- t(apply(x.el, 1, get.var.v3.1, x.nb = x.data, 
                s = s))
            okay = okay.2
            depths <- apply(x.el, 1, function(x) round(mean(depth.sorted.each[as.logical(x)])))
            temp <- cbind(temp, depths, okay)
            var.mean <- rbind(var.mean, temp)
        }
        colnames(var.mean) <- c("mean", "var", "NumberOfPoints", 
            "step", "depth", "okay")
        var.mean.list[[i]] <- var.mean
        names(var.mean.list)[i] <- paste(names(data)[i], ".var.mean", 
            sep = "")
    }
    return(list(var.mean.all = var.mean.list, mag.out = mag.out))
}
