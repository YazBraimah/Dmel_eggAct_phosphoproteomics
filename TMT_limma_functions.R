read.peptides <- function(dat, cha){
    output <- NULL
    
    dat$Sequence <- as.character(dat$Sequence)
    dat$Protein.Group.Accessions <- as.character(dat$Protein.Group.Accessions)
    dat$Quan.Usage <- as.character(dat$Quan.Usage)
    dat$Quan.Info <- as.character(dat$Quan.Info)
    dat$Isolation.Interference <- as.numeric(as.character(dat$Isolation.Interference))
    
    dat <- subset(dat, Isolation.Interference>=30)  
    dat <- subset(dat, Quan.Usage=="Used")
    dat <- subset(dat, Protein.Group.Accessions!="")
    dat <- subset(dat, !apply(dat[cha], 1, f <- function(x) any(is.na(x))))  
}


quantify.proteins <- function(dat, cha){
    e.function <- function(x, seq) tapply(x, seq, median)
    output <- NULL
    
    dat$Sequence <- toupper(dat$Sequence) # Capital letters
    accessions <- as.character(unique(dat$Protein.Group.Accessions))
    n.proteins <- length(accessions)
    n.cha <- length(cha)
    
    for(k in 1:n.proteins){
        id <- accessions[k]
        sdat <- subset(dat, Protein.Group.Accessions==id)[c("Sequence", cha)]
        sdat[cha] <- log2(sdat[cha])
        sdat[cha] <- sdat[cha] - apply(sdat[cha], 1, median)
        pdat <- sdat[, -1]
        n.spectra <- ifelse(is.integer(dim(pdat)), nrow(pdat), 1)
        temp <- apply(sdat[,-1], 2, e.function,seq=sdat[, 1])          
        n.peptides <- ifelse(is.integer(dim(temp)), nrow(temp), 1)    
        if(is.integer(dim(pdat))) pdat <- apply(pdat, 2, median)
        pdat <- c(pdat, n.peptides=n.peptides, n.spectra=n.spectra)
        output <- rbind(output, pdat)
    }
    output[,1:n.cha] <- sweep(output[,1:n.cha],1,apply(output[,1:n.cha],1,median))
    output[,1:n.cha] <- sweep(output[,1:n.cha],2,apply(output[,1:n.cha],2,median))
    output[,1:n.cha] <- sweep(output[,1:n.cha],1,apply(output[,1:n.cha],1,median))
    output[,1:n.cha] <- round(output[,1:n.cha],3)
    row.names(output) <- accessions
    output <- as.data.frame(output)
    return(output)
}


eb.fit <- function(dat, design){
    n <- dim(dat)[1]
    fit <- lmFit(dat, design)
    fit.eb <- eBayes(fit)
    log2FC <- fit.eb$coefficients[, 2]
    df.r <- fit.eb$df.residual
    df.0 <- rep(fit.eb$df.prior, n)
    s2.0 <- rep(fit.eb$s2.prior, n)
    s2 <- (fit.eb$sigma)^2
    s2.post <- fit.eb$s2.post
    t.ord <- fit.eb$coefficients[, 2]/fit.eb$sigma/fit.eb$stdev.unscaled[, 2]
    t.mod <- fit.eb$t[, 2]
    p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
    p.mod <- fit.eb$p.value[, 2]
    q.ord <- qvalue(p.ord)$q
    q.mod <- qvalue(p.mod)$q
    results.eb <- data.frame(log2FC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
    results.eb <- results.eb[order(results.eb$p.mod), ]
    return(results.eb)
}

eb.fit.byCoef <- function(x, ebayes, thresh = 0){
    log2FC <- ebayes$coefficients[,x]
  t.ord <- ebayes$coefficients[,x]/ebayes$sigma/ebayes$stdev.unscaled[,x]
  t.mod <- ebayes$t[,x]
  p.ord <- 2*pt(-abs(t.ord), ebayes$df.residual)
  p.mod <- ebayes$p.value[,x]
  q.ord <- qvalue(p.ord)$q
  q.mod <- qvalue(p.mod)$q
  result <- data.frame(log2FC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod)
  result <- result[order(result$p.mod), ]
result$sig = ifelse(result$q.mod < 0.05 & (result$log2FC > thresh | result$log2FC < -thresh), "YES", "NO")
result$Accession = rownames(result)
result$comparison= colnames(ebayes)[x]
    return(result)
}

eb.fit.mult <- function(dat, design){
    n <- dim(dat)[1]
    fit <- lmFit(dat, design)
    fit.eb <- eBayes(fit)
    log2FC <- fit.eb$coef[, "tr2"]
    df.0 <- rep(fit.eb$df.prior, n)
    df.r <- fit.eb$df.residual
    s2.0 <- rep(fit.eb$s2.prior, n)
    s2 <- (fit.eb$sigma)^2
    s2.post <- fit.eb$s2.post
    t.ord <- fit.eb$coef[, "tr2"]/fit.eb$sigma/fit.eb$stdev.unscaled[, "tr2"]
    t.mod <- fit.eb$t[, "tr2"]
    p.ord <- 2*pt(-abs(t.ord), fit.eb$df.residual)
    p.mod <- fit.eb$p.value[, "tr2"]
    q.ord <- q.mod <- rep(NA,n)
    ids <- which(!is.na(p.ord))
    k <- 0
    q1 <- qvalue(p.ord[!is.na(p.ord)])$q
    q2 <- qvalue(p.mod[!is.na(p.mod)])$q
    for(i in ids){ 
        k <- k+1
        q.ord[i] <- q1[k] 
        q.mod[i] <- q2[k]
    }
    results.eb.mult <- data.frame(log2FC, t.ord, t.mod, p.ord, p.mod, q.ord, q.mod, df.r, df.0, s2.0, s2, s2.post)
    results.eb.mult <- results.eb.mult[order(results.eb.mult$p.mod), ]
    return(results.eb.mult)
}

## Modifications of functions to compare groups of lists 
## (from http://stackoverflow.com/questions/23559371/how-to-get-the-list-of-items-in-venn-diagram-in-r)
Intersect <- function (x) {  
    # Multiple set version of intersect
    # x is a list
    if (length(x) == 1) {
        unlist(x)
    } else if (length(x) == 2) {
        intersect(x[[1]], x[[2]])
    } else if (length(x) > 2){
        intersect(x[[1]], Intersect(x[-1]))
    }
}
#
Union <- function (x) {  
    # Multiple set version of union
    # x is a list
    if (length(x) == 1) {
        unlist(x)
    } else if (length(x) == 2) {
        union(x[[1]], x[[2]])
    } else if (length(x) > 2) {
        union(x[[1]], Union(x[-1]))
    }
}
#
Setdiff <- function (x, y) {
    # Remove the union of the y's from the common x's. 
    # x and y are lists of characters.
    xx <- Intersect(x)
    yy <- Union(y)
    setdiff(xx, yy)
}

#

calcIntErr <- function(data){
    tmpDF = select(data, contains("Abundance Ratio:"))
colnames(tmpDF) = gsub("Abundance Ratio: ", "", colnames(tmpDF))
colnames(tmpDF) = gsub("_", " ", colnames(tmpDF))
colnames(tmpDF) = gsub(", Egg", "", colnames(tmpDF))
colnames(tmpDF) = gsub(", Oocyte", "", colnames(tmpDF))
colnames(tmpDF) = gsub(" WT", "", colnames(tmpDF))
colnames(tmpDF) = gsub(" CanB2", "", colnames(tmpDF))
colnames(tmpDF) = gsub(" CnAact", "", colnames(tmpDF))
colnames(tmpDF) = gsub(" )", ")", colnames(tmpDF))
tmpDF = data.frame(first = abs((log2(tmpDF$`(129N) / (126)`))-(log2(tmpDF$`(129C) / (127C)`))),
          second = abs((log2(tmpDF$`(129C) / (126)`))-(log2(tmpDF$`(129N) / (127C)`))),
          third = abs((log2(tmpDF$`(130N) / (127N)`))-(log2(tmpDF$`(130C) / (128N)`))),
          fourth = abs((log2(tmpDF$`(130N) / (127N)`))-(log2(tmpDF$`(130C) / (128N)`))),
          fifth = abs((log2(tmpDF$`(130N) / (127N)`))-(log2(tmpDF$`(131) / (128C)`))),
          sixth = abs((log2(tmpDF$`(130N) / (127N)`))-(log2(tmpDF$`(131) / (128N)`))),
          seventh = abs((log2(tmpDF$`(130N) / (128N)`))-(log2(tmpDF$`(130C) / (127N)`))),
          eighth = abs((log2(tmpDF$`(130N) / (128N)`))-(log2(tmpDF$`(130C) / (127N)`))),
          ninth = abs((log2(tmpDF$`(130N) / (128N)`))-(log2(tmpDF$`(130C) / (128C)`))),
          tenth = abs((log2(tmpDF$`(130N) / (128N)`))-(log2(tmpDF$`(131) / (127N)`))),
          eleventh = abs((log2(tmpDF$`(130N) / (128N)`))-(log2(tmpDF$`(131) / (128C)`))),
          twelvth = abs((log2(tmpDF$`(130N) / (128C)`))-(log2(tmpDF$`(130C) / (127N)`))),
          thirteenth = abs((log2(tmpDF$`(130N) / (128C)`))-(log2(tmpDF$`(130C) / (128N)`))),
          fourteenth = abs((log2(tmpDF$`(130N) / (128C)`))-(log2(tmpDF$`(131) / (127N)`))),
          fifteenth = abs((log2(tmpDF$`(130N) / (128C)`))-(log2(tmpDF$`(131) / (128N)`))))
newDF = as.matrix(data.frame(quantile(tmpDF$first, .95, na.rm = TRUE), 
                  quantile(tmpDF$second, .95, na.rm = TRUE),
                  quantile(tmpDF$third, .95, na.rm = TRUE),
                  quantile(tmpDF$fourth, .95, na.rm = TRUE),
                  quantile(tmpDF$fifth, .95, na.rm = TRUE),
                  quantile(tmpDF$sixth, .95, na.rm = TRUE),
                  quantile(tmpDF$seventh, .95, na.rm = TRUE),
                  quantile(tmpDF$eighth, .95, na.rm = TRUE),
                  quantile(tmpDF$ninth, .95, na.rm = TRUE),
                  quantile(tmpDF$tenth, .95, na.rm = TRUE),
                  quantile(tmpDF$eleventh, .95, na.rm = TRUE),
                  quantile(tmpDF$twelvth, .95, na.rm = TRUE),
                  quantile(tmpDF$thirteenth, .95, na.rm = TRUE),
                  quantile(tmpDF$fourteenth, .95, na.rm = TRUE),
                  quantile(tmpDF$fifteenth, .95, na.rm = TRUE)))
myMedian = apply(newDF, 1, median)
return(myMedian)
    }

## Miscellaneous operators
'%!in%' <- function(x,y)!('%in%'(x,y))


PlotPhosphoPeps <- function(acc){
  text <- paste("  Gene name:", subset(canB_phosphoproteins_msa.annots, Accession == acc)$`Gene Symbol`,"\n",
    "FlyBase ID:", subset(canB_phosphoproteins_msa.annots, Accession == acc)$`Ensembl Gene ID`,"\n",
    "Description:", subset(canB_phosphoproteins_msa.annots, Accession == acc)$`Description`,"\n",
    "MF:", subset(canB_phosphoproteins_msa.annots, Accession == acc)$`Molecular Function`,"\n",
    "BP:", subset(canB_phosphoproteins_msa.annots, Accession == acc)$`Biological Process`,"\n",
    "CC:", subset(canB_phosphoproteins_msa.annots, Accession == acc)$`Cellular Component`,"\n",
     sep = "  ")
  text.p <- ggparagraph(text, face = "italic", size = 14)

  tmpPhosphoPepMSA = subset(canB_phosphopeptides.msa.se, Accession == acc)
  tmpPhosphoPepMSA.se = tmpPhosphoPepMSA %>% 
                separate(sample, c("tissue", "treatment"), "_")
  tmpPhosphoPepMSA.se$tissue = factor(tmpPhosphoPepMSA.se$tissue, levels = c("Oocyte", "Egg"))
  tmpPhosphoPepMSA.se$treatment = factor(tmpPhosphoPepMSA.se$treatment, levels = c("WT", "CanB2"))
  
  PhosphoPepMSA.gg = ggplot(tmpPhosphoPepMSA.se, aes(tissue, abundance, fill = tissue)) + 
              geom_bar(position=position_dodge(), stat="identity") + 
              geom_errorbar(aes(ymin=abundance-se, ymax=abundance+se), width=.2, position=position_dodge(.9)) + 
              facet_grid(treatment~Modification) + 
              theme_bw() + 
              theme(legend.position = "none") +
              scale_fill_manual(values = wes_palette("Royal1")) +
              labs(title = "Phosphopeptide (MSA)")
  
  tmpPhosphoPepNL = subset(canB_phosphopeptides.nl.se, Accession == acc)
  tmpPhosphoPepNL.se = tmpPhosphoPepNL %>% 
               separate(sample, c("tissue", "treatment"), "_")
  tmpPhosphoPepNL.se$tissue = factor(tmpPhosphoPepNL.se$tissue, levels = c("Oocyte", "Egg"))
  tmpPhosphoPepNL.se$treatment = factor(tmpPhosphoPepNL.se$treatment, levels = c("WT", "CanB2"))
  
  PhosphoPepNL.gg = ggplot(tmpPhosphoPepNL.se, aes(tissue, abundance, fill = tissue)) + 
              geom_bar(position=position_dodge(), stat="identity") + 
              geom_errorbar(aes(ymin=abundance-se, ymax=abundance+se), width=.2, position=position_dodge(.9)) + 
              facet_grid(treatment~Modification) + 
              theme_bw() + 
              theme(legend.position = "none") +
              scale_fill_manual(values = wes_palette("Royal1")) +
              labs(title = "Phosphopeptide (NL)")

    tmpProtein = subset(proteins.se, Accession == acc)
  tmpProtein.se = tmpProtein %>% 
            separate(sample, c("tissue", "treatment"), "_")
  tmpProtein.se$tissue = factor(tmpProtein.se$tissue, levels = c("Oocyte", "Egg"))
  tmpProtein.se$treatment = factor(tmpProtein.se$treatment, levels = c("WT", "CanB2"))

  wholeProt.gg = ggplot(tmpProtein.se, aes(tissue, abundance, fill = tissue)) + 
              geom_bar(position=position_dodge(), stat="identity") + 
              geom_errorbar(aes(ymin=abundance-se, ymax=abundance+se), width=.2, position=position_dodge(.9)) + 
              facet_grid(.~treatment) + 
              theme_bw() + 
              theme(legend.position = "none") +
              scale_fill_manual(values = wes_palette("Royal1")) +
              labs(title = "Whole protein")

    p <- ggarrange(text.p, wholeProt.gg, PhosphoPepMSA.gg, PhosphoPepNL.gg)
    return(p)
}

corr_eqn <- function(x,y, digits = 2) {
  corr_coef <- round(cor(x, y), digits = digits)
  paste("italic(r) == ", corr_coef)
}