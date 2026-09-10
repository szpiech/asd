
# File: lodestarPlots.R
# Date: 10 August 2025
# Author: T. Quinn Smith
# Principal Investigator: Dr. Zachary A. Szpiech
# Purpose: Create plots from LODESTAR analysis.

# I am keeping this code redundant so it can be copied and modified.

# NOTE: I remove all dropped blocks and blocks that did not converge in cMDS.

library(ggplot2, quietly = TRUE, warn.conflicts = FALSE)
library(jsonlite)
library(dplyr, quietly = TRUE, warn.conflicts = FALSE)

exit <- function() { invokeRestart("abort"); }

# Print the usage statement for the script.
printUsage <- function() {
    cat("lodestarPlots.R v1.0\n");
    cat("----------------------\n\n");
    cat("Written by T. Quinn Smith\n");
    cat("Principal Investigator: Zachary A. Szpiech\n");
    cat("The Pennsylvania State University\n\n");
    cat("Usage: Rscript lodestarPlots.R [option] <blocks.json> (<pops.txt>)\n");
    cat("       <block.json> is the JSON file produced by LODESTAR.\n");
    cat("       <pops.txt> is optional. Labels sample by population.\n");
    cat("Option:\n");
    cat("-------\n");
    cat("mds b i j              Plot component j v. component i of the b'th block.\n");
    cat("target i j             Plot the centered and normalized components j v. i  of\n");
    cat("                           points Procrustes is performed against.\n");
    cat("global i j             Plot the centered and normalized components j v. i  of\n");
    cat("                           the global relatedness plot.\n");
    cat("axis i                 Plot the i'th component along the genome for each sample.\n");
    cat("tvals                  Plot the t-statistic along the genome. Ignores <pops.txt>.\n");
    cat("pvals                  Plot the p-values along the genome. Ignores <pops.txt>.\n");
    cat("var                    Plot captured variance along the genome. Ignores <pops.txt>.\n");
    cat("print b                Prints the coordinates of the b'th block. Ignores <pops.txt>.\n");
    cat("\n");
}

# Print the coordinates of the block.
printBlock <- function(JSON, w) {
    points <- as.data.frame(JSON$Blocks$X[JSON$Blocks["BlockNumber"] == w]);
    write.table(points, sep = "\t", row.names = FALSE, col.names = FALSE);
}

# Plot the global points.
global <- function(JSON, popsFile, i , j) {
    filename = paste("global_", i, "_", j, ".png", sep = "");
    if (popsFile != "") {
        points <- as.data.frame(JSON$GlobalX);
        labels <- read.delim(popsFile, header = FALSE)[,1];
        df <- data.frame(x = points[,as.integer(i)], y = points[,as.integer(j)], Populations = labels);
        p <- ggplot(df, aes(x = x, y = y, color = Populations)) +
            labs(Populations = "Populations") +
            geom_point() + 
            labs(x = paste("Axis ", i), y = paste("Axis ", j), title="Global") +
            theme(text = element_text(size = 12), plot.background = element_rect(fill = "white") );
        ggsave(filename, plot = p, width = 6, height = 4, dpi = 300);
    } else {
        points <- as.data.frame(JSON$GlobalX);
        df <- data.frame(x = points[,as.integer(i)], y = points[,as.integer(j)]);
        p <- ggplot(df, aes(x = x, y = y)) +
            geom_point() + 
            labs(x = paste("Axis ", i), y = paste("Axis ", j), title="Global") +
            theme(text = element_text(size = 12), plot.background = element_rect(fill = "white"));
        ggsave(filename, plot = p, width = 6, height = 4, dpi = 300);
    }
}

# Plot the target points.
target <- function(JSON, popsFile, i , j) {
    if (is.null(JSON$Y)) {
        print("No user defined coordinates were supplied. Exiting!\n");
        exit();
    }
    filename = paste("target_", i, "_", j, ".png", sep = "");
    if (popsFile != "") {
        points <- as.data.frame(JSON$Y);
        labels <- read.delim(popsFile, header = FALSE)[,1];
        df <- data.frame(x = points[,as.integer(i)], y = points[,as.integer(j)], Populations = labels);
        p <- ggplot(df, aes(x = x, y = y, color = Populations)) +
            labs(Populations = "Populations") +
            geom_point() + 
            labs(x = paste("Axis ", i), y = paste("Axis ", j), title="Y") +
            theme(text = element_text(size = 12), plot.background = element_rect(fill = "white") );
        ggsave(filename, plot = p, width = 6, height = 4, dpi = 300);
    } else {
        points <- as.data.frame(JSON$Y);
        df <- data.frame(x = points[,as.integer(i)], y = points[,as.integer(j)]);
        p <- ggplot(df, aes(x = x, y = y)) +
            geom_point() + 
            labs(x = paste("Axis ", i), y = paste("Axis ", j), title="Y") +
            theme(text = element_text(size = 12), plot.background = element_rect(fill = "white"));
        ggsave(filename, plot = p, width = 6, height = 4, dpi = 300);
    }
}

# Plot an axis's value at each window for all samples.
axis <- function(JSON, popsFile, i) {
    filename = paste("axis_", i, ".png", sep = "");

    numSamples = nrow(JSON$GlobalX);
    chroms = unlist(lapply(JSON$Blocks["Chromosome"], function (x) rep(x, numSamples)));
    bp = unlist(lapply(as.numeric(JSON$Blocks["StartCoordinate"][,1]), function (x) rep(x, numSamples)));
    components = unlist(lapply(JSON$Blocks$X, function(x) x[,as.integer(i)]));

    # If we have population membership, then label the points.
    if (popsFile != "") {
        labels <- read.delim(popsFile, header = FALSE)[,1];
        labels <- unlist(rep(labels, length(JSON$Blocks$X)));
        data <- data.frame(
            CHR = chroms,
            BP = bp,
            COMP = components,
            Populations = labels
        );
        colnames(data) <- c("CHR", "BP", "COMP", "Populations");
        g <- ggplot(data, aes(x = BP, y = COMP, color = as.factor(Populations))) +
            geom_point( aes(color=as.factor(Populations)), alpha=0.8, size=1.3) +
            labs(x = "Chromosome", y = paste("Axis ", i)) +
            theme_bw() +
            labs(color = "Populations") +
            theme(
                panel.border = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank()
            )
        ggsave(filename=filename, plot=g, width = 10, height = 5, dpi = 300);
    } else {
        data <- data.frame(
            CHR = chroms,
            BP = bp,
            COMP = components
        );
        colnames(data) <- c("CHR", "BP", "COMP");
        g <- ggplot(data, aes(x = BP, y = COMP)) +
            geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
            # Alternate two colours over however many chromosomes are present, rather than
        #  a hard-coded 44.
        scale_color_manual(values = rep(c("red", "blue"), length.out = length(unique(don$CHR)) )) +
            labs(x = "Chromosome", y = paste("Axis ", i)) +
            theme_bw() +
            theme( 
                text = element_text(size = 12),
                legend.position="none",
                panel.border = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank()
            )
        ggsave(filename=filename, plot=g, width = 10, height = 5, dpi = 300);
    }
}

# Plot the captured variance along the genome.
var <- function(JSON) {
    # Drop invalid blocks.
    filename = "variance_captured.png";
    data <- data.frame(
        CHR = JSON$Blocks["Chromosome"],  
        BP = as.numeric(JSON$Blocks["StartCoordinate"][,1]),  
        V = as.numeric(JSON$Blocks["VarianceCaptured"][,1]),
        SNP = "."
    );
    colnames(data) <- c("CHR", "BP", "V", "SNP");
    # NOTE: this was as.numeric(gsub("chr", "", CHR)), which turns any non-numeric
    #  chromosome name -- scaffolds, X, Y, and Drosophila-style 2L/3R -- into NA, and
    #  those blocks then vanish from the plot without a warning. Map names to positions
    #  in order of first appearance instead, and keep the name for the axis labels.
    data <- data %>%
        mutate(CHRNAME = as.character(CHR),
               CHR = match(CHRNAME, unique(CHRNAME)));
    # https://r-graph-gallery.com/101_Manhattan_plot.html
    don <- data %>% 
        group_by(CHR) %>% 
        summarise(chr_len=max(BP)) %>% 
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>%
        left_join(data, ., by=c("CHR"="CHR")) %>%
        arrange(CHR, BP) %>%
        mutate( BPcum=BP+tot);
    axisdf = don %>%
        group_by(CHR) %>%
        summarize(center=( max(BPcum) + min(BPcum) ) / 2, CHRNAME=first(CHRNAME) );
    if (length(unique(data$CHR)) > 1) {
        geopoint = geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3);
        xlab = scale_x_continuous("Chromosome", label = axisdf$CHRNAME, breaks= axisdf$center);
    } else {
        geopoint = geom_point();
        xlab = scale_x_continuous("Chromosome Position");
    }
    plot <- ggplot(don, aes(x=BPcum, y=V)) +
        geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
        # Alternate two colours over however many chromosomes are present, rather than
        #  a hard-coded 44.
        scale_color_manual(values = rep(c("red", "blue"), length.out = length(unique(don$CHR)) )) +
        scale_y_continuous("Proportion of Variance", expand = c(0, 0), limits = c(0, 1) ) + 
        xlab + geopoint +
        theme_bw() +
        theme( 
            text = element_text(size = 12),
            legend.position="none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
        ) + ggtitle("Variance Captured");
    ggsave(filename, plot, width = 10, height = 6, dpi = 300);
}

# Plot the pvals along the genome.
pvals <- function(JSON) {
    # Drop invalid blocks.
    filename = "pvalues.png";
    data <- data.frame(
        CHR = JSON$Blocks["Chromosome"],  
        BP = as.numeric(JSON$Blocks["StartCoordinate"][,1]),  
        P = -log10(as.numeric(JSON$Blocks["ProcrustesPValue"][,1])),
        SNP = "."
    );
    colnames(data) <- c("CHR", "BP", "P", "SNP");
    # NOTE: this was as.numeric(gsub("chr", "", CHR)), which turns any non-numeric
    #  chromosome name -- scaffolds, X, Y, and Drosophila-style 2L/3R -- into NA, and
    #  those blocks then vanish from the plot without a warning. Map names to positions
    #  in order of first appearance instead, and keep the name for the axis labels.
    data <- data %>%
        mutate(CHRNAME = as.character(CHR),
               CHR = match(CHRNAME, unique(CHRNAME)));
    # https://r-graph-gallery.com/101_Manhattan_plot.html
    don <- data %>% 
        group_by(CHR) %>% 
        summarise(chr_len=max(BP)) %>% 
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>%
        left_join(data, ., by=c("CHR"="CHR")) %>%
        arrange(CHR, BP) %>%
        mutate( BPcum=BP+tot);
    axisdf = don %>%
        group_by(CHR) %>%
        summarize(center=( max(BPcum) + min(BPcum) ) / 2, CHRNAME=first(CHRNAME) );
    if (length(unique(data$CHR)) > 1) {
        geopoint = geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3);
        xlab = scale_x_continuous("Chromosome", label = axisdf$CHRNAME, breaks= axisdf$center);
    } else {
        geopoint = geom_point();
        xlab = scale_x_continuous("Chromosome Position");
    }
    plot <- ggplot(don, aes(x=BPcum, y=P)) +
        geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
        # Alternate two colours over however many chromosomes are present, rather than
        #  a hard-coded 44.
        scale_color_manual(values = rep(c("red", "blue"), length.out = length(unique(don$CHR)) )) +
        scale_y_continuous("-log(p)", expand = c(0, 0)) + 
        xlab + geopoint +
        theme_bw() +
        theme( 
            text = element_text(size = 12),
            legend.position="none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
        ) + ggtitle("P Values");
    ggsave(filename, plot, width = 10, height = 6, dpi = 300);
}

# Plot the t-statistic along the genome.
tvals <- function(JSON) {
    # Drop invalid blocks.
    filename = "procrustes_statistic.png";
    data <- data.frame(
        CHR = JSON$Blocks["Chromosome"],  
        BP = as.numeric(JSON$Blocks["StartCoordinate"][,1]),  
        T = as.numeric(JSON$Blocks["ProcrustesStatistic"][,1]),
        SNP = "."
    );
    colnames(data) <- c("CHR", "BP", "t", "SNP");
    # NOTE: this was as.numeric(gsub("chr", "", CHR)), which turns any non-numeric
    #  chromosome name -- scaffolds, X, Y, and Drosophila-style 2L/3R -- into NA, and
    #  those blocks then vanish from the plot without a warning. Map names to positions
    #  in order of first appearance instead, and keep the name for the axis labels.
    data <- data %>%
        mutate(CHRNAME = as.character(CHR),
               CHR = match(CHRNAME, unique(CHRNAME)));
    # https://r-graph-gallery.com/101_Manhattan_plot.html
    don <- data %>% 
        group_by(CHR) %>% 
        summarise(chr_len=max(BP)) %>% 
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        select(-chr_len) %>%
        left_join(data, ., by=c("CHR"="CHR")) %>%
        arrange(CHR, BP) %>%
        mutate( BPcum=BP+tot);
    axisdf = don %>%
        group_by(CHR) %>%
        summarize(center=( max(BPcum) + min(BPcum) ) / 2, CHRNAME=first(CHRNAME) );
    if (length(unique(data$CHR)) > 1) {
        geopoint = geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3);
        xlab = scale_x_continuous("Chromosome", label = axisdf$CHRNAME, breaks= axisdf$center);
    } else {
        geopoint = geom_point();
        xlab = scale_x_continuous("Chromosome Position");
    }
    plot <- ggplot(don, aes(x=BPcum, y=t)) +
        geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
        # Alternate two colours over however many chromosomes are present, rather than
        #  a hard-coded 44.
        scale_color_manual(values = rep(c("red", "blue"), length.out = length(unique(don$CHR)) )) +
        scale_y_continuous("Procrustes t", expand = c(0, 0), limits = c(0, 1) ) + 
        xlab + geopoint +
        theme_bw() +
        theme( 
            text = element_text(size = 12),
            legend.position="none",
            panel.border = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()
        ) + ggtitle("Procrustes t");
    ggsave(filename, plot, width = 10, height = 6, dpi = 300);
}

# Plot MDS for the i and j axes for the given window.
mds <- function(JSON, popsFile, w, i, j) {
    filename = paste("cmds_", w, "_", i, "_", j, ".png", sep = "");
    if (popsFile != "") {
        points <- as.data.frame(JSON$Blocks$X[JSON$Blocks["BlockNumber"] == w])
        labels <- read.delim(popsFile, header = FALSE)[,1];
        df <- data.frame(x = points[,as.integer(i)], y = points[,as.integer(j)], Populations = labels);
        p <- ggplot(df, aes(x = x, y = y, color = Populations)) +
            labs(Populations = "Populations") +
            geom_point() + 
            labs(x = paste("Axis ", i), y = paste("Axis ", j), title="cMDS") +
            theme(text = element_text(size = 12), plot.background = element_rect(fill = "white") );
        ggsave(filename, plot = p, width = 6, height = 4, dpi = 300);
    } else {
        points <- as.data.frame(JSON$Blocks$X[JSON$Blocks["BlockNumber"] == w])
        df <- data.frame(x = points[,as.integer(i)], y = points[,as.integer(j)]);
        p <- ggplot(df, aes(x = x, y = y)) +
            geom_point() + 
            labs(x = paste("Axis ", i), y = paste("Axis ", j), title="cMDS") +
            theme(text = element_text(size = 12), plot.background = element_rect(fill = "white"));
        ggsave(filename, plot = p, width = 6, height = 4, dpi = 300);
    }
}

# Warning messages and run a give command.
#   No error checking.
cmd <- function(cmd, blocksFile, popsFile, args) {
    JSON = fromJSON(blocksFile);
    # Blocks whose cMDS did not converge now carry null (read as NA) rather than the -1
    #  sentinel, and dropped blocks are no longer emitted at all. Filter on NA, and keep
    #  the -1 test so a JSON file written by v1.0 still loads.
    JSON$Blocks = JSON$Blocks[!is.na(JSON$Blocks["ProcrustesStatistic"]) &
                              JSON$Blocks["ProcrustesStatistic"] != -1 &
                              JSON$Blocks["Chromosome"] != "GLOBAL",]
    switch(cmd,
        mds={
            mds(JSON, popsFile, args[1], args[2], args[3]);
        },
        tvals={
            tvals(JSON);
        },
        pvals={
            pvals(JSON);
        },
        axis={
            axis(JSON, popsFile, args[1]);
        },
        print={
            printBlock(JSON, args[1]);
        },
        var={
            var(JSON);
        },
        target={
            target(JSON, popsFile, args[1], args[2]);
        },
        global={
            global(JSON, popsFile, args[1], args[2]);
        },
        {
            return;
        }
    )
}

# Parse input arguments.

args <- commandArgs(trailingOnly = TRUE);

if (length(args) < 2) {
    printUsage();
    exit();
}

if (args[1] != "bh" && args[1] != "qvals" && args[1] != "target" && args[1] != "var" && args[1] != "print" && args[1] != "mds" && args[1] != "tvals" && args[1] != "axis" && args[1] != "pvals" && args[1] != "global") {
    cat("Unrecognized command! Exiting!\n");
    exit();
}

jsonFile = "";
popsFile = "";

if (substring(args[length(args)], nchar(args[length(args)]) - 3, nchar(args[length(args)])) == ".txt") {
    popsFile = args[length(args)];
    if (!file.exists(popsFile)) {
        cat("<pops.txt> does not exist! Exiting!\n");
        exit();
    }
    jsonFile = args[length(args) - 1];
    if (!file.exists(jsonFile)) {
        cat("<blocks.json> does not exist! Exiting!\n");
        exit();
    }
    cmd(args[1], jsonFile, popsFile, args[seq(2, length(args) - 2)]);
} else if (substr(args[length(args)], nchar(args[length(args)]) - 4, nchar(args[length(args)])) == ".json") {
    jsonFile = args[length(args)];
    if (!file.exists(jsonFile)) {
        cat("<blocks.json> does not exist! Exiting!\n");
        exit();
    }
    cmd(args[1], jsonFile, popsFile, args[seq(2, length(args) - 1)]);
} else {
    cat("<blocks.json> not provided! Exiting!\n");
    exit();
}