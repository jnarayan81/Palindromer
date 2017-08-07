#!/usr/bin/env Rscript

# check for the argparse package installed
if (!('argparse' %in% installed.packages()[, 'Package'])) {
  install.packages('argparse', repos='http://cran.rstudio.com/')
}

suppressPackageStartupMessages(library('argparse'))

parse.args <- function() {
    parser <- ArgumentParser()
    parser$add_argument('-s', '--size', type='integer', default=640,
                        help='width and height of dotplot figure in px')
    parser$add_argument('-r', '--resolution', type='integer', default=72,
                        help='resoultion of dotplot figure in dpi')
    parser$add_argument('-q', '--query', default='query',
                        help='query label')
    parser$add_argument('-t', '--target', default='target',
                        help='target label')
    parser$add_argument('--title', default='', help='dotplot title')
    parser$add_argument('input', help='rdotplot file')
    parser$add_argument('-o', '--output', help='output figure file')
    return(parser$parse_args())
}

args <- parse.args()
if (is.null(args$output)) {
  args$output <- gsub('\\..*', '.png', args$input)
}
x <- read.table(args$input, as.is=TRUE)
png(args$output, width=args$size, height=args$size, res=args$resolution)
plot(x, type='l', xlab=args$target, ylab=args$query, main=args$title, col="red")
dev.off()
