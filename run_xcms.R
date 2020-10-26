library(xcms)
library(magrittr)
library(optparse)

parser <- OptionParser()
parser <- add_option(parser, c("--ppm"), type="integer", default=15, help="Centwave Parts-Per-Million", metavar="number")
parser <- add_option(parser, c("--pwlower"), type="integer", default=15, help="Centwave lower bound for peakwidth", metavar="number")
parser <- add_option(parser, c("--pwupper"), type="integer", default=80, help="Centwave upper bound for peakwidth", metavar="number")
parser <- add_option(parser, c("--snthresh"), type="integer", default=5, help="Centwave snthresh", metavar="number")
parser <- add_option(parser, c("--noise"), type="integer", default=1000, help="Centwave noise", metavar="number")
parser <- add_option(parser, c("--prefilterlower"), type="integer", default=3, help="Centwave lower bound for prefilter", metavar="number")
parser <- add_option(parser, c("--prefilterupper"), type="integer", default=500, help="Centwave upper bound for prefilter", metavar="number")

args <- parse_args(parser, positional_arguments=TRUE)
output = args$args[1]
files = args$args[-1]

cwp <- CentWaveParam(
                        ppm = args$options$ppm, 
                        peakwidth = c(args$options$pwlower, args$options$pwupper), 
                        snthresh = args$options$snthresh, 
                        noise = args$options$noise, 
                        prefilter = c(args$options$prefilterlower, args$options$prefilterupper)
                    )

for (filename in files) {
  raw_data <- readMSData(files = filename,
                         pdata = new("NAnnotatedDataFrame",
                                     data.frame(sample_name=tools::file_path_sans_ext(filename))), msLevel.=1,
                         mode = "onDisk")
  xdata <- findChromPeaks(raw_data, param = cwp)
  cp = chromPeaks(xdata)
  if(length(cp) >= 1){
    data_out = data.frame(1:nrow(cp), cp[,'mz'], cp[,'rt'],
                            cp[,'rtmin'], cp[,'rtmax'],
                            cp[,'maxo'], cp[,'into'],
                            cp[,'mzmin'], cp[,'mzmax'])
    data_out = data_out[order(data_out[,2]),]
    data_out[,1] = 1:nrow(cp)
  }else{
    data_out = data.frame(vector(), double(), double(),
                            double(), double(),
                            double(), double(),
                            double(), double())
  }
  bname = basename(filename)
  colnames(data_out) = c('row ID', 'row m/z', 'row retention time', paste(bname, 'Peak RT start'),
                         paste(bname, 'Peak RT end'), paste(bname, 'Peak height'),
                         paste(bname, 'Peak area'), paste(bname, 'Peak m/z min'),
                         paste(bname, 'Peak m/z max'))
  write.csv(data_out, file=paste0(file.path(output, basename(tools::file_path_sans_ext(filename))), '_xcms_box.csv'), row.names=FALSE)
}
