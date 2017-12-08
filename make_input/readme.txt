This jar file is used to create contact matrices in sparse matrix format. Contact matrices are normalized by ICE normalization.
This program is implemented for data from this GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM455133

Run the program as follow:
> java -jar MakeInput.jar chr_id resolution output_file input_files

chr_id: chromosome id (e.g 1,2,3 ...)
resolution: size of bins
output_file: an output file
input_files: a list of input files (separated by space)

Example:
> java -jar MakeInput.jar 1 1000000 chr1_1mb.txt input\GSE18199_RAW\GSM455137_30305AAXX.1.maq.hic.summary.binned.txt input\GSE18199_RAW\GSM455138_30305AAXX.2.maq.hic.summary.binned.txt

