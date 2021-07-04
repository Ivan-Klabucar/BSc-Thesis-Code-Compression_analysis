# Compression_analysis
Tool for generating csv files of phred quality score distributions of FASTQ reads,
and for calculating the mean error of phred quality compression used in the Biosoup module.

Example usage:
./compression_analyzer -t -f <path to csv file to fill with quality distribution> <path to FASTQ file>
The -t option signals the analyzer to calculate mean compression error on dataset and print it to standard output. 
