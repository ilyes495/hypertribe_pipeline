rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz ../genome_data/
cd ../genome_data/
gzip -d hg38.fa.gz
