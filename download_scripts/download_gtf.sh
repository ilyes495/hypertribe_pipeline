rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ensGene.gtf.gz ../genome_data/
cd ../genome_data/
gzip -d hg38.ensGene.gtf.gz
mv hg38.ensGene.gtf hg38.gtf
