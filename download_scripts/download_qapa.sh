# Download QAPQ repository
git clone https://github.com/morrislab/qapa.git ../software/qapa
cd ../software/qapa
conda env create -f environment.yml

# Create data directory
cd ../../
mkdir -p genome_data/qapa

# Gene Ensembl
mysql \
--user=anonymous --host=martdb.ensembl.org --port=5316 -A ensembl_mart_88 \
-e "select stable_id_1023 as 'Gene stable ID', stable_id_1066 as 'Transcript stable ID', \
biotype_1020 as 'Gene type', biotype_1064 as 'Transcript type', \
display_label_1074 as 'Gene name' from hsapiens_gene_ensembl__transcript__main" \
> genome_data/qapa/gene_ensembl.txt

# Gene Gencode
mysql \
--user=genome --host=genome-mysql.cse.ucsc.edu -A \
-e "select * from wgEncodeGencodeBasicV26" hg38 \
> genome_data/qapa/gene_gencode.txt

# PolyA Gencode
 mysql \
 --user=genome --host=genome-mysql.cse.ucsc.edu -A \
 -e "select chrom, txStart, txEnd, name2, score, strand \
 from wgEncodeGencodePolyaV26 where name2 = 'polyA_site'" -N hg38 \
 > genome_data/qapa/polya_gencode.bed
 
 # PolyA Polysite
wget --no-check-certificate https://polyasite.unibas.ch/download/atlas/2.0/GRCh38.96/atlas.clusters.2.0.GRCh38.96.bed.gz -P genome_data/qapa/
cd genome_data/qapa/
gzip -d atlas.clusters.2.0.GRCh38.96.bed.gz
mv atlas.clusters.2.0.GRCh38.96.bed polya_polyasite.bed
