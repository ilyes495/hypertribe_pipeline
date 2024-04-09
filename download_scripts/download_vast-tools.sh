# Download vast-tools repository
git clone https://github.com/vastgroup/vast-tools.git ../software/vast-tools
tmp=`pwd`
echo "export PATH=${pwd}/../software/vast-tools" >> ~/.bashrc

# Create data directory
mkdir ../genome_data/vast-tools

# Vastdb
wget http://vastdb.crg.eu/libs/vastdb.hs2.23.06.20.tar.gz -P ../genome_data/vast-tools/
cd ../genome_data/vast-tools/
tar xzvf vastdb.hs2.23.06.20.tar.gz
rm vastdb.hs2.23.06.20.tar.gz
