imname=$1

header=$(head -n 1 ${imname}-sources.txt)
tail -n +2 -q ${imname}-sources.txt >> sources.tmp
tail -n +2 -q SRC*/*sources.txt >> sources.tmp
awk '{sub(/[^,]*/,"");sub(/,/,"")} 1' sources.tmp \
| awk '{printf("s0c%d,%s\n", NR-1, $1)}' \
| tail -n +1 > ${imname}-sources-all.txt
sed -i "1s/^/$header \n/" ${imname}-sources-all.txt
rm sources.tmp

tail -n +2 -q SRC*/*sources.txt >> outliers.tmp
awk '{sub(/[^,]*/,"");sub(/,/,"")} 1' outliers.tmp \
| awk '{printf("s0c%d,%s\n", NR-1, $1)}' \
| tail -n +1 > ${imname}-outliers.txt
sed -i "1s/^/$header \n/" ${imname}-outliers.txt
rm outliers.tmp