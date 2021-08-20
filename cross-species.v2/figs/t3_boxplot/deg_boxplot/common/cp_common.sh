cat common_up_deg_t3.tsv common_dn_deg_t3.tsv > tmp

dos2unix tmp

cat tmp | while read id
do
	cp ../${id}_*.png .
done

wait
rm tmp
