for i in {0..47}
do
sed 's/ //g' plotData$i.txt > new$i.txt
tail -r new$i.txt > plotData$i.txt
rm new$i.txt
done
