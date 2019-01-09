mkdir 1
mkdir 2
mkdir 4
mkdir 8
mkdir 16
mkdir 32
mkdir 48



touch 1/plot0.txt
echo "\"prog\",\"comp\",\"io\",\"tree\",\"comm\",\"gath\"" >> 1/plot0.txt
sed -n '301, 350p' plotData0.txt >> 1/plot0.txt
for i in {0..1}
do 
touch 2/plot$i.txt
echo "\"prog\",\"comp\",\"io\",\"tree\",\"comm\",\"gath\"" >> 2/plot$i.txt
sed -n '251, 300p' plotData$i.txt >> 2/plot$i.txt
done
for i in {0..3}
do 
touch 4/plot$i.txt
echo "\"prog\",\"comp\",\"io\",\"tree\",\"comm\",\"gath\"" >> 4/plot$i.txt
sed -n '201, 250p' plotData$i.txt >> 4/plot$i.txt
done
for i in {0..7}
do 
touch 8/plot$i.txt
echo "\"prog\",\"comp\",\"io\",\"tree\",\"comm\",\"gath\"" >> 8/plot$i.txt
sed -n '151, 200p' plotData$i.txt >> 8/plot$i.txt
done
for i in {0..15}
do 
touch 16/plot$i.txt
echo "\"prog\",\"comp\",\"io\",\"tree\",\"comm\",\"gath\"" >> 16/plot$i.txt
sed -n '101, 150p' plotData$i.txt >> 16/plot$i.txt
done
for i in {0..31}
do 
touch 32/plot$i.txt
echo "\"prog\",\"comp\",\"io\",\"tree\",\"comm\",\"gath\"" >> 32/plot$i.txt
sed -n '51, 100p' plotData$i.txt >> 32/plot$i.txt
done
for i in {0..47}
do 
touch 48/plot$i.txt
echo "\"prog\",\"comp\",\"io\",\"tree\",\"comm\",\"gath\"" >> 48/plot$i.txt
sed -n '1, 50p' plotData$i.txt >> 48/plot$i.txt
done

