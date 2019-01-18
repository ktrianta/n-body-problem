

for i in {0..47}
do 
cp plot$i.txt plotTest$i.txt
echo "\"prog\",\"comp\",\"io\",\"tree\",\"comm\",\"gath\"" >> 48/plot$i.txt
sed -n '2, 22p' plot$i.txt >> plotTest$i.txt
sed -n '14, 22p' plot$i.txt >> plotTest$i.txt
mv plotTest$i.txt plot$i.txt
done

