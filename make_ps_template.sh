echo '"\n \' > template.tps
while read line
do 
    echo "$line\\n \\" >> template.tps
done < template.eps 
echo '\n "' >> template.tps
