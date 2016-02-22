rm lambda.txt
for file in ./snapmats/*.txt
do ./fennel -v -e "$file" "MM" $1 10000
done
