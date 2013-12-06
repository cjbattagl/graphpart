rm lambda.txt
for file in ./snapmats/*.mtx
do ./fennel -v -e "$file" "MM" $1
done
