for file in ../../data/snapmats/*.mtx
do ./fennel -v -e "$file" "MM" 8
done
