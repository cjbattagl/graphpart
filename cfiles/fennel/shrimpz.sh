for file in ../../data/snapmats/*.mtx
do ./fennel -v -e "$file" "MM"
done
