fsuffix="/*.png"
fsuff2=".gif"
for f in exp_*/fittedimages; do
  if [ -d "$f" ]; then
     #echo "$f""$fsuffix"
     arr1=(${f//// })
     arr1len=${arr1}
     arrind=${arr1len-1}
     dataname="${arr1[$arrind]}"
     gifname="$f""/""$dataname""$fsuff2"
     echo "Making ""$gifname"
     convert "$f""$fsuffix" "$gifname"
  fi
done


for f in exp_*/Inputimages; do
  if [ -d "$f" ]; then
     #echo "$f""$fsuffix"
     arr1=(${f//// })
     arr1len=${arr1}
     arrind=${arr1len-1}
     dataname="${arr1[$arrind]}"
     gifname="$f""/""$dataname""$fsuff2"
     echo "Making ""$gifname"
     convert "$f""$fsuffix" "$gifname"
  fi
done

