fsuffix="/*.png"
fsuff2=".gif"

BASEPATH=${PWD##*/}

while getopts "h?vf:" opt; do
    case "$opt" in
    h|\?)
        show_help
        exit 0
        ;;
    f)  BASEPATH=$OPTARG
        ;;
    esac
done
echo $BASEPATH/exp_*/fittedimages
for f in $BASEPATH/exp_*/fittedimages; do
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


for f in $BASEPATH/exp_*/Inputimages; do
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

