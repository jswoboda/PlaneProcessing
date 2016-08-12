#!/bin/bash

ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done
echo "$ALLVAR"

if [ "$ALLVAR" = "" ]; then
    echo expression evaluated as true
fi
#grep $ALLVAR *.py
