#!/bin/bash

# Check if two arguments are provided
if [ $# -ne 3 ]
then
    echo "Invalid number of arguments supplied. Please provide three arguments."
    exit 1
fi

# Check if first argument is an integer
if ! [[ $1 =~ ^[0-9]+$ ]]
then
    echo "Invalid first argument. Please provide an integer."
    exit 1
fi

# Create temp directory if it doesn't exist
mkdir -p temp

# Perform the replacements and save to a new file
sed -e "s|\(\"pheno\": \".*/FCsTopo/\).*\.parquet\"|\1pconns.parquet\"|" \
    -e "s|\(\"out\": \".*/FCs/\).*\"|\1pconns.${2}.FE.${1}.csv\"|" \
    -e "s|\(\"Method\": \).*|\1\"${2}\",|" \
    ${3}/feExample.json > "${3}/temp/${2}.FE.${1}.json"
