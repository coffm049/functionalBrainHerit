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
sed -e "s|\(\"pheno\": \".*/pconns/\).*\.parquet\"|\1part.${1}.parquet\"|" \
    -e "s|\(\"out\": \".*/FCs/\).*\"|\1pconns.${2}.RE.${1}.csv\"|" \
    -e "s|\(\"Method\": \).*|\1\"${2}\",|" \
    ${3}/reExample.json > "${3}/temp/${2}.RE.${1}.json"