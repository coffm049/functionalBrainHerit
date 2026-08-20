#!/bin/bash

source ~/software/MASH/.venv/bin/activate

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

# transpose here
python flipParquets.py --file /panfs/jay/groups/31/rando149/coffm049/ABCD/Workflow/02_Phenotypes/FCsTopo/data/probaConns/part.${1}.parquet

# Perform the replacements and save to a new file
sed -e "s|\(\"pheno\": \".*/probaConns/\).*\.parquet\"|\1part.${1}_trans.parquet\"|" \
    -e "s|\(\"out\": \".*/FCs/probaConns/\).*\"|\1pconns.${2}.FE.${1}.csv\"|" \
    -e "s|\(\"Method\": \).*|\1\"${2}\",|" \
    ${3}/feExample.json > "${3}/temp/${2}.FE.${1}.json"
