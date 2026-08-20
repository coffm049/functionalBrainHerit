import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser(description="Read a file path using --file")
    parser.add_argument('--file', type=str, required=True, help="Path to the input file")
    args = parser.parse_args()
    # alter filepath in some way for output file
    # or maybe mv original file into .orig space then act on it 
    print("Reading: ")
    print(args.file)
    tempfile = args.file.split(".")
    pre = tempfile.pop(-1)
    newfile = ".".join(tempfile) + "_trans" + "." + pre
    
    temp = pd.read_parquet(args.file, engine = "pyarrow")
    print("Writing: ")
    print(newfile)
    temp.T.to_parquet(newfile, engine = "pyarrow")


if __name__ == "__main__" :
    main()
