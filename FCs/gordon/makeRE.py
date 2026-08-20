import json
import argparse


def update_mpheno(iteration, input_json_file = "reExample2.json", estMethod = "AdjHE"):
    
    new_mpheno = [f"o{i}" for i in range(iteration * 208, (iteration + 1) * 208 - 1)]

    with open(input_json_file, 'r') as f_in:
        data = json.load(f_in)
        data['mpheno'] = new_mpheno
        data["out"] = "/panfs/jay/groups/31/rando149/coffm049/ABCD/Results/03_heritability/Topography/FCs/pconns.%s.RE.%s.csv" % (estMethod, iteration)
        data["Method"] = estMethod
    with open(output_json_file, 'w') as f_out:
        json.dump(data, f_out)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Update mpheno field in JSON file')
    parser.add_argument('-i', '--iteration', type=int, required=True, help='Iteration number')
    parser.add_argument('-m', '--estMethod', required=True, help='Estimation method, one of AdjHE, GCTA, Twin')
    args = parser.parse_args()

    update_mpheno(args.iteration, input_json_file = "reExample2.json", estMethod = args.estMethod)
