import os
import pandas as pd
import argparse
from jsonapi_client import Session as APISession
import requests
from pathlib import PurePath as pp
import sourmash
from Bio import SeqIO
import time

def parse_arguments():
    """
    This function is there to capture the arguments in the command line and to parse them to use them correctly.

    Please refer to the usage how to use this tool.

    """

    parser = argparse.ArgumentParser(
        prog="mgnify_search",
        description="queue MAGs against the MGnify db",
        usage="mgnify_search fasta_file_path output_path",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=True,
    )

    parser.add_argument(
        "fasta_file_path", 
        type=str, 
        help="enter the path where all the fasta files are stored to get queued",
        default="./"
    )

    parser.add_argument(
        'output_path',
        type=str,
        help='enter the output path where the result should be stored',
        default="./"
    )

    parser.add_argument("--version", action="version", version="1.0.0")

    parser.print_usage = parser.print_help

    args = parser.parse_args()

    return args

def check_path(output_path):
    if not os.path.exists(output_path):
        os.makedirs(output_path)
        print(f"[INFO] Created folder: {output_path}")
    else:
        print(f"[INFO] Folder already exists: {output_path}")

def load_fasta(fasta_file_path, output_path):
    check_path(output_path + '/sign')
    for mag in os.listdir(fasta_file_path):
        sketch = sourmash.MinHash(n=0, ksize=31, scaled=1000)

        for index, record in enumerate(SeqIO.parse(os.path.join(fasta_file_path, mag), 'fasta')):
            sketch.add_sequence(str(record.seq))

        signature = sourmash.SourmashSignature(sketch, name=record.name)
        with open(os.path.join(output_path + '/sign', pp(mag).stem + '.sig'), 'wt') as fp:
            sourmash.save_signatures([signature], fp)

def load_catalogues():
    catalogue_endpoint = "genome-catalogues"
    with APISession("https://www.ebi.ac.uk/metagenomics/api/v1") as mgnify:
        catalogues = map(lambda r: r.json, mgnify.iterate(catalogue_endpoint))
        catalogues = pd.json_normalize(catalogues)

    return catalogues

def load_to_mgnify(sig_path, catalogue_ids):
    endpoint = 'https://www.ebi.ac.uk/metagenomics/api/v1/genomes-search/gather'

    sign = [open(os.path.join(sig_path, sig), 'rb') for sig in os.listdir(sig_path)]
    sketch_uploads = [('file_uploaded', signature) for signature in sign]

    submitted_job = requests.post(endpoint, data={'mag_catalogues': catalogue_ids}, files=sketch_uploads).json()

    map(lambda fp: fp.close(), sign)

    job_done = False
    while not job_done:
        print('[INFO] Checking status...')

        query_result = None
        while not query_result:
            query_result = requests.get(submitted_job['data']['status_URL'])
            print('[INFO] Still waiting for jobs to complete. Current status of jobs')
            print('[INFO] Will check again in 15 seconds')
            time.sleep(15) 
        
        queries_status = {sig['job_id']: sig['status'] for sig in query_result.json()['data']['signatures']}
        job_done = all(map(lambda q: q == 'SUCCESS', queries_status.values()))
    
    print('Job done!')

    return pd.json_normalize(query_result.json()['data']['signatures'])

def get_taxonomy_of_mgnify_mag(match_row):
    mgyg_accession = match_row['result.match']
    with APISession("https://www.ebi.ac.uk/metagenomics/api/v1") as mgnify:
        genome_document = mgnify.get('genomes', mgyg_accession)
        return genome_document.resource.taxon_lineage

def save_taxa(matches, output_path):
    with open(os.path.join(output_path, 'taxa_matches.tsv'), 'w') as tsv:
        tsv.write(f'filename\tmatch id\tbest match lineage\n')
        for row, match in matches.iterrows():
            tsv.write(f"{match['filename']}\t{match['result.match']}\t{match['best_match_taxonomy']}\n")


def main():
    args = parse_arguments()

    check_path(args.output_path)

    print(f'[INFO] Load fasta to sourmash data')
    load_fasta(args.fasta_file_path, args.output_path)
    print(f'[INFO] Finish loading fasta to sourmash data')

    print(f'[INFO] Download catalogues')
    catalogues = load_catalogues()
    print(f'[INFO] Finish downloading')

    print(f'[INFO] Start API calls to MGnify to get loaded MAG data')
    result_df = load_to_mgnify(args.output_path + '/sign', list(catalogues['id']))
    print(f'[INFO] Finish all API calls')


    matches = result_df.dropna(subset=['result.match'])
    matches.to_csv(os.path.join(args.output_path, 'matches.tsv'), sep='\t', index=False)

    matches = matches.copy()
    matches['best_match_taxonomy'] = matches.apply(get_taxonomy_of_mgnify_mag, axis=1)

    save_taxa(matches, args.output_path)

    result = result_df.groupby('filename').apply(
    lambda query: query['result.matches'].sum() == 0)

    result.to_csv(os.path.join(args.output_path, "novelty.tsv"), header=True, sep="\t")

    result_df_copy = result_df[result_df['result.matches'] > 0].copy()

    result_df_copy.to_csv(os.path.join(args.output_path, 'novel_matches.tsv'), sep='\t', index=False)

if __name__ == "__main__":
    main()