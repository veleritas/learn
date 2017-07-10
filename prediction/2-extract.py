import argparse
import json
import threading
import concurrent.futures
import bz2
import csv
import time
import multiprocessing as mp

import pandas
import py2neo

from tqdm import tqdm


# Override the default py2neo timeout
py2neo.packages.httpstream.http.socket_timeout = 1e8

def main():

    def generate_parameters(max_elems=None):
        """Generate compound, disease, metapath combinations"""
        for n, metapath_dict in enumerate(metapaths):
            metapath = metapath_dict['abbreviation']
            query = metapath_dict['dwpc_query']

            for pair_info in pairs:
                if max_elems is not None and n == max_elems:
                    break
                yield {
                    'neo': neo,
                    'hetnet': 'hetio-ind',
                    'chemical_id': pair_info.chemical_id,
                    'disease_id': pair_info.disease_id,
                    'metapath': metapath,
                    'query': query,
                    'w': 0.4,
                }


    def compute_dwpc(neo, hetnet, query, metapath, chemical_id, disease_id, w):
        """Execute the neo4j query and write results to file"""

        # moved into main so that it can access the writer_lock variable

        start = time.time()
        results = neo.run(query, source=chemical_id, target=disease_id, w=w)

        # py2neo 3 uses cursors now, which consumes data differently
        all_records = results.data()
        assert len(all_records) == 1

        record = all_records[0]

        seconds = '{0:.4g}'.format(time.time() - start)
        row = (
            hetnet, chemical_id, disease_id, metapath,
            record['PC'], w, '{0:.6g}'.format(record['DWPC']), seconds
        )

        with writer_lock:
            writer.writerow(row)

#-------------------------------------------------------------------------------

    parser = argparse.ArgumentParser(description="Extract metapaths (prediction)")
    parser.add_argument(
        '--workers', type=int, default=mp.cpu_count(),
        choices=range(1, mp.cpu_count()+1)
    )
    args = parser.parse_args()


    with open('../all-features/servers.json') as read_file:
        instances = json.load(read_file)

    for instance in instances:
        if instance['name'] == 'rephetio-v2.0':
            neo = py2neo.Graph(
                host="localhost", http_port=instance['port']+1,
                bolt_port=instance['port'], bolt=True
            )

            break

#-------------------------------------------------------------------------------

    print("Number of workers: {}".format(args.workers))

    with open('features/metapaths.json') as read_file:
        metapaths = json.load(read_file)

    metapaths.sort(key=lambda x: x['join_complexities'][0])
    print("Number of metapaths: {}".format(len(metapaths)))

#-------------------------------------------------------------------------------

    pair_df = pandas.read_table('features/compound-disease-pairs.tsv.bz2')
    pairs = list(pair_df.itertuples())

#-------------------------------------------------------------------------------

    # Total number of queries
    total_queries = len(metapaths) * len(pair_df)

#-------------------------------------------------------------------------------

    # ## Execute queries

    # Parameters
    max_elems = None

    # Prepare writer
    path = 'features/dwpc.tsv.bz2'
    write_file = bz2.open(path, 'wt')
    writer = csv.writer(write_file, delimiter='\t')
    writer.writerow([
        'hetnet', 'chemical_id', 'disease_id', 'metapath', 'PC',
        'w', 'DWPC', 'seconds'
    ])

    # Create ThreadPoolExecutor
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=args.workers)
    writer_lock = threading.Lock()

    # Submit jobs
    for params in tqdm(generate_parameters(max_elems), total = total_queries):

        while executor._work_queue.qsize() > 10000:
            time.sleep(1)
        executor.submit(compute_dwpc, **params)

    # Shutdown and close
    executor.shutdown()
    write_file.close()

if __name__ == "__main__":
    main()
