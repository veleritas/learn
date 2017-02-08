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

def generate_parameters(max_elems=None, parts=None, metapaths=None):
    """Generate compound, disease, metapath combinations"""
    for n, metapath_dict in enumerate(metapaths):
        metapath = metapath_dict['abbreviation']
        query = metapath_dict['dwpc_query']
        for part_info in parts:
            if max_elems is not None and n == max_elems:
                break
            yield {
                'neo': part_info.neo,
                'hetnet': part_info.hetnet,
                'compound_id': part_info.compound_id,
                'disease_id': part_info.disease_id,
                'metapath': metapath,
                'query': query,
                'w': 0.4,
            }


def compute_dwpc(neo, hetnet, query, metapath, compound_id, disease_id, w):
    """Execute the neo4j query and write results to file"""
    start = time.time()
    results = neo.run(query, source=compound_id, target=disease_id, w=w)
    record = results.one
    seconds = '{0:.4g}'.format(time.time() - start)
    row = (
        hetnet, compound_id, disease_id, metapath,
        record['PC'], w, '{0:.6g}'.format(record['DWPC']), seconds
    )

    with writer_lock:
        writer.writerow(row)


def main():
    with open('servers.json') as read_file:
        instances = json.load(read_file)

    name_to_neo = dict()
    for instance in instances:
        neo = py2neo.Graph(
            host="localhost", http_port=instance['port']+1,
            bolt_port=instance['port'], bolt=True
        )

        name_to_neo[instance['name']] = neo

#-------------------------------------------------------------------------------

    with open('data/metapaths.json') as read_file:
        metapaths = json.load(read_file)

    metapaths.sort(key=lambda x: x['join_complexities'][0])
    print("Number of metapaths: {}".format(len(metapaths)))

#-------------------------------------------------------------------------------

    part_df = pandas.read_table('data/partitions.tsv')
    part_df['neo'] = part_df.hetnet.map(name_to_neo)
    parts = list(part_df.itertuples())

#-------------------------------------------------------------------------------

    # Total number of queries
    total_queries = len(metapaths) * len(part_df)

#-------------------------------------------------------------------------------

    # ## Execute queries

    # Parameters
    workers = 2
#    workers = mp.cpu_count()
    max_elems = None

    # Prepare writer
    path = 'data/dwpc.tsv.bz2'
    write_file = bz2.open(path, 'wt')
    writer = csv.writer(write_file, delimiter='\t')
    writer.writerow([
        'hetnet', 'compound_id', 'disease_id', 'metapath', 'PC',
        'w', 'DWPC', 'seconds'
    ])

    # Create ThreadPoolExecutor
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=workers)
    writer_lock = threading.Lock()

    # Submit jobs
    for params in tqdm(
        generate_parameters(max_elems, parts, metapaths),
        total = total_queries
    ):

        while executor._work_queue.qsize() > 10000:
            time.sleep(1)
        executor.submit(compute_dwpc, **params)

    # Shutdown and close
    executor.shutdown()
    write_file.close()

if __name__ == "__main__":
    main()
