""" This is a modified version of simple.py script bundled with Read Until API"""
import argparse
import logging
import traceback
import time
import numpy
import os
import h5py
import glob
import dyss

def _get_parser():
    parser = argparse.ArgumentParser('Dyss -- a tiny program for selective sequencing on MinION')
    parser.add_argument('--min_chunk_size', type=int, default=5000,
                        help='Minimum read chunk size to receive.')
    parser.add_argument('--num_scouts', default=14, type=int,
                        help='number of scouts. Default is 14')
    parser.add_argument('--num_packs', default=3, type=int,
                        help='number of packs. Default is 3')
    parser.add_argument('--reference', required=True,
                        help='reference seqence to be amplified. Currently reference size is bounded by 100Kbp.')
    parser.add_argument('--model', required=True,
                        help='model file.')
    parser.add_argument('--param', required=True,
                        help='training data.')
    parser.add_argument('--power', default=9, type=int,
                        help='chunking power. Integer type. Default is 9.')
    parser.add_argument('--referencesize', default=200000, type=int,
                        help='Reference size(Event num = 2*bp). Default is 200000')
    parser.add_argument('--test', required=True,
                        help='directly to generate test data. When test flag is on, this argument is required.')
    return parser


class MockRead(object):
    """
    Mock read for debug
    """
    def __init__(self,
                 read_number,
                 raw_data):
        self.number = read_number
        self.raw_data = raw_data.tostring()

def read_folder(path, batch_size=10, chunk_size=3000):
    """
    debug function
    """
    search_query = path + "*.fast5"
    print(search_query)
    queries = glob.glob(search_query)
    channel = 100
    read = 0
    result = []
    print("{}".format(len(queries)))
    for file in queries:
        if os.path.isfile(file):
            with h5py.File(file, 'r') as f:
                try:
                    read_num = [x for x in f['/Raw/Reads'].keys()][0]
                    id = f['/Raw/Reads/' + read_num].attrs['read_id'].decode()
                    signals = f['/Raw/Reads/' + read_num + '/Signal'].value[0:chunk_size]
                except:
                    continue
                else:
                    result.append((channel, MockRead(read, signals)))
                    channel += 1
                    read += 1
                    if len(result) >= batch_size:
                        break
    return result


def mock_running(classifer, path, batch_size=10, chunk_size=3000):
    """
    debug function
    :param dyss: an instance of Dyss object.
    :param path: a path to test dataset direclty. That directly should contain .fast5 file.
    """
    read_batch = read_folder(path, batch_size=batch_size, chunk_size=chunk_size)
    t0 = time.time()
    queries = [("mock_read_id",
                channel, read.number,
                numpy.fromstring(read.raw_data, dtype=numpy.int16).tolist()
               )for (channel, read) in read_batch]
    result = classifer.batch_classify(queries)
    if result != None:
        for (status, channel, number, id) in result:
            print(channel, number, status)
    t1 = time.time()
    print('process{} reads in {} sec.'.format(len(queries), t1-t0))

def main():
    args = _get_parser().parse_args()
    logger = logging.getLogger('Manager')
    classifier = dyss.Dyss(num_scouts=args.num_scouts,
                           num_packs=args.num_packs,
                           reference=args.reference,
                           model=args.model,
                           param=args.param,
                           power=args.power,
                           referencesize=args.referencesize)
    print(classifier.classifier)
    mock_running(classifier, args.test, 500, args.min_chunk_size)
    classifier.free()

if __name__=="__main__":
    main()
