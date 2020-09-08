""" This is a modified version of simple.py script bundled with Read Until API"""
import argparse
import logging
import sys
import traceback
import time
import numpy
import read_until
import cffi
import os
import h5py
import glob
import concurrent.futures
import dyss

def _get_parser():
    parser = argparse.ArgumentParser('Dyss -- a tiny program for selective sequencing on MinION')
    parser.add_argument('--port', type=int, default=8000,
                        help='MinKNOW server port.')
    parser.add_argument('--analysis_delay', type=int, default=1,
                        help='Period to wait before starting analysis.')
    parser.add_argument('--run_time', type=int, default=900,
                        help='Period to run the analysis.')
    parser.add_argument('--min_chunk_size', type=int, default=3500,
                        help='Minimum read chunk size to receive.')
    parser.add_argument('--control_group', default=2, type=int,
                        help='Inverse proportion of channels in control group.')
    parser.add_argument('--batch_size', default=30, type=int,
                        help='Inverse proportion of channels in control group.')
    parser.add_argument(
        '--debug', help="Print all debugging information",
        action="store_const", dest="log_level",
        const=logging.DEBUG, default=logging.WARNING,
    )
    parser.add_argument(
        '--verbose', help="Print verbose messaging.",
        action="store_const", dest="log_level",
        const=logging.INFO,
    )
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
    parser.add_argument('--referencesize', default=400000, type=int,
                        help='Reference size(Event num = 2*bp). Default is 400000')
    return parser


def signal_based_analysis(client, classifier, batch_size=30, delay=1, throttle=0.5, control_group=16):
    """A tiny analysis function based on raw signal comparison.
    
    :param client: an instance of a `ReadUntilClient` object.
    :param batch_size: number of reads to pull from `client` at a time.
    :param delay: number of seconds to wait before starting analysis.
    :param throttle: minimum interval between requests to `client`.
    :param dyss: an instance of Dyss object constructed by libdyss.construct_dyss().
    :param debug: flag whether or not output every query into stdout.
    """
    logger = logging.getLogger('Dyss')
    logger.warn(
        'Initialization of Dyss classification'
        'When debug and verbose flag is on, it generates literaly all inputs.'
        'If you want to apply this application to real sequencing experiment,'
        'it is highly recommended to turn these flags off.'
    )
    # we sleep a little simply to ensure the client has started initialised
    logger.info('Starting analysis of reads in {}s.'.format(delay))
    time.sleep(delay)
    
    while client.is_running:
        # If thre are too many queries, reject them all.
        if client.queue_length > 300:
            read_batch = client.get_read_chunks(batch_size = client.queue_length, last = True)
            for (channel, read) in read_batch:
                read.raw_data = read_until.NullRaw
                if channel % control_group != 0:
                    client.unblock_read(channel, read.number)
                client.stop_receiving_read(channel, read.number)
        t0 = time.time()
        # Then, running usual classification step
        read_batch = client.get_read_chunks(batch_size=batch_size, last=True)
        # convert the read data into a numpy array of correct type
        queries = [(read.id, channel,read.number, 
                    numpy.fromstring(read.raw_data, client.signal_dtype).tolist()
        ) for (channel,read) in read_batch if channel % control_group != 0]
        querylen = len(queries)
        # clear the raw reads from allocated memory
        for (channel,read) in read_batch:
            read.raw_data = read_until.NullRaw
            if channel % control_group == 0:
                client.stop_receiving_read(channel, read.number)
        result = classifier.batch_classify(queries)
        if result is not None:
            for (status,channel,number,id) in result:
                if status == 0:
                    # The i th read doesn't seems to be in the target region. Reject it.
                    client.unblock_read(channel, number)
                    client.stop_receiving_read(channel, number)
                    logger.info('Rejected {} {} {}'.format(id,channel,number))
                elif status  == 1:
                    # The i th read seems to be in the target region.
                    client.stop_receiving_read(channel, number)
                    logger.info('Accepted {} {} {}'.format(id,channel,number))
                else:
                    logger.info('Chunked {} {} {}'.format(id,channel, number))
                    # else, the i th read didn't have enough signal. Keep going.
        # limit the rate at which we make requests            
        t1 = time.time()
        if t0 + throttle > t1:
            time.sleep(throttle + t0 - t1)
        logger.info('process {} reads in {}.'.format(querylen, t1-t0))
        
    logger.info('Finished analysis of reads.')
    classifier.free()



def main():
    args = _get_parser().parse_args()
    logging.basicConfig(format='[%(asctime)s - %(name)s] %(message)s',
                        datefmt='%H:%M:%S', level=args.log_level)
    logger = logging.getLogger('Manager')
    classifier = dyss.Dyss(num_scouts=args.num_scouts,
                     num_packs=args.num_packs,
                     reference=args.reference,
                     model=args.model,
                     param=args.param,
                     power=args.power,
                     referencesize=args.referencesize)

    read_until_client = read_until.ReadUntilClient(
        mk_port=args.port, one_chunk=False, filter_strands=True
    )
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = list()
        futures.append(executor.submit(
            read_until_client.run, runner_kwargs={
                'run_time':args.run_time, 'min_chunk_size':args.min_chunk_size
            }
        ))
        futures.append(executor.submit(
            signal_based_analysis, read_until_client, classifier,
            batch_size=args.batch_size, delay=args.analysis_delay,control_group=args.control_group
        ))
        for f in concurrent.futures.as_completed(futures):
            if f.exception() is not None:
                logger.warning(f.exception())
                
if __name__=="__main__":
    main()
