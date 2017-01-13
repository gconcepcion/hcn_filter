#!/usr/bin/env python
'''
High Copy Number filter for pacbio raw input data
'''
import os
import sys
import glob
import subprocess
import tempfile
from functools import partial
import logging
import argparse
import collections
import multiprocessing as mp


#from signal import signal, SIGPIPE, SIG_DFL
# signal(SIGPIPE,SIG_DFL)

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


def setup_log(alog, level=logging.INFO, file_name=None, log_filter=None,
              str_formatter='[%(levelname)s] %(asctime)-15s ' \
			    '[%(name)s %(funcName)s %(lineno)d] %(message)s'):
    """Core Util to setup log handler

    :param alog: a log instance
    :param level: (int) Level of logging debug
    :param file_name: (str, None) if None, stdout is used, str write to file
    :param log_filter: (LogFilter, None)
    :param str_formatter: (str) log formatting str
    """
    alog.setLevel(logging.DEBUG)
    if file_name is None:
        handler = logging.StreamHandler(sys.stdout)
    else:
        handler = logging.FileHandler(file_name)
    formatter = logging.Formatter(str_formatter)
    handler.setFormatter(formatter)
    handler.setLevel(level)
    if log_filter:
        handler.addFilter(log_filter)
    alog.addHandler(handler)


def run_HPC_REPmask(dbpath):
    """
    Generate a set of daligner commands to compare the diagonal
    :param dbpath:
    :return:
    """

    prefix = os.path.basename(dbpath).lstrip('.db')
    cmd = ['HPC.REPmask',
           '-f{n}'.format(n=prefix),
           '-c10',
           '-g1',
           dbpath]

    log.info("Running HPC.REPmask: %s", cmd)
    subprocess.call(cmd)
    return


def get_input_quantity(database):
    """ parse DBstats output for total quantity of raw input bases"""
    cmd = ['DBstats -nu {d}'.format(d=database)]
    log.info("Calculating input quantity: %s", cmd)
    proc = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, _ = proc.communicate()

    basepairs = int([i for i in out.splitlines() if "base pairs" in i]
		    [0].split()[0].replace(',', ''))

    return basepairs


def run_DBsplit(database, genome_size):
    """Split the database up into blocks based on input coverage and expected genome_size"""
    input_quantity = get_input_quantity(database)
    coverage = float(input_quantity) / float(genome_size * 1000000)
    log.info("Approximate starting coverage given Genome Size: %s is %s", genome_size, coverage)
    block_size = int(float(input_quantity / coverage) / 1000000)
    cmd = ["DBsplit", "-f", "-s{c}".format(c=block_size), database]
    log.info("Running DBsplit: %s", cmd)
    subprocess.call(cmd)
    return


def run_fasta2DB(fasta_path):
    """ Convert fasta data to dazzler DB
    :param fasta_path: path to fasta file to be filtered
    :return:
    """
    prefix = os.path.basename(fasta_path).rstrip('.fasta')
    dbpath = os.path.join(os.getcwd(), '{p}.db'.format(p=prefix))
    if not os.path.exists(dbpath):
        cmd = ['fasta2DB', prefix, fasta_path]
        log.info("Running fasta2DB: %s", cmd)
        subprocess.call(cmd)
    else:
        log.info("%s already exists, skipping fasta2DB", dbpath)
    return dbpath


def run_cmds(commands, slots, ncpus, queue, local):
    """

    :param commands: list of commands to run
    :param slots: number of cores needed per task
    :param ncpus: Total number of cluster slots to request
    :param queue: what SGE queue to run on - default is 'default'
    :param local: run locally - currently not fully implemented
    :return:
    """
    log.info("Running %s grid jobs", commands)
    processes = set()
    max_processes = ncpus / slots

    if local is True:
        ##not working
        for cmd in commands[1:]:
            command = cmd.split()
            processes.add(subprocess.Popen([command, 'name']))
            if len(processes) >= max_processes:
                os.wait()
                processes.difference_update(
                    [p for p in processes if p.poll() is not None])
    else:
        log.info("Starting grid pool of size %s", max_processes)
        partial_qsub = partial(qsub_cmd, queue=queue)
        pool = mp.Pool(max_processes)
        pool.map(partial_qsub, commands[1:])
        pool.close()
    return


def qsub_cmd(cmd, queue):
    """
    Generate a qsub command for each task
    :param cmd: command to run
    :param queue: what queue to run on
    :return:
    """

    bash = []
    bash.append(cmd + '\n')
    if 'daligner' in cmd:
        slots = 4
    else:
        slots = 1
    with tempfile.NamedTemporaryFile(dir=os.getcwd()) as bashscript:
        bashscript.write("\n".join(bash))
        bashscript.flush()

        cmd = "qsub -sync y -q {q} -pe smp {s} -cwd -j y -V -S " \
		      "/bin/bash {x}".format(q=queue,
                                     s=slots,
                                     x=bashscript.name)
        log.info("Submitting: %s", cmd)
        subprocess.call(cmd, shell=True)

    return


def run_overlaps(ncpus, local, queue):
    """Run daligner overlap jobs

    :param ncpus: total number of cluster slots to request
    :param local: (bool) run locally or not
    :param queue: what queue to run on
    :return:
    """

    overlap_cli = glob.glob('*.OVL')
    slots = 4
    if not len(overlap_cli) > 1:

        with open(overlap_cli[0], 'r') as fin:
            daligner_cmds = fin.read().splitlines()

        run_cmds(daligner_cmds, slots, ncpus, queue, local)

    else:
        log.info('more than one overlap script')

    return


def run_LAcheck(ncpus, local, queue):
    """ run LAcheck on *.las files"""

    lacheck_cli = glob.glob('*.CHECK.OPT')
    slots = 1
    if not len(lacheck_cli) > 1:
        with open(lacheck_cli[0], 'r') as fin:
            lacheck_cmds = fin.read().splitlines()
        run_cmds(lacheck_cmds, slots, ncpus, queue, local)
    else:
        log.info('more than one lacheck script')
    return


def run_REPmask(ncpus, local, queue):
    """
    Run REPmask
    :param ncpus:
    :param local:
    :param queue:
    :return:
    """

    lamask_cli = glob.glob('*.MASK')
    slots = 1
    with open(lamask_cli[0], 'r') as fin:
        lamask_cmds = fin.read().splitlines()

    log.info("Running REPmask jobs")
    run_cmds(lamask_cmds, slots, ncpus, queue, local)
    return


def run_Catrack(tmpdir):
    """Run Catrack to merge tracks
    :param tmpdir:
    :return:
    """
    database = glob.glob(os.path.join(tmpdir, '*.db'))[0]

    cmd = ['Catrack {d} rep1'.format(d=database)]
    log.info("Running Catrack")
    subprocess.call(cmd, shell=True)
    return


def dump_readids(tmpdir, percent=0.8):
    """ Dump readids from database
    :param tmpdir:
    :param percent:
    :return:
    """
    database = glob.glob(os.path.join(tmpdir, '*.db'))[0]
    cmd = 'DBdump -mrep1 -r -h {d}'.format(d=database)
    proc = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, _ = proc.communicate()
    readlist = dump_to_readlist(out)
    readids = []
    for read in readlist:
        if (float(read.coverage) / float(read.length)) >= percent:
            readids.append(read.readname.rstrip())
    log.info("Excluding %s high copy number reads", len(readids))

    return readids


def dump_to_readlist(dump):
    '''
    read 'P' and 'C' lines  from LAdump
    returns alignment within same ctg
    '''
    lasAlign = collections.namedtuple(
        'lasAlign', 'readname coverage gap egap length')
    readlist = []
    for line in dump.splitlines():
        fin = line.strip().split()
        if fin[0] == 'R':
            # f[1] is read index
            readID = int(fin[1])
        elif fin[0] == 'H':
            movie = fin[2]
        elif fin[0] == 'L':
            holenum = fin[1]
            fr = int(fin[2])
            to = int(fin[3])
        elif fin[0] == 'T0':
            nitv = int(fin[1])
            intervals = [(int(fin[2 * i]), int(fin[2 * i + 1]))
                         for i in range(1, nitv + 1)]
            name = movie + '/' + holenum + '/' + str(fr) + '_' + str(to)
            length = to - fr
            cov = sum(itv[1] - itv[0] for itv in intervals)
            if nitv != 0:
                prev_itv = intervals[0]
                gap = 0
                for itv in intervals[1:]:
                    gap = max(itv[0] - prev_itv[1], gap)
                    prev_itv = itv
                edge_gap = max(intervals[0][0], len - intervals[-1][1])
            else:
                gap = length
                edge_gap = length
            readlist.append(lasAlign(readname=name, coverage=cov,
                                     gap=gap, egap=edge_gap, length=l))
    return readlist


def filter_fasta(fasta_path, readids):
    """Read input fasta and write output by filtering HCN reads
    """
    output_path = fasta_path.replace('.fasta', '_filtered.fasta')
    excluded_ids = fasta_path.replace('.fasta', '_excluded.ids')
    log.info("Writing filtered output: %s", output_path)
    log.info("Writing excluded readids: %s", excluded_ids)

    with open(excluded_ids, 'w') as excl:
        excl.write("\n".join(readids))

    with open(output_path, 'w+') as output:
        with open(fasta_path, 'r') as fin:
            for line in fin:
                if line.startswith('>'):
                    _, readid = line.split('>')
                    if readid.rstrip() not in readids:
                        out = True
                        output.write(line)
                    else:
                        out = False
                else:
                    if out is True:
                        output.write(line)


def main():
    """Main function to drive the script"""

    parser = get_parser()
    args = parser.parse_args()
    fasta_path = os.path.abspath(args.fastafile)
    genome_size = args.genome_size
    ncpus = args.ncpus
    local = args.local
    queue = args.queue
    job_root = os.getcwd()
    prefix = os.path.basename(fasta_path).rstrip('.fasta')
    tmpdir = os.path.join(job_root, 'tmp_{p}'.format(p=prefix))

    if args.debug:
        setup_log(log, level=logging.INFO)
    else:
        log.addHandler(logging.NullHandler())

    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    os.chdir(tmpdir)

    dbpath = run_fasta2DB(fasta_path)
    run_DBsplit(dbpath, genome_size)
    run_HPC_REPmask(dbpath)
    run_overlaps(ncpus, local, queue)
    run_LAcheck(ncpus, local, queue)
    run_REPmask(ncpus, local, queue)
    run_Catrack(tmpdir)
    read_ids = dump_readids(tmpdir)
    filter_fasta(fasta_path, read_ids)
    return


def get_parser():
    """Return an argparse instance"""
    parser = argparse.ArgumentParser()
    parser.add_argument('fastafile', action='store',
                        help="Path to fasta file to filter.")
    parser.add_argument('--genome_size', dest='genome_size', action='store',
                        type=int, default=10, help="Genome Size (in megabasases)")
    parser.add_argument('--ncpus', dest='ncpus', action='store',
                        type=int, default=12, help="# of cores to use. Daligner uses 4 by default")
    parser.add_argument('--queue', dest='queue', action='store', default='default',
                        help="If not running locally, the queue you want to run in")
    parser.add_argument('--local', dest='local', action='store_true',
                        help="Run locally, or on SGE")
    parser.add_argument('--debug', dest='debug', action='store_true',
                        help="Turn on verbose logging mode")
    return parser


if __name__ == "__main__":
    sys.exit(main())
