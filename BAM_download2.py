import os.path
import re
import sys
from subprocess import Popen
from subprocess import PIPE
from multiprocessing import Pool
delete_bams = False

def format_gtdownload(url):
    cmd = ""
    cmd += "ionice -c2 /pod/pstore/groups/brookslab/bin/GeneTorrent/GeneTorrent-3.8.7/cghub/bin/gtdownload -c "

    gtrepo = re.search("(gtrepo\S+?.com)", url).group(1)

    if gtrepo == "gtrepo-osdc-tcga.annailabs.com":
        cmd += "/scratch/jakutagawa/RNA-seq/keys/AKUTAGAWA_bionimbus_pcawg_tcga_180209.pem "
    else:
        cmd += ""

#    cmd += "-v "
    cmd += "%s" % url
    return cmd

def remove_bams(batch_file):
    #print ("Removing bams from file %s" % batch_file
    p = Popen("awk '{print $2}' %s | xargs -n 1 rm " % batch_file, shell=True)
    p.wait()
    p = Popen("awk '{print $2\".bai\"}' %s | xargs -n 1 rm " % batch_file, shell=True)
    p.wait()
    return

def run_md5(num,batch):
    p = Popen("python do-md5.py %s" % (batch), shell=True)
    p.wait()
    return


def runCMD(cmd):

    p = Popen(cmd, shell=True)
    p.wait()

def run_commands(threads,cmd_list):
    waitlist = []
    uuid_list = []
    jb_proc = None
    batch = 1

    p = Pool(threads)
    

    for i, _ in enumerate(p.imap_unordered(runCMD, cmd_list), 1):
        sys.stderr.write('Downloading bams...\rdone {0:%}'.format(i/len(cmd_list)))
    sys.stderr.write('\n')

    p.close()



urls = sys.argv[1]
infile = open(urls,'r')
download_commands = []

for line in infile:
    uuid = ((line.rstrip()).split("/"))[-1]
    download_commands.append(format_gtdownload(line.rstrip()))

run_commands(10,download_commands)
