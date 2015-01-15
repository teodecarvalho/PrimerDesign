# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Design of Primers for Functional Genes

# <markdowncell>

# Author: Teotonio Soares de Carvalho

# <markdowncell>

# Please contact me at teotonio@msu.edu if you have any suggestions or questions.

# <headingcell level=2>

# Setup for Parallel Processing

# <codecell>

# Requires an IPcluster running with at least one cluster.
# An error will be raised if no cluster is found.
try:
    import os
    from IPython.parallel import Client
    rc = Client()
    cwd = os.getcwd() 
    # To make sure that the works are operating on the same work dir.
    rc[:].execute('import os;os.chdir(%s)'%cwd)
    dview = rc[:]
    dview.block = False
    lview = rc.load_balanced_view()
except:
    raise Exception("Please, start an IPython cluster to proceed.")

# <codecell>

def wait_on(ar, verbose = True):
    """
    Tracks the progress of a task running on a IPcluster. Downloaded from the internet.
    """
    from datetime import datetime
    N = len(ar.msg_ids)
    rc = ar._client
    submitted = rc.metadata[ar.msg_ids[0]]['submitted']
    while not ar.ready():
        ar.wait(1)
        msgs = [msg_id not in rc.outstanding for msg_id in ar.msg_ids]
        progress = sum(msgs)
        dt = (datetime.now()-submitted).total_seconds()
        if verbose:
            clear_output()
            print "%3i/%3i tasks finished after %4i s" % (progress, N, dt),
            sys.stdout.flush()
    if verbose:
        print
        print "done"

# <headingcell level=3>

# Import modules

# <codecell>

import time, os, shutil, math, sys, copy, random, time 
import regex    # Make sure you have installed this module. Use: pip install regex
import pandas   # Make sure you have installed this module as well and its dependencies. 
import numpy as np  # If you installed pandas successfully, this package is already installed
from selenium import webdriver # Selenium and webdriver are installed separately
from selenium.webdriver.support.select import Select
from Bio.Seq import Seq  # This also have to be installed: pip install Biopython
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from IPython.core.display import clear_output
from itertools import combinations, permutations, izip, product, izip_longest
import datetime as dt
from functools import partial
import primer3
from IPython.core.display import clear_output
import commands
with dview.sync_imports():
    from functools import partial
    import numpy
    import primer3
    from Bio import SeqIO, Seq
    from string import Template
    import regex
    import os
    from Bio.Seq import Seq
    import pandas
    from itertools import combinations, izip, product, permutations

# <headingcell level=2>

# Download Archaea amoA Sequences from Fungene

# <codecell>

def load_driver(gene_url):
    """
    Loads the fungene web page in chrome and assign it to the global variable "driver".
    I could not use the module ghost, which is faster, because the repository and the 
    analysis page in fungene are in separate tabs.
    Requires webdriver, chrome, and selenium.
    """
    global driver
    if "driver" in dir(): 
        driver.quit()
    driver = webdriver.Chrome()
    driver.get(gene_url);
    driver.find_element_by_partial_link_text("Display Options").click()
    seq = Select(driver.find_element_by_id("seqsPerPage"))
    seq.select_by_value("2000")
    driver.find_element_by_id("displayCmd").submit()

# <codecell>

def filter_by_score(score):
    """
    Selects the minimum score as 400 in the repository page.
    """
    form = driver.find_element_by_name("hmmId")
    driver.find_element_by_link_text("Show/Hide filter options").click()
    min_score = driver.find_element_by_name("min_bits")
    min_score.send_keys(str(score))
    form.submit()

# <codecell>

def find_last_page():
    """
    Detects the total number of pages using regular expression.
    """
    pages = regex.findall(r"page=([0-9]*)", driver.page_source)
    if not pages:
        # Added to correct an error when there is only one page
        last_page = 1
    else:
        last_page = max([int(page) for page in pages])
    return last_page

# <codecell>

def load_page_repository(page_number, gene_url):
    """
    Loads a page, specified by "page_number", for the archaeal amoA
    from the fungene repository and select all the sequences in it.
    Arguments:
        - page_number: an integer specifying the desired page.
    """
    url = gene_url + "&page=%d"%page_number
    driver.get(url)
    driver.find_element_by_link_text("Select Entire Page").click()
    get_nucl2prot_accession()

# <codecell>

def remove_fungene(download_path):
    """
    Removes files previously download from fungene in the default
    download directory used by chrome.
    Arguments:
        - download_path: a string indicating  the default download folder.
    """
    for file in os.listdir(download_path):
        if regex.match(r"fungene.*?aligned_(nucleotide)*(protein)*_seqs", file):
            os.remove(download_path + file)

# <codecell>

def get_nucl2prot_accession():
    """
    Extracts the respective acc for protein and nucleotide for each sequence
    and saves the results in a file "./data/nucleotide2protein" where each
    line is as follow:
    acc for protein | acc for nucleotide
    """
    reg_exp = regex.compile(r"gpprotdata.jsp\?seqAccno=([0-9A-Z]+).+?"
                             "gbnucdata.jsp\?seqAccno=([0-9A-Z]+)", 
                             regex.DOTALL|regex.MULTILINE|regex.VERBOSE)
    accession = reg_exp.findall(driver.page_source)
    accession = "\n".join(["|".join(pair) for pair in accession])
    with open("./data/download/nucleotide2protein", "a") as handle:
        handle.write(accession + "\n")

# <codecell>

def download_hmm(default_download_path):
    """
    Download the hmm for AOA and copy it to ./data/
    Arguments:
        - download path: string indicating the default download folder
    """
    for file in os.listdir(default_download_path):
        if ".hmm" in file:
            os.remove(default_download_path + file)
    driver.find_element_by_link_text("(download HMM)").click()
    time.sleep(1)
    while True:
        files = os.listdir(default_download_path)
        file = [file for file in files if ".hmm" in file]
        if file:
            file = file[0]
            break
    time.sleep(2)
    shutil.copy(default_download_path + file, "./data/download/" + file)

# <codecell>

def switch_to_analysis_window():
    """
    Switchs from repository windows to the analysis windows.
    Fails if the analysis windows is not already open.
    """
    driver.switch_to.window(driver.window_handles[1])

# <codecell>

def switch_to_repository_window():
    """
    Switches back to the repository windows.
    """
    driver.switch_to.window(driver.window_handles[0])

# <codecell>

def deselect_all_sequences():
    """
    Deselect sequences already downloaded in the repository page.
    """
    try:
        driver.find_element_by_link_text("Deselect All Sequences").click()
    except:
        pass

# <codecell>

def download_sequences(default_download_path, count):
    """
    Initiate the "Begin Analysis" link for the fungene repository to
    download the sequences. Downloads the unaligned nucleotides.
    """
    driver.find_element_by_link_text("Begin Analysis").click()
    switch_to_analysis_window()
    for seq_type in ["Nucleotide", "Protein"]:
        driver.find_element_by_id("download_%s_seqs"%seq_type).click()
        # Uncheck the aligned option if checked. 
        aligned_option = driver.find_element_by_id("aligned1")
        if aligned_option.is_selected(): 
            aligned_option.click()
        # Download the sequences
        driver.find_element_by_name("download").click()
        # If the internet is slow it is better to increase this time
        time.sleep(2) 
        move_file_to_data(default_download_path, count, seq_type)
    switch_to_repository_window()
    deselect_all_sequences()

# <codecell>

def move_file_to_data(default_download_path, count, seq_type):
    exp = r"fungene.*?aligned_%s_seqs"%seq_type.lower()
    files = os.listdir(default_download_path)
    file = [file for file in files if regex.match(exp, file)][0]
    with open(default_download_path + file, "r") as handle:
        fasta = handle.read()
    os.remove(default_download_path + file)
    file_name = "./data/download/%s_%d"%(seq_type.lower(), count)
    with open(file_name, "w") as handle:
        handle.write(fasta)

# <codecell>

def gather_all_fasta(protein = False):
    """
    Gather all download nucleotide unaligned sequences in one fasta
    file at "./data/arch_amoa_all.fasta"
    Arguments:
        - download path: string indicating the default download folder
    """
    if protein:
        seq_type = "protein"
    else:
        seq_type = "nucleotide"
    with open("./data/download/all_%s"%seq_type, "w") as handle:
        for file in os.listdir("./data/download/"):
            if regex.match(r"%s_[0-9]+"%seq_type, file):
                with open("./data/download/" + file, "r") as handle_split:
                    handle.write(handle_split.read())
                

# <codecell>

def main_download(hmm_id, default_download_path, score):
    """
    Downloads the proteins, nucleotides, hmm file and the correspondence between acession numbers from proteins
    and nucleotides.
    Arguments:
        -hmm_id: an integer giving the fungene hmm_id for a particular gene. You can get it by clicking in the gene link
                 in the fungene database. In the url that appears in the address bar you will see the hmm_id.
        -default_download_path: string giving the default path for download for chrome.
        -score: an integer giving the minimum hmm score for the sequences. Be careful with this parameter as it is very
                gene dependent. That is why there is no default value.
    This function uses chromedriver to automate the download. While the browser is working to download the files you should 
    not interact with it or unexpected results may arise. It usually takes about 2 minutes to complete the download on a 
    good network. If you notice that any of the download fails, remove the data folder, restart the notebook and run again.
    """
    # Remove previous files to avoid errors
    if os.path.isdir("./data/"):
        shutil.rmtree("./data/")
    os.makedirs("./data/")
    os.makedirs("./data/download/")
    gene_url = "http://fungene.cme.msu.edu/hmm_details.spr?hmm_id=%d"%hmm_id
    remove_fungene(default_download_path)
    load_driver(gene_url)
    download_hmm(default_download_path)
    filter_by_score(score)
    last_page = find_last_page()
    deselect_all_sequences()
    count = 1
    for page_number in range(1, last_page + 1):
        load_page_repository(page_number, gene_url)
        # Download the sequences when repository page is multiple 
        # of 5 or the last page
        if not page_number % 5 or page_number == last_page:
            download_sequences(default_download_path, count)
            time.sleep(1)
            count += 1
    gather_all_fasta()
    gather_all_fasta(protein = True)

# <headingcell level=2>

# Align Protein

# <codecell>

def read_fasta(filename):
    """
    Reads a fasta file and return it as a list of records.
    Arguments:
        - filename: string indicating the name of the file 
                    including the path.
    Returns a list of records formatted by BioPython.
    """
    with open(filename, "rU") as handle:
        records = list(SeqIO.parse(handle, "fasta"))
    return records

def write_fasta(filename, records_list):
    """
    Takes a list of records (BioPython) and writes it to filename.
    """
    with open(filename, "w") as handle:
        fasta_writer = SeqIO.FastaIO.FastaWriter(handle)
        fasta_writer.write_file(records_list)

def make_dict_records(fasta_file_name):
    """
    Returns a dict of records where the keys are the accession
    number and the values are the records (Biopython)
    """
    records = read_fasta(fasta_file_name)
    records_dict = {record.name:record for record in records}
    return records_dict

# <codecell>

def align_protein():
    """
    Takes the unaligned proteins in ./data/download and align them using the 
    hmm profile.
    """
    file = [file for file in os.listdir("./data/download/") if ".hmm" in file][0]
    cmd = ("hmmalign " 
           "./data/download/%s "
           "./data/download/all_protein "
           "> ./data/download/aligned_prot ")%file
    print "Aligning sequences"
    process = os.system(cmd)
    if process:
        print cmd
        raise
    print "Reading Alignment"
    alignment = AlignIO.read(open("./data/download/aligned_prot"), "stockholm")
    print "Writing Alignment"
    write_fasta("./data/download/aligned_prot", alignment)
    sys.stdout.flush()

# <headingcell level=2>

# Align Nucleotides Using Proteins

# <codecell>

def get_nuc2prot():
    """
    Returns a dict of nucleotide accessions numbers as keys and 
    protein acession numbers as values.
    """
    nuc2prot_acc = {}
    with open("./data/download/nucleotide2protein", "r") as handle:
        line = handle.readline()
        while line:
            prot, nuc = line.split("|")
            nuc2prot_acc[nuc[:-1]] = prot
            line = handle.readline()
    return nuc2prot_acc

# <codecell>

def align_nucleotides():
    """
    Takes the protein alignment positions and align the nucleotide codons
    based on that.
    Raises an error if the frame is not corrected.
    """
    nucleotides = make_dict_records("./data/download/all_nucleotide")
    proteins = make_dict_records("./data/download/aligned_prot")
    nuc2prot_acc = get_nuc2prot()
    aligned_nucleotides = []
    for nuc_acc, nucleotide in nucleotides.iteritems():
        protein = proteins[nuc2prot_acc[nuc_acc]]
        prot_unaligned = regex.sub(r"[-\.]", "", str(protein.seq))
        prot_transl = str(nucleotide.seq.translate())
        #Checking if frames are right
        nuc_seq = str(nucleotide.seq)
        protein_seq = str(protein.seq)
        codons = [nuc_seq[i: i+3] for i in xrange(0, len(nuc_seq), 3)]
        aligned_seq = []
        for position in protein_seq:
            if position == "." or position == "-":
                aligned_seq += ["---"]
            else:
                if codons:  
                    codon = codons.pop(0)
                    aligned_seq += codon
                else:
                    aligned_seq += ["---"]
        nucleotide.seq = Seq("".join(aligned_seq))
        aligned_nucleotides += [nucleotide]
        if not prot_transl.upper()[1:-1] in prot_unaligned.upper():
            print "Incorrect frame!"
            raise
    write_fasta("./data/download/aligned_nucleotides", aligned_nucleotides)

# <headingcell level=2>

# Trimm Alignment

# <codecell>

def count_invalid_pos(file_name = None, records = None):
    """
    Takes an aligment and counts the number of gaps "." or "-".
    Arguments:
        - file_name: if given, the function reads the fasta file.
        - records: a list of records from an aligment (optional)
    Only one of these two arguments must be provided.
    """
    if not records:
        records = read_fasta(file_name)
    nseq = len(records)
    seqs = [str(record.seq.upper()) for record in records]
    align_pos = izip(*seqs)
    count_pos = {}
    for count, bases in enumerate(align_pos):
        count_pos[count] = len(regex.findall(r"[-\.]", "".join(bases)))
    count = nseq - pandas.Series(count_pos)
    return (count, nseq)

# <codecell>

def trimm_columns(records, prop_non_gap = .001, column_info_prop_discard = .3):
    """
    Removes columns of the aligment in both ends where the proportion of
    gaps is higher than 30%. It DOESN'T remove positions in the middle of the
    aligment.
    """
    count, nseq = count_invalid_pos(records = records)
    drop_threshold = prop_non_gap * nseq
    richer_pos = count[count / float(nseq) > column_info_prop_discard]
    start = richer_pos.index[0]
    end = richer_pos.index[-1]
    count_richer_pos = count[start:end + 1]
    positions_to_drop = list(count_richer_pos[count_richer_pos <= drop_threshold].index)
    positions_to_keep = [i for i in range(start, end + 1)]
    trimmed_records = []
    for record in records:
        seq = str(record.seq)
        seq = "".join([seq[i] for i in positions_to_keep])
        record.seq = Seq(seq)
        trimmed_records += [record]
    return trimmed_records

# <codecell>

def dealign_seq(in_file_name = "./data/trimmed_align",
                out_file_name = "./data/unaligned_trimmed"):
    """
    Removes gaps from sequences.
    """
    records = read_fasta(in_file_name)
    for record in records:
        record.seq = Seq(regex.sub(r"[\.-]", "", str(record.seq))).upper()
        record.description = ""
    write_fasta(out_file_name, records)

# <codecell>

def trimm_records(max_prop_gap_ends = .1, column_info_prop_discard = .3):
    """
    Removes positions of the alignment with less than a specified threshold of information (bases) and
    discards sequences (after trimming the positions) with more a threshold of its length as gaps at any
    of its ends.
    Arguments:
        -column_info_prop_discard: minimal allowed proportion of information (bases) in the columns at both ends
                                   of the aligment. The purporse is to delete columns in both ends for which information
                                   is not available for most sequences.
        -max_prop_gap_ends: maximal proportion of the length a sequence (after trimming columns as above) that consists of
                            continuous gaps in any of the ends of the sequence. The purpose is to remove sequences are too
                            short and cannot be used to test the primers.
    """
    records = read_fasta("./data/download/aligned_nucleotides")
    trimmed_columns = trimm_columns(records, column_info_prop_discard = column_info_prop_discard)
    npos = len(trimmed_columns[0].seq)
    max_gap_ends = max_prop_gap_ends * npos
    records_to_drop = 0
    records_to_keep = []
    dropped = 0
    for count, record in enumerate(trimmed_columns):
        seq = str(record.seq)
        gaps = regex.search("(?P<start>^[-\.]*).+?(?P<end>[-\.]*$)", seq)
        if (len(gaps["start"]) >= max_gap_ends) or (len(gaps["end"]) >= max_gap_ends):
            dropped += 1
            records_to_drop += 1
        else:
            records_to_keep += [record]
    write_fasta("./data/trimmed_align", records_to_keep)
    print "%d sequences had more gaps in their ends than the specified threshold and were deleted."%dropped

# <headingcell level=2>

# Remove redundancy at 100% Similarity

# <codecell>

def cluster_sequences(input_file, output_file, similarity, word_size):
    """
    Cluster all sequences and write two files:
        - ./data/arch_amoa_repr.clstr with information for each cluster
        - ./data/arch_amoa_repr a fasta file with representatives (not used)
    """
    if os.path.isfile(output_file):
        os.remove(output_file)
    cmd = Template(
     "cd-hit-est -i $input_file "   # Input File
                "-o $output_file "  # Output File
                "-c $similarity "   # Threshold similarity
                "-n $word_size "    # Word size (see manual)
                "-T 0 "             # Use all processors
                "-M 0")             # Use all memory
    cmd =  cmd.substitute(input_file = input_file,
                          output_file = output_file,
                          similarity = similarity,
                          word_size = word_size)
    print cmd
    process = os.system(cmd)
    if process:
        print cmd
        raise

# <headingcell level=2>

# Create fasta file for each cluster

# <markdowncell>

# This code was reused from previous versions. That is the reason for the fancy Arch_Group class

# <codecell>

def make_split_groups_out(file_name):
    """
    Breaks the output of cd-hit and return it as a list of
    strings were each element correspond to a group.
    """
    with open(file_name, "r") as handle:
        groups = handle.read()
    groups_split = groups.split(">Cluster")[1:]
    return groups_split

# <codecell>

class Arch_Group:
    """
    Used to parse the results from cd-hist. An instance of Arch_Group is a group of
    nucleotide sequences that are 97% similar. 
    The instance have the following properties:
        - name: the name of the representative sequence
        - representative(deprecated): same as above
        - n_members: number of sequences in group
        - members: a list of strings with members accession number
        - nuc_members: list of integers giving the number of nuc for 
                       each member.
        - n_invalid(deprecated): number of non ACTG characters in the representative seq
        - seq(deprecated): the seq of the representative.
    Note:
        The representative is not used in the subsequent analysis. Instead of using
        a representative sequence, I use the consensus in the next steps. But the 
        groups will be identified by the representative acc number.
    """
    def __init__(self,group, dict_fasta):
        self.get_representative(group)
        self.name = self.representative
        self.get_members(group)
        self.get_nuc_numbers(group)
        self.n_members = len(self.members)
        self.seq = dict_fasta[self.name].seq
        self.get_n_invalid_bases_representative()
    
    def __repr__(self):
        return "Group: %s, %d sequence(s)"%(self.name, self.n_members)
    
    def get_representative(self, group):
        repr_id = regex.findall(r"\>([A-Z|0-9]*?[_0-9]*)\.\.\.\ \*", group)[0]
        self.representative = repr_id
    
    def get_members(self, group):
        seqs_id = regex.findall(r"\>([A-Z|0-9]*?[_0-9]*)\.\.\.", group)
        self.members = seqs_id
            
    def get_nuc_numbers(self, group):
        nuc_numbers = regex.findall(r"([0-9]{1,4})nt,", group)
        nuc_numbers = [int(num) for num in nuc_numbers]
        self.nuc_numbers = nuc_numbers
    
    def get_sequences(self, records_dict):
        self.sequences = [records_dict[id] for id in self.sequences_ids]
        
    def get_n_invalid_bases_representative(self):
        n_invalid = len(regex.findall(r"[^atgc]", str(self.seq)))
        self.n_invalid = n_invalid
    
    def get_consensus_record(self):
        record = read_fasta("./data/consensus/%s.fasta"%self.name)
        self.consensus_record = record
    def get_members_record_list(self):
        records = read_fasta("./data/groups/%s.fasta"%self.name)
        self.members_record_list = records

# <codecell>

def make_dict_groups(group_file_name = "./data/cluster_97.clstr",
                     fasta_file_name = "./data/contigs_100"):
    """
    Parses the arch_amoa_repr.clstr and returns a dict
    where the key is the name of the group and the values
    are instances of the Arch_Group
    """
    dict_fasta = make_dict_records(fasta_file_name)
    split_groups_out = make_split_groups_out(group_file_name)
    groups = [Arch_Group(group, dict_fasta) for group in split_groups_out]
    groups = {group.name:group for group in groups}
    return groups

# <codecell>

def make_fasta_for_groups(dict_groups, dict_fasta, path_name):
    """
    Creates one fasta file for each groups containing the sequences
    for that group.
    Arguments:
        - dict_groups: a dict of groups as returned by make_dict_groups 
        - dict_fasta: a dict of sequences
        - path_name: the path where the fasta file for the grouped must be saved
    """
    if os.path.isdir(path_name):
        shutil.rmtree(path_name)
    os.makedirs(path_name)
    for group_name, group in dict_groups.iteritems():
        records = [dict_fasta[member] for member in group.members]
        write_fasta(path_name + group_name, records)

# <codecell>

def main_make_fasta_from_groups(group_name, fasta_file_name, path_name):
    groups = make_dict_groups(group_name, fasta_file_name)
    dict_fasta = make_dict_records(fasta_file_name)
    make_fasta_for_groups(groups, dict_fasta, path_name)

# <headingcell level=2>

# Step 05. Make Consensus Sequences

# <codecell>

def make_consensus(group_name, aligned_path, consensus_path, plurality):
    """
    Find the consensus sequence for each aligned fasta file
    in ./data/aligned/
    Arguments:
        - group_name: string indicating the name of the group
        - aligned_path: string indicating the path where the aligned sequences are.
        - consensus_path: string indicating where the consensus are.        
        - plurality: the minimal proportion of agreement that must be at a given position for
                     the consensus to receive the most frequent base at that position.

    Writes the consensus sequences to ./data/consensus
    """
    n_seq = len(read_fasta(aligned_path + group_name))
    min_agreement = int(plurality * n_seq)
    cmd = Template("cons "                         # From emboss
           "$aligned_path$group_name "      # Input fasta
           "$consensus_path$group_name "    # Output fasta
           "-name $group_name "  # The consensus sequence is named with the group name
           "-plurality $min_agreement " # See description of the function above
          )
    cmd = cmd.substitute(aligned_path = aligned_path,
                         consensus_path = consensus_path,
                         group_name = group_name,
                         min_agreement = min_agreement)
    process = os.system(cmd)
    if process:
        raise RuntimeError('program {} failed!'.format(cmd))

# <codecell>

def gather_all_consensus(consensus_path, consensus_name):
    """
    Arguments:
        Gather the individual consensus file in a unique file.
        - consensus_path: the path where the consensus file are.
        - consensus_name: the name of the resulting consensus file.
    """
    records = ""
    for file in os.listdir(consensus_path):
        with open(consensus_path + file, "r") as handle:
            records += handle.read()
    with open(consensus_name, "w") as handle:
        handle.write(records)

# <codecell>

def main_parallel_consensus(consensus_path, aligned_path, 
                            consensus_name, fasta_file_name, 
                            group_file_name,
                            plurality):
    """
    Run the previous functions in parallel.
    Arguments:
        - consensus_path: string indicating the path where the individual consensus are.
        - aligned_path: the path were the individual aligments are.
        - consensus_name: the name of the consensus file with all consensus sequences.
        - fasta_file_name: the original sequences from which the consensus were made.
        - group_file_name: the file indicating groups from cd-hit
        - plurality: the minimal proportion of agreement that must be at a given position for
                     the consensus to receive the most frequent base at that position.
    """
    if os.path.isdir(consensus_path):
        shutil.rmtree(consensus_path)
    os.makedirs(consensus_path)
    groups = make_dict_groups(group_file_name = group_file_name,
                              fasta_file_name = fasta_file_name)
    cmd_template = Template("cp $aligned_path$group_name $consensus_path$group_name")
    for group in groups.values():
        if group.n_members == 1:
            cmd = cmd_template.substitute(consensus_path = consensus_path, 
                                         aligned_path = aligned_path,
                                         group_name = group.name)
            os.system(cmd)
    # Clusters are only sent to work if the number of members > 1
    clusters = [group_name for group_name, group in groups.iteritems()\
                           if group.n_members > 1]
    # Import os in all engines (nodes)
    dview.execute("import os")
    # Send the make_consensus function to all engines
    dview.push({"make_consensus":make_consensus,
                "read_fasta":read_fasta})
    # Map/Reduce
    task = lview.map(partial(make_consensus,
                             aligned_path = aligned_path, 
                             consensus_path = consensus_path, 
                             plurality = plurality),
                             clusters)
    wait_on(task)
    if not task.successful():
        raise Exception("Consensus failed!")
    gather_all_consensus(consensus_path = consensus_path, 
                         consensus_name = consensus_name)

# <headingcell level=3>

# Classify Sequences

# <codecell>

def get_seqs_id_from_pester():
    """ 
    Extracts relevant information about sequences (name, full_name, acc,
    Taxon_Level1 and Seq_Len) as given by Pester et al 2012 in the arb
    database provided as supplementary material.
    This function only applies for the analysis of Archaea amoA
    """
    # Read data from the NDS file exported from arb
    seq_all = pandas.read_table("./Pester Consensus Data/pester.nds", 
                             header=None)
    # Rename the columns
    seq_all.columns = ["name", "full_name", "taxonomy", "acc", 
                       "Taxon_Level_1", "Taxon_Level_2", "Taxon_Level_3", 
                       "Seq_Len", "Habitat", "unknown"]
    # Select rows from sequences in the tree.
    selected_rows = seq_all.Taxon_Level_1.notnull()
    # Select relevant columns for subsequent analysis
    selected_columns = ["name", "Taxon_Level_1", "Taxon_Level_2", "Taxon_Level_3"]
    seq_valid = seq_all.ix[selected_rows, selected_columns]
    return seq_valid

# <headingcell level=3>

# Screen Oligos

# <codecell>

def enumerate_oligos(starts, kmer_sizes, seqs_all, look_ahead = 50):
    """
    Enumerate all possible oligos (kmers) with sizes kmer_sizes from
    the aligment.
    """
    unique = {}
    for start in starts:
        # To account for gaps, I choose a region much bigger than kmer_size
        seqs = [seq[start:start + look_ahead].replace("-", "") for seq in seqs_all]
        for kmer_size in kmer_sizes:
            oligos = [seq[:kmer_size] for seq in seqs]
            for count, oligo in enumerate(oligos):
                if oligo:
                    if not oligo in unique:
                        unique[oligo] = [count]
                    else:
                        unique[oligo] += [count]
    unique_set = {oligo:set(accs) for oligo, accs in unique.iteritems()}
    return unique_set

# <codecell>

def filter_oligos(oligos_dict, primer_conc,
                  hairpin_tm_max, homo_tm_max,
                  tm_max, tm_min, min_occurrence, no_3_T, 
                  no_poly_3_GC, no_poly_run, max_degen,
                  mv_conc, dv_conc, rev):
    """
    Removes oligos that don't have desirable properties.
    """
    # Rescale concentration of primers
    primer_conc_resc = primer_conc/float(max_degen)
    # Remove oligos that occur less than min_occurrence
    if type(oligos_dict.values()[0]) == int:
        oligos = [oligo for oligo, value in oligos_dict.iteritems() if \
                                               value >= min_occurrence]
    else:
        oligos = [oligo for oligo, value in oligos_dict.iteritems() if \
                                         len(value) >= min_occurrence]
    # It is necessary to use the reverse of the complement for the reverse primers
    if rev:
        oligos = [str(Seq(oligo).reverse_complement()) for oligo in oligos]
    # Apply empirical rules
    oligos = apply_rules(oligos, rev, no_3_T, no_poly_3_GC)

    # Remove kmers with hairpin tm > threshold
    calc_hairpin = partial(primer3.calcHairpinTm,
                           dna_conc = primer_conc_resc,
                           mv_conc = mv_conc, 
                           dv_conc = dv_conc)
    hairpin = map(calc_hairpin, oligos)
    oligos = [oligo for oligo, hp_tm in zip(oligos, hairpin) if hp_tm <= hairpin_tm_max]

    # Remove kmers with homodimers tm > threshold
    calc_homo_tm = partial(primer3.calcHomodimerTm,    
                           dna_conc = primer_conc_resc,
                           mv_conc = mv_conc, 
                           dv_conc = dv_conc)
    homo_tm = map(calc_homo_tm, oligos)
    oligos = [oligo for oligo, homo_tm in zip(oligos, homo_tm) if homo_tm <= homo_tm_max]

    # Remove kmers with poly runs
    if no_poly_run:
        polys = ["AAAA", "TTTT", "GGGG", "CCCC"]
        find_poly = lambda x: not any([True for poly in polys if poly in x])
        oligos = filter(find_poly, oligos)
    
    # Remove primers with tm above threshold
    calc_tm = partial(primer3.calcTm,                           
                      dna_conc = primer_conc_resc,
                      mv_conc = mv_conc, 
                      dv_conc = dv_conc)
    tms = map(calc_tm, oligos)
    for oligo, tm in zip(oligos, tms):
        if tm < tm_min or tm > tm_max:
            oligos.remove(oligo)
        else:
            pass
    # For compatibility with the dictionaries, I will return the rev primers for their original 
    if rev:
        oligos = [str(Seq(oligo).reverse_complement()) for oligo in oligos]
    return oligos

# <codecell>

def apply_rules(oligos, rev, no_3_T, no_poly_3_GC):
    """
    Apply some empirical rules for primer design.
    """
    # I invert the oligo when it is the reverse so that I can treat the 3 terminal equally
    if no_poly_3_GC:
        find_poly_GC = lambda x: not any([True for poly in ["GGG", "CCC"] if poly in x[-3:]])
        oligos = filter(find_poly_GC, oligos)
    if no_3_T:
        oligos = filter(lambda x: not x[-1] == "T", oligos)
    return oligos

# <codecell>

def discard_redundant(oligo_series, max_degen, rev, 
                      primer_conc, verbose,  min_diff,
                      mv_conc, dv_conc):
    """
    From reduntant oligos (oligos that detects exactly the same sequences) select those who have
    highest melting temperature and keep the remaining in a dictionary for further reuse if the 
    representant is discarded for any reason.
    """
    # Rescale primer concentration
    primer_conc_resc = primer_conc/float(max_degen)
    
    to_remove = []
    if rev:
        oligos = [str(Seq(oligo).reverse_complement()) for oligo in oligo_series.index]
    else:
        oligos = [oligo for oligo in oligo_series.index]
    redundant = {}
    calc_tm = partial(primer3.calcTm,                           
                  dna_conc = primer_conc_resc,
                  mv_conc = mv_conc, 
                  dv_conc = dv_conc)
    tms = map(calc_tm, oligos)
    tms = {oligo:tm for tm, oligo in zip(tms, oligo_series.index)}
    for count, base_primer in enumerate(oligo_series.index):
        if verbose:
            clear_output()
            print "Iteration %d of %d" % (count, len(oligo_series))
            sys.stdout.flush()
        if base_primer in to_remove:
            continue
        for primer in oligo_series.index:
            if primer == base_primer:
                continue
            if primer in to_remove:
                continue
            union =  oligo_series[base_primer].union(oligo_series[primer])
            diff = len(union) - len(oligo_series[base_primer])
            size_diff = len(oligo_series[base_primer]) - len(oligo_series[primer])
            if numpy.abs(diff) <= min_diff and size_diff <= min_diff:
                # If two oligos match the same sequences, keep the one with highest Tm
                if tms[primer] < tms[base_primer]:
                    to_remove += [primer]
                    # Keep a record of the redundant oligos because they may be reused if the
                    # one chose here is discarded later
                    if not base_primer in redundant:
                        redundant[base_primer] = [primer]
                    else:
                        redundant[base_primer] += [primer]
                    # If the primer that was removed contains other redundant primers
                    # add them to the dict as well in the key corresponding to the selected primer
                    if primer in redundant:
                        redundant[base_primer] += redundant[primer]
                else:
                    to_remove += [base_primer]
                    if not primer in redundant:
                        redundant[primer] = [base_primer]
                    else:
                        redundant[primer] += [base_primer]
                    if base_primer in redundant:
                        redundant[primer] += redundant[base_primer]
                    break
    to_keep = [index for index in oligo_series.index if not index in to_remove]
    return {"oligo_series":oligo_series[to_keep], "redundant":redundant}

# <codecell>

def find_valid_positions(seqs_all, max_gap_prop):
    """
    Find positions in the aligment that are not mostly gaps.
    """
    positions = zip(*seqs_all)
    count_gaps = lambda x: len([base for base in x if base == "-"])
    gaps = map(count_gaps, positions)
    max_gaps = len(seqs_all) * max_gap_prop
    pos = range(len(positions))
    valid_pos = filter(lambda x: x[1] <= max_gaps, zip(pos, gaps))
    valid_pos = [p[0] for p in valid_pos]
    return valid_pos

# <codecell>

def enumerate_positions_for_screen(fasta_file, kmer_sizes, step, 
                                   max_gap_prop):
    """
    Calculates the positions where oligos should be enumerated based on the
    step parameter.
    """
    # Select valid columns
    records = read_fasta(fasta_file)
    n_seq = float(len(records))
    seqs_all = [str(record.seq) for record in records]
    valid_pos = find_valid_positions(seqs_all, max_gap_prop)
    
    # Define positions were oligos will be enumerated
    last_pos = valid_pos[ :-min(kmer_sizes)]
    
    pos = [p for p in valid_pos if p < last_pos]
    starts = [pos[i] for i in range(0, len(pos), step)]
    return (starts, seqs_all, n_seq)

# <codecell>

def enumerate_kmers_screen(start, starts, seqs_all, n_seq, kmer_sizes, 
                           hairpin_tm_max, homo_tm_max,
                           tm_max, tm_min, min_occurrence, no_3_T, 
                           no_poly_3_GC, max_degen, no_poly_run,
                           primer_conc, min_diff, mv_conc, dv_conc, look_ahead):
    """
    Apply the functions above to enumerate all possible oligos that match some criteria
    specified in its arguments along a given alignemnt. This is used for screening positions
    in the alignment that might be useful for designing primers.
    """
    # Rescale primer concentration
    primer_conc_resc = primer_conc/float(max_degen)
    # After the middle of the alignment, primers will be tested as reverse
    rev = False
    if start > starts[len(starts) / 2]:
        rev = True
    unique_kmers = enumerate_oligos(starts = [start], 
                                    kmer_sizes = kmer_sizes, 
                                    seqs_all = seqs_all,
                                    look_ahead = look_ahead)
    
    kmers = filter_oligos(oligos_dict = unique_kmers, 
                          rev = rev, 
                          primer_conc = primer_conc,
                          hairpin_tm_max = hairpin_tm_max, 
                          homo_tm_max = homo_tm_max, 
                          tm_max = tm_max, 
                          tm_min = tm_min, 
                          min_occurrence = min_occurrence, 
                          no_3_T = no_3_T, 
                          no_poly_3_GC = no_poly_3_GC, 
                          max_degen = max_degen, 
                          no_poly_run = no_poly_run,
                          mv_conc = mv_conc, 
                          dv_conc = dv_conc)

    unique_kmers = {i:unique_kmers[i] for i in kmers}
    kmers_series = pandas.Series(unique_kmers)
    cover = kmers_series.map(len).sort(inplace = False, ascending = False)
    kmers_series = kmers_series[cover.index]
    not_redundant = discard_redundant(kmers_series, 
                                      max_degen = max_degen, 
                                      rev = rev, 
                                      primer_conc = primer_conc, 
                                      verbose = False,
                                      min_diff = min_diff,
                                      mv_conc = mv_conc,
                                      dv_conc = dv_conc)
    
    kmers_series = not_redundant["oligo_series"]
    unique_kmers = {kmer:unique_kmers[kmer] for kmer in \
                                            kmers_series.index}
    # Selected the n (max_degen) best oligos
    sequences_detected = kmers_series.map(len)
    sequences_detected.sort(ascending = False)
    best_primers = sequences_detected.index[:max_degen]
    detected = set()
    for oligo in best_primers:
        detected = detected.union(kmers_series[oligo])
    coverage = len(detected) / n_seq
    unique_kmers = {key:unique_kmers[key] for key in best_primers}
    
    calc_tm = partial(primer3.calcTm,                           
                      dna_conc = primer_conc_resc,
                      mv_conc = mv_conc, 
                      dv_conc = dv_conc)
    tm = map(calc_tm, unique_kmers.keys())
    if len(tm):
        tms_median = numpy.median(tm)
        tms_10 = numpy.percentile(tm, 10)
        tms_90 = numpy.percentile(tm, 90)
    else:
        tms_median, tms_10, tms_90 = (0, 0, 0)
    # After these filtering, calculate coverage and richness of the most
    # abundant oligos
    richness = len(unique_kmers)
    return ((start, {"Richness":richness, 
                    "Coverage":coverage, 
                    "Median_Tm":tms_median,
                    "Tm_10":tms_10,
                    "Tm_90":tms_90}), 
             unique_kmers, tm)

# <codecell>

def screen_oligos(fasta_file, kmer_sizes, hairpin_tm_max = 35, homo_tm_max = 35,
                  tm_max = 65, tm_min = 50, min_occurrence = 5, no_3_T = True, 
                  no_poly_3_GC = True, max_degen = 60, no_poly_run = True,
                  step = 3, primer_conc = 200, max_gap_prop = .1, 
                  min_diff = 0, mv_conc = 50, dv_conc = 1.5, look_ahead = 50):
    """
    Apply the function enumerate_kmers_screen in parallel for enumerating all possible oligos that match some criteria
    specified in its arguments along a given alignment. This is used for screening positions
    in the alignment that might be useful for designing primers.
    Arguments:
        -fasta_file: a string given the name of the fasta file with aligned sequences to be used for enumeration of oligos;
        -kmer_sizes: a list of integers with the desired size of oligos;
        -hairpin_tm_max: maximal hairpin melting temperature allowed (in Celsius degrees);
        -homo_tm_max: maximal homodimer melting temperature allowed (in Celsius degrees);
        -tm_max: maximal melting temperature allowed for a oligo (in Celsius degrees);
        -tm_min: minimal melting temperature allowed for a oligo (in Celsius degrees);
        -min_occurrence: integer. Minimal allowed occurrence of a oligo along all sequences;
        -no_3_T: boolean. Should oligos with a T in the 3' end be discarded?
        -no_poly_3_GC: boolean. Should oligos with three G's of C's in the 3' end be discarded?
        -max_degen: the maximal number of subprimers desired. This will also be used to rescale the oligo concentration.
        -no_poly_run: boolean. Should oligos with four or more runs of the same bases be discarded?
        -step: distance between positions in the aligment from which primers should be enumerated. A step of 1 implies that
               all positions will be used.
        -primer_conc: total primer concentration in nM. This concentration will be rescaled automatically by the max_degen.
        -max_gap_prop: float. Maximal proportion of gaps allowed for any position where oligos will be enumerated.
        -min_diff: minimal difference of sequences detected for a oligo to be considered redundant. This parameter should be
                   kept at its default value unless you have strong reasons to change it.
        -mv_conc: monovalent ions concentration in mM.
        -dv_conc: divalent ions conentration in mM.
    """
    starts, seqs_all, n_seq = enumerate_positions_for_screen(\
                                    fasta_file,
                                    kmer_sizes = kmer_sizes, 
                                    max_gap_prop = max_gap_prop,
                                    step = step)
    kwargs = {"starts":starts, 
            "seqs_all":seqs_all, 
            "n_seq":n_seq, 
            "kmer_sizes":kmer_sizes, 
            "hairpin_tm_max":hairpin_tm_max, 
            "homo_tm_max":homo_tm_max,
            "tm_max":tm_max, 
            "tm_min":tm_min, 
            "min_occurrence":min_occurrence, 
            "no_3_T":no_3_T, 
            "no_poly_3_GC":no_poly_3_GC, 
            "max_degen":max_degen, 
            "no_poly_run":no_poly_run,
            "primer_conc":primer_conc,
            "min_diff":min_diff,
            "mv_conc":mv_conc,
            "dv_conc":dv_conc,
            "look_ahead":look_ahead}
    p = partial(enumerate_kmers_screen, **kwargs)
    dview.push({"enumerate_kmers_screen":enumerate_kmers_screen,
               "discard_redundant":discard_redundant,
               "filter_oligos":filter_oligos,
               "enumerate_oligos":enumerate_oligos,
               "apply_rules":apply_rules})
    task = lview.map(p, starts, chunksize = 20)
    wait_on(task)
    pos = {res[0][0]:res[0][1] for res in task.result}
    unique_kmers = [(res[0][0], res[1]) for res in task.result]
    data = pandas.DataFrame(pos).T
    data["Pos"] = data.index
    return {"data":data, "unique_kmers":unique_kmers}

# <codecell>

def combine_positions(positions):
    all = []
    for pos in positions:
        all += unique_sets[pos]
    return set(all)

# <headingcell level=3>

# Enumerate oligos using specified positions

# <markdowncell>

# Some of the functions used in this section were defined in the previous section

# <codecell>

def unite_pairs(pair):
    set_pair = fwd_unique[pair[0]].intersection(rev_unique[pair[1]])
    return (pair, set_pair)

# <codecell>

def enumerate_pairs(fwd_unique, rev_unique):
    """
    Combine fwd and rev primers as pairs and return it as a dict where the key is the pair itself
    and the values are the set of sequences detected by the pair.
    """
    pairs_oligos = list(product(fwd_unique.index, rev_unique.index))
    dview.push({"fwd_unique":fwd_unique,
                "rev_unique":rev_unique,
                "unite_pairs":unite_pairs})
    if not len(pairs_oligos):
        print  "No primers were found!"
        return
    task = lview.map(unite_pairs, pairs_oligos, chunksize = 1000)
    wait_on(task)
    pairs = {pair:set for pair, set in task.result}
    return pairs

# <codecell>

def filter_pair(fwd, rev, max_delta, max_tm_ht, dv_conc, mv_conc, primer_conc):
    """
    Tests whether or not a pair o primer is compatible
    """
    rev_comp = str(Seq(rev).reverse_complement())
    fwd_tm = primer3.calcTm(fwd, mv_conc = mv_conc,
                            dv_conc = dv_conc,
                            dna_conc = primer_conc)
    rev_tm = primer3.calcTm(rev_comp, mv_conc = mv_conc,
                            dv_conc = dv_conc,
                            dna_conc = primer_conc)
    is_delta_high = np.abs(fwd_tm - rev_tm) > max_delta
    ht_tm = primer3.calcHeterodimerTm(fwd, rev_comp, 
                                      mv_conc = mv_conc,
                                      dv_conc = dv_conc,
                                      dna_conc = primer_conc)
    is_heterodimer = ht_tm > max_tm_ht
    return is_delta_high or is_heterodimer

# <codecell>

def update_redundancy(redundancy, new_key, rev = False):
    oligo_key = "rev" if rev else "fwd"
    for key, value in redundancy[oligo_key].iteritems():
        if new_key in value:
            redundancy[oligo_key][new_key] = redundancy[oligo_key][key]
    return redundancy

# <codecell>

def test_heterodimer(fwd, rev, included_fwd, included_rev_comp, max_tm_ht, mv_conc, dv_conc, primer_conc):
    if not len(included_rev_comp) and not len(included_fwd):
        return False
    rev_comp = str(Seq(rev).reverse_complement())
    for oligo in included_rev_comp + included_fwd:
        for candidate in [rev_comp, fwd]:
            ht_tm = primer3.calcHeterodimerTm(oligo, 
                                              candidate, 
                                              mv_conc = mv_conc,
                                              dv_conc = dv_conc,
                                              dna_conc = primer_conc)
            if ht_tm >= max_tm_ht:
                return True

# <codecell>

def select_pairs(fwd_dict, rev_dict, redundancy_dict, pairs, max_degen = 100, max_delta = 5, 
                 max_tm_ht = 35, min_increase = 5, primer_conc=200, mv_conc = 50, dv_conc = 1.5):
    """
    Select primers as pairs based on their coverage and compatibility.
    Arguments:
        -fwd_dict: a dict of forward primers as returned by the function enumerate_primers;
        -rev_dict: a dict of reverse primers as returned by the function enumerate_primers;
        -redundancy_dict: a dict of primers grouped by their redundancy as returned by enumerate_primers;
        -pairs: a dict with pairs of primers as keys and the sequences detected by the pair as values. This
                dict is created by the function enumerate_primers;
        -max_degen: maximal number of subprimers to be kept;
        -max_delta: maximal absolute difference in melting temperature between primers in a pair. 
        -max_tm_ht: maximal heterodimer melting temperature;
        -min_increase: minimal number of new sequences detected for a candidate pair to be selected;
    """
    best = set()
    included_fwd = []
    included_rev = []
    included_rev_comp = []
    pairs_series = pandas.Series(pairs)
    cover = pairs_series.map(len).sort(inplace = False, ascending = False)
    pairs_series = pairs_series[cover.index]
    # Rescale concentration of primers
    primer_conc_resc = primer_conc/float(max_degen)
    # Increase the min_increase for the first iterations
    current_min_increase = min_increase + 100
    while True:
        sys.stdin.readline() # Just to display the results in real time.
        increase = 0
        to_del = []
        reject_pair = False
        fwd_degen_reached = len(included_fwd) > max_degen
        rev_degen_reached = len(included_rev) > max_degen
        for fwd, rev in pairs_series.index:
            # If the number of sequences detected by the pair is smaller than min_increase, 
            # delete it and go to next iteration.
            if len(pairs[(fwd, rev)]) < min_increase:
                to_del += [(fwd, rev)]
                continue
            will_degen_exceed_fwd =  (fwd not in included_fwd) and fwd_degen_reached 
            will_degen_exceed_rev =  (rev not in included_rev) and rev_degen_reached
            is_not_compatible = filter_pair(fwd = fwd, rev = rev, max_delta = max_delta, 
                                            max_tm_ht = max_tm_ht, dv_conc = dv_conc, 
                                            mv_conc = mv_conc, primer_conc = primer_conc_resc)
            if is_not_compatible:
                # When a pair is not compatible, try to look for redundant oligos that are compatible
                new_pair = reuse_redundant(fwd = fwd, rev = rev, redundancy_dict = redundancy_dict, 
                                           max_delta = max_delta, max_tm_ht = max_tm_ht,
                                           dv_conc = dv_conc, mv_conc = mv_conc, 
                                           primer_conc = primer_conc_resc)
                if new_pair:
                    new_fwd, new_rev = new_pair
                    pairs[(new_fwd, new_rev)] = pairs[(fwd, rev)]
                    fwd, rev = new_fwd, new_rev
                    is_not_compatible = False
            # If a pair is incompatible and no substitute could be found or 
            # if the degeration was reached, delete the pair and go to next iteration.
            if will_degen_exceed_fwd or will_degen_exceed_rev or is_not_compatible:
                to_del += [(fwd, rev)]
                continue
            # If both primers were already included in previous pairs, add the detected sequences 
            # to the set of detected sequences, delete the pair, and go to next iteration.
            if fwd in included_fwd and rev in included_rev:
                best = best.union(pairs[(fwd, rev)])
                to_del += [(fwd, rev)]
                continue
            is_heterodimer = test_heterodimer(fwd = fwd, rev = rev, 
                                              included_fwd = included_fwd, 
                                              included_rev_comp = included_rev_comp, 
                                              max_tm_ht = max_tm_ht, 
                                              mv_conc = mv_conc, dv_conc = dv_conc, 
                                              primer_conc = primer_conc_resc)
            if is_heterodimer:
                to_del += [(fwd, rev)]
                continue
            # If a pair survived the previous conditions, test it.
            union = pairs[fwd, rev].union(best)
            increase = len(union) - len(best)
            # If the pair increases the coverage, add it to the included oligos list
            # delete the pair and stop this internal loop
            if increase >= current_min_increase:
                best = best.union(pairs[(fwd, rev)])
                if not fwd in included_fwd:
                    included_fwd += [fwd]
                if not rev in included_rev:
                    included_rev += [rev]
                    included_rev_comp = [str(Seq(rev_i).reverse_complement()) for rev_i in included_rev_comp]
                to_del += [(fwd, rev)]
                break
            elif increase < min_increase:
                to_del += [(fwd, rev)]
        # If no pair gave a high enough increase, reduce the current_min_increase by 20
        # until reach the min_increase specified by the user
        if increase < current_min_increase:
            current_min_increase -= 20
            if current_min_increase < min_increase:
                break
        if fwd_degen_reached and rev_degen_reached:
            break
        if to_del:
            to_keep = [idx for idx in pairs_series.index if not idx in to_del]
            pairs = {key:pairs[key] for key in to_keep}
            pairs_series = pairs_series[to_keep]
        clear_output()
        print "Total Coverage: %d" % (len(best))
        print "Forward Degeneration: %d" % (len(included_fwd))
        print "Reverse Degeneration: %d" % (len(included_rev))
        print "Remaing pairs: %d" % (len(pairs))
        sys.stdout.flush()
    return {"fwd":included_fwd,
        "rev":included_rev,
        "covered":best}

# <codecell>

def enumerate_primers(target_file_name,
                      fwd_starts,
                      rev_starts,
                      hairpin_tm_max = 30, 
                      primer_conc = 200, 
                      homo_tm_max = 30,
                      kmer_sizes = [18, 19, 20, 21, 23, 24, 25, 26, 27, 28],
                      tm_min = 55,
                      tm_max = 60,
                      min_occurrence = 10,
                      no_3_T = True,
                      no_poly_3_GC = True,
                      no_poly_run = True,
                      max_degen = 60,
                      mv_conc=50, 
                      dv_conc=1.5,
                      look_ahead = 50):
    """
    Enumerates forward and reverse primers from two regions of an alignment and filters them according to user-defined criteria.
    Arguments:
        -target_file_name: a string given the name of the fasta file with aligned sequences to be used for enumeration of oligos;
        -fwd_starts: a list of integers giving the starting positions for enumerating forward primers;
        -rev_starts: a list of integers giving the starting positions for enumerating reverse primers;
        -hairpin_max_tm: maximal hairpin melting temperature allowed (in Celsius degrees);
        -primer_conc: total primer concentration in nM. This concentration will be rescaled automatically by the max_degen.
        -homo_tm_max: maximal homodimer melting temperature allowed (in Celsius degrees);
        -kmer_sizes: a list of integers with the desired size of oligos;
        -tm_min: minimal melting temperature allowed for a oligo (in Celsius degrees);
        -tm_max: maximal melting temperature allowed for a oligo (in Celsius degrees);
        -min_occurrence: integer. Minimal allowed occurrence of a oligo along all sequences;
        -no_3_T: boolean. Should oligos with a T in the 3' end be discarded?
        -no_poly_3_GC: boolean. Should oligos with three G's of C's in the 3' end be discarded?
        -max_degen: the maximal number of subprimers desired. This will also be used to rescale the oligo concentration.
        -no_poly_run: boolean. Should oligos with four or more runs of the same bases be discarded?
        -step: distance between positions in the aligment from which primers should be enumerated. A step of 1 implies that
               all positions will be used.
        -mv_conc: monovalent ions concentration in mM.
        -dv_conc: divalent ions concentration in mM.
    """
    records = read_fasta(target_file_name)
    seqs_all = [str(record.seq) for record in records]
    # Enumerate oligos
    fwd_unique = enumerate_oligos(starts = fwd_starts,
                                  kmer_sizes = kmer_sizes, 
                                  seqs_all = seqs_all,
                                  look_ahead = look_ahead)
    rev_unique = enumerate_oligos(starts = rev_starts,
                                  kmer_sizes = kmer_sizes, 
                                  seqs_all = seqs_all,
                                  look_ahead = look_ahead)
    # Filter oligos
    fwd_oligos = filter_oligos(fwd_unique, 
                               rev = False, 
                               hairpin_tm_max = hairpin_tm_max, 
                               homo_tm_max = homo_tm_max, 
                               tm_max = tm_max, 
                               tm_min = tm_min,
                               min_occurrence = min_occurrence, 
                               primer_conc = primer_conc,
                               no_3_T = no_3_T, 
                               no_poly_3_GC = no_poly_3_GC,
                               no_poly_run = no_poly_run,
                               max_degen = max_degen,
                               mv_conc = mv_conc, 
                               dv_conc = dv_conc)
    rev_oligos = filter_oligos(rev_unique, 
                               rev = True, 
                               hairpin_tm_max = hairpin_tm_max, 
                               homo_tm_max = homo_tm_max, 
                               tm_max = tm_max, 
                               tm_min = tm_min,
                               min_occurrence = min_occurrence, 
                               primer_conc = primer_conc,
                               no_3_T = no_3_T, 
                               no_poly_3_GC = no_poly_3_GC, 
                               no_poly_run = no_poly_run,
                               max_degen = max_degen,
                               mv_conc = mv_conc, 
                               dv_conc = dv_conc)
    # Remove redundancy
    fwd_unique = pandas.Series({oligo:fwd_unique[oligo] for oligo in fwd_oligos})
    fwd_reduced = discard_redundant(fwd_unique,
                                    max_degen = max_degen, 
                                    rev = False, 
                                    primer_conc = primer_conc, 
                                    verbose = False,
                                    min_diff = 0,
                                    mv_conc = mv_conc,
                                    dv_conc = dv_conc)
    fwd_unique = fwd_reduced["oligo_series"]
    fwd_cover = fwd_unique.map(len).sort(inplace = False, ascending = False)
    rev_unique = pandas.Series({oligo:rev_unique[oligo] for oligo in rev_oligos})
    rev_reduced = discard_redundant(rev_unique,
                                    max_degen = max_degen, 
                                    rev = True, 
                                    primer_conc = primer_conc, 
                                    verbose = False,
                                    min_diff = 0,
                                    mv_conc = mv_conc,
                                    dv_conc = dv_conc)
    redundancy = {"rev":rev_reduced["redundant"], "fwd":fwd_reduced["redundant"]}
    # Calculate coverage
    rev_unique = rev_reduced["oligo_series"]
    rev_cover = rev_unique.map(len).sort(inplace = False, ascending = False)
    # Make all possible combinations of primers as pairs
    pairs = enumerate_pairs(fwd_unique, rev_unique)
    if not pairs:
        return
    all_pairs = set()
    for pair in pairs.values():
        all_pairs = all_pairs.union(pair)
    all_fwd = set()
    for fwd in fwd_unique.values:
        all_fwd = all_fwd.union(fwd)
    all_rev = set()
    for rev in rev_unique.values:
        all_rev = all_rev.union(rev)
    seqs_detected = all_fwd.intersection(all_rev)
    # Report results
    print "Coverage of all fwd: %d " % len(all_fwd)
    print "Coverage of all rev: %d " % len(all_rev)
    print "Joint coverage of fwd and rev: %d" % len(all_fwd.intersection(all_rev))
    print "Max possible coverage: %d" % len(all_pairs)
    print "Number of Foward Oligos: %d" % len(fwd_unique)
    print "Number of Reverse Oligos: %d" % len(rev_unique)
    return (fwd_unique, rev_unique, redundancy, pairs, seqs_detected)

# <codecell>

def increase_degeneracy(included_oligos, redundancy, oligo_series, min_increase):
    """
    After all pairs of primers have been exhausted, this function tries to add new primers independently
    if the degeneracy is below the maximum specified by the user. This is done in the hope that these primers
    will pair with other primers when they have one or more mismatch with the template sequence.
    Arguments:
        -included_oligos: oligos already included as returned by the function select pairs;
        -redundancy: the dict of oligo groups as returned by enumerate_primers;
        -oligo_series: a pandas series of oligos detected by each primer as returned by enumerate_primers;
        -min_increase: minimal number of new sequences detected for a new primer to be detected.
    """
    # Some oligos in the included oligos are not in the key of the dict redundancy. When that is the case
    # I have to look for it in the values of this dict. The purpose is to have a set of all detected sequences
    found = []
    for oligo in included_oligos:
        if oligo in oligo_series.index:
            accs = oligo_series[oligo]
        else:
            for value in redundancy.values():
                if oligo in value:
                    for oligo in value:
                        if oligo in oligo_series.index:
                            accs = oligo_series[key]
                            break
                        else:
                            pass
                    break
        for acc in accs:
            found += [acc]
    found = set(found)
    not_included = list(set(oligo_series.index) - set(included_oligos))
    while not_included and \
          len(included_oligos) < max_degen:
        unions = pandas.Series()
        for oligo in not_included:
            unions[oligo] = len(found.union(oligo_series[oligo]))
        unions = unions - len(found)
        increase = unions.max()
        best_oligo = unions.idxmax()
        if increase > min_increase:
            found = found.union(oligo_series[best_oligo])
            included_oligos += [best_oligo]
            not_included.remove(best_oligo)
        else:
            break
        "Perfect matches are %d." % len(found)
    return included_oligos

# <codecell>

def reuse_redundant(fwd, rev, redundancy_dict, max_delta, max_tm_ht,
                                           dv_conc, mv_conc, primer_conc):
    if fwd in redundancy_dict["fwd"]:
        fwd_list = redundancy_dict["fwd"][fwd]
    else:
        fwd_list = [fwd]
    if rev in redundancy_dict["rev"]:
        rev_list = redundancy_dict["rev"][rev]
    else:
        rev_list = [rev]
    for fwd in fwd_list:
        for rev in rev_list:
            is_compatible = not filter_pair(fwd = fwd, rev = rev, max_delta = max_delta, 
                                            max_tm_ht = max_tm_ht, dv_conc = dv_conc, 
                                            mv_conc = mv_conc, primer_conc = primer_conc)
            if is_compatible:
                return (fwd, rev)

# <codecell>

def test_heterodimers_post(fwds, revs_comp, primer_conc, mv_conc, dv_conc, max_ht_tm):
    """
    Tests if there are hereterodimers among the primers returned by select_pairs.
    Arguments:
        -best: the object returned by select_pairs;
        -max_degen: maximal degeneration allowed. Used to rescaled primer concentration.
        -mv_conc: monovalent ions concentration in mM;
        -dv_conc: divalent ions concentration in mM.
    """
    # Rescale primer concentration
    degeneracy = max(len(fwds), len(revs_comp))
    primer_conc_resc = primer_conc/float(degeneracy)
    pairs = list(product(fwds, revs_comp))
    calc_ht_tm = partial(primer3.calcHeterodimerTm, 
                         dna_conc = primer_conc_resc,
                         mv_conc = mv_conc, 
                         dv_conc = dv_conc)
    hetero_tms = map(lambda x: calc_ht_tm(x[0], x[1]) > max_ht_tm, pairs)
    if any(hetero_tms):
        print "Heterodimers found!" # I will improve this later
        print [pair for count, pair in enumerate(pairs) if hetero_tms[count]]
    else:
        print "No Heterodimers found!"

# <headingcell level=3>

# Test Primers

# <codecell>

def search_versions_oligos(oligo, substitution, database):
    """
    Not used
    """
    exp = "(%s){s<=%d}" % (oligo, substitution)
    def search_oligo(seq, exp = exp):
        match = regex.search(exp, seq[1])
        if match:
            return match.groups()
    dview.push({"search_oligo":search_oligo,
                "exp":exp})
    task = lview.map(search_oligo, database, chunksize = 500)
    wait_on(task)
    found = [s for s in task.result if s]

# <codecell>

def search_oligos(fwds, revs, substitution, database, return_results = False, return_coverage = False, verbose = True):
    """
    Searches oligos in a list of sequences using regular expression.
    Arguments:
        -fwds: a list of strings giving the forward primers to be searched.
        -revs: a list of strings giving the reverse primers to be searched.
        -substitution: an integer giving the maximum number of mismatches allowerd.
        -database: a list of strings giving the sequences to be used as template.
    Fails if the template and the primers are not in the same orientation.
    """
    fwd_exp = "|".join(["(%s){s<=%d}"% (primer, substitution) for primer in fwds])
    rev_exp = "|".join(["(%s){s<=%d}"% (primer, substitution) for primer in revs])
    fwd_exp = fwd_exp.replace("I", "[ACTG]")
    rev_exp = rev_exp.replace("I", "[ACTG]")
    def search_pair(seq, fwd_exp = fwd_exp, rev_exp = rev_exp):
        fwd_match = regex.search(fwd_exp, seq[1])
        rev_match = regex.search(rev_exp, seq[1])
        if fwd_match and rev_match:
            return (fwd_match.groups(), rev_match.groups())
    dview.push({"search_pair":search_pair,
                "fwd_exp":fwd_exp,
                "rev_exp":rev_exp})
    task = lview.map(search_pair, database, chunksize = 500)
    wait_on(task, verbose = verbose)
    found = [s for s in task.result if s]
    coverage = len(found)
    if return_coverage:
        return coverage
    elif return_results:
        return task.result
    else:
        print "Coverage is %d out of %d sequences" % (coverage, len(database))

# <codecell>

def make_mfe_cmd(primers, 
                 database,
                 output,
                 mfe_exec,
                 ppc,
                 min_tm,
                 oligo_conc,
                 mv_conc,
                 dv_conc):
    """
    Make the MFEprimer command to test the primers' coverage and/or specificity.
    """
    cmd = "%s " % mfe_exec +\
          "-i %s " % primers +\
          "--oligo_conc=%f " % oligo_conc +\
          "-d %s " % database +\
          "--mono_conc=%f " % mv_conc +\
          "--diva_conc=%f " % dv_conc +\
          "--tm_start=%f " % min_tm +\
          "--ppc %d " % ppc +\
          "--tab " +\
          "-o %s " % output
    return cmd

# <codecell>

def run_mfe_parallel(database, 
                     fwd_list, 
                     rev_list,
                     output,
                     min_tm = 40,
                     ppc = 30,
                     mfe_exec = "./MFEprimer/MFEprimer.py",
                     oligo_conc = 200,
                     mv_conc = 50,
                     dv_conc = 1.5):
    """
    Runs MFEprimer in parallel to test the specificity and coverage of a set of primers.
    Arguments:
        -database: the name of the fasta file where the unaligned sequences are.
        -fwd_list: a list of strings giving the fwd primers to be tested.
        -rev_list: a list of strings giving the rev primers to be tested.
        -output: the name of the output file including the path.
        -min_tm: minimal melting temperature for a primer to detected a target sequence.
        -ppc: see MFEprimer manual or just keep it as it is.
        -mfe_exec: the path for the MFEprimer python file.
        -oligo_conc: primer concentration. It should be manually rescaled if degenerate primers
                     of multiplex is being used.
        -mv_conc: monovalent ions concentration in mM.
        -dv_conc: divalent ions concentration in mM.
    """
    # There was a conflict with some other product object, so I decided to import it here
    from itertools import product
    if os.path.isdir("./data/.temp_mfe"):
        shutil.rmtree("./data/.temp_mfe")
    os.makedirs("./data/.temp_mfe")
    os.makedirs("./data/.temp_mfe/primers")
    os.makedirs("./data/.temp_mfe/results")
    def write_primer_pair(pair, count_id):
        with open("./data/.temp_mfe/primers/pair_%d" % count_id, "w") as handle:
            handle.write(">pair_%d_%s_fp\n%s\n>pair_%d_%s_rp\n%s\n" % \
                        (count_id, pair[0], pair[0], count_id, pair[1], pair[1]))
    primers_pairs = product(fwd_list, rev_list)
    pair_dict = {}
    for count_id, pair in enumerate(primers_pairs):
        pair_dict["pair_%d"%count_id] = pair # Used to index the oligo_conc dict
        write_primer_pair(pair, count_id)
    primers = os.listdir("./data/.temp_mfe/primers/")
    cmds = []
    for count, primer in enumerate(primers):
        # This is for using a dictionary of oligos concentration
        if type(oligo_conc) == dict:
            pair = pair_dict[primer]
            if pair[0] in oligo_conc:
                curr_oligo_conc = oligo_conc[pair[0]]
            elif pair[1] in oligo_conc:
                curr_oligo_conc = oligo_conc[pair[1]]
            else:
                raise Exception("Oligo concentration invalid! Are the oligos in the correct strand?")
        else:
            try:
                curr_oligo_conc = float(oligo_conc)
            except TypeError:
                raise Exception("Oligo_conc must be either a number or a dictionary")
        cmd = make_mfe_cmd(primers = "./data/.temp_mfe/primers/%s" % primer,
             database = database,
             output = "./data/.temp_mfe/results/result_%d" % count,
             min_tm = min_tm,
             ppc = ppc,
             mfe_exec = mfe_exec,
             oligo_conc = curr_oligo_conc,
             mv_conc = mv_conc,
             dv_conc = dv_conc)
        cmds += [cmd]
    run_cmd = lambda cmd: os.system(cmd)
    dview.push({"run_cmd":run_cmd})
    task = lview.map(run_cmd, cmds)
    wait_on(task)
    cmd = ("cd ./data/.temp_mfe/results/;"
           "cat $(ls) > all_results;"
           "awk 'NR==1{print $0} !/AmpID/ {print $0}' all_results > primers_out.txt;"
           "cd ../../../;"
           "cp ./data/.temp_mfe/results/primers_out.txt %s")%output
    os.system(cmd)
    #shutil.rmtree("./data/.temp_mfe")
    return task

# <codecell>

#***************************** In development

# <codecell>

def run_tntblast(nproc, 
                 fwd_list,
                 rev_list,
                 output_name, 
                 database_name, 
                 min_tm,
                 max_tm, 
                 primer_conc, 
                 mv_conc,
                 dntp_conc = 0.8, 
                 dv_conc = 1.5,
                 plex = False, 
                 clamp = None, 
                 rescale_conc = False,
                 target_strand =  "plus", 
                 lighter_output = True):
    """
    Tests the primers set using thermonuclotide blast.
    Arguments:
        -nproc: an integer giving the number of processors to be used by mpi
        -query_name: a string with the name of the file containing the primers
        -output_name: a string giving the name of the output file
        -database_name: a string giving the name of the fasta file to be used as input
        -min_tm: an integer with the minimal melting temperature
        -max_tm: an integer with the maximal melting temperature
        -primer_conc: a string with the primer concentration formatted as float.

    """
    # To avoid conflict with pylab
    from itertools import product 
    pairs = product(fwd_list, rev_list)
    with open("./data/.tnt_primers", "w") as handle:
        for count, pair in enumerate(pairs):
            handle.write("pair_%d\t%s\t%s\n" % (count, pair[0], pair[1]))
    try:
        mv_corrected = mv_conc + 120 * (dv_conc - dntp_conc)**.5
    except ValueError:
        mv_corrected = mv_conc
    cmd = "mpirun -np %d "%nproc +\
          "tntblast -i %s "% './data/.tnt_primers' +\
          "-s %fe-3 "%mv_corrected +\
          "-o %s "%output_name +\
          "-d %s "%database_name +\
          "-e %d "%min_tm +\
          "-x %d "%max_tm +\
          "-t %fe-9 "%primer_conc +\
          "--target-strand=%s "%target_strand
    if plex:
        cmd += " --plex=F "
    if clamp:
        cmd += "--primer-clamp=%d "%clamp
    if not rescale_conc:
        cmd += "--rescale-ct=F "
    if lighter_output:
        cmd += " -a F -M F " 
    print cmd
    process = os.system(cmd)
    return process

# <codecell>

#********************* End of development section

# <codecell>

def process_data(fwd_data_name, rev_data_name, delta_tm = 5, return_raw = False):
    """
    Combines the results for fwd and rev primer in a unique dataframe and keeps only the best
    match for each sequence detected.
    Arguments:
        -fwd_data_name: the name of the data file with the results of MFEprimer for foward primers.
        -rev_data_name: the name of the data file with the results of MFEprimer for reverse primers.
        -delta_tm: maximal allowed absolute difference between fwd and reverse melting temperature in Celsius degrees.
    """
    data_fwd = pandas.read_csv(fwd_data_name, sep = "\t")
    data_rev = pandas.read_csv(rev_data_name, sep = "\t")
    data_fwd = data_fwd[["FpID", "HitID", "FpTm", "BindingStart"]]
    data_rev = data_rev[["RpID", "HitID", "RpTm", "BindingStop"]]
    data = pandas.merge(data_rev, data_fwd, on="HitID", how="outer")
    data = data.dropna()
    data["DeltaTm"] = np.abs(data.FpTm - data.RpTm)
    data["AmpLen"] = data.BindingStop - data.BindingStart
    data = data.ix[data.DeltaTm <= delta_tm, :]
    data["fwd_primer"] = data.FpID.map(lambda x: x.split("_")[2])
    data["rev_primer"] = data.RpID.map(lambda x: x.split("_")[2])
    data["Lowest_Tm"] = data.apply(lambda row: min(row["FpTm"], row["RpTm"]), axis = 1)
    if return_raw:
        return data
    grouped = data.groupby(["HitID"], as_index = False)
    data = grouped.apply(lambda group: group.ix[group.Lowest_Tm.idxmax()])
    fwds_tm = data.groupby("fwd_primer").FpTm.max()
    revs_tm = data.groupby("rev_primer").RpTm.max()
    def calculate_max_diff(row, fwds_tm, revs_tm):
        fwd = row["fwd_primer"]
        rev = row["rev_primer"]
        fwd_tm_max = fwds_tm[fwd]
        rev_tm_max = revs_tm[rev]
        fwd_diff = np.abs(row["FpTm"] - fwd_tm_max)
        rev_diff = np.abs(row["RpTm"] - rev_tm_max)
        return max(fwd_diff, rev_diff)
    data["Diff"] = data.apply(lambda x:calculate_max_diff(x, fwds_tm, revs_tm), axis = 1)
    return data

# <headingcell level=3>

# Add inosine to primers

# <codecell>

def calc_tm_general(oligos, 
                    complements, 
                    oligo_conc, 
                    mv_conc, 
                    dv_conc, 
                    nn_method = "all97",
                    salt_method = "san04",
                    in_file = "./data/oligos_for_melting",
                    verbose = False):
    if not complements:
        complements = [str(Seq(oligo).complement()) for oligo in oligos]
        complements = [complement.replace("I", "A") for complement in complements]
    assert type(oligo_conc) in [float, int], "Oligo concentration must be numeric."
    assert type(mv_conc) in [float, int], "Ion concentration must be numeric."
    assert type(dv_conc) in [float, int], "Ion concentration must be numeric."
    with open(in_file, "w") as handle:
        for oligo, complement in zip(oligos, complements):
            handle.write("A%sA T%sT\n" % (oligo, complement))
    melting_path = os.path.abspath("./MELTING/executable/melting-batch") 
    file_path = os.path.abspath(in_file)
    cmd = Template(
    "$melting_path "
    "-H dnadna "
    "-nn $nn_method "
    "-ion $salt_method "
    "-P ${oligo_conc}e-9 "
    "-E Na=${mv_conc}e-3:Mg=${dv_conc}e-3 "
    "$file_path")
    cmd = cmd.substitute(oligo = oligo,
                         salt_method = salt_method,
                         nn_method = nn_method,
                         mv_conc = mv_conc,
                         dv_conc = dv_conc,
                         oligo_conc = oligo_conc,
                         melting_path = melting_path,
                         file_path = file_path)
    if verbose:
        print cmd
    os.system(cmd)
    out_file = in_file + ".results.csv"
    # To remove a ^M character that is preventing the file to be read
    cmd = "awk '!/Delta/{print $1, $2, $3, $4, $5}' %s > %s" % \
          (out_file, out_file + ".modified")
    os.system(cmd)
    results = pandas.read_table(out_file + ".modified", sep = "[\s\t]", header = None)
    results.columns = ["Oligo", "Match", "DeltaH",
                       "DeltaS", "Tm"]
    #os.remove(in_file)
    #os.remove(out_file)
    #os.remove(out_file + ".modified")
    #results["Oligo"] = results["Oligo"].map(lambda x: x[1:-1])
    #results["Match"] = results["Match"].map(lambda x: x[1:-1])
    # To make it comparable with the calculations using primer3
    results["Tm"] = results["Tm"] - 1
    return results[["Oligo", "Match", "Tm"]]

# <codecell>

def find_match_mfe(data, HitID, oligo, binding, rev = False):
    if rev:
        start = binding - len(oligo)
        stop = binding
    else:
        start = binding -1
        stop = binding + len(oligo) - 1
    match = data[HitID][int(start):int(stop)]
    if rev:
        match = str(match[::-1])
    else:
        match = str(Seq(match).complement())
    return match

# <codecell>

def find_all_matches_from_mfe(data, mfe_database, fwd_dummy, min_tm):
    records = make_dict_records("./data/mfe/Nitrososphaera_dummy.genomic")
    records = {key:str(value.seq) for key, value in records.iteritems()}
    msg = ("Database does not contain dummy oligos. Please provide the same database"
       "used in run_mfe_parallel function.")
    assert fwd_dummy in records.values()[0], msg
    data_valid = data.ix[(data.FpTm > min_tm) & (data.RpTm > min_tm), :]
    for oligo, binding, match in [("fwd_primer", "BindingStart", "MatchFwd"),
                                  ("rev_primer", "BindingStop", "MatchRev")]:
        unique = set(zip(data_valid["HitID"], data_valid[oligo], data_valid[binding]))
        unique = list(unique)
        rev = oligo == "rev_primer"
        matches = map(lambda x: find_match_mfe(data = records, 
                                               HitID = x[0], 
                                               oligo = x[1], 
                                               binding = x[2],
                                               rev = rev),
                      unique)
        data_unique = pandas.DataFrame(unique, columns = ["HitID", oligo, binding])
        data_unique[match] = matches
        data_valid = pandas.merge(data_valid, data_unique, on = ["HitID", oligo, binding])
    return data_valid

# <codecell>

def correct_tm(data, total_oligo_conc, mv_conc = 50, dv_conc = 1.5, min_tm = 40, verbose = False):
    data_valid = data.ix[(data.FpTm > min_tm) & (data.RpTm > min_tm), :]
    for column_primer, column_match in [("fwd_primer", "MatchFwd"), 
                                        ("rev_primer", "MatchRev")]:
        unique_comb = set(zip(data_valid[column_primer], data_valid[column_match]))
        oligos = [unique[0] for unique in unique_comb]
        complements = [unique[1] for unique in unique_comb]
        oligo_conc_rescaled = total_oligo_conc / float(len(set(oligos)))
        tm = calc_tm_general(oligos = oligos,
                             complements = complements,
                             oligo_conc = oligo_conc_rescaled,
                             mv_conc = mv_conc,
                             dv_conc = dv_conc,
                             verbose = verbose)
        rename_dict = {"Oligo":column_primer, 
                       "Match":column_match, 
                       "Tm":("tm_" + column_primer[:3])}
        tm = tm.rename(columns = rename_dict)
        data_valid = pandas.merge(data_valid, tm, on = [column_primer, column_match])
    return data_valid

# <codecell>

def find_low_agreement_pos(seqs, min_agreement):
    """
    Find positions of most frequent mismatches in a given primer.
    It is used in the function add_inosine.
    """
    pos = zip(*seqs)
    comp = {}
    for count, p in enumerate(pos):
        comp_i = pandas.Series(p).value_counts() / len(p)
        max_agreement = comp_i.max()
        comp[count] = max_agreement
    comp = pandas.Series(comp)
    comp = comp[comp < min_agreement]
    return comp

# <codecell>

def filter_low_agreement_positions(pos, oligo, allowed_bases, 
                                   min_distance, max_inos_add):
    """
    Discard regions of high variability in the primer if they are not suitable for inosine addition.
    It is used in the function add_inosine.
    """
    pos.sort()
    # Discard positions whose base is not allowed to be substituted by inosine
    good_pos = [p for p in pos.index if oligo[p] in allowed_bases]
    included_pos = []
    for count, p in enumerate(good_pos):
        if count == 0:
            included_pos += [p]
            continue
        distances = [np.abs(p - p_inc) for p_inc in included_pos]
        close = [d for d in distances if d < min_distance]
        if close:
            continue
        else:
            included_pos += [p]
        if len(included_pos) >= max_inos_add:
            break
    return included_pos

# <codecell>

def add_inosine(oligos, records, min_agreement = .85, allowed_bases = "ATG", 
                min_distance = 5, max_inos_add = 3):
    """
    Adds inosine to the primers.
    Arguments:
        -oligos: a list of primers where the inosine should be added.
        -records: a list of records as returned by the function read_fasta.
        -min_agreement: float. Highest proportion allowed of the most abundant base in a given position 
                        for making it a candidate to inosine addition.
        -allowed_bases: string giving the bases that can be substituted by inosine in the original primer.
        -min_distance: minimal spacing in number of bases that should be between inosines.
        -max_inos_add: maximal number of inosines per primer.
    Currently there is no way to test for heterodimers, homodimers of hairpins in primers with inosine using
    primer3, so the user should test if the inosine addition is causing this kind of problem using tools
    like bioanalyzer.
    """
    oligos_degen = []
    for oligo in oligos:
        exp = "(%s){s<=3}"%oligo
        seqs = [regex.findall(exp, str(record.seq)) for record in records]
        seqs = [s[0] for s in seqs if s]
        low_agree = find_low_agreement_pos(seqs, min_agreement)
        included_pos = filter_low_agreement_positions(pos = low_agree, 
                                                      oligo = oligo, 
                                                      allowed_bases = allowed_bases, 
                                                      min_distance = min_distance,
                                                      max_inos_add = max_inos_add)
        oligo = list(oligo)
        if included_pos:
            for p in included_pos:
              oligo[p] = "I"  
        oligos_degen += ["".join(oligo)]
    return oligos_degen

# <codecell>

def remove_redundant_after_inosine_addition(data, min_increase):
    """
    Removes oligos that don't increase coverage after inosine addition.
    Arguments:
        -data: a pandas dataframe with sequences detected by the oligo as returned by measure_coverage_oligos
        -min_increase: minimal increase of coverage for a oligo to be kept.
    """
    cover = data.sum().sort(inplace = False, ascending = False)
    to_del = []
    for primer in cover.index:
        if primer == cover.index[0]:
            previously_detected = data[primer]
            continue
        union = data[primer] | previously_detected
        increase = union.sum() - previously_detected.sum()
        if increase < min_increase:
            to_del += [primer]
    return to_del

# <codecell>

def measure_coverage_oligo(oligos, database, substitution, rev = False):
    """
    Measure the coverage of oligos. Used to discard redundancy.
    Arguments:
        -oligos: a list of oligos to search.
        -database: the template sequences to search.
        -substitution: maximal number of substitutions (I highly recommend to keep it as 0).
        -rev: is the oligo a reverse primer?
    """
    covers = {oligo:None for oligo in oligos}
    for oligo in oligos:
        if rev:
            fwds = ["I"]
            revs = [oligo]
        else:
            fwds = [oligo]
            revs = ["I"]
        covers[oligo] = search_oligos(fwds = fwds, 
                                    revs = revs,
                                    substitution = substitution, 
                                    database = database, 
                                    return_results = True)
    covers = pandas.DataFrame(covers)
    covers = covers.notnull().astype(int)
    return covers

# <codecell>

def distribute_conc_by_cover(data, total_conc):
    """
    Divides the total concentration of the primers proportionally to the number of sequences detected by each primer.
    Arguments:
        -data: a pandas dataframe with sequences detected by the oligo as returned by remove_redundant_after_inosine_addition
        -total_conc: total primer concentration in nM
    """
    detected = data.sum(axis = 1) > 0
    data = data.ix[detected, :]
    effec_conc = 200 / float(len(data.index))
    ind_conc = (float(1) / data.sum(axis = 1)) * effec_conc
    conc = data.apply(lambda x: x*ind_conc, axis = 0).sum()
    return dict(conc)

# <codecell>

def make_complement_conc_rev(conc_rev):
    new_dict = {}
    for primer, conc in conc_rev.iteritems():
        comp = str(Seq(primer).reverse_complement())
        new_dict[comp] = conc
    return new_dict

# <codecell>

def prune_inosine(fwd_degen_list,
                  fwd_degen_dict,
                  rev_degen_list,
                  rev_degen_dict,
                  orig_coverage,
                  min_increase,
                  rev = False):
    """
    Removes inosines that don't increase coverage in relation to sequences already detected by other primers.
    Arguments:
        -fwd_degen_list: a list of forward degenerate primers
        -rev_degen_list: a list of reverse degenerate primers
        -fwd_degen_dict: a dictionary of forward degenerate primers
        -rev_degen_dict: a dictionary of reverse degenerate primers
        -orig_coverage: initial coverage.
        -min_increase: minimal decrease in number of sequences detected to keep a inosine.
        -rev: is the oligo a reverse primer?
    """
    if rev:
        oligo_list = rev_degen_list
        oligo_dict = rev_degen_dict
    else:
        oligo_list = fwd_degen_list
        oligo_dict = fwd_degen_dict
    for pos, oligo in enumerate(oligo_list):
        inosines = [p for p, b in enumerate(oligo) if b == "I"]
        orig_oligo = oligo_dict[oligo]
        for inos_pos in inosines:
            oligo_wo_1_inos = list(oligo)
            oligo_wo_1_inos[inos_pos] = orig_oligo[inos_pos]
            oligo_wo_1_inos = "".join(oligo_wo_1_inos)
            oligo_list[pos] = oligo_wo_1_inos
            if rev:
                rev_list_for_test = oligo_list
                fwd_list_for_test = fwd_degen_list
            else:
                rev_list_for_test = rev_degen_list
                fwd_list_for_test = oligo_list
            mod_coverage = search_oligos(fwds = fwd_list_for_test, 
                                         revs = rev_list_for_test, 
                                         substitution = 0, 
                                         database = database,
                                         return_coverage = True, 
                                         verbose = False)
            if orig_coverage - mod_coverage >= min_increase:
                oligo_list[pos] = oligo # Return to previous value
            else:
                oligo = oligo_list[pos] # Update the reference for next iterations
                orig_coverage = mod_coverage
        clear_output()
        msg_oligo = "reverse" if rev else "forward"
        print "Prunning %s primer %d of %d." % (msg_oligo, pos+1, len(oligo_list))
        sys.stdout.flush()
    return oligo_list

# <codecell>

def prune_inosine_fwd_rev(rev_degen_list,
                          fwd_degen_list,
                          fwd_degen_dict,
                          rev_degen_dict,
                          database,
                          min_increase = 20):
    """
    Applies the function prune_inosine to remove inosines that don't increase coverage in relation to 
    sequences already detected by other primers.
    Arguments:
        -fwd_degen_list: a list of forward degenerate primers
        -rev_degen_list: a list of reverse degenerate primers
        -fwd_degen_dict: a dictionary of forward degenerate primers
        -rev_degen_dict: a dictionary of reverse degenerate primers
        -database: a list of sequences to be searched against the oligos.
        -min_increase: minimal decrease in number of sequences detected to keep a inosine.
    """
    orig_coverage = search_oligos(fwds = fwd_degen_list, 
                                  revs = rev_degen_list, 
                                  substitution = 0, database = database,
                                  return_coverage = True,
                                  verbose = False)
    fwd_degen_list = prune_inosine(fwd_degen_list = fwd_degen_list,
                                   fwd_degen_dict = fwd_degen_dict,
                                   rev_degen_list = rev_degen_list,
                                   rev_degen_dict = rev_degen_dict,
                                   orig_coverage = orig_coverage,
                                   min_increase = min_increase,
                                   rev = False)
    orig_coverage = search_oligos(fwds = fwd_degen_list, 
                                  revs = rev_degen_list, 
                                  substitution = 0, database = database,
                                  return_coverage = True,
                                  verbose = False)
    rev_degen_list = prune_inosine(fwd_degen_list = fwd_degen_list,
                                   fwd_degen_dict = fwd_degen_dict,
                                   rev_degen_list = rev_degen_list,
                                   rev_degen_dict = rev_degen_dict,
                                   orig_coverage = orig_coverage,
                                   min_increase = min_increase,
                                   rev = True)
    return {"fwd_degen":fwd_degen_list, 
            "rev_degen":rev_degen_list}
    

# <codecell>

def classify_seqs(HitID, records_class):
    try:
        seq_class = records_class[HitID]
    except KeyError:
        seq_class = "unknown"
    return seq_class

# <headingcell level=3>

# Discard redundant oligos

# <codecell>

def remove_redundant_one(data_raw_filtered, min_increase, rev = False):
    """
    Removes sequences that are redundant based on thermodynamic simulations.
    Arguments:
        -data_raw_filtered: raw data from a MFEprimer run.
        -min_increase: minimal increase in coverage to keep a primer.
        -rev: is the oligo a reverse primer?
    """
    primer_pos = "rev_primer" if rev else "fwd_primer"
    grouped = data_raw_filtered.groupby([primer_pos])
    covered = grouped.apply(lambda x: set(x["HitID"].unique()))
    coverage = covered.map(len).sort(inplace = False)
    for primer in coverage.index:
        original_coverage = len(data_raw_filtered.HitID.unique())
        data_wo_one = data_raw_filtered.ix[data_raw_filtered[primer_pos] != primer, :]
        mod_coverage = len(data_wo_one.HitID.unique())
        if (original_coverage - mod_coverage) <= min_increase:
            data_raw_filtered = data_wo_one
    return data_raw_filtered

# <codecell>

def remove_redundant_all(data_raw, min_increase):
    """
    Removes sequences that are redundant based on thermodynamic simulations.
    Arguments:
        -data_raw: raw data from MFEprimer (see example)
        -min_increase: minimal increase in coverage to keep a primer.
    """
    min_increase = 60
    grouped = data_raw.groupby(["fwd_primer", "rev_primer"])
    grouped_fwd = data_raw.groupby(["fwd_primer"])
    grouped_rev = data_raw.groupby(["rev_primer"])
    mean_fwd = grouped_fwd.FpTm.max().mean()
    mean_rev = grouped_rev.RpTm.max().mean()
    c1 = data_raw.FpTm.map(lambda x: numpy.abs(x - mean_fwd) <= 4) # Is Tm 4C above or below median Tm?
    c2 = data_raw.RpTm.map(lambda x: numpy.abs(x - mean_rev) <= 4)
    c3 = data_raw.DeltaTm <= 5 # Is delta Tm > 5?
    data_raw_filtered = data_raw.ix[c1 & c2 & c3, :]
    data_filtered = remove_redundant_one(data_raw_filtered, min_increase, rev = False)
    data_filtered = remove_redundant_one(data_filtered, min_increase, rev = True)
    return data_filtered

# <codecell>

def substitute_inosine(oligo):
    iupac = {
    "A":"A",
    "C":"C",
    "T":"T",
    "G":"G",
    "R":"AG",
    "Y":"CT",
    "S":"GC",
    "W":"AT",
    "K":"GT",
    "M":"AC",
    "B":"CGT",
    "D":"AGT",
    "H":"ACT",
    "V":"ACG",
    "N":"ACTG",
    "I":"ACTG",
    }
    oligo = list(oligo)
    sub_primers = [iupac[base] for base in oligo]
    sub_primers = list(product(*sub_primers))
    sub_primers = ["".join(sub_primer) for sub_primer in sub_primers]
    return sub_primers

# <codecell>

def make_degenerate(oligos):
    iupac = { 
     'A': 'A',
     'AC': 'M',
     'ACG': 'V',
     'ACT': 'H',
     'ACTG': 'N',
     'AG': 'R',
     'AGT': 'D',
     'AT': 'W',
     'C': 'C',
     'CGT': 'B',
     'CT': 'Y',
     'G': 'G',
     'GC': 'S',
     'GT': 'K',
     'T': 'T'
     }
    for key in iupac.keys():
        permuted_keys = permutations(key)
        for perm_key in permuted_keys:
            iupac["".join(perm_key)] = iupac[key]
    bases = [set(column) for column in zip(*oligos)]
    bases = ["".join(pos) for pos in bases]
    print "|".join(bases)
    bases = [iupac[pos] for pos in bases]
    return "".join(bases)

# <codecell>

def verify_dimers_degen_primers(fwd_primer_list, 
                                rev_primer_list, 
                                oligo_conc,
                                mv_conc = 50, 
                                dv_conc = 1.5,
                                max_dimer_tm = 35):
    # Use highest concentration to calculate Tm (conservative)
    
    oligos = fwd_primer_list + rev_primer_list
    expanded_oligos = []
    for oligo in oligos:
        expanded_oligos += substitute_inosine(oligo)
    comb_oligos = combinations(expanded_oligos, 2)
    found = False
    for pair in comb_oligos:
        tm = primer3.calcHeterodimerTm(pair[0], 
                                       pair[1],
                                       dna_conc = oligo_conc,
                                       dv_conc = dv_conc,
                                       mv_conc = mv_conc)
        if tm > max_dimer_tm:
            print "Dimer found!"
            found = True
            print pair
            print "Melting temperature: %f" % tm
    if not found:
        print "No dimer found."

# <codecell>

def verify_hairpins(oligo_list,
                    oligo_conc,
                    max_hairpin_tm = 35, 
                    dv_conc = 1.5, 
                    mv_conc = 1.5):
    found = False
    for oligo in oligo_list:
        tm = primer3.calcHairpinTm(oligo, 
                                   dna_conc = oligo_conc,
                                   dv_conc = dv_conc,
                                   mv_conc = mv_conc)
        if tm > max_hairpin_tm:
            print "Hairpin found!"
            found = True
            print oligo
            print "Melting temperature: %f" % tm
    if not found:
        print "No hairpin found."

# <headingcell level=3>

# Tools for alignment visualization

# <codecell>

def deduplicate_for_visualization(aligned_fasta_file):
    print "Removing duplicated sequences."
    cluster_sequences(input_file = "./data/ordered/unaligned.fasta",
                      output_file = "./data/ordered/unaligned_not_amb.fasta",
                      similarity = 1,
                      word_size = 10)
    not_amb = make_dict_records("./data/ordered/unaligned_not_amb.fasta")
    aligned_amb = make_dict_records(aligned_fasta_file)
    aligned_not_amb = [aligned_amb[seq_id] for seq_id in not_amb.keys()]
    write_fasta("./data/ordered/aligned_not_amb.fasta", aligned_not_amb)

# <codecell>

def find_start_end(seq1, seq2):
    start_1 = regex.search(r"^[-\.]*",  seq1)
    start_2 = regex.search(r"^[-\.]*",  seq2)
    end_1 = regex.search(r"[-\.]*$",  seq1)
    end_2 = regex.search(r"[-\.]*$",  seq2)
    start = max(start_1.end(), start_2.end())
    end = min(end_1.start(), end_2.start())
    return (start, end)

# <codecell>

def calculate_distance(seq1, seq2):
    start, end = find_start_end(seq1, seq2)
    distance = sum(a != b for a,b in zip(seq1[start:end], seq2[start:end])) / float(end - start)
    return distance

# <codecell>

def write_ordered_fasta(ref_id,
                        out_fasta, similarity_threshold,
                        original_aligned_file,
                        fasta_file = "./data/ordered/aligned_not_amb.fasta"):
    records = make_dict_records(fasta_file)
    # The ref_id could be removed after deduplication.
    if ref_id in records.keys():
        ref_seq = str(records[ref_id].seq)
    else:
        original_records = make_dict_records(original_aligned_file)
        try:
            ref_seq = str(original_records[ref_id].seq)
        except KeyError:
            raise Exception("%s is not in the provided fasta file!" % ref_seq)
    targets = [str(record.seq) for record in records.values()]
    dists = map(lambda target: calculate_distance(ref_seq, target), targets)
    dists = {seq_id:dist for seq_id, dist in zip(records.keys(), dists)}
    dists = pandas.Series(dists)
    dists.sort(inplace = True)
    dists = dists[dists <= (1 - similarity_threshold)]
    records_sorted = []
    for count, seq_id in enumerate(dists.index):
        records_sorted += [records[seq_id]]
        records_sorted[count].id += "_%.2f%%"% ((1 - dists[seq_id]) * 100)
    write_fasta(out_fasta, records_sorted)
    return records_sorted

# <codecell>

def calculate_composition(records_ordered, reference_seq_id):
    print "Calculating composition."
    seqs = [str(record.seq) for record in records_ordered]
    pos = zip(*seqs)
    compositions = []
    for column in pos:
        composition = {"A":0, "C":0, "T":0, "G":0, "Others":0}
        for base in column:
            try:
                composition[base.upper()] += 1
            except KeyError:
                composition["Others"] += 1
        compositions += [composition]
    data = pandas.DataFrame(compositions) / len(seqs)
    data["Seq_ID"] = reference_seq_id
    data = data.reset_index()
    if not os.path.isdir("./data/ordered/compositions"):
        os.makedirs("./data/ordered/compositions")
    data.to_csv("./data/ordered/compositions/%s.csv" % reference_seq_id)

# <codecell>

def concatenate_compositions():
    files = os.listdir("./data/ordered/compositions/")
    cur_file = files.pop()
    data = pandas.read_csv("./data/ordered/compositions/%s" % cur_file)
    while files:
        cur_file = files.pop()
        cur_data = pandas.read_csv("./data/ordered/compositions/%s" % cur_file)
        data = pandas.concat([data, cur_data])
    data.to_csv("./data/ordered/compositions.csv", index = False)

# <codecell>

def main_visualization(reference_seq_ids, 
                       aligned_fasta_file, 
                       out_fasta, 
                       similarity_threshold,
                       composition_name,
                       classes = None,
                       deduplicate = True):
    """
    This function reorders the provided aligned sequences in a fasta file according to their 
    similarity to a reference sequence in the file. It also plot the base composition for the 
    resulting ordered fasta file.
    Arguments:
        -reference_seq_ids: must be a list with one or more references. Note the reference is the 
                            first "word" that follows tha > sign in a fasta file. For example, a
                            record where the sequence description is as follow:
                            >Bradyrhizobium_3 japonicum
                            The reference for this sequence would be "Bradyrhizobium_3".
                            Also note that this argument must be a list, meaning that the reference(s)
                            must be enclosed by square brackets [] as in the example.
        -aligned_fasta_file: a string giving the name (including the path) of the aligned fasta file
                             to be ordered. This function takes only one fasta file at a time.
        -out_fasta: the prefix of the name of the output file. The reference sequence ID will be added
                    to this prefix to make the name of the file. Make sure you include the path.
        -similarity_threshold: them minimum similarity to the reference. Sequences will be discarded 
                               if they are less similar then the threshold. The distance used is the
                               Hamming distance proportional to sequence size. Pairwise deletion is 
                               used to deal with missing data.
        -composition_name: the name of the pdf file with the composition plot, including the path.
        -classes: an optional argument in case you want to give a meaningful title for each composition
                  plot instead of the reference_id. Note that it must be a list and the order of the 
                  elements is correspondent to the order of the elements in reference_seq_ids.
        -deduplicate: a boolean (True or False) argument indicating whether or not the fasta file
                      should be deduplicated.
    """
    if not classes:
        classes = reference_seq_ids
    if type(reference_seq_ids) != list:
        raise Exception("reference_seq_ids must be given as a list!")
    if os.path.isdir("./data/ordered"):
        shutil.rmtree("./data/ordered")
    records_test = make_dict_records(aligned_fasta_file)
    for seq_id in reference_seq_ids:
        assert seq_id in records_test.keys(), "%s not in the provided fasta file!" % seq_id
    os.makedirs("./data/ordered")
    dealign_seq(in_file_name = aligned_fasta_file, 
                out_file_name = "./data/ordered/unaligned.fasta")
    # Create the "./data/ordered/aligned_not_amb.fasta"
    if deduplicate:
        deduplicate_for_visualization(aligned_fasta_file = aligned_fasta_file)
    else:
        shutil.copyfile(aligned_fasta_file, 
                        "./data/ordered/aligned_not_amb.fasta")
    out_aligned_files = []
    for pos, reference_seq_id in enumerate(reference_seq_ids):
        records_ordered = write_ordered_fasta(ref_id = reference_seq_id,
                                              out_fasta = out_fasta + reference_seq_id, 
                                              original_aligned_file = aligned_fasta_file,
                                              similarity_threshold = similarity_threshold)
        calculate_composition(records_ordered, classes[pos])
        out_aligned_files += [out_fasta + reference_seq_id]
    concatenate_compositions()
    if not ".pdf"  in composition_name:
        composition_name = composition_name + ".pdf"
    make_r_script(composition_name)
    process = os.system("Rscript ./data/ordered/plot_composition.R")
    if process:
        print "It was not possible to plot the composition."
    else:
        shutil.rmtree("./data/ordered")
        print "Done!"
        print "The pdf file with compositions was saved as: %s." % composition_name
        print "The ordered fasta file(s) were/was saved as follows:"
        for file_out in out_aligned_files:
            print file_out

# <codecell>

def make_r_script(composition_name):
    script = ("""
    library(plyr)
    library(reshape2)
    library(ggplot2)

    orig_data <- read.csv("./data/ordered/compositions.csv", stringsAsFactors = FALSE)

    data <- melt(orig_data[, -1], 
                 id.vars = c("index", "Seq_ID"), 
                 value.name = "Freq",
                 variable.name = "Base")

    data$Freq <- as.numeric(data$Freq)

    width = length(unique(data$index)) / 8
    max_index = max(data$index)

    data = ddply(data, c("index", "Seq_ID"), transform, Order = order(Freq, decreasing = TRUE))
    data = ddply(data, c("index", "Seq_ID"), function(x) x[x$Order, ])
    data = ddply(data, c("index", "Seq_ID"), transform, FreqSum = cumsum(Freq))
    data$Base <- as.character(data$Base)
    data$Base[data$Base == "Others"] <- "-"
    data$Seq_ID <- as.factor(data$Seq_ID)
    add_space = function(x, n) paste(rep(paste(rep(" ", 100), collapse = ""), n), x, collapse = "")
    n = max(data$index) / 30
    for(level in levels(data$Seq_ID)){
      levels(data$Seq_ID)[levels(data$Seq_ID) == level] <- add_space(level, n)
    }
    
    
    theme_set(theme_minimal(9))
    graph = ggplot(data) +
      aes(x = index, fill = reorder(Base, Freq), y = Freq) +
      facet_wrap(~ Seq_ID, ncol = 1) +
      geom_bar(stat = "identity", 
               colour = "black", 
               size = .2, 
               show_guide = FALSE,
               width = 1) +
      scale_fill_manual(values = c("C" = "red", 
                                   "A" = "blue", 
                                   "T" = "green" , 
                                   "-" = "grey", 
                                   "G" = "orange")) +
      geom_text(aes(y = FreqSum - Freq / 2, 
                    label = Base,
                    size = Freq), 
                show_guide = FALSE) +
      scale_size_continuous(range = c(0.5, 4)) +
      labs(x = "Position", 
           y = "Frequency", 
           fill = "Bases") +
      scale_x_continuous(expand = c(0, 0), 
                         breaks = seq(0, max_index),
                         labels = seq(1, max_index + 1)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            legend.position = "left",
            strip.text.x = element_text(size=8))


    ggsave("%s", width = width, height = 1.5 * length(unique(data$Seq_ID)), limitsize = FALSE)
    """) % composition_name
    with open("./data/ordered/plot_composition.R", "w") as handle:
        handle.write(script)

# <codecell>

def create_test_train_datasets(original_data_algn, original_data_not_algn, prop_test = .3):
    records_algn = make_dict_records(original_data_algn)
    records_not_algn = make_dict_records(original_data_not_algn)
    n_seq = len(records_algn)
    test_size = int(prop_test * n_seq)
    sample_test = random.sample(xrange(n_seq), test_size)
    sample_train = set(xrange(n_seq)) - set(sample_test)
    acc_test = [records_algn.keys()[i] for i in sample_test]
    acc_train = [records_algn.keys()[i] for i in sample_train]
    test_records_algn = [records_algn[i] for i in acc_test]
    test_records_not_algn = [records_not_algn[i] for i in acc_test]
    train_records_algn = [records_algn[i] for i in acc_train]
    train_records_not_algn = [records_not_algn[i] for i in acc_train]
    write_fasta("./data/aligned_train", train_records_algn)
    write_fasta("./data/unaligned_train", train_records_not_algn)
    write_fasta("./data/aligned_test", test_records_algn)
    write_fasta("./data/unaligned_test", test_records_not_algn)

# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


# <codecell>


