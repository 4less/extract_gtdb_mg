import sys
from typing import Optional
import argparse
import os
import subprocess
from lib.Scan.PfamScan import PfamScan

hmmsearch_cmd = 'hmmsearch'
prodigal_cmd = 'prodigal'


###

def has_tool(cmd: str):
    try:
        env = os.environ.copy()
        proc = subprocess.Popen(['which', cmd], stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, env=env, encoding='utf-8')

        output, error = proc.communicate()
        return len(output) > 0

    except:
        return False


##########################################################################
### Run commands ###
##########################################################################

def run_prodigal(genome_file: str, tmp_folder: str, closed_ends=False):
    closed_ends = '-c' if closed_ends else ''

    genome_basename = os.path.basename(genome_file)
    genome_basename = genome_basename.rsplit('.', 1)[0]

    aa_out = os.path.join(tmp_folder, "{}.faa".format(genome_basename))
    nt_out = os.path.join(tmp_folder, "{}.fna".format(genome_basename))

    args = ['prodigal', '-m', '-p', 'single', '-q', '-g', '11',
            '-a', aa_out, '-d', nt_out, '-i', genome_file]

    print(' '.join(args))

    p = subprocess.Popen(args, stdout=subprocess.PIPE, encoding='utf-8')
    stdout, stderr = p.communicate()

    print(stdout)
    print(stderr)

    if p.returncode != 0:
        print('Non-zero exit code returned when running prodigal: {stdout}')
        return False, None, None

    return True, aa_out, nt_out


def run_pfam_search(genes_aa: str, hmm_dir: str, tmp_folder: str, cpu: int = 1):
    genes_basename = os.path.basename(genes_aa)
    genes_basename = genes_basename.rsplit('.', 1)[0]

    pfam_output = os.path.join(tmp_folder, "{}_pfam.hits".format(genes_basename))

    pfam_scan = PfamScan(cpu=cpu, fasta=genes_aa, dir=hmm_dir)
    pfam_scan.search()
    pfam_scan.write_results(pfam_output, None, None, None, None)

    return True, pfam_output


def run_tigrfam_search(genes_aa: str, hmm_file: str, tmp_folder: str, cpu: int = 1):
    genes_basename = os.path.basename(genes_aa)
    genes_basename = genes_basename.rsplit('.', 1)[0]

    tigr_output = os.path.join(tmp_folder, "{}_tigr".format(genes_basename))
    tigr_output_hits = os.path.join(tmp_folder, "{}_tigr.hits".format(genes_basename))

    args = [hmmsearch_cmd, '-o', tigr_output, '--tblout', tigr_output_hits,
            '--noali', '--notextw', '--cut_nc', '--cpu', str(cpu), hmm_file, genes_aa]

    p = subprocess.Popen(args, stdout=subprocess.PIPE, encoding='utf-8')
    stdout, stderr = p.communicate()

    # print(stdout)
    # print(stderr)

    if p.returncode != 0:
        print('Non-zero exit code returned when running prodigal: {stdout}')
        return False, None, None

    return True, tigr_output_hits


##########################################################################
### Extract commands ###
##########################################################################

class Hit(object):
    """More significant hits have a higher bit-score and lower e-value.

    Running sorted() on a list of Hit objects will yield a list in order of
    least significant to most significant hits.
    """

    def __init__(self, gene_id: str, hmm_id: str, e_val: float, bit_score: float):
        """Store information about this hit."""
        self.gene_id = gene_id
        self.hmm_id = hmm_id
        self.e_val = e_val
        self.bit_score = bit_score

    def __repr__(self) -> str:
        return f'{self.gene_id} {self.hmm_id} ({self.e_val}/{self.bit_score})'

    def __eq__(self, other) -> bool:
        """Returns True if self is equal to other."""
        return isinstance(other, Hit) and other.gene_id == self.gene_id \
               and other.hmm_id == self.hmm_id and other.e_val == self.e_val \
               and other.bit_score == self.bit_score

    def __lt__(self, other) -> bool:
        """Is self less significant than other?"""
        if self.bit_score < other.bit_score:
            return True
        elif self.bit_score == other.bit_score:
            if self.e_val > other.e_val:
                return True
            elif self.e_val == other.e_val:
                if self.hmm_id > other.hmm_id:
                    return True
                elif self.hmm_id == other.hmm_id:
                    return self.gene_id > other.gene_id
        return False

    def __gt__(self, other) -> bool:
        """Is self more significant than other?"""
        raise NotImplemented
        # if self.bit_score > other.bit_score:
        #     return True
        # elif self.bit_score == other.bit_score:
        #     if self.e_val < other.e_val:
        #         return True
        #     elif self.e_val == other.e_val:
        #         if self.hmm_id < other.hmm_id:
        #             return True
        #         elif self.hmm_id == other.hmm_id:
        #             return self.gene_id < other.gene_id
        # return False

    def __hash__(self) -> int:
        return hash(f'{self.gene_id}_{self.hmm_id}_{self.e_val}_{self.bit_score}')

    def hmm_str(self) -> str:
        """Report the e-value and bit-score for the hmm."""
        return f'{self.hmm_id},{self.e_val},{self.bit_score}'


def add_hit(hit_dict, gene_id: str, hmm_id: str, e_val: float, bit_score: float):
    """Store the most significant HMM hit for each gene."""
    if gene_id not in hit_dict:
        hit_dict[gene_id] = dict()

    # Check if this hit already exists.
    new_hit = Hit(gene_id, hmm_id, e_val, bit_score)
    if hmm_id in hit_dict[gene_id]:
        # Store the new hit if it's more significant.
        if hit_dict[gene_id][hmm_id] < new_hit:
            hit_dict[gene_id][hmm_id] = new_hit
    else:
        hit_dict[gene_id][hmm_id] = new_hit


def get_top_hit(hit_dict, gene_id: str) -> Optional[Hit]:
    """Returns the most significant hit for a given gene id or None."""
    if gene_id not in hit_dict or len(hit_dict[gene_id]) == 0:
        return None
    return sorted(hit_dict[gene_id].values(), reverse=True)[0]


def contains_gene_id(hit_dict, gene_id: str) -> bool:
    """Returns True if the gene_id is found in the top hit file."""
    return gene_id in hit_dict


def add_hit_tigrfam(hit_dict, gene_id: str, hmm_id: str, e_val: float, bit_score: float):
    """Only store the most significant HMM hit per gene as per GTDBNCBI."""
    if contains_gene_id(hit_dict, gene_id):
        # The new hit was more significant, remove any other hits.
        if get_top_hit(hit_dict, gene_id) < Hit(gene_id, hmm_id, e_val, bit_score):
            hit_dict[gene_id] = dict()
            add_hit(hit_dict, gene_id, hmm_id, e_val, bit_score)
    else:
        add_hit(hit_dict, gene_id, hmm_id, e_val, bit_score)


def print_tophits(hit_dict):
    header = ['Gene Id', 'Top hits (Family id,e-value,bitscore)']
    print('\t'.join(header))

    for gene_id, hits in sorted(hit_dict.items()):
        out_hits = list()
        for cur_hit in sorted(hits.values(), reverse=True):
            out_hits.append(cur_hit.hmm_str())
        concat_hits = ';'.join(out_hits)
        # fh.write(f'{gene_id}\t{concat_hits}\n')
        print(f'{gene_id}\t{concat_hits}')


def write_tophits(hit_dict, output: str):
    header = ['Gene Id', 'Top hits (Family id,e-value,bitscore)']
    with open(output, 'w') as outfile:
        outfile.write('\t'.join(header) + '\n')

        for gene_id, hits in sorted(hit_dict.items()):
            out_hits = list()
            for cur_hit in sorted(hits.values(), reverse=True):
                out_hits.append(cur_hit.hmm_str())
            concat_hits = ';'.join(out_hits)
            # fh.write(f'{gene_id}\t{concat_hits}\n')
            outfile.write(f'{gene_id}\t{concat_hits}\n')


def tophit_pfam(pfam_file: str):
    hit_dict = dict()
    with open(pfam_file, 'r') as fh_pfam:
        for line in fh_pfam:
            if line[0] == '#' or not line.strip():
                continue

            line_split = line.split()
            gene_id = line_split[0]
            hmm_id = line_split[5]
            evalue = float(line_split[12])
            bitscore = float(line_split[11])
            add_hit(hit_dict, gene_id, hmm_id, evalue, bitscore)

    return hit_dict


def tophit_tigr(tigrfam_file: str):
    hit_dict = dict()
    with open(tigrfam_file, 'r') as fh_tigrfam:
        for line in fh_tigrfam:
            if line[0] == '#':
                continue

            line_split = line.split()
            gene_id = line_split[0]
            hmm_id = line_split[3]
            evalue = float(line_split[4])
            bitscore = float(line_split[5])

            # add_hit(hit_dict ,gene_id, hmm_id, evalue, bitscore)
            add_hit_tigrfam(hit_dict, gene_id, hmm_id, evalue, bitscore)

    return hit_dict


class Record:
    def __init__(self, header: str):
        self.header = header
        self.sequence = ""

    def append(self, sequence: str):
        self.sequence += sequence


def ReadFasta(file_path: str):
    header = ""
    sequence = ""
    with open(file_path, 'r') as file:
        for line in file:
            line = line.rstrip()

            if line.startswith('>'):
                header = line
                if sequence and header:
                    record = Record(header)
                    record.sequence = sequence
                    yield record
                sequence = ""
            else:
                sequence += line
    return None


class Hit(object):
    """More significant hits have a higher bit-score and lower e-value.

    Running sorted() on a list of Hit objects will yield a list in order of
    least significant to most significant hits.
    """

    def __init__(self, gene_id: str, hmm_id: str, e_val: float, bit_score: float):
        """Store information about this hit."""
        self.gene_id = gene_id
        self.hmm_id = hmm_id
        self.e_val = e_val
        self.bit_score = bit_score

    def __repr__(self) -> str:
        return f'{self.gene_id} {self.hmm_id} ({self.e_val}/{self.bit_score})'

    def __eq__(self, other) -> bool:
        """Returns True if self is equal to other."""
        return isinstance(other, Hit) and other.gene_id == self.gene_id \
               and other.hmm_id == self.hmm_id and other.e_val == self.e_val \
               and other.bit_score == self.bit_score

    def __lt__(self, other) -> bool:
        """Is self less significant than other?"""
        if self.bit_score < other.bit_score:
            return True
        elif self.bit_score == other.bit_score:
            if self.e_val > other.e_val:
                return True
            elif self.e_val == other.e_val:
                if self.hmm_id > other.hmm_id:
                    return True
                elif self.hmm_id == other.hmm_id:
                    return self.gene_id > other.gene_id
        return False

    def __gt__(self, other) -> bool:
        """Is self more significant than other?"""
        raise NotImplemented
        # if self.bit_score > other.bit_score:
        #     return True
        # elif self.bit_score == other.bit_score:
        #     if self.e_val < other.e_val:
        #         return True
        #     elif self.e_val == other.e_val:
        #         if self.hmm_id < other.hmm_id:
        #             return True
        #         elif self.hmm_id == other.hmm_id:
        #             return self.gene_id < other.gene_id
        # return False

    def __hash__(self) -> int:
        return hash(f'{self.gene_id}_{self.hmm_id}_{self.e_val}_{self.bit_score}')

    def hmm_str(self) -> str:
        """Report the e-value and bit-score for the hmm."""
        return f'{self.hmm_id},{self.e_val},{self.bit_score}'


def get_header_to_mg(top_hits_pfam_filepath: str, top_hits_tigr_filepath: str):
    hit_dictionary = dict()
    for file_path in (top_hits_pfam_filepath, top_hits_tigr_filepath):
        with open(file_path, 'r') as file:
            for line in file:
                line = line.rstrip()

                tokens = line.split('\t')

                gene_id = tokens[0]

                hmm_id, dummy1, dummy2 = tokens[1].split(',')

                if hmm_id not in hit_dictionary:
                    hit_dictionary[hmm_id] = list()
                hit_dictionary[hmm_id].append(gene_id)

        extract_list = dict()

        for marker_name, gene_id_list in hit_dictionary.items():
            if len(gene_id_list) == 1:
                extract_list[gene_id_list[0]] = marker_name

    return extract_list

def get_header_to_mg_meta(top_hits_pfam_filepath: str, top_hits_tigr_filepath: str):
    hit_dictionary = dict()
    for file_path in (top_hits_pfam_filepath, top_hits_tigr_filepath):
        with open(file_path, 'r') as file:
            for line in file:
                line = line.rstrip()

                tokens = line.split('\t')

                gene_id = tokens[0]

                hmm_id, dummy1, dummy2 = tokens[1].split(',')

                if hmm_id not in hit_dictionary:
                    hit_dictionary[hmm_id] = list()
                hit_dictionary[hmm_id].append(gene_id)

    return hit_dictionary

def write_markers(output_path: str, header2mg, sequence_file: str):
    with open(output_path, 'w') as output:
        for record in ReadFasta(sequence_file):
            gene_name = record.header[1:].split(' # ')[0]
            # print("{} -> {}".format(gene_name, (gene_name in extract_list)))
            if gene_name in header2mg:
                record.header = ">" + header2mg[gene_name]
                output.write(record.header + '\n')
                output.write(record.sequence + '\n')


def make_dir(path: str):
    if os.path.isdir(path):
        print("Directory {} already exists.".format(path))
        return
    try:
        os.makedirs(path)
    except OSError:
        print("Creation of the directory %s failed" % path)
    else:
        print("Successfully created the directory %s" % path)


def rm_dir(path: str):
    try:
        for f in os.listdir(path):
            os.remove(os.path.join(path, f))
        os.rmdir(path)
    except OSError:
        print("Deletion of the directory %s failed" % path)
    else:
        print("Successfully deleted the directory %s" % path)


##########################################################################
### Set up variables ###
##########################################################################

os.path.dirname(os.path.realpath(__file__))
script_file = os.path.realpath(os.path.join(os.getcwd(), sys.argv[0]))
script_dir = os.path.dirname(script_file)

# Relative paths to hmm files
PFAM_FOLDER = os.path.join(script_dir, 'hmm/pfam')
TIGRFAM_FILE = os.path.join(script_dir, 'hmm/tigrfam/tigrfam.hmm')

if not has_tool(hmmsearch_cmd):
    print("Need tool {}".format(hmmsearch_cmd))
    exit()
if not has_tool(prodigal_cmd):
    print("Need tool {}".format(prodigal_cmd))
    exit()

# Minimal script for extracting PFAM MARKERS
parser = argparse.ArgumentParser()

# Group 1
parser.add_argument("-a", "--aa", type=str,
                    help="One or more files with genes (Protein) (Mutually exclusive to (--genomes))")
parser.add_argument("-n", "--nt", type=str,
                    help="One or more files with genes (DNA) (Mutually exclusive to (--genomes) and requires --aa)")

# Group 2
parser.add_argument("-g", "--genome", type=str,
                    help="Genome file(s) (Mutually exclusive to (--aa/--nt))")

parser.add_argument("-o", "--output_dir", type=str, required=True,
                    help="Marker_gene output folder")
parser.add_argument("-p", "--tmp", type=str, required=True,
                    help="Temporary folder")
parser.add_argument("-k", "--keep_tmp", action="store_true",
                    help="Keep temporary files")
parser.add_argument("-x", "--overwrite_tmp", action="store_true",
                    help="Ignore if tmp dir already exists")
parser.add_argument("-t", "--threads", type=int,
                    help="Number of threads to use.")
parser.add_argument("-l", "--list_only", action="store_true",
                    help="Do not output marker genes as sequence. Output a single tab-delimited file with sequence header in column 1 and marker gene name in column 2. The output file is located in the specified output folder with the filename 'marker_genes.tsv'")

args = parser.parse_args()

from_genomes: bool = True

if (args.aa or args.nt) and args.genome:
    print("options (--aa|--nt) and --genomes are mutually exclusive.")
    exit()

if not args.genome and not args.aa:
    print("Either --genome or --aa must be specified")
    exit()

input_len = 0
if args.aa:
    from_genomes = False
    input_aa = args.aa.split(',')
    input_nt = None
    if args.nt:
        input_nt = args.nt.split(',')

        if len(input_aa) != len(input_nt):
            print("If --aa and --nt are specified, both must be of same length")
            exit()
    input_len = len(input_aa)

if args.genome:
    from_genomes = True
    input_genomes = args.genome.split(',')
    input_len = len(input_genomes)

output_dir = args.output_dir

threads = 1

if args.threads: threads = args.threads

list_only = False
list_output = None
marker_genes_file = "marker_genes.tsv"
marker_genes_meta_file = "marker_genes_meta.tsv"
if args.list_only:
    list_only = True
    list_output = open(os.path.join(output_dir, marker_genes_file), 'w')
    list_meta_output = open(os.path.join(output_dir, marker_genes_meta_file), 'w')

##########################################################################
### Run ###
##########################################################################

tmp_dir = args.tmp

if os.path.isdir(tmp_dir) and not args.overwrite_tmp:
    print("Tmp dir {} already exists. Please delete dir first".format(tmp_dir))
    exit()

make_dir(tmp_dir)

for index in range(input_len):
    aa_file = None
    nt_file = None

    if from_genomes:
        input_genome = input_genomes[index]

        print("file: {}".format(input_genome))

        input_basename = os.path.basename(input_genome)
        input_basename = input_basename.rsplit('.', 1)[0]

        success, aa_file, nt_file = run_prodigal(input_genome, tmp_dir)
    else:
        aa_file = input_aa[index]
        nt_file = input_nt[index] if input_nt else None

        print("aa file: {}".format(aa_file))
        print("nt file: {}".format(nt_file))

        input_basename = os.path.basename(aa_file)
        input_basename = input_basename.rsplit('.', 1)[0]

    pfam_success, pfam_hits = run_pfam_search(aa_file, PFAM_FOLDER, tmp_dir, cpu=threads)
    tigr_success, tigr_hits = run_tigrfam_search(aa_file, TIGRFAM_FILE, tmp_dir, cpu=threads)

    # print("{} -> {}".format(pfam_success, pfam_hits))
    # print("{} -> {}".format(tigr_success, tigr_hits))

    dict_tigr = tophit_tigr(tigr_hits)
    dict_pfam = tophit_pfam(pfam_hits)

    tigr_tophit_file = os.path.join(tmp_dir, '{}_tigr_tophit.tsv'.format(input_basename))
    pfam_tophit_file = os.path.join(tmp_dir, '{}_pfam_tophit.tsv'.format(input_basename))

    write_tophits(dict_tigr, tigr_tophit_file)
    write_tophits(dict_pfam, pfam_tophit_file)

    header2mg = get_header_to_mg(pfam_tophit_file, tigr_tophit_file)
    header2mg_meta = get_header_to_mg_meta(pfam_tophit_file, tigr_tophit_file)



    # print(header2mg)

    mg_basename = os.path.join(output_dir, input_basename)

    if list_only:
        for header, mg in header2mg.items():
            list_output.write("{}\t{}\n".format(header, mg))
        for header, mg in header2mg_meta.items():
            for ele in mg:
                list_meta_output.write("{}\t{}\n".format(header, ele))

    else:
        write_markers("{}.faa".format(mg_basename), header2mg, aa_file)
        print("{}.faa".format(mg_basename))
        if nt_file:
            write_markers("{}.fna".format(mg_basename), header2mg, nt_file)
            print("{}.fna".format(mg_basename))

if list_output:
    list_output.close()
    list_meta_output.close()


if not args.keep_tmp:
    rm_dir(tmp_dir)
