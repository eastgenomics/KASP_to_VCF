"""
To convert SNP genotypes in KASP csv format into vcf format as specified by GEL.
See docs/SNP_Spec_for_Positive_Sample_Identification_v1.0.docx for details of vcf format

Matt Garner 190621
matthew.garner@addenbrookes.nhs.uk

"""

import os
import csv
import pprint
import copy
from datetime import datetime
import sys
import hashlib
import urllib.request

def QC():
    # TO DO
    return True


class Vcf(object):

    def __init__(self, sample):
        self.fileformat =  "VCFv4.2"
        self.sample     = sample
        self.source     = "KASP_GEL_WGS_v0.1"
        self.reference  = "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        self.contigs    = self._get_SNP_contigs()
        self.filter     = ['ID=PASS,Description="All filters passed"']
        self.format     = ['ID=GT,Number=1,Type=String,Description="Genotype"']
        self.column_headers = ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",self.sample.sample_id]
        self.header     = self._build_header()
        self.sample.sort_snps()

    def _get_SNP_contigs(self):
        snp_contig_ids = set([snp.chrom for snp in self.sample.snps])
        snp_contig_ids = sorted(snp_contig_ids, key=lambda x: int(x[3:]))   # Won't work if we have non-integer chroms, ok for GEL 24
        
        ref_contigs = []
        with open("/mnt/storage/projects/GEL_WGS_SNP_genotyping/data/genome/GRCh38_full_analysis_set_plus_decoy_hla.contignames.fa") as fh:
            for line in fh:
                if line.startswith(">"):
                    ref_contigs.append(line.strip()[1:])
        
        ref_contig_ids = {ref_contig.split(" ")[0]:ref_contig for ref_contig in ref_contigs}
        
        vcf_contigs = []
        for ref_contig_id in ref_contig_ids:
            if ref_contig_id in snp_contig_ids:
                vcf_contigs.append(ref_contig_ids[ref_contig_id])
        
        return sorted(vcf_contigs, key=lambda x: int(x.split("  ")[0][3:]))

    def _build_header(self):
        header_lines = []

        header_kvps={"fileformat":self.fileformat,
                     "fileDate": self._get_datetime(),
                     "source": self.source,
                     "reference": self.reference,
                     }
        
        header_order = ["fileformat",
                        "fileDate",
                        "source",
                        "reference",
                        ]
        
        for header_field in header_order:
            header_lines.append("##%s=%s\n" % (header_field, header_kvps[header_field]))

        # Contigs are pulled from reference file
        for contig in self.contigs:
            formatted_contig = self._format_contig(contig)
            header_lines.append("##contig=<%s>\n" % formatted_contig )            

        header_lines.append("#"+"\t".join(self.column_headers)+"\n")
      
        self.header_lines = header_lines

    def _format_contig(self, contig):
        contig_string = ""
        
        fields = contig.strip().split("  ")
        contig_id = fields[0]
        
        kvps = {"ID": contig_id}
        
        key_order = ["ID"]

        for field in fields[1:]:
            key, value = field.split(":")
            kvps[key] = value
            key_order.append(key)

        first=True
        for key in key_order:
            str_to_append = "%s=%s" % (key, kvps[key])
            if not first:
                str_to_append = ","+str_to_append
            first=False
            contig_string+=str_to_append

        return contig_string

    def write(self, output_filepath):

        temp_out_file = output_filepath+".tmp"

        with open(temp_out_file, "w") as out_vcf_fh:
            for line in self.header_lines:
                out_vcf_fh.write(line)

            for snp in self.sample.snps:
                snp.set_vcf_record()
                out_vcf_fh.write(snp.vcf_record + "\n")

        # To safeguard against partial vcfs
        os.rename(temp_out_file, output_filepath)
        
    def _get_datetime(self):
        return datetime.today().strftime('%Y%m%d')

   
class SNP(object):

    def __init__(self, rsid, chrom=None, pos=-1, ref=None, alt=None, filter="PASS", qual="1", info="."):
        self.rsid   = rsid
        self.chrom  = chrom
        self.pos    = pos
        self.ref    = ref
        self.alt    = alt
        self.calls  = None
        self.format = "GT"
        self.GT     = []
        self.vcf_record = None
        self.filter = filter
        self.info = info
        self.qual = qual

    def __hash__(self):
        return hash((self.rsid))

    def __eq__(self, other):
        if not isinstance(other, type(self)): return NotImplemented
        return self.rsid == other.rsid

    def set_genotype(self, call):

        mapping = {"A":"T",
                   "T":"A",
                   "G":"C",
                   "C":"G",
                   "?":"?"}

        self.alleles = call.split(":")
        self.alleles_strand = 1

        # Figure out which strand the calls are on
        for allele in self.alleles:
            if allele not in [self.ref, self.alt]:
                self.alleles_strand = 0
                break

        # Convert the calls to the plus strand
        if self.alleles_strand == 0:
            self.alleles = [mapping[allele] for allele in self.alleles]
            alleles_strand = 1
       
        # Sort alleles into ref, alt
        sort_order = {self.ref:0,self.alt:1,"?":2}
        self.alleles = sorted(self.alleles, key=lambda x: sort_order[x])

        # Now record ref as 0, alt as 1, no call as .
        for allele in self.alleles:
            if allele == self.ref:
                self.GT.append("0")
            elif allele == self.alt:
                self.GT.append("1")
            elif allele == "?":
                self.GT.extend([".","."])

    def set_vcf_record(self):
        vcf_record = "\t".join([self.chrom, self.pos, self.rsid, self.ref, self.alt, self.qual, self.filter, self.info, self.format, "/".join(self.GT)])
        self.vcf_record = vcf_record


class Sample(object):
    """
    A sample has a list of SNPs, and a method to generate a vcf containing these SNPs

    Sample(sample_id)
    """


    def __init__(self, sample_id):
        self.sample_id = sample_id
        assert sample_id.startswith("GM"), "sample_ids should start with 'GM'"
        self.snps = []

    def add_snp(self, snp):
        """
        Add a SNP to the list of SNPs assocaited with this sample
        """
        self.snps.append(snp)

    def __hash__(self):
        return hash((self.sample_id))

    def __eq__(self, other):
        if not isinstance(other, type(self)): return NotImplemented
        return self.sample_id == other.sample_id

    def sort_snps(self):
        """
        Sort list of SNPs associated with this Sample by chrom then pos.

        May not handle non-integer chroms (MT, X, Y) correctly but doesn't need to right now
        """
        self.snps = sorted(self.snps, key=lambda x: (int(x.chrom[3:]), int(x.pos)))

    def generate_vcf(self):
        self.vcf = Vcf(self)


def get_data_dict(csv_reader):
    section = None                          # Header is always the first section
    sections = ["Statistics","SNPs","Data"] # These sections follow the header

    data = {s:[] for s in sections} # We will extract the data into this dict

    section_column_headers = [] # Each section has different column headers. The headers for the current section are stored here
        
    for row in csv_reader:
        
        # Skip empty rows
        if not any(row):
            continue

        # Track which section of data we are in
        if row[0] in sections:
            section = row[0]
            section_column_headers = []
            #print("Getting '%s' section..." % section)
            continue

        # Skip header
        if not section:
            continue

        # Capture the section column headers
        if not section_column_headers:
            section_column_headers = row
            continue

        # Capture data rows from each section
        row_data = {}
        for index, column in enumerate(section_column_headers):
            if not column:  # Skip empty columns
                continue
            row_data[column] = row[index]

        data[section].append(row_data)

    return data

def get_sample(sample_id, samples):
    """
    Find the sample object in a list of samples with a matching sample_id
    """
    return next((sample for sample in samples if sample.sample_id == sample_id), None)

def get_snp(rsid, snps):
    """
    Find the snp object in a list of snps with a matching rsid
    """
    return next((snp for snp in snps if snp.rsid == rsid), None)


def check_reference_genome():
    project_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    reference_genome_fa_filepath = os.path.join(project_dir, "data/genome/GRCh38_full_analysis_set_plus_decoy_hla.fa")
    
    # Download the file if it doesn't exist already. Too big for github.
    if not os.path.exists(reference_genome_fa_filepath):
        print("Error: Reference genome missing")
        print("\n1. Download reference genome from: ")
        print("http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa")
        print("\n2. Save to:")
        print(reference_genome_fa_filepath)
        print("\n3. Try again\n")

    # Check reference genome integrity
    print("Calc md5...")
    with open(reference_genome_fa_filepath, "rb") as f:
        md5 = hashlib.md5()
        while True:
            data = f.read(8192)
            if not data:
                break
            md5.update(data)
    
        md5sum = md5.hexdigest()
        assert md5sum == "64b32de2fc934679c16e83a2bc072064", "Reference genome md5 incorrect. Aborting"


def main(csv_filepath):

    check_reference_genome()
    exit()
    
    project_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    target_snp_filepath = os.path.join(project_dir, "data/SNPs/GRCh38_coords_ref_alt_min_chr")
    samples = []
    snps = []

    # These SNPs act as a template for the genotyped SNPs associated with each Sample
    with open(target_snp_filepath) as tgt_snps_fh:
        for line in tgt_snps_fh:
            fields = line.strip().split("\t")
            chrom, pos, rsid, ref, alt = fields
            if not SNP(rsid) in snps:
                snps.append(SNP(rsid, chrom, pos, ref, alt))


    # Now we process the KASP data
    with open(csv_filepath) as csv_fh:
        csv_reader = csv.reader(csv_fh)
        data = get_data_dict(csv_reader)

        for entry in data["Data"]:

            sample_id = entry["SubjectID"]

            # NTC does not need a vcf submitted to GEL
            if sample_id == "NTC":
                continue

            # Make a new Sample if this is a new sample_id
            if not Sample(sample_id) in samples:
                samples.append(Sample(sample_id))
            
            # Get the sample object for this sample
            sample = get_sample(sample_id, samples)

            # Get the SNP template object, make a copy to associate with this sample
            rsid = entry["SNPID"]
            call = entry["Call"]
            snp = copy.deepcopy(get_snp(rsid, snps))
            
            # Add genotype info to the SNP
            snp.set_genotype(call)

            # Add the SNP to the sample
            sample.add_snp(snp)

    print("Generating vcfs...")
    csv_dir = os.path.dirname(os.path.abspath(csv_filepath))
    for sample in samples:
        output_filepath = os.path.join(csv_dir, sample.sample_id + ".vcf")
        if os.path.exists(output_filepath):
            print("{output_filepath} already exists, skipping".format(output_filepath=output_filepath))
            continue
        
        sample.generate_vcf()
        sample.vcf.write(output_filepath)

    done_flag_file = os.path.join(csv_dir, "{runfolder}.done".format(runfolder=os.path.basename(csv_dir)))
    open(done_flag_file, 'a').close()


if __name__ == "__main__":
    main(sys.argv[1])