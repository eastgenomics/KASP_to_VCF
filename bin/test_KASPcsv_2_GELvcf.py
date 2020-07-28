"""
Unit tests for GEL_WGS_SNP_genotyping KASPcsv_2_GELvcf.py

Matt Garner 190627
matthew.garner@addenbrookes.nhs.uk

"""

import KASPcsv_2_GELvcf
import datetime

class Test_SNP(object):


    def test_init_empty(self):
        rsid = "rs1886176"

        snp = KASPcsv_2_GELvcf.SNP(rsid)
        assert snp.rsid == rsid
        assert snp.chrom == None
        assert snp.pos == -1
        assert snp.ref == None
        assert snp.alt == None


    def test_init_populated(self):
        rsid = "rs1886176"
        chrom = "13"
        pos = "20141662 "
        ref = "C"
        alt = "A"

        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)
        assert snp.rsid == rsid
        assert snp.chrom == chrom
        assert snp.pos == pos
        assert snp.ref == ref
        assert snp.alt == alt


    def test_set_genotype_het(self):
        rsid = "rs1886176"
        chrom = "13"
        pos = "20141662 "
        ref = "C"
        alt = "A"
        call = "C:A"

        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)

        snp.set_genotype(call)
        assert snp.GT == ["0","1"]


    def test_set_genotype_homalt(self):
        rsid = "rs1886176"
        chrom = "13"
        pos = "20141662 "
        ref = "C"
        alt = "A"
        call = "A:A"

        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)

        snp.set_genotype(call)
        assert snp.GT == ["1","1"]


    def test_set_genotype_homref(self):
        rsid = "rs1886176"
        chrom = "13"
        pos = "20141662 "
        ref = "C"
        alt = "A"
        call = "C:C"

        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)

        snp.set_genotype(call)
        assert snp.GT == ["0","0"]


    def test_set_genotype_nocall(self):
        rsid = "rs1886176"
        chrom = "13"
        pos = "20141662 "
        ref = "C"
        alt = "A"
        call = "?"

        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)

        snp.set_genotype(call)
        assert snp.GT == [".","."]


    def test_set_genotype_het_minus_strand(self):
        rsid = "rs1886176"
        chrom = "13"
        pos = "20141662 "
        ref = "C"
        alt = "A"
        call = "G:T"

        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)

        snp.set_genotype(call)
        assert snp.GT == ["0","1"]


    def test_set_genotype_het_alt_ref_minus_strand(self):
        rsid = "rs1886176"
        chrom = "13"
        pos = "20141662 "
        ref = "C"
        alt = "A"
        call = "T:G"

        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)

        snp.set_genotype(call)
        assert snp.GT == ["0","1"]


    def test_set_vcf_record(self):
        rsid = "rs1886176"
        chrom = "13"
        pos = "20141662"
        ref = "C"
        alt = "A"
        call = "C:A"
        qual= "1"

        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt, qual=qual)
        snp.GT = ["0","1"]
        snp.set_vcf_record()
        
        assert snp.vcf_record == "13\t20141662\trs1886176\tC\tA\t1\tPASS\t.\tGT\t0/1"


    def test_get_sample(self):
        samples = []

        for i in range(1,10):
            sample = KASPcsv_2_GELvcf.Sample("GM99."+str(i))
            samples.append(sample)

        sample_4 = KASPcsv_2_GELvcf.get_sample("GM99.4", samples)
        assert sample_4 == KASPcsv_2_GELvcf.Sample("GM99.4")


    def test_get_snp(self):
        snps = []

        for i in range(1,10):
            snp = KASPcsv_2_GELvcf.SNP("rs%s" % str(i))
            snps.append(snp)

        snp_4 = KASPcsv_2_GELvcf.get_snp("rs4", snps)
        assert snp_4 == KASPcsv_2_GELvcf.SNP("rs4")


class Test_Sample(object):


    def test_init(self):
        sample_id = "GM99.99999"
        sample = KASPcsv_2_GELvcf.Sample(sample_id)
        assert sample.sample_id == sample_id


    def test_add_snp(self):
        sample_id = "GM99.99999"
        sample = KASPcsv_2_GELvcf.Sample(sample_id)

        rsid = "rs1886176"
        chrom = "chr13"
        pos = "20141662 "
        ref = "C"
        alt = "A"
        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)        

        assert len(sample.snps) == 0    # No snps initially
        
        sample.add_snp(snp)

        assert len(sample.snps) == 1    # 1 snp added

        assert sample.snps[0] == snp    # The correct snp is added!


    def test_sort_snps(self):

        # Generate some SNPs, sort, check list is sorted by chrom then pos
        sample_id = "GM99.99999"
        sample = KASPcsv_2_GELvcf.Sample(sample_id)

        #chr2    227032260       rs10203363      C       T
        #chr2    168932506       rs497692        T       C
        #chr3    4362083 rs2819561       A       G
        #chr1    179551371       rs1410592       G       A
        #chr1    67395837        rs2229546       C       A

        rsid = "rs10203363"
        chrom = "chr2"
        pos = "227032260"
        ref = "C"
        alt = "T"
        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)
        sample.add_snp(snp)
        
        rsid = "rs497692"
        chrom = "chr2"
        pos = "168932506"
        ref = "T"
        alt = "C"
        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)
        sample.add_snp(snp)

        rsid = "rs2819561"
        chrom = "chr3"
        pos = "4362083"
        ref = "A"
        alt = "G"
        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)
        sample.add_snp(snp)

        rsid = "rs1410592"
        chrom = "chr1"
        pos = "179551371"
        ref = "G"
        alt = "A"
        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)
        sample.add_snp(snp)

        rsid = "rs2229546"
        chrom = "chr1"
        pos = "67395837"
        ref = "C"
        alt = "A"
        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)
        sample.add_snp(snp)

        sample.sort_snps()

        chrom_order = [snp.chrom for snp in sample.snps]
        pos_order = [snp.pos for snp in sample.snps]

        assert chrom_order == ["chr1",
                               "chr1",
                               "chr2",
                               "chr2",
                               "chr3",
                               ]
        assert pos_order   == ["67395837",
                               "179551371",
                               "168932506",
                               "227032260",
                               "4362083",
                               ]


    def test_generate_write_vcf(self):

        # Generate some SNPs, sort, check list is sorted by chrom then pos
        sample_id = "GM99.99999"
        sample = KASPcsv_2_GELvcf.Sample(sample_id)

        #chr2    227032260       rs10203363      C       T
        #chr2    168932506       rs497692        T       C
        #chr3    4362083 rs2819561       A       G
        #chr1    179551371       rs1410592       G       A
        #chr1    67395837        rs2229546       C       A

        rsid = "rs10203363"
        chrom = "chr2"
        pos = "227032260"
        ref = "C"
        alt = "T"
        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)
        snp.set_genotype("C:T")
        sample.add_snp(snp)
        
        rsid = "rs497692"
        chrom = "chr2"
        pos = "168932506"
        ref = "T"
        alt = "C"
        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)
        snp.set_genotype("T:C")
        sample.add_snp(snp)

        rsid = "rs2819561"
        chrom = "chr3"
        pos = "4362083"
        ref = "A"
        alt = "G"
        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)
        snp.set_genotype("A:A")
        sample.add_snp(snp)

        rsid = "rs1410592"
        chrom = "chr1"
        pos = "179551371"
        ref = "G"
        alt = "A"
        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)
        snp.set_genotype("A:A")
        sample.add_snp(snp)

        rsid = "rs2229546"
        chrom = "chr1"
        pos = "67395837"
        ref = "C"
        alt = "A"
        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt)
        snp.set_genotype("C:A")
        sample.add_snp(snp)

        sample.sort_snps()

        sample.generate_vcf()
        output_vcf = "test_generate_write_vcf.vcf"
        sample.vcf.write(output_vcf)
        
        expected_vcf = \
"""##fileformat=VCFv4.2
##fileDate={datetime}
##source=KASP_to_VCF_v1.0
##reference=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
##contig=<ID=chr1,AC=CM000663.2,gi=568336023,LN=248956422,rl=Chromosome,M5=6aef897c3d6ff0c78aff06ac189178dd,AS=GRCh38>
##contig=<ID=chr2,AC=CM000664.2,gi=568336022,LN=242193529,rl=Chromosome,M5=f98db672eb0993dcfdabafe2a882905c,AS=GRCh38>
##contig=<ID=chr3,AC=CM000665.2,gi=568336021,LN=198295559,rl=Chromosome,M5=76635a41ea913a405ded820447d067b0,AS=GRCh38>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tGM99.99999
chr1\t67395837\trs2229546\tC\tA\t1\tPASS\t.\tGT\t0/1
chr1\t179551371\trs1410592\tG\tA\t1\tPASS\t.\tGT\t1/1
chr2\t168932506\trs497692\tT\tC\t1\tPASS\t.\tGT\t0/1
chr2\t227032260\trs10203363\tC\tT\t1\tPASS\t.\tGT\t0/1
chr3\t4362083\trs2819561\tA\tG\t1\tPASS\t.\tGT\t0/0
""".format(datetime=datetime.datetime.today().strftime('%Y%m%d'))
        
        with open(output_vcf) as out_vcf_fh:
            output_vcf_text = out_vcf_fh.read()
            assert expected_vcf == output_vcf_text


class Test_VCF(object):


    def test_init(self):
        sample_id = "GM99.99999"
        sample = KASPcsv_2_GELvcf.Sample(sample_id)

        rsid = "rs1886176"
        chrom = "chr13"
        pos = "20141662 "
        ref = "C"
        alt = "A"
        snp = KASPcsv_2_GELvcf.SNP(rsid, chrom=chrom, pos=pos, ref=ref, alt=alt) 
        snp.set_genotype("C:A")
        sample.add_snp(snp)

        vcf = KASPcsv_2_GELvcf.Vcf(sample)
        assert vcf.fileformat == "VCFv4.2"
        assert vcf.sample == sample
        assert vcf.source == "KASP_to_VCF_v1.0"
        assert vcf.reference == "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"
        # contigs in separate test
        assert vcf.filter == ['ID=PASS,Description="All filters passed"']
        assert vcf.format == ['ID=GT,Number=1,Type=String,Description="Genotype"']
        assert vcf.column_headers == ["CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","GM99.99999"]
        # header in separate test
