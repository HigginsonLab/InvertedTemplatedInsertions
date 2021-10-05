from cyvcf2 import VCF, Writer
from pyfaidx import Fasta
from argparse import ArgumentParser
from Bio.Seq import Seq

#Main function to extract iTINS and uTINS from insertion in vcf file
def main():
    #Set up variables and initialize TINS counts
    args = get_args()
    genome = Fasta(args.genome)
    vcf = VCF(args.vcf)
    direct = 0
    rev_com = 0
    tins = 0
    #For every entry in vcf file
    for v in vcf:
        #For every individual insertion
        for ins in v.ALT:
            #If greater than 5 bp, check for iTINS and uTINS
            if len(ins) >= 5:
                #Get 50 bp on either side of insertion site
                ref_nucl = get_reference_nucleotide(v.CHROM, v.start, genome)
                ins_seq = Seq(ins)
                ins_rev = Seq(ins).reverse_complement()
                #If reverse complement template found --> iTINS
                if str(ins_rev) in ref_nucl:
                    rev_com += 1
                    tins += 1
                    #Return length of iTINS
                    print("iTINS Insertion Length:", len(ins))
                    index = 0
                    position = []
                    #Return position of iTINS, closest to the insertion site
                    while index < len(ref_nucl):
                        index = ref_nucl.find(str(ins_rev), index)
                        if index == -1:
                            break
                        position.append(index)
                        index += 1
                    closest = min(position, key=lambda x:abs(x-50))
                    print("iTINS Position:", closest)
                #If direct template found --> uTINS
                elif str(ins_seq) in ref_nucl:
                    direct += 1
                    tins += 1
                    #Return length of uTINS
                    print("uTINS Insertion Length:", len(ins))
                    index = 0
                    position = []
                    #Return position of uTINS, closest to the insertion site
                    while index < len(ref_nucl):
                        index = ref_nucl.find(str(ins_seq), index)
                        if index == -1:
                            break
                        position.append(index)
                        index += 1
                    closest = min(position, key=lambda x:abs(x-50))
                    print("uTINS Position:", closest)
    vcf.close()
    print(direct)
    print(rev_com)
    print(tins)

#Find the reference sequence, 50 bp on either side of the insertion
def get_reference_nucleotide(chrom, position, fai):
    return fai[chrom][position - 50:position + 50].seq

#Get the reference genome and vcf file to check
def get_args():
    parser = ArgumentParser(description="")
    parser.add_argument("-g", "--genome", help="genome fasta to match", required=True)
    parser.add_argument("-v", "--vcf", help="input vcf file", required=True)
    return parser.parse_args()


if __name__ == '__main__':
    main()