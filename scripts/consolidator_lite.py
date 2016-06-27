import StringIO
import sys
import csv
import os
from optparse import OptionParser
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

if __name__ == '__main__':
    usage = "\nGiven a DE file with merged ortholog mapping, 2 fasta files each for coding features (protein and nucleotide) used to make the mapping file,\nthis program will select the set of DE genes, BLAST their nucleotide/protein sequences for orthologous genes,\n and produce data tables representing blast results and identify uniquely DE genes between the two input organisms. \nUsage: %prog arg"
    parser = OptionParser(usage)

    parser.add_option("-1", "--de1", dest="de1_file", help="File representing differentially expressed genes for first comparison condition. Output from DESeq2.")
    parser.add_option("-x", "--fasta1", dest="fasta1_file", help="File containing protein coding features for organism 1 in FASTA format.")
    parser.add_option("-y", "--fasta2", dest="fasta2_file", help="File containing protein coding features for organism 2 in FASTA format.")
    parser.add_option("-a", "--nucl1", dest="nucl1_file", help="File containing nucleotide coding features for organism 1 in FASTA format.")
    parser.add_option("-b", "--nucl2", dest="nucl2_file", help="File containing nucleotide coding features for organism 2 in FASTA format.")
    parser.add_option("-o", "--output", dest="output", help="Prefix for output files.")
	
    (options, args) = parser.parse_args()
	
    if len(args) != 0  or not options.de1_file or not options.fasta1_file or not options.fasta2_file or not options.nucl1_file or not options.nucl2_file or not options.output:
        parser.error("Incorrect number of arguments.\n\t\tUse -h to get more information")

    prot1_features = []      #Biopython parse of protein fasta1_file.
    prot2_features = []      #Biopython parse of protein fasta2_file.
    nucl1_features = []      #Biopython parse of nucleotide fasta1_file.
    nucl2_features = []      #Biopython parse of nucleotide fasta2_file.
    genes = []               #Expression gene list for condition 1.
    seqDict = {}             #Dictionary of sequence pairs for condition 1.

    #Output file locations.
    output_dir = options.output + "_tmp"
    prot_outfile = output_dir + "/" + options.output + "_protData.dat"
    nuc_outfile = output_dir + "/" + options.output + "_nucData.dat"

    #Prepare fasta files.
    try:
        print "Reading FASTA files..."
        prot1_features = list(SeqIO.parse(options.fasta1_file, "fasta"))
        prot2_features = list(SeqIO.parse(options.fasta2_file, "fasta"))
        nucl1_features = list(SeqIO.parse(options.nucl1_file, "fasta"))
        nucl2_features = list(SeqIO.parse(options.nucl2_file, "fasta"))
    except Exception,e:
        print "Error reading FASTA files. Terminating execution."
        print str(e)
        sys.exit()
	
    #Parse out geneID's from descriptions and save them as part of the biopython object for nucleotide files.
    print "\tParsing nucleotide geneIDs for strain 1..."
    for i in range(len(nucl1_features)):
        try:
            nucl1_features[i].geneID = nucl1_features[i].description[nucl1_features[i].description.index("gene=")+5:nucl1_features[i].description.index(']')]    
        #If there is no 'gene=...' in the description, mark it and continue.
        except ValueError:
            nucl1_features[i].geneID = -1

    print "\tParsing nucleotide geneIDs for strain 2..."
    for i in range(len(nucl2_features)):
        try:
            nucl2_features[i].geneID = nucl2_features[i].description[nucl2_features[i].description.index("gene=")+5:nucl2_features[i].description.index(']')]    
        #If there is no 'gene=...' in the description, mark it and continue.
        except ValueError:
            nucl2_features[i].geneID = -1

    #Parse out geneID's from descriptions and save them as part of the biopython object for protein files.
    print "\tParsing protein geneIDs for strain 1..."
    for i in range(len(prot1_features)):
        try:
            #Parse and save the GeneID.
            prot1_features[i].geneID = prot1_features[i].description[prot1_features[i].description.index("gene=")+5:prot1_features[i].description.index(']')]
        #If there is no 'gene=...' in the description, mark it and continue.
        except ValueError:
            prot1_features[i].geneID = -1

    print "\tParsing protein geneIDs for strain 2..."
    for i in range(len(prot2_features)):
        try:
            #Parse and save the GeneID.
            prot2_features[i].geneID = prot2_features[i].description[prot2_features[i].description.index("gene=")+5:prot2_features[i].description.index(']')]
         #If there is no 'gene=...' in the description, mark it and continue.
        except ValueError:
            prot2_features[i].geneID = -1

    #Prepare differential expression data file.
    try:
        print "Reading DE file..."
        rows = []
        with open(options.de1_file, 'r') as f:
            reader = csv.reader(f)
            next(reader) #Skip the headers line.
            rows = list(reader)
        for i in range(len(rows)):
            if rows[i][9] == "NA":
                continue
            if float(rows[i][9]) < 0.05:
                genes.append([rows[i][1], rows[i][2], 1])
            if float(rows[i][9]) >= 0.05:
                genes.append([rows[i][1], rows[i][2], 0])
        print "\t" + str(len(genes)) + " genes identified."

    except Exception,e:
        print "Error processing differential expression files. Terminating execution."
        print genes
        print str(e)
        sys.exit()

    print "Identifying sequences for gene pairs."
    for i in range(len(genes)):
        nuc1_seq = -1
        nuc2_seq = -1
        prot1_seq = -1
        prot2_seq = -1

        #Finding genes[i][0] in nucl1 file.
        for j in range(len(nucl1_features)):
            if nucl1_features[j].geneID.upper() == genes[i][0].upper():
                nuc1_seq = nucl1_features[j].seq
                break
        #Finding genes[i][1] in nucl2 file.
        for j in range(len(nucl2_features)):
            if nucl2_features[j].geneID.upper() == genes[i][1].upper():
                nuc2_seq = nucl2_features[j].seq
                break
        #Finding genes[i][0] in prot1 file.
        for j in range(len(prot1_features)):
            if prot1_features[j].geneID.upper() == genes[i][0].upper():
                prot1_seq = prot1_features[j].seq
                break
        #Finding genes[i][1] in prot2 file.
        for j in range(len(prot2_features)):
            if prot2_features[j].geneID.upper() == genes[i][1].upper():
                prot2_seq = prot2_features[j].seq
                break

        if nuc1_seq == -1 or nuc2_seq == -1 or prot1_seq == -1 or prot2_seq == -1:
            print "Not all sequences found for " + str(genes[i])
            print "\tNuc1_Seq: " + str(nuc1_seq)
            print "\tNuc2_Seq: " + str(nuc2_seq)
            print "\tProt1_Seq: " + str(prot1_seq)
            print "\tProt2_Seq: " + str(prot2_seq)
            sys.exit()
        else:
            seqDict[genes[i][0]] = [SeqRecord(nuc1_seq,id=genes[i][0]),SeqRecord(nuc2_seq,id=genes[i][1]),SeqRecord(prot1_seq,id=genes[i][0]), SeqRecord(prot2_seq,id=genes[i][1]), genes[i][2]]

    #Create a directory for our output files if it doesn't already exist.
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    #Blast our sequences.
    try:
        protfile = open(prot_outfile,"w")
        protfile.write("Gene\tIdentity\tScore\tLength\tMismatches\tDE\n")
        nucfile = open(nuc_outfile,"w")
        nucfile.write("Gene\tIdentity\tScore\tLength\tMismatches\tDE\n")

        for gene, seqSet in seqDict.iteritems():
            isDE = seqSet[4]

            #Nucleotide
            nuc1 = seqSet[0]
            nuc2 = seqSet[1]
            nuc1_file = output_dir + "/" + gene + "_nuc1.fasta"
            nuc2_file = output_dir + "/" + gene + "_nuc2.fasta"
            SeqIO.write(nuc1, nuc1_file, "fasta")
            SeqIO.write(nuc2, nuc2_file, "fasta")

            #Protein
            prot1 = seqSet[2]
            prot2 = seqSet[3]
            prot1_file = output_dir + "/" + gene + "_seq1.fasta"
            prot2_file = output_dir + "/" + gene + "_seq2.fasta"
            SeqIO.write(prot1, prot1_file, "fasta")
            SeqIO.write(prot2, prot2_file, "fasta")

            #Run blast...
            blast_nuc = NcbiblastnCommandline(query=nuc1_file, subject=nuc2_file, outfmt=5)()[0]
            blast_nuc_record = NCBIXML.read(StringIO.StringIO(blast_nuc))
            blast_prot = NcbiblastpCommandline(query=prot1_file, subject=prot2_file, outfmt=5)()[0]
            blast_prot_record = NCBIXML.read(StringIO.StringIO(blast_prot))

            for alignment in blast_nuc_record.alignments:
                for hsp in alignment.hsps:
                    if isDE == 1:
                        nucfile.write(gene + "\t" + str(hsp.identities) + "\t" + str(hsp.score) + "\t" + str(hsp.align_length) + "\t" + str(hsp.align_length - hsp.identities)  + "\t" + "1\n")
                    else:
                        nucfile.write(gene + "\t" + str(hsp.identities) + "\t" + str(hsp.score) + "\t" + str(hsp.align_length) + "\t" + str(hsp.align_length - hsp.identities) + "\t" + "0\n")
 
            for alignment in blast_prot_record.alignments:
                for hsp in alignment.hsps:
                    if isDE == 1:
                        protfile.write(gene + "\t" + str(hsp.identities) + "\t" + str(hsp.score) + "\t" + str(hsp.align_length) + "\t" + str(hsp.align_length - hsp.identities)  + "\t" + "1\n")
                    else:
                        protfile.write(gene + "\t" + str(hsp.identities) + "\t" + str(hsp.score) + "\t" + str(hsp.align_length) + "\t" + str(hsp.align_length - hsp.identities) + "\t" + "0\n")
    except Exception,e:
        print "Error processing differential expression BLAST output. Terminating execution."
        print str(e)
        sys.exit()
    finally:
        protfile.close()
        nucfile.close()

    print "Done."