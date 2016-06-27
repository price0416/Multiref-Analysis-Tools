#!/usr/bin/python
import sys, os
import StringIO
from optparse import OptionParser

if __name__ == '__main__':

    usage = "\n\nGiven reciprocal BLAST output for two organisms, this script will identify the set of non-ambiguous contigs/genes.\nUsage: %prog arg"
    parser = OptionParser(usage)

    parser.add_option("-c", "--contigs", dest="contig_file", help="Blast output (tab format) for contigs blasted against reference genome features database.")
    parser.add_option("-f", "--features", dest="features_file", help="Blast output (tab format) for reference genome features blasted against contigs database.")
    parser.add_option("-g", "--fasta", dest="fasta_file", help="Fasta file of reference genome features (complete) as downloaded from NCBI.")
    parser.add_option("-o", "--output", dest="output", help="Output file destination.")


    (options, args) = parser.parse_args()
    if len(args) != 0 or not options.contig_file or not options.features_file or not options.fasta_file or not options.output:
        parser.error("incorrect number of arguments.\n\t\tUse -h to get more information")

    contig_string = ""
    feature_string = ""
    fasta_string = ""
    contig_dict = {}
    feature_dict = {}
    
    #Open files and read to strings.
    try:
        contigs = open(options.contig_file, "r")
        print "Opened " + options.contig_file
        features = open(options.features_file, "r")
        print "Opened " + options.features_file
        fasta = open(options.fasta_file, "r")
        print "Opened " + options.fasta_file
        
        try:
            contig_string = contigs.read()
            print "\tRead " + options.contig_file + "."
            feature_string = features.read()
            print "\tRead " + options.features_file + "."
            fasta_string = fasta.read()
            print "\tRead " + options.fasta_file + "."
        finally:
            contigs.close()
            features.close()
            fasta.close()
            print "Closed all files."            
    except IOError, e:
        print "IO Error. Halting execution."
        print e[0], e[1]
        
    #Split each file into lists.
    li_contigs  = contig_string.split("\n")[:-1]
    li_features = feature_string.split("\n")[:-1]
    li_fasta    = fasta_string.split(">")[:-1]

    #Split contigs rows into individual items.
    for i in range(len(li_contigs)):
        li_contigs[i] = li_contigs[i].split("\t")
        
    for i in range(len(li_features)):
        li_features[i] = li_features[i].split("\t") 

    #Column values for BLAST input:
    #0:  query name
    #1:  subject name
    #2:  percent identities
    #3:  aligned length
    #4:  number of mismatched positions
    #5:  number of gap positions
    #6:  query sequence start
    #7:  query sequence end
    #8:  subject sequence start
    #9:  subject sequence end
    #10: e-value
    #11: bit score
    
    contig_one2one = {}
    feature_one2one = {}
    master_dictionary = {}
    min_identity = 98.0
    max_mismatch = 1
    min_length   = 200
	
    #A li_contigs row now looks like this: 
    #['TR|845|c4_g7_i2|', 'NC_004459.3_gene_3035', '99.27', '1509', '11', '0', '703', '2211', '1509', '1', '0.0', ' 2726']
    for i in range(len(li_contigs)):
        #Trim up the ID's so they match in both input files. Will look like this: 845|c4_g7_i2
        li_contigs[i][0] = li_contigs[i][0].upper()
        if li_contigs[i][0][-1] == "|":
            li_contigs[i][0] = li_contigs[i][0][:len(li_contigs[i][0])-1]
			
    #A li_features row now looks like this:
    #['lcl|NC_004459.3_gene_1', 'tr|845|c4_g7_i2', '99.10', '1002', '9', '0', '1', '1002', '3253', '2252', '0.0', ' 1801']
    for i in range(len(li_features)):
        splitPos = li_features[i][0].index("|")
        li_features[i][0] = li_features[i][0][splitPos+1:]
        li_features[i][1] = li_features[i][1].upper()

    foundIt=0
    #Go through li_contigs and start building up a dictionary.  Keys are the contig_id's.  Values will be a list of rows representing possible gene matches.
    for i in range(len(li_contigs)):
        if float(li_contigs[i][2]) < min_identity or int(li_contigs[i][4]) > max_mismatch or int(li_contigs[i][3]) < min_length:
            continue
        if li_contigs[i][0] not in master_dictionary:
            curRow = [[li_contigs[i][1], li_contigs[i][0][:2] + li_contigs[i][0][3:],li_contigs[i][2], li_contigs[i][3], li_contigs[i][6], li_contigs[i][7], li_contigs[i][8], li_contigs[i][9]]]
            master_dictionary[li_contigs[i][0]] = curRow
        else:
            curItems = master_dictionary.get(li_contigs[i][0])
            curRow = [li_contigs[i][1], li_contigs[i][0][:2] + li_contigs[i][0][3:],li_contigs[i][2],li_contigs[i][3], li_contigs[i][6], li_contigs[i][7], li_contigs[i][8], li_contigs[i][9]]
            curItems.append(curRow)
            master_dictionary[li_contigs[i][0]] = curItems



    print str(len(master_dictionary)) + " keys in master dictionary after contig scan."

    #Now go through the li_features and add any matches we missed.
    for i in range(len(li_features)):
        if float(li_features[i][2]) < min_identity or int(li_features[i][4]) > max_mismatch or int(li_features[i][3]) < min_length:
            continue
        if li_features[i][1] not in master_dictionary:
            curRow = [li_features[i][0], li_features[i][1][:2] + li_features[i][1][3:],li_features[i][2],li_features[i][3], li_features[i][6], li_features[i][7], li_features[i][8], li_features[i][9]]
            master_dictionary[li_features[i][1]] = curRow
        else:
            curItems = master_dictionary.get(li_features[i][1])
            exists = 0
            for j in range(len(curItems)):
                if curItems[j][0] == li_features[i][1]:
                    exists = 1
                    break
            if exists == 0:
                curRow = [li_features[i][0], li_features[i][1][:2] + li_features[i][1][3:],li_features[i][2], li_features[i][3], li_features[i][6], li_features[i][7], li_features[i][8], li_features[i][9]]
                curItems.append(curRow)
                master_dictionary[li_features[i][1]] = curItems


    print str(len(master_dictionary)) + " keys in master dictionary after feature scan."
    print master_dictionary
    matches = {}

    failures = 0
    #For matches with multiple possible gene matches, select the best match by identity and length.
    for key, rows in master_dictionary.items():
        #If this is str it means there were no matches meeting the cutoff criteria, so skip this key.
        if type(master_dictionary[key][0]) is str:
            failures = failures + 1
            continue

        print "----------------------------"
        print "Key: " + str(key)
        if len(rows) == 1:
            matches[key] = rows[0]
            print str(key) + " has only a single match:\n\t" + str(rows[0]) + "\n"
            continue
        else:
            print str(key) + " has " + str(len(rows)) + " possible matches:" 
            len_set = []   #The set of lengths by value for all gene matches.
            ident_set = [] #The set of identities by value for all gene matches.
            len_dups = 0   #boolean for if duplicate lengths exist.
            ident_dups = 0 #boolean for if duplicate identity values exist.

            #Build up the set of possible identities and lengths.
            for item in rows:
                print str(item)		
                ident_set.append(float(item[2]))
                len_set.append(float(item[3]))

            #Check the identity and length sets for duplicate entries.
            if len(ident_set) != len(set(ident_set)):
                ident_dups = 1
            if len(len_set) != len(set(len_set)):
                len_dups = 1

            print "\t\tident_set: " + str(ident_set)
            print "\t\tlen_set: " + str(len_set)
            print "\t\tident_dups: " + str(ident_dups) + "\tlen_dups: " + str(len_dups)

            #There are no duplicate identities or duplicate length gene matches.
            if ident_dups == 0:
                best_match = max(ident_set)
                best_match_index = ident_set.index(best_match)
                selected_row = rows[best_match_index]
                print "\t\t\tMatch identified: " + str(selected_row)
                matches[key] = selected_row
                continue

            #There are duplicate identities, but lengths are different for all gene matches.
            elif ident_dups == 1:
                dup_index = []
                max_ident = max(ident_set)
                curIndex = 0
                for pos in ident_set:
                    if pos == max_ident:
                        dup_index.append(curIndex)
                    curIndex = curIndex + 1
                print "\t\tIndexes with duplicates at maximum identity: " + str(dup_index)
                
                #If the highest identity occurs only once, select that one.
                if len(dup_index) == 1:
                    print "\t\t\tMatch identified: " + str(rows[dup_index[0]])
                    matches[key] = rows[dup_index[0]]

                #Otherwise, pick the one with the higher length.
                else:
                    cur_len_set = []
                    for pos in dup_index:
                        cur_len_set.append(len_set[pos])
                    best_match = max(cur_len_set)
                    #Blast results are ordered by length, so the first item in the length set is always the longest.
                    #the list.index() function returns the first value that matches from the list, so index will automatically
                    #identify the match with the highest length here.
                    best_match_index = len_set.index(best_match)
                    selected_row = rows[best_match_index]
                    print "\t\t\tMatch identified: " + str(selected_row)
                    matches[key] = selected_row
                    continue
            
    print "\n\n"
    print str(len(master_dictionary)) + " items in master dictionary."
    print str(len(matches)) + " items in final match set."
    print str(failures) + " items failed to meet matching criteria."
	
    #Sanity check.  The number of failures should be the starting dictionary size minus the ending dictionary size.
    if (len(master_dictionary) - len(matches)) != failures:
        print "Length mismatch error.  Number of failed matches doesnt line up with initial and final dictionary lengths."
        sys.exit()

    print "Matching complete."
    print matches	

    #Convert the match dictionary to a list to simplify writing the outfile.
    uniqueList = []
    for key,value in matches.items():
        if len(value) == 0:
            continue
			
        keyName = key
        contigMatch = value[0]
        contigIdLong = value[1]
        contigLen = value[3]
        contigStart = value[4]
        contigEnd = value[5]
        featureStart = value[6]
        featureEnd = value[7]
        uniqueList.append([keyName,contigMatch,contigIdLong,contigLen,contigStart,contigEnd,featureStart,featureEnd])

    #Now, go through the fasta file and parse out the locus_tags and locations for these genes.
    #Each item in li_fasta currently looks like: lcl|NC_004459.3_gene_1 [locus_tag=VV1_RS00005] [location=complement(1..1002)] ...
    print "Finding locus tags and location information for unique reciprocal matches..."
    for i in range(len(li_fasta)):
        if not "|" in li_fasta[i]:
            continue
        else:
            splitPos = li_fasta[i].index("|")
            li_fasta[i] = li_fasta[i][splitPos+1:]
            splitPos = li_fasta[i].index("[")
            fasta_featureID = li_fasta[i][:splitPos-1]
            fasta_fetureID = fasta_featureID.strip()
            for j in range(len(uniqueList)):
                if uniqueList[j][1] == fasta_featureID:
                    locus_tag = li_fasta[i].index("=")
                    stopPos = li_fasta[i].index("]")
                    uniqueList[j].append(li_fasta[i][locus_tag+1:stopPos])
                    locationStart = li_fasta[i].index("n=")  #This is the end of the 'location=' part.
                    locationEnd = li_fasta[i].rindex("]")
                    location = li_fasta[i][locationStart+2:locationEnd]
                    strand = "+"
                    if li_fasta[i][locationStart+2] == 'c':    #If its complement its '-' strand, else its '+':
                        strand = "-"
                    uniqueList[j].append(location)
                    uniqueList[j].append(strand)
					
    print "Done." 
    print uniqueList
	
    #Turns out there are some case sensitive names in the denovo assemblies.  So reformat those to match whats in the bam files.
    for i in range(len(uniqueList)):
        splitPoint = uniqueList[i][2].index("|")
        uniqueList[i][2] = uniqueList[i][2][:splitPoint] + uniqueList[i][2][splitPoint:].lower()
	
    #Write the output to a file.
    #uniqueList looks like: [contigKey,contigMatch,contigIdLong,contigLen,contigStart,contigEnd,featureStart,featureEnd, fasta_locusID, feature_location]
    print "Writing output..."
    try: 
        outfile = open(options.output, "w")
        outfile.write("ContigID\tFeatureID\tContigIdLong\tLength\tContigStart\tContigEnd\tFeatureStart\tFeatureEnd\tLocusID\tLocation\tStrand\n")
        print "Opened " + options.output + " for writing."
        for i in range(len(uniqueList)):
            for item in uniqueList[i][:-1]:
                outfile.write(item + "\t")
            outfile.write(uniqueList[i][-1] + "\n")
    finally:
        outfile.close()
    print "Process complete."	

	#MODIFY THIS SO THE OUPUT IS THE SAME AS VERSION1, NOW THERE ARE SOME FORMATTING ISSUES:
	
#SHOULD BE: 
#734|c1_g3_i1    NC_004459.3_gene_2836   TR734|c1_g3_i1  314     1       314     364     51      VV1_RS14210
#BUT IS
#TR|579|C0_G7_I1 NC_004459.3_gene_17     TR579|C0_G7_I1  1860    102     1961    1860    1       VV1_RS00085 
