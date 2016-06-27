#!/usr/bin/python
import sys, os
import StringIO
from optparse import OptionParser

if __name__ == '__main__':

    usage = "\n\nMerges two tab delimeted files based on a mapping file.\nUsage: %prog arg"
    parser = OptionParser(usage)

    parser.add_option("-1", "--CMCP6", dest="file1", help="CMCP6_YJ016 input file.")
    parser.add_option("-2", "--YJ016", dest="file2", help="YJ016_CMCP6 input file.")
    parser.add_option("-m", "--mapFile", dest="map_file", help="File with map of values.")
    parser.add_option("-o", "--output", dest="output", help="Output file destination.")

    (options, args) = parser.parse_args()
    if len(args) != 0 or not options.file1 or not options.file2 or not options.map_file or not options.output:
        parser.error("Incorrect number of arguments.\n\t\tUse -h to get more information")

    f1_string = ""
    f2_string = ""
    map_string = ""

    #Open files and read to strings.
    try:
        file1 = open(options.file1, "r")
        print "Opened " + options.file1
        file2 = open(options.file2, "r")
        print "Opened " + options.file2
        map = open(options.map_file, "r")
        print "Opened " + options.map_file
        
        try:
            f1_string = file1.read()
            print "\tRead " + options.file1 + "."
            f2_string = file2.read()
            print "\tRead " + options.file2 + "."
            map_string = map.read()
            print "\tRead " + options.map_file + "."
        finally:
            file1.close()
            file2.close()
            map.close()
            print "Closed all files."            
    except IOError, e:
        print "IO Error. Halting execution."
        print e[0], e[1]

    #Split each file into lists.
    li_f1  = f1_string.split("\n")[:-1]
    li_f2  = f2_string.split("\n")[:-1]
    li_map = map_string.split("\n")[:-1]

    map_dict = {}
    CMCP6_dict = {}
    YJ016_dict = {}

    #Split rows into individual items, make a mapping dictionary.      
    for i in range(len(li_map)):
        li_map[i] = li_map[i].split("\t")
        map_dict[li_map[i][0]] = li_map[i][1]

    for i in range(len(li_f1)):
        li_f1[i] = li_f1[i].split("\t")
        CMCP6_dict[li_f1[i][0]] = li_f1[i]

    for i in range(len(li_f2)):
        li_f2[i] = li_f2[i].split("\t")
        YJ016_dict[li_f2[i][0]] = li_f2[i]

    #Whenever we find mapable positions, write a line to the output file.
    try:
        outfile = open(options.output, "w")
        #outfile.write("CMCP6_YJ016_HS1\tCMCP6_YJ016_HS2\tCMCP6_YJ016_ASW1\tCMCP6_YJ016_ASW2\tYJ016_CMCP6_HS1\tYJ016_CMCP6_HS2\tYJ016_CMCP6_ASW1\tYJ016_CMCP6_ASW2\n")
	outfile.write("IAI_ID\tK12_ID\tIAI1_batch1\tIAI1_batch2\tIAI1_chemostat1\tIAI1_chemostat2\tIAI1_starvation1\tIAI1_starvation2\tK12_IAI1ref_batch1\tK12_IAI1ref_batch2\tK12_IAI1ref_chemostat1\tK12_IAI1ref_chemostat2\tK12_batch1\tK12_batch2\tK12_chemostat1\tK12_chemostat2\tK12_starvation1\tK12_starvation2\tIAI1_K12ref_batch1\tIAI1_K12ref_batch2\tIAI1_K12ref_chemostat1\tIAI1_K12ref_chemostat2\tIAI1_K12ref_starvation1\tIAI1_K12ref_starvation2\n")
        for i in range(len(li_f2)):
            if li_f2[i][0] in map_dict:
                CMCP6_row = CMCP6_dict.get(map_dict.get(li_f2[i][0]))
                YJ016_row = YJ016_dict.get(li_f2[i][0])
                outfile.write(CMCP6_row[0] + "\t" +	YJ016_row[0] + "\t" + CMCP6_row[1] + "\t" + CMCP6_row[2] + "\t" + CMCP6_row[3] + "\t" + CMCP6_row[4] + "\t" + CMCP6_row[5] + "\t" + CMCP6_row[6] + "\t" + CMCP6_row[7] + "\t" + CMCP6_row[8] + "\t" + CMCP6_row[9] + "\t" + CMCP6_row[10] + "\t" + CMCP6_row[11] + "\t" + CMCP6_row[12] + "\t" + YJ016_row[1] + "\t" + YJ016_row[2] + "\t" + YJ016_row[3] + "\t" + YJ016_row[4] + "\t" + YJ016_row[5] + "\t" + YJ016_row[6] + "\t" + YJ016_row[7] + "\t" + YJ016_row[8]  + "\t" + YJ016_row[9] + "\t" + YJ016_row[10] + "\t" + YJ016_row[11] + "\t" + YJ016_row[12] +  "\n")
    finally:
        outfile.close()
    print "Process complete."
