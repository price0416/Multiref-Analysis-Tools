#!/usr/bin/python
import sys, os
import StringIO
from optparse import OptionParser

if __name__ == '__main__':

    usage = "\n\nReplace values for a specified column in a file with corresponding values from a second file when matching data exist in both files.\nUsage: %prog arg"
    parser = OptionParser(usage)

    parser.add_option("-t", "--targetFile", dest="target_file", help="File with values to be replaced.")
    parser.add_option("-m", "--mapFile", dest="map_file", help="File with map of values. Col1 = value to match, Col2 = replacement value.")
    parser.add_option("-c", "--targetColumn", dest="target_col", help="Which column in the target file to scan. (Integer, 0 based).")
    parser.add_option("-o", "--output", dest="output", help="Output file destination.")

    (options, args) = parser.parse_args()
    if len(args) != 0 or not options.target_file or not options.map_file or not options.target_col or not options.output:
        parser.error("Incorrect number of arguments.\n\t\tUse -h to get more information")

    target_string = ""
    map_string = ""
    fasta_string = ""
    map_dict = {}
    target_col = int(options.target_col)
    
    #Open files and read to strings.
    try:
        target = open(options.target_file, "r")
        print "Opened " + options.target_file
        map = open(options.map_file, "r")
        print "Opened " + options.map_file
        
        try:
            target_string = target.read()
            print "\tRead " + options.target_file + "."
            map_string = map.read()
            print "\tRead " + options.map_file + "."
        finally:
            target.close()
            map.close()
            print "Closed all files."            
    except IOError, e:
        print "IO Error. Halting execution."
        print e[0], e[1]


    #Split each file into lists.
    li_target  = target_string.split("\n")[:-1]
    li_map = map_string.split("\n")[:-1]

    #Split rows into individual items, make a mapping dictionary.      
    for i in range(len(li_map)):
        li_map[i] = li_map[i].split("\t")
        map_dict[li_map[i][0]] = li_map[i][1]
    #print map_dict

    #Iterate file, scanning target column and replacing values from our mapping dictionary.
    for i in range(len(li_target)):
        li_target[i] = li_target[i].split("\t")
        print li_target[i][target_col]
        if li_target[i][target_col] in map_dict:
            li_target[i][target_col] = map_dict.get(li_target[i][target_col])

    #Write the new file to the output file destination.
    print "Writing output..."
    try:
        outfile = open(options.output, "w")
        print "Opened " + options.output + " for writing."
        for i in range(len(li_target)):
            for item in li_target[i][:-1]:
                outfile.write(item + "\t")
            outfile.write(li_target[i][-1] + "\n")
    finally:
        outfile.close()
    print "Process complete."