#!/usr/bin/env python

print("Initializing...")

import argparse
import re

convert_nuc = { # Dictionary of degenerate nucleotides
    "U":["T","t"], # Uracil # also include U and for anywhere there is a t
    "W":["A","T","a","t"], # Weak
    "S":["C","G","c","g"], # Strong
    "M":["A","C","a","c"], # Amino
    "K":["G","T","g","t"], # Keto
    "R":["A","G","a","g"], # Purine
    "Y":["C","T","c","t"], # Pyrimidine
    "B":["C","G","T","c","g","t"], # Not A
    "D":["A","G","T","a","g","t"], # Not C
    "H":["A","C","T","a","c","t"], # Not G
    "V":["A","C","G","a","c","g"], # Not T
    "N":["A","C","G","T","a","c","g","t"], # Unknown
    "A":["A","a"], # A
    "T":["T","t"], # T
    "C":["C","c"], # C
    "G":["G","g"]  # G
}

def get_args():
    '''Argparse code to allow for in-line command interface.'''

    parser = argparse.ArgumentParser(description=
    """This script takes in a FASTA file and list of motifs and returns a SVG image showing the locations of the motifs in each sequence""")
    parser.add_argument('-f','--file',  action='store', nargs='?', type=str, 
                    required=True, help='Name of fasta file')
    parser.add_argument('-m','--motif',  action='store', nargs='?', type=str, 
                    required=True, help='Name of file containing motifs')

    return parser.parse_args()

def readMotifs(file2:str)->list:
    '''This function takes the motif file and returns a list of the motifs'''

    # Create an empty list
    list_motifs = list()

    # Open the file and iterate through each line
    with open(file2, "r") as fh:
        for line in fh: 
            line = line.rstrip("\n") # remove newline character
            line = line.upper() # make all characters uppercase
            list_motifs.append(line) # add to the list

    return(list_motifs)

def iterate_motif(seq:str, motif:str, l:int, n:int, counter:int):
    '''Recursion function to iterate through the sequence and motif for matches. Returns true if complete match.'''
    if seq[l+counter] in convert_nuc[motif[counter]]: # if the nucleotide is in the degenerate dictionary (e.g. C or c is in Y)
        counter += 1 # counter for nucleotide base 
        
        if counter < n: # if the counter is less than the length of the motif, keep iterating through
            return iterate_motif(seq, motif, l, n, counter)
        else: # counter >= n
            return True

    else: # if the nucleotide is not in the degenerate dictionary aka "not a match"
        return False

def findFeatures(gene:str, genes:dict, seq:str, list_motifs:list) -> dict:
    '''This function the "genes" data structure and adds the length, exon postions, and motif positions for the specified gene'''

    # Return the length of the sequence
    genes[gene]["length"] = len(seq)

    # Return the bp numbers for the exon(s)
    exon_position = list()

    for l in range(len(seq)):
        if seq[l].isupper(): # exons are capitalized
            exon_position.append(l+1) # add the bp to the list
    
    genes[gene]["exon"] = findRanges(exon_position) # return left and right bp for the exon(s)

    # Return the leftmost bp of each motif in the seq
    motif_position = dict()

    for motif in list_motifs: # check each motif
        n = len(motif)
        motif_position[motif] = list()

        for l in range(len(seq)-n+1): # iterate through the sequence base by base to check for motif match
            counter = 0
            if iterate_motif(seq, motif, l, n, counter): # iterate_motif check for subsequent bp matches
                motif_position[motif].append(l+1) # if complete match, add bp to list

    genes[gene]["motif"] = motif_position

    return genes

def findRanges(numbers:list)->list:
    '''Takes in a list of numbers, removes continuous sequences of numbers besides the first and last in the sequence.
    It is used to return the left and right most position of the exon.'''

    ranges = list()
    a = b = numbers[0] # a and b are equal to the 1st number in the list

    for num in numbers[1:]: # for subsequent numbers
        if num == b + 1: # if the subsequent number is equal to the previous number plus one
            b = num # set b as that new number
        else:
            ranges.append(a if a==b else (a,b)) # for discontinuous sequence, add a to list
            a = b = num # set a and b to the new number
    ranges.append(a if a==b else (a,b)) # append a and b to list as a list 

    return ranges

def drawImage(genes:dict, file_name:str) -> None:
    '''This function draws the SVG from the "genes" data structure. Motifs are marked the leftmost position with a colored circle.
    Exons are denoted by a black box and intron by a grey line.'''
    # import libraries 
    import math
    import cairo
    from IPython.display import SVG, display, Image

    a = len(genes) # number of seuqences

    width, height = 1000, 200*(a+1) # all sequences are less than 1000 bases, number of sequences times 200

    # create the coordinates to display your graphic, designate output
    surface = cairo.SVGSurface(file_name+".svg", width, height)

    # create the coordiantes you will be drawing on (like a transparency) - you can create a transformation
    context = cairo.Context(surface)

    # color scheme for the motifs
    color = [
        (255/255, 0/255, 0/255), # red
        (255/255, 165/255, 0/255), # orange    
        (255/255, 255/255, 0/255), # yellow
        (0/255, 128/255, 0/255), # green
        (0/255, 0/255, 255/255), # blue
        (127/255, 0/255, 255/255) # violet               
    ]

    b = 1 # counter for the sequences

    for seq in genes:
        context.set_source_rgb(0, 0, 0)
        context.move_to(5, 15+200*b)
        context.show_text(seq)

        # draw a line for the intron
        context.set_source_rgb(128/255, 128/255, 128/255) # gray line
        context.set_line_width(4) 
        context.move_to(0, 98+200*b) 
        context.line_to(genes[seq]['length'], 98+200*b) # for the length of the exon
        context.stroke() 

        # draw a rectange for each exon
        for exon in genes[seq]['exon']:
            context.set_source_rgb(0, 0, 0) # black box
            context.rectangle(exon[0], 88+200*b, exon[1]-exon[0], 20)
            context.fill()

        c = 0 # motif counter
        d = len(genes[seq]['motif']) # number of motifs ####<-------------- HERE

        for key in genes[seq]['motif']:
            l = len(genes[seq]['motif'][key])

            for i in genes[seq]['motif'][key]:
                context.set_source_rgb(color[c][0], color[c][1], color[c][2])
                context.set_line_width(1)
                context.move_to(i, 98+200*b)
                context.line_to(i, 25+150*c/d+200*b)
                context.stroke()

                context.arc(i, 23+150*c/d+200*b, 3, 0, 2*math.pi)
                context.fill_preserve()
                context.stroke()

            c += 1
        b += 1

    ## set up legend template
    context.set_source_rgb(0, 0, 0)
    context.move_to(5, 15)
    context.show_text("Legend")

    d = 5 # number of motifs

    # Exon
    context.set_source_rgb(0, 0, 0)
    context.rectangle(5, 25, 5, 5)
    context.fill()

    context.set_font_size(8.0)
    context.move_to(15, 30)
    context.show_text("Exon")

    # Intron
    context.set_source_rgb(128/255, 128/255, 128/255)
    context.set_line_width(2) # for the length of the exon
    context.move_to(5, 28+(200-15)*1/(d+2)) 
    context.line_to(10, 28+(200-15)*1/(d+2))
    context.stroke()

    context.set_font_size(8.0)
    context.move_to(15, 30+(200-15)*1/(d+2))
    context.show_text("Intron")

    # Motifs
    c = 0

    for key in genes[seq]['motif']:
        context.set_source_rgb(color[c][0], color[c][1], color[c][2])
        context.move_to(7, 26+(200-15)*(c+2)/(d+2))
        context.arc(7, 26+(200-15)*(c+2)/(d+2), 2, 0, 2*math.pi)
        context.fill_preserve()
        context.stroke()

        context.set_font_size(8.0)
        context.move_to(15, 30+(200-15)*(c+2)/(d+2))
        context.show_text(key)

        c += 1

    surface.finish()

    return None

def main(file1:str, list_motifs:list)->None:
    '''Main function that reads in FASTA file and parses out sequence and header information.
    Calls other functions as needed.'''
    # Pull file name
    file_name = re.search(r"(.*)(\.fasta)", file1)

    # intialize a empty string and dictionary
    seq = "" # to store sequence information
    genes = dict() # to store gene name and feature information

    with open(file1, "r") as fh:
        for line in fh:
            if line[0] == ">": # pull header information
                if seq != "": # if there is a sequence stored from previous loop, call feature finding function
                    genes = findFeatures(gene, genes, seq, list_motifs)
                    seq = "" # reset string to be empty

                # reset with new gene
                gene = re.findall(r">([a-zA-Z]*)", line)[0]
                genes[gene] = dict()
            else: 
                # remove newline character and add line to sequence string
                line = line.rstrip("\n")
                seq += line

    # account for that last case
    if seq != "":
        genes = findFeatures(gene, genes, seq, list_motifs)
        seq = ""

    print("Drawing...")
    drawImage(genes, file_name.group(1)) # draw the features, output to same file name as FASTA input

    print("Complete.")

if __name__ == "__main__":
    args = get_args()
    file1, file2 = args.file, args.motif

    list_motifs = readMotifs(file2)
    main(file1, list_motifs)

