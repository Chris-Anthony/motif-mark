#!/usr/bin/env python

print("Initializing...")

import argparse
import re
import math
import cairo
from IPython.display import SVG, display, Image

# Global constants
# LEFT_MARGIN = # left margin
# HEIGHT_GENE_GROUP = # height of gene groups
# Y_OFFSET_GENE = # y-offset

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

convert_nuc = { # Dictionary of degenerate nucleotides <- change to argparse
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

color = [
        (255/255, 0/255, 0/255), # red
        (255/255, 165/255, 0/255), # orange    
        (255/255, 255/255, 0/255), # yellow
        (0/255, 128/255, 0/255), # green
        (0/255, 0/255, 255/255), # blue
        (127/255, 0/255, 255/255) # violet               
]

class Gene:
    def __init__(self, seq, the_rank):
        '''Gene object contains the FASTA sequence and its length'''
    
        self.start = 0
        self.length = len(seq)
        self.rank = the_rank

    ## methods
    def draw(self, context):
        context.set_source_rgb(0, 0, 0)
        context.set_line_width(4) 
        context.move_to(0, 98+200*self.rank) 
        context.line_to(self.length, 98+200*self.rank)
        context.stroke() 

class Exon:
    def __init__(self, seq, the_rank):
        '''Exon object contains leftmost starting bp relative to FASTA sequence and length of gene'''
        exon_position = list()

        for l in range(len(seq)):
            if seq[l].isupper(): # exons are capitalized
                exon_position.append(l+1) # add the bp to the list
    
        # change list of exon positions to a range
        exon_ranges = list()
        a = b = exon_position[0] # a and b are equal to the 1st number in the list

        for num in exon_position[1:]: # for subsequent exon_position
            if num == b + 1: # if the subsequent number is equal to the previous number plus one
                b = num # set b as that new number
            else:
                exon_ranges.append(a if a==b else (a,b)) # for discontinuous sequence, add a to list
                a = b = num # set a and b to the new number
        exon_ranges.append(a if a==b else (a,b)) # append a and b to list as a list 

        self.start = exon_ranges
        self.rank = the_rank

    ## methods
    def draw(self, context):
        for exon in self.start:
            context.set_source_rgb(0, 0, 0) # black box
            context.rectangle(exon[0], 88+200*self.rank, exon[1]-exon[0], 20)
            context.fill()

class Motif:
    def __init__(self, seq, list_motifs, the_rank):
        '''Motif object contains the leftmost bp positions of motifs relative to FASTA sequence and their length'''

        motif_position = dict()

        for motif in list_motifs: # check each motif
            n = len(motif)
            motif_position[motif] = list()

            for l in range(len(seq)-n+1): # iterate through the sequence base by base to check for motif match
                counter = 0
                if iterate_motif(seq, motif, l, n, counter): # iterate_motif check for subsequent bp matches
                    motif_position[motif].append(l+1) # if complete match, add bp to list

        self.start = motif_position
        self.rank = the_rank
    
    ## methods
    def draw(self, context):
        c = 0 # motif counter
        d = len(self.start)

        for key in self.start:
            l = len(self.start[key])

            for i in self.start[key]:
                context.set_source_rgb(color[c][0], color[c][1], color[c][2])
                context.rectangle(i, 88+200*self.rank-12.5*c/d, len(key), 25)
                context.fill()

                # context.set_line_width(1)
                # context.move_to(i, 98+200*self.rank)
                # context.line_to(i, 25+150*c/d+200*self.rank)
                # context.stroke()

                # context.arc(i, 23+150*c/d+200*self.rank, 3, 0, 2*math.pi)
                # context.fill_preserve()
                # context.stroke()

            c += 1

class FastaHeader:
    def __init__(self, the_header, the_rank):
        '''FastaHeader object contains Gene name and location relative to FASTA sequence'''

        self.fastaheader = re.findall(r">([a-zA-Z]*)", the_header)[0]
        self.rank = the_rank

    ## methods
    def draw(self, context):
        context.set_source_rgb(0, 0, 0)
        context.move_to(5, 15+200*self.rank)
        context.show_text(self.fastaheader)

class GeneGroup:
    def __init__(self, header, sequence, motifs, rank):
        self.fastaheader = FastaHeader(header, rank)
        self.gene = Gene(sequence, rank)
        self.exon = Exon(sequence, rank) # list of exons
        self.motif = Motif(sequence, motifs, rank) # list of motifs

    ## methods
    def draw(self, context):
        self.fastaheader.draw(context)
        self.gene.draw(context)
        self.exon.draw(context)
        self.motif.draw(context)

def main(file1:str, list_motifs:list)->None:
    '''Main function that reads in FASTA file and parses out sequence and header information.
    Calls other functions as needed.'''
    # Pull file name
    file_name = re.search(r"(.*)(\.fasta)", file1)

    # code to start drawing
    print("Drawing...")
    width, height = 1000, 200*(10)
    surface = cairo.SVGSurface(file_name.group(1)+".svg", width, height)
    context = cairo.Context(surface)

    # intialize a empty string and dictionary
    seq = "" # to store sequence information
    rank = 1 # height rank for proper drawing of genes

    with open(file1, "r") as fh:
        for line in fh:
            if line[0] == ">": # pull header information
                if seq != "": # if there is a sequence stored from previous loop, call feature finding function
                    gene = GeneGroup(header, seq, list_motifs, rank)
                    gene.draw(context)
                    rank += 1
                    seq = "" # reset string to be empty

                # reset with new gene
                header = line
            else: 
                # remove newline character and add line to sequence string
                line = line.rstrip("\n")
                seq += line

    # account for that last case
    if seq != "":
        gene = GeneGroup(header, seq, list_motifs, rank)
        gene.draw(context)
        rank += 1
        seq = ""

    # code to finish drawing
    ## set up legend template
    context.set_source_rgb(0, 0, 0)
    context.move_to(5, 15)
    context.show_text("Legend")

    d = len(list_motifs)

    # Exon
    context.set_source_rgb(0, 0, 0)
    context.rectangle(5, 25, 5, 5)
    context.fill()

    context.set_font_size(8.0)
    context.move_to(15, 30)
    context.show_text("Exon")

    # Intron
    context.set_source_rgb(0, 0, 0)
    context.set_line_width(2) # for the length of the exon
    context.move_to(5, 28+(200-15)*1/(d+2)) 
    context.line_to(10, 28+(200-15)*1/(d+2))
    context.stroke()

    context.set_font_size(8.0)
    context.move_to(15, 30+(200-15)*1/(d+2))
    context.show_text("Intron")

    # Motifs
    c = 0

    for key in list_motifs:
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

    print("Complete.")

if __name__ == "__main__":
    args = get_args()
    file1, file2 = args.file, args.motif

    list_motifs = readMotifs(file2)
    main(file1, list_motifs)



