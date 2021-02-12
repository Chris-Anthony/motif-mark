# Motif-Mark
#### version: v2.0
#### updated: 02/12/2021
#### created: 02/04/2021
#### author: Chris Chua

## Summary

This script takes a FASTA file and a list of motifs (.txt) and generates a .SVG image showing the motifs locations on each sequence. 
The name of the generated .SVG file is the same as the file name of the FASTA file. 
Currently, this script can only handle sequences up to 1kb in length and up to five motifs. 

## Set-Up

Ensure that you have the pycairo package and python v3.7.6+ installed. 
Download motif-mark in the location of your fasta files.

## Running Motif-Mark

./motif-mark-v2.py -f Figure_1.fasta -m Fig_1_motifs.txt 

## Commands

-f	--file	Name of fasta file

-m	--motif	Name of motif file
