#!/usr/bin/env python

#########################################################################################
#Imports
from cairo import SVGSurface, Context, Matrix
from Bio import Seq
from itertools import product

import argparse
import re
import sys
import os
import cairo 
import math
#########################################################################################
def get_arguments():
    parser = argparse.ArgumentParser(description="This is a program that takes two arguments: a FASTA file and a simple text file with a list of up to 10 motifs.  Returns an PNG and SVG files (one per fasta intry) illustrating where the Exon, introns, and motifs are in each fasta entry")
    parser.add_argument("-f","--fasta_file", help="-f <path><file>, Requires a fasta file, assumes sequences are in the correct orientation",required=True,type=str)
    parser.add_argument("-m","--motif_file", help="-m <path><file>, Requires a simple text file, each line is a motif",required=True,type=str)

    return parser.parse_args()

args = get_arguments()
fasta_file = args.fasta_file #fasta file
motif_file = args.motif_file #motif

#########################################################################################
#Global Variables

lines={}
org_list=[]
motif_list=[]

#fasta="/Users/mandiedriskill/Desktop/Spliceing_Assignment/INSR.fasta"
#motifs="/Users/mandiedriskill/Desktop/Spliceing_Assignment/test_motif.txt"
#########################################################################################
#dictionaries

color_dict={"Blue":(0,0,205),"Cyan":(0,250,25),"Steel":(0,0,0), "Red":(255,0,0), "purple":(50,0,50), "yellow":(255,255,0),"Pink":(50,.5,.7), "Orange":(250,.4,.1), "Green":(0,.9,.3),"Dark Blue":(.094,.094,.42)}
colors=[*color_dict]
#########################################################################################
#fuctions
def get_ambiguous(seq):
    d = Seq.IUPAC.IUPACData.ambiguous_dna_values
    ra = []
    for i in product(*[d[j] for j in seq]):
        ra.append("".join(i))
    return ra

#########################################################################################
#Exon and Motif Plotting

with open(fasta_file, "r") as fa, open(motif_file, "r") as mo:
	for mo_line in mo:
		mo_line=mo_line.strip("\n")
		motif_cap=mo_line.upper()
		motif_cap=motif_cap.replace("U", "T")
		org_list.append(mo_line)
		exon_motifs=get_ambiguous(motif_cap)
		motif_list.append(exon_motifs)
        

	for line in fa:
		line = line.strip()
		if line.startswith(">"):
			id = line
			lines[id] = " "
		else:
			lines[id] += line 

for header, seq in lines.items():
	parts=header.split()
	name=parts[1]
	line_length=len(seq)
	pattern_exon='[A-Z]+'
	pattern_intron='[a-z]+'
	exon=re.findall(pattern_exon, seq)
	lengths=([len(i) for i in exon])
	exons=[(m.start(0)) for m in re.finditer(pattern_exon, seq)]
    
	width=line_length+100
	height=500
	surface = cairo.SVGSurface(name+".svg", width, height)
	context = cairo.Context(surface)
    
	context.set_line_width(5)
	context.move_to(0,75) #(X,Y)
	context.line_to((line_length), 75)#Match above
	context.stroke()
	pat = cairo.LinearGradient(0.0, 0.0, 0.0, 1.0)
	pat.add_color_stop_rgba(0, .65, .65, .65, 1)

    

	for exon_starts,exon_lengths in zip(exons,lengths):
		context.rectangle(int(exon_starts), 25, int(exon_lengths), 100)
		context.set_source(pat)
		context.fill()
		context.stroke()
        
	upper_seq=seq.upper()

	ctr=0
	for motif in motif_list:
		ctr_color=colors[ctr]
		motif_color=color_dict[ctr_color]
        
		pat1 = cairo.LinearGradient(0.0, 0.0, 0.0, 1)
		pat1.add_color_stop_rgba(1,motif_color[0],motif_color[1],motif_color[2], .5)
    
		legend_ctr=ctr*25
		context.rectangle(0,(180+legend_ctr),15,15)
		context.set_source(pat1)
		context.fill()
		context.select_font_face("pyrus", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
		context.set_font_size(20)

		context.move_to(20, 195+ legend_ctr)
		context.set_source_rgb(0, 0, 0)
		context.show_text(org_list[ctr])
		context.set_source(pat1)
		context.fill()
		context.stroke()
        
		ctr+=1
		for item in motif:            
			pattern_motif=item
			motif_start=[(m.start(0)) for m in re.finditer(pattern_motif, upper_seq)]            
			motif_lengths=re.findall(pattern_motif, upper_seq)
			motif_length=([len(i) for i in motif_lengths])
			for mo_start,mo_length in zip(motif_start,motif_length):
				context.rectangle(int(mo_start), 25, int(mo_length), 100)

			context.set_source(pat1)
			context.fill()
			context.stroke()
            
	context.rectangle(0,150,15,15)
	context.set_source(pat)
	context.fill()
	context.select_font_face("pyrus", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
	context.set_font_size(20)
	context.move_to(20, 165)
	context.set_source_rgb(0, 0, 0)
	context.show_text("Exon")

	context.move_to(100, 160)
	context.line_to(80,160)
	context.stroke()
	context.move_to(110, 165)
	context.show_text("Intron")
	
	surface.write_to_png(name+".png")
    

        
        

    

    
    
    
    
    
    



