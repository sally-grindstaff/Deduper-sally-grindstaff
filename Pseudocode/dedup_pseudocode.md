## Reference-based removal of PCR duplicates: Pseudocode  
### Sally Grindstaff  
  
**The problem:** 

read in UMI file (STL96.txt) and put UMI sequences in a set, close UMI file

open coordinate-sorted input SAM file for reading
open output SAM file for writing
create read_info_dict (empty dictionary)
loop through input file line by line
  output all lines that start with '@' (header section)
  split line by tabs, which will store each column entry as a item in a list
  call adjust_start_pos function on indexes 3 and 5 of list (4th and 6th columns of SAM entry) and assign output to a variable called adj_start_pos
  call get_strand function on index 1 (bitwise flag) of list
  call make_str function on adj_start_pos, strand, and index 2 (chromosome) of list and assign output to a variable called pos
  call get_UMI function on index 0 (QNAME) of list and assign output to variable called umi
  check whether umi is already a key in the dictionary
    if yes, check whether pos is a value of that key
      if yes, go to beginning of loop
      if no, add pos as value under that key and write line to output file
    if no, add umi as a key to the dictionary and add pos as value of that key and write line to output file
 close input and output files
  
 functions:
 
 def adjust_start_pos(start_pos: str, CIGAR: str) -> str:
    # takes two strings (start position and CIGAR string) and returns string of adjusted start position, which is the position where the read starts actually aligning
    convert start_pos to int
    use regex to check if there is soft clipping at the beginning, which will be represented by a number and an S at the beginning of the CIGAR string. ([0-9]+)S.*
      if there is soft clipping (if regex search returns a number), then convert returned number to int and subtract from start_pos. return string of this number
      if no soft clipping (if regex search returns null), then return string of start_pos
      
adjust_start_pos(222, 71M)
  returns 222
adjust_start_pos(222, 2S69M)
  returns 220
   
def get_strand(flag: str) -> str:
    # takes bitwise flag and returns a string that indicates the strand (either + or -)
    convert flag string into int and check whether 16th bit is set
      if yes, return '-'
      if no, return '+'
   
get_strand(0)
  returns +
get_strand(16)
  returns -
  
 def make_str(adj_start_pos: str, strand: str, chrom: str) -> str:
    # takes adjusted start position, strand, and chromosome; extracts bitwise flag from qname; returns string containing adjusted start pos, strand, and chrom
    concatenate the three strings into one string
    
make_str(222, +, 2)
  returns 222+2
    
def get_UMI(qname:str) -> str:
    # takes QNAME and extracts UMI
    split QNAME string by ':' and return the last item in the list
    
get_UMI(NS10050:001100:AATTCCGG)
  returns AATTCCGG
