#!/usr/bin/env python

'''
Assignment 1 
Bi623
August 18, 2016

Tabular Format for Blastn/blastp columns 
 1.	 qseqid	 query (e.g., gene) sequence id
 2.	 sseqid	 subject (e.g., reference genome) sequence id
 3.	 pident	 percentage of identical matches
 4.	 length	 alignment length
 5.	 mismatch	 number of mismatches
 6.	 gapopen	 number of gap openings
 7.	 qstart	 start of alignment in query
 8.	 qend	 end of alignment in query
 9.	 sstart	 start of alignment in subject
 10.	 send	 end of alignment in subject
 11.	 evalue	 expect value
 12.	 bitscore	 bit score
'''


#read in the query results of blastp 

#file1 = "/Users/jsmith/Documents/bi623/160818_class4_PS1/PS1/stickbyplat_test.txt"
#file2 = "/Users/jsmith/Documents/bi623/160818_class4_PS1/PS1/ortholog_test.txt"
file1 = "/home12/jsmith16/bi623/160818_PS1/blastp_results/stickbyplat_10hits.txt"
file2 = "/home12/jsmith16/bi623/160818_PS1/blastp_results/platbystick_10hits.txt" 
blastquery1 = open(file1, "r")
blastquery2 =  open(file2, "r")

#initialize and empty dictionary to hold the queries, hits, and evalues
#this is for stickleback(Gasterosteus_aculeatus = GA) blasted against platyfish
Ga_top_hits = {}


#for loop to iterate through the blastp output (in tabular format). 
for entry in blastquery1:
    entry = entry.strip('\r') #strip return characters
    entry = entry.strip('\n') #strip new line characters
    entry = entry.split('\t') #split each line into an array to select for gene_id based on a tab delimiter  
    query_id = entry[0] #query_id is the Ga gene
    evalue = float(entry[10]) #evalue of sequence match
    subject_id = entry[1] #matched gene from the reference DB, platyfish in this case
    if query_id not in Ga_top_hits: #add the query_id to the dictionary, top_hits, the first time it is seen
        Ga_top_hits[query_id] = [subject_id, evalue] #value of the key, query_id, is an array with two values, subject_id and evalue
    if query_id in Ga_top_hits:
        vals_top_hits = Ga_top_hits.get(query_id) #this is the current values in the dictionary top_hits for each query_id 
        top_hit_eval = float(vals_top_hits[1]) #select the 2nd position in the current values, which corresponds to the evalue
        if evalue < top_hit_eval: #if the evalue is less than the current evalue for a gene
            Ga_top_hits[query_id] = [subject_id, evalue] #update the dictionary with the hit/match gene and its evalue
             

'''
#check that the output is correct 
for key, val in Ga_top_hits.items():
    print(key, val)
'''

#initialize and empty dictionary to hold the queries, hits, and evalues
#this is for platyfish( Xiphophorus_maculatus= Xm) blasted against stickleback
Xm_top_hits = {}

#for loop to iterate through the blastp output (in tabular format). 
for entry in blastquery2:
    entry = entry.strip('\r')
    entry = entry.strip('\n')
    entry = entry.split('\t') #split each line into an array to select for gene_id based on a tab delimiter  
    query_id = entry[0]
    evalue = float(entry[10])
    subject_id = entry[1]
    if query_id not in Xm_top_hits: #add the query_id to the dictionary, top_hits, the first time it is seen
        Xm_top_hits[query_id] = [subject_id, evalue] #value of the key, query_id, is an array with two values, subject_id and evalue
    if query_id in Xm_top_hits:
        vals_top_hits = Xm_top_hits.get(query_id) #this is the current values in the dictionary top_hits for each query_id 
        top_hit_eval = float(vals_top_hits[1]) #select the 2nd position in the current values, which corresponds to the evalue
        if evalue < top_hit_eval:
            Xm_top_hits[query_id] = [subject_id, evalue] #update the dictionary for this gene with lower evalue

'''
for key, val in Xm_top_hits.items():
    print(key, val)          
'''

#initialize an empty array
Ga_hits = []

#for loop to list the items of the Ga_top_hits dictionary.
for key,val in Ga_top_hits.items():
    Ga_gene = key   #gene is the key
    Xm_hit = val[0] #select position 1 of the values in the dictionary, this is the hit/match gene
    Ga_hits.append([Ga_gene, Xm_hit]) #append the Ga gene and the platyfish hit to an array

#sort the array by gene names
Ga_hits.sort()

#initialize an empty array
Xm_hits = []

#for loop to list the items in the Xm_top_hits dictionary
for key, val in Xm_top_hits.items():
    Xm_gene = key   #the gene is the key 
    Ga_hit = val[0] #select position 1 of the values in the dictionary, this is the hit/match
    Xm_hits.append([Ga_hit, Xm_gene]) #append the gene name and stickleback hit/match to the array
    								 	#ensure that Ga_hit is in column 1 of the array to that
    								 	#it matches the Ga_gene in column 1 of Ga_hits array
    								 	# and Xm_hit is in column 2, so Xm_gene is in column 2

#sort the array by gene names
Xm_hits.sort()

#the number of subject-query hits for each is the length of the top_hits arrays 
Ga_query_subject_hits = len(Ga_top_hits)
Xm_query_subject_hits = len(Xm_top_hits)

#initialize an empty array for the reciprocal best hits
orthologs = []

#for loop to compare the best hits for each species to find Reciprocal Best Hits (RPH) 
#RBH are used to identify orthologs between two species, platyfish and stickleback in this case
for entry in Ga_hits:
    if entry in Xm_hits:
        orthologs.append(entry) #if orthologs are found in the array, then append them to another empty array called orthologs

#number of orthologs is the length of the orthologs array 
num_orthologs = len(orthologs)

#create a file write the orthologs list using append mode. 
orthologsfilePlatyfish = open("orthologsFound.txt", "a")

#write the orthologs to the file
for ortholog in orthologs:
    ortholog = '\t'.join(ortholog) #use join function to remove the brackets and commas from each entry
    orthologsfilePlatyfish.write('%s \n' % ortholog)    

#close the file 
orthologsfilePlatyfish.close()

#create another file for the results of the number of orthologs using append mode.
resultfile1 = open("resultsfile.txt", "a")

#write the number of query-subject matches and the total number of orthologs to a result file
resultfile1.write("The number of subject query matches in Stickleback is %s \n" % (Ga_query_subject_hits))
resultfile1.write("The number of subject query matches in Platyfish is %s \n" % (Xm_query_subject_hits))
resultfile1.write("The number of orthologs between Stickleback and Platyfish is %s \n" % (num_orthologs)) 

#close the file
resultfile1.close()

#print the results to the screen to provide information and indicate the program is complete.
print("The number of subject query matches in Stickleback is", Ga_query_subject_hits,)
print("The number of subject query matches in Platyfish is", Xm_query_subject_hits)
print("The number of orthologs between Stickleback and Platyfish is", num_orthologs)

        
        
        
        
        
        
        
        
        
        
        
        
        
        