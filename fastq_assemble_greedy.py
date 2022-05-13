from Bio import SeqIO


def calculateOverlap(read1, read2):
    #By Ryan Byrne
    #Calculates the overlap between two sequences
    
    #Start at last char of first string and first char of second string
    i = -1
    j = 0
    count = 0
    
    #Iterates to the minimum length of one of the strings
    minLin = min(len(read1), len(read2))
    while ( j < minLin):
        #For window size of 1 test first and last character
        if (j == 0):
            if (read1[i] == read2[j]):
                count = 1
        #For all other window sizes check if characters in the windows of both strings are
        #the same
        else:
            if (read1[i:] == read2[0:j + 1]):
                count = j + 1
        i -= 1
        j += 1
    return count

def buildGraph(Reads):
    #By Ryan Byrne
    #Builds Overlap Graph from a list of reads
    
    Graph = {}
    
    #Go through reads and build graph based on the overlap between two reads
    for read in Reads:
        Edges = {}
        #Compare each read to all others
        for otherreads in Reads:
            if (otherreads != read):
                overlap = calculateOverlap(read, otherreads)
                #If overlap exists add to graph
                if (overlap > 0):
                    Edges[otherreads] = overlap
        Graph[read] = Edges
    return Graph

def mergeReads(read1, read2, overlap):
    #By Ryan Byrne
    #Merges two read strings together based on overlap
    
    mergeString = ''
    
    #Sets output strings to all chars in first read
    for char in read1:
        mergeString += char
    
    #Sets the output strings of all strings after overlap in the second read
    for i in range(overlap, len(read2)):
        mergeString += read2[i]
    return mergeString

def calculateContig(Reads, Graph):
    #By Ryan Byrne
    #Calculates contigs of reads based on overlap graph
    
    #Variables for determining max overlap
    maxEdge = -1
    mergeString = ''
    maxEdgeVert1 = ''
    maxEdgeVert2 = ''
    
    #If the list only has one item there is no more merging so return only contig
    if (len(Reads) == 1):
        return Reads
    
    #Iterate through graph to find the largest overlap
    for vertice in Graph:
        for edge in Graph[vertice]:
            weight = Graph[vertice][edge]
            
            #If overlap is larger than others set vertix as max
            if (weight > maxEdge):
                maxEdge = weight
                mergeString = mergeReads(vertice, edge, weight)
                maxEdgeVert1 = vertice
                maxEdgeVert2 = edge
            
            #If the overlap equals largest overlap if the merged string is
            #alphabetically ahead use that string as max
            elif (weight == maxEdge):
                testString = mergeReads(vertice, edge, weight)
                if (testString < mergeString):
                    mergeString = testString
                    maxEdgeVert1 = vertice
                    maxEdgeVert2 = edge
    
    #If there is no edges no more merging can be down so output list
    if (maxEdge == -1):
        return Reads
    
    #Remove from the reads that were merged into a merged string
    Reads.remove(maxEdgeVert1)
    Reads.remove(maxEdgeVert2)
    
    #Add merged string
    Reads.append(mergeString)
    
    #Build new graph based on merged reads list
    Overlap_Graph = buildGraph(Reads)
    
    #Recursively Call function with new reads list and graph
    contig = calculateContig(Reads, Overlap_Graph)
    return contig

def fastq_assemble_greedy(Reads):
    #By Ryan Byrne
    #Assignment io function
    
    try:
        
        #Try to open file if input is a file
        blast_records = SeqIO.parse(Reads, 'fastq')
        
        #Store sequences in fastq file to list
        Seqs = []
        for record in blast_records:
            Seqs.append(str(record.seq))
        
        #Build inital graph
        Overlap_Graph = buildGraph(Seqs)
        
        #Calculate contig of file
        contig = calculateContig(Seqs, Overlap_Graph)
        
        return contig
    except:
        
        #Build inital graph
        Overlap_Graph = buildGraph(Reads)
        
        #Calculate contig of Reads
        contig = calculateContig(Reads, Overlap_Graph)
        
        return contig
                        
                    
        


