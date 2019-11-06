# Kmer2SNP: 

**Kmer2SNP** is a fast and accurate variants calling tool for next-generation sequencing data. It's a reference-free and alignment free method.
 
# Downloading Kmer2SNP

To download **Kmer2SNP**, you have to install the following software
<pre><code>
(1) DSK: 
https://github.com/GATB/dsk
(2) findGSE: 
https://github.com/schneebergerlab/findGSE/blob/master/INSTALL
(3) Python package networkx
https://networkx.github.io/documentation/networkx-1.11/install.html

  
Then clone the **Kmer2SNP** repository to your machine.
<pre><code> git clone https://github.com/yanboANU/Kmer2SNP.git </code></pre>

After you download Kmer2SNP, you need to change path in the runkmercalling.sh.

# Input Format

The input of **Kmer2SNP**  can be fasta, fastq, either gzipped or not. 

# Ouput Format

See Kmer2SNP/example, k_31_pair.snp stores isolated SNP kmer pair, k_31_pair.non and k_31_pair.non.sep store non-isolated SNP kmer pair. 
k_31_pair.non.sep put one SNP in the middle position and length of kmer pair is 31+31 if you set k=31. 
k_31_pair.non store the result after merge two non-isolated SNP kmer pair. 

# Example Usage

<pre><code>  sh runkmercalling.sh k cov single_chr22_30x_data.fq  </code></pre>

K is the size of k-mer, cov is reads coverage.


