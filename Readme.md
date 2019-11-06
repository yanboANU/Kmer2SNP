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

# Input Format

The input of **Kmer2SNP**  can be fasta, fastq, either gzipped or not. 

# Example Usage

<pre><code>  sh runkmercalling.sh k cov single_chr22_30x_data.fq  </code></pre>

K is the size of k-mer, cov is reads coverage.

# References
