# salmap
This script is for mapping and counting RNA-seq reads to transcriptome data.  
RNA-seqのリードをトランスクリプトーム配列にマッピングし、リード数をカウントするためのスクリプトです。

# Download
```bash:
git clone https://github.com/KodaiGibu/salmap.git
```

# Environment Construction
```bash:
mamba create -n salmap python=3.12 -y
conda activate salmap

mamba install -c bioconda seqfu==1.22.3 -y
mamba install -c bioconda trimmomatic==0.39 -y
mamba install -c bioconda fastqc==0.12.1 -y
mamba install -c bioconda salmon==1.10.3 -y
```

# Usege
```bash:
chmod +x salmap.py #add permisson
python salmap.py -d /path/to/rawdata_dir -o /path/to/out_put_dir/ -r /path/to/rna.fna -1 AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG 
```

# Options
`-d`　Path of the directory containing rawdata  
`-r`　Fasta file of the reference transcriptome sequence  
`-1`　Adapter sequence information to be removed #1   
`-2`　Adapter sequence information to be removed #2   
`-o`　Path of the directory to output the output file  





