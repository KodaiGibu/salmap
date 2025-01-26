# salmap
This script is for mapping and counting RNA-seq reads to transcriptome data.
RNA-seqのリードをトランスクリプトーム配列にマッピングし、リード数をカウントするためのスクリプトです。

# Download
```bash:
git clone https://github.com/KodaiGibu/rnamap.git
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
`-d`　rawdataのあるディレクトリ  
`-r`　リファレンスとなるトランスクリプトーム配列のfastaファイル   
`-1`　除去するアダプター配列情報 #1  
`-2`　除去するアダプター配列情報 #2  
`-o`　アウトプットを出力するディレクトリ
