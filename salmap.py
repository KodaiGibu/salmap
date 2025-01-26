import os
import subprocess
import argparse
import shutil
import glob
from datetime import datetime

# 引数のパーサーを設定
parser = argparse.ArgumentParser(description='RNA-seq Analysis with salmon')
parser.add_argument('-raw_dir', '-d', type=str, help='生データのディレクトリ')
parser.add_argument('-output', '-o', type=str, help='アウトプットするディレクトリ')
parser.add_argument('-ref', '-r', type=str, help='リファレンスとなるトランスクリプトームfasta')
parser.add_argument('-primer_1', '-1', type=str, help='アダプター配列1')
parser.add_argument('-primer_2', '-2', type=str, help='アダプター配列2')
args = parser.parse_args()

class AnalysisPipeline:
    def __init__(self, raw_dir, output, ref, primer_1, primer_2):
        self.raw_dir = raw_dir
        self.output = output
        self.ref = ref
        self.primer_1 = primer_1
        self.primer_2 = primer_2
        self.id_list = []
    
    def get_started_datetime(self): # 現在の日付と時刻を取得
        start = datetime.now()

    def setup_output_directory(self):
        os.makedirs(self.output, exist_ok=True)
        os.chdir(self.output)

    def generate_statistics_and_manifest(self):
        subprocess.run(f"seqfu stats {self.raw_dir}/*gz --csv > rawdata_stat.csv", shell=True)
        subprocess.run(f"seqfu metadata --format manifest {self.raw_dir}/ > manifest.tsv", shell=True)
        with open("manifest.tsv", "r") as infile, open("temp.tsv", "w") as outfile:
            next(infile)  # Skip the first line
            for line in infile:
                outfile.write(line)
        with open("temp.tsv", "r") as file:
            self.id_list = [line.strip().split('\t')[0] for line in file]

    def create_adapters_file(self):
        with open("adapters.fa", "w") as file:
            file.write(">adapter1/1\n")
            file.write(f"{self.primer_1}\n")
            file.write(">adapter1/2\n")
            file.write(f"{self.primer_2}\n")

    def run_trimmomatic(self):
        with open("temp.tsv", "r") as file:
            lines = file.readlines()
        for line in lines:
            ID, raw1, raw2 = line.strip().split('\t')
            subprocess.run(f"trimmomatic PE -threads 4 -phred33 -trimlog {ID}_TrimLog.txt {raw1} {raw2} "
                           f"{ID}_R1_trimmed.fastq.gz {ID}_R1_unpaired.fastq.gz {ID}_R2_trimmed.fastq.gz {ID}_R2_unpaired.fastq.gz "
                           f"ILLUMINACLIP:adapters.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:36", shell=True)

    def run_fastqc(self):
        os.makedirs("fastqc", exist_ok=True)
        for id in self.id_list:
            subprocess.run(f"fastqc -t 10 --nogroup -o fastqc -f fastq {id}_R1_trimmed.fastq.gz {id}_R2_trimmed.fastq.gz", shell=True)

    def run_salmon(self):
        os.makedirs("salmon_out", exist_ok=True)
        subprocess.run(f"salmon index -t {self.ref} -i ./salmon_out/transcripts_index_salmon -k 31", shell=True)
        os.makedirs("salmon_result", exist_ok=True)
        for id in self.id_list:
            result_dir = f"{id}_exp_salmon"
            os.makedirs(f"salmon_result/{result_dir}", exist_ok=True)
            subprocess.run(f"salmon quant -i ./salmon_out/transcripts_index_salmon -p 4 --validateMappings -l A "
                           f"-1 {id}_R1_trimmed.fastq.gz -2 {id}_R2_trimmed.fastq.gz -o salmon_result/{result_dir}", shell=True)

    def extract_mapping_rate(self):
        with open("MappingRate.txt", "w") as outfile:
            for id in self.id_list:
                with open(f"salmon_result/{id}_exp_salmon/logs/salmon_quant.log", "r") as infile:
                    for line in infile:
                        if "Mapping rate" in line:
                            rate = line.split('= ')[1].strip()
                            outfile.write(f"{rate}\n")
        with open("MappingRate.txt", "r") as file:
            mapping_rate = [line.strip() for line in file]
        with open("MappingRate.tsv", "w") as file:
            for i, id in enumerate(self.id_list):
                file.write(f"{id}\t{mapping_rate[i]}\n")

    def get_version(self):
        seqfu_version = subprocess.check_output(['seqfu', '-v']).decode('utf-8').strip()
        trimmomatic_version = subprocess.check_output(['trimmomatic', '-version']).decode('utf-8').strip()
        fastqc_version = subprocess.check_output(['fastqc', '-v']).decode('utf-8').strip()
        salmon_version = subprocess.check_output(['salmon', '-v']).decode('utf-8').strip()

    def output_parameters(self):
        with open("reports.txt", "w") as file:
            file.write(f"Date and time the run started: {start().strftime('%Y-%m-%d %H:%M:%S')}\n")
            file.write(f"Date and time when the run ended: {datetime.end().strftime('%Y-%m-%d %H:%M:%S')}\n")
            file.write(f"Comand: salmap.py -d {self.raw_dir} -o {self.ref} -1 {self.primer_1} -2 {self.primer_2}\n")
            file.write(f"Reference: {self.ref}\n")
            file.write(f"adapter1: {self.primer_1}\n")
            file.write(f"adapter2: {self.primer_2}\n")
            file.write(f"\n")
            file.write(f"#version\n")
            file.write(f"seqfu {seqfu_version}\n")
            file.write(f"trimmomatic {trimmomatic_version}\n")
            file.write(f"{fastqc_version}\n")
            file.write(f"{salmon_version}\n")
            file.write(f"\n")
            file.write(f"salmap version 1.0.0\n")
            file.write(f"Thank you for using!\n")

    def organize_files(self):
        os.makedirs('salmap_out', exist_ok=True)
        os.makedirs('tmp', exist_ok=True)
        os.makedirs('trimmed', exist_ok=True)
        shutil.move('rawdata_stat.csv', 'salmap_out/rawdata_stat.csv')
        shutil.move('reports.txt', 'salmap_out/reports.txt')
        shutil.move('MappingRate.tsv', 'salmap_out/MappingRate.tsv')
        shutil.move('salmon_result', 'salmap_out/salmon_result')
        shutil.move('fastqc', 'salmap_out/fastqc')
        for file in glob.glob('*trimmed.fastq.gz'):
            shutil.move(file, 'trimmed/' + os.path.basename(file))
        for file in glob.glob('*_TrimLog.txt'):
            shutil.move(file, 'trimmed/' + os.path.basename(file))
        shutil.move('trimmed', 'tmp/trimmed')
        shutil.move('salmon_out', 'tmp/salmon_out')
        os.remove('MappingRate.txt')
        os.remove('adapters.fa')
        os.remove('manifest.tsv')
        for file in glob.glob('*_unpaired.fastq.gz'):
            os.remove(file)
        os.remove('temp.tsv')

    def run_pipeline(self):
        self.get_started_datetime()
        self.setup_output_directory()
        self.generate_statistics_and_manifest()
        self.create_adapters_file()
        self.run_trimmomatic()
        self.run_fastqc()
        self.run_salmon()
        self.extract_mapping_rate()
        self.get_version()
        self.output_parameters()
        self.organize_files()

if __name__ == "__main__":
    pipeline = AnalysisPipeline(args.raw_dir, args.output, args.ref, args.primer_1, args.primer_2)
    pipeline.run_pipeline()
