##@title Analyze your protein
from google.colab import drive
drive.mount('/content/drive')
!cp /content/drive/MyDrive/WuLab/esm1b_t33_650M_UR50S-contact-regression.pt ./

import os
from google.colab import files
import datetime
import re

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

########## input
INPUT = "RPL3L"#@param ["RPL3L", "MYC"] {allow-input: true}

#@markdown - Input format: one raw protein sequence; space allowed
#@markdown - Example: copy & paste a multi-line sequence from a FASTA file (without the header)
#@markdown - To run: click `Runtime` -> `Run all` in the menu bar, or click the triangle play/run button on the left

MUTATION = None #@param {allow-input: true}
#@markdown - Input format: seperate each mutation using commas. Each mutation should contain the wild type amino acid, its position in the sequence, and the mutated amino acid e.g."M1A, E2C"
#@markdown - Specify the sequence that contains the mutations in the input above

INDEL = None #@param {allow-input: true}
#@markdown - Input format: seperate using commas the mutated sequence with indel mutations and the index of the mutation
#@markdown - Specify the wild type sequence in the input above

seq = INPUT
mut = MUTATION
indel = INDEL

if seq == "RPL3L":
  seq = "MSHRKFSAPRHGHLGFLPHKRSHRHRGKVKTWPRDDPSQPVHLTAFLGYKAGMTHTLREVHRPGLKISKREEVEAVTIVETPPLVVVGVVGYVATPRGLRSFKTIFAEHLSDECRRRFYKDWHKSKKKAFTKACKRWRDTDGKKQLQKDFAAMKKYCKVIRVIVHTQMKLLPFRQKKAHIMEIQLNGGTVAEKVAWAQARLEKQVPVHSVFSQSEVIDVIAVTKGRGVKGVTSRWHTKKLPRKTHKGLRKVACIGAWHPARVGCSIARAGQKGYHHRTELNKKIFRIGRGPHMEDGKLVKNNASTSYDVTAKSITPLGGFPHYGEVNNDFVMLKGCIAGTKKRVITLRKSLLVHHSRQAVENIELKFIDTTSKFGHGRFQTAQEKRAFMGPQKKHLEKETPETSGDL"
elif seq == "MYC":
  seq = "MDFFRVVENQQPPATMPLNVSFTNRNYDLDYDSVQPYFYCDEEENFYQQQQQSELQPPAPSEDIWKKFELLPTPPLSPSRRSGLCSPSYVAVTPFSLRGDNDGGGGSFSTADQLEMVTELLGGDMVNQSFICDPDDETFIKNIIIQDCMWSGFSAAAKLVSEKLASYQAARKDSGSPNPARGHSVCSTSSLYLQDLSAAASECIDPSVVFPYPLNDSSSPKSCASQDSSAFSPSSDSLLSSTESSPQGSPEPLVLHEETPPTTSSDSEEEQEDEEEIDVVSVEKRQAPGKRSESGSPSAGGHSKPPHSPLVLKRCHVSTHQHNYAAPPSTRKDYPAAKRVKLDSVRVLRQISNNRKCTSPRSSDTEENVKRRTHNVLERQRRNELKRSFFALRDQIPELENNEKAPKVVILKKATAYILSVQAEEQKLISEEDLLRKRREQLKHKLEQLRNSCA"
else: # user input
  # clean up sequence: upper case, remove space
  seq = seq.upper().replace(' ','')
  # if contains non aa letters:
  if not all(char in 'ACDEFGHIKLMNPQRSTVWY' for char in seq):
    print("\n\n")
    print('\n'+ bcolors.BOLD +bcolors.FAIL + "WARNING: Your sequence contains letters other than ACDEFGHIKLMNPQRSTVWY!"+bcolors.ENDC)
    L0  = len(seq)
    seq = re.sub('[^ACDEFGHIKLMNPQRSTVWY]+', '', seq)
    L1 = len(seq)
    print(L0-L1,'non-aa letters removed!'+bcolors.ENDC)
    exit()

if mut:
  mut = mut.replace(' ','')
if indel:
  indel = indel.replace(' ','')

######### options

# set model

MODEL = "esm1b_t33_650M_UR50S" #@param ["esm1v_t33_650M_UR90S_1", "esm1b_t33_650M_UR50S"]

# remove files from a previous run
if os.path.exists("ESMScan-all-mutants.txt"):
  print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+': Removing files from a previous run')
  !rm ESMScan-* res.zip run.sh

if not os.path.exists("ESM-Scan"):
  print("\n")
  print('\n\n'+ bcolors.BOLD +bcolors.OKBLUE + "Installing packages"  +bcolors.ENDC)
  print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
  !pip install biopython
  !pip install fair-esm
  !git clone https://github.com/katherineda/ESM-Scan.git
  !cd /content
  print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

model_location="/content/"+MODEL+".pt"
if not os.path.exists(model_location ):
  print('\n\n'+ bcolors.BOLD +bcolors.OKBLUE + "Downloading pre-trained ESM model"  +bcolors.ENDC)
  if MODEL == "esm1b_t33_650M_UR50S":
    !wget https://dl.fbaipublicfiles.com/fair-esm/models/esm1b_t33_650M_UR50S.pt
  else:
    !wget https://dl.fbaipublicfiles.com/fair-esm/models/esm1v_t33_650M_UR90S_1.pt

print('\n\n'+ bcolors.BOLD +bcolors.OKBLUE + "Running saturation mutagenesis"  +bcolors.ENDC)

if type(mut) == str:
  cmd="python /content/ESM-Scan/esmscan.py --model-location "+model_location+" --sequence "+seq + " --dms-mutation " + mut
elif type(indel) == str:
  cmd="python /content/ESM-Scan/esmscan.py --model-location "+model_location+" --sequence "+seq + " --dms-indel " + indel + " --scoring-strategy " + "indel"
else:
  cmd="python /content/ESM-Scan/esmscan.py --model-location "+model_location+" --sequence "+seq


print(cmd)

with open("run.sh",'w') as f:
  f.write(cmd+'\n')

!chmod +x /content/run.sh
!/content/run.sh

'''
import subprocess

proc = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)

(out, err) = proc.communicate()
print("Screen output:", out)
print("Screen error:", err)
'''
#os.system(cmd)

print('\n\n'+ bcolors.BOLD +bcolors.OKBLUE + "Downloading results"  +bcolors.ENDC)

if os.path.exists('ESMScan-res-in-matrix.csv'):
  os.system(f'zip res.zip *.pdf *.csv')
  files.download(f"res.zip")
  print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+': Done! Please see results in res.zip')
elif type(mut) == str or type(indel) == str:
  os.system(f'zip res.zip *.pdf *.csv')
  files.download(f"res.zip")
  print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+': Done! Please see results in res.zip')
else:
  print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+': No output files generated')
