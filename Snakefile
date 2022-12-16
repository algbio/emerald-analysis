configfile: 'parameters.yaml'

from datetime import datetime
import urllib
import os

print("Preprocessing config", flush = True)

USERNAME = os.environ.get("USER")

# Paths
DB_FILE = config["db_file"].replace("<username>", USERNAME)
DIAMOND_PATH = config["diamond"].replace("<username>", USERNAME)
RAXML_PATH = config["raxml"].replace("<username>", USERNAME)
MUSCLE_PATH = config["muscle"].replace("<username>", USERNAME)
HMMER_PATH = config["hmmer"].replace("<username>", USERNAME)
SAFETY_PATH = config["safety"].replace("<username>", USERNAME)
MMSEQS_PATH = config["mmseqs"].replace("<username>", USERNAME)
#TMSCORE_PATH = config["tmscore"].replace("<username>", USERNAME)
WORK_DIR = config["wrkdir"].replace("<username>", USERNAME)
DATA_DIR = config["datadir"].replace("<username>", USERNAME)
TEMP_DIR = config["tempdir"].replace("<username>", USERNAME)
STRIDE = config["stride"].replace("<username>", USERNAME)
DSSP = config["dssp"].replace("<username>", USERNAME)
PDB_DIR = config["pdbdir"].replace("<username>", USERNAME)
FAMILY_DB = config["family_db"].replace("<username>", USERNAME)

# Parameters
MIN_IDENTITY = config["min_identity"]
CLUSTER_ALGO = config["clustering_algorithm"]
CLUSTER_SENSITIVITY = config["sensitivity"]
FAMILY_PROFILE = FAMILY_DB.replace("seed", "hmm")
REF_CRITERION = config["ref_criterion"]
UNIFORM = (f"--uniform --bin_size {config['bin_size']} --bin_width {config['bin_width']}") if config["uniform"] else ""
ALPHA_SAFE = config["alpha_safe"]
DELTA_SAFE = config["delta_safe"]
BENCHMARK = config["benchmark"]
# turns 0.75f into ".a75"
ALPHA_SAFE_S = "a" + f"{ALPHA_SAFE:.2f}".split(".")[-1]

if CLUSTER_SENSITIVITY == "auto":
    min_id = int(MIN_IDENTITY)
    if min_id <= 25:
        CLUSTER_SENSITIVITY = "--ultra-sensitive"
    elif min_id <= 40:
        CLUSTER_SENSITIVITY = "--very-sensitive"
    elif min_id <= 90:
        CLUSTER_SENSITIVITY = "--sensitive"
    elif 90 < min_id:
        CLUSTER_SENSITIVITY = "--fast"
    else:
        assert False, "Error choosing the sensitivity"


# removes path and file extension from filename
def clean_filename(filename):
    return filename.split("/")[-1].split(".")[0]

DB_FILENAME = clean_filename(DB_FILE)
WORK_DIR = os.path.join(WORK_DIR, f"{DB_FILENAME}.{MIN_IDENTITY}.{CLUSTER_ALGO}")
CLUSTER_IDS, = glob_wildcards(os.path.join(WORK_DIR, "fasta", "{id}.fasta"))


def get_timestamp():
    return datetime.today().strftime('%Y-%m-%dT%H:%M')

print(f"WORK DIR: {WORK_DIR}")
print(f"DATA DIR: {DATA_DIR}")
print(f"TEMP DIR: {TEMP_DIR}")
LOG_DIR = os.path.join("logs", get_timestamp(), "%j.out")
print(f"LOG DIR: {LOG_DIR}")
print("Finished config preprocessing", flush = True)



rule all:
    input:
        #expand(os.path.join(WORK_DIR, f"safety.{ALPHA_SAFE_S}.d{DELTA_SAFE}", "{id}.out"), id=CLUSTER_IDS),
        expand(os.path.join(WORK_DIR, "phmmer", "{id}.out"), id=CLUSTER_IDS),
        expand(os.path.join(WORK_DIR, "hmmsearch", "{id}.out"), id=CLUSTER_IDS),
        expand(os.path.join(WORK_DIR, "hmmscan", "{id}.out"), id=CLUSTER_IDS),
        # expand(os.path.join(WORK_DIR, "id", "{id}.out"), id=CLUSTER_IDS),

def get_ref_reqs():
    if REF_CRITERION == "--identity" or REF_CRITERION == "--highlow":
        return "id"
    elif REF_CRITERION == "--similarity":
        return "hmmsearch"
    return ""

rule change_ref:
    input:
        db = DB_FILE,
        path = WORK_DIR,
        req = WORK_DIR + "/" + get_ref_reqs()
    threads: workflow.cores
    resources:
        logdir = LOG_DIR,
        mem_mb = 2048,
        time = "02:00:00",
        queue = "short"
    output:
        temp(os.path.join(WORK_DIR, "validate_ref.txt"))
    shell:
        "python3 scripts/change_ref.py {input.db} {input.path} {REF_CRITERION} | tee {output}"



def estimate_mem(filename,  wildcards, attempt):
    with open(filename, "r") as f:
        f.readline()
        t = len(f.readline())
        return int(min(t * t * 0.0015 + 128, 14000)) * attempt

def estimate_cpu_time(filename, wildcards, attempt):
    with open(filename, "r") as f:
        f.readline()
        k = len(f.readline())
        t = min(k * k / 1000 * 0.3 + 5, 500) * attempt
        return f"{int(t/60):02}:{int(t%60):02}:00"

def get_ref(wildcards, attempt):
    path = os.path.join(WORK_DIR, "refs", wildcards.id + ".ref.fasta")
    with open(path, "r") as f:
        line = f.readline().split(" ")[0][1:]
        if "|" in line:
            return line.split("|")[1]
        return line.split(":")[1]
         

rule safe:
    input:
        expand(os.path.join(WORK_DIR, f"safety.{ALPHA_SAFE_S}.d{DELTA_SAFE}", "{id}.out"), id=CLUSTER_IDS)

rule _safe:
    input:
        fasta = ancient(os.path.join(WORK_DIR, "fasta", "{id}.fasta"))
    output:
        out = os.path.join(WORK_DIR, f"safety.{ALPHA_SAFE_S}.d{DELTA_SAFE}", "{id}.out"),
        benchmark = os.path.join(WORK_DIR, f"safety.{ALPHA_SAFE_S}", f"benchmark", "{id}.out") if BENCHMARK else []
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt, input: estimate_mem(str(input), wildcards, attempt),
        time = lambda wildcards, attempt, input: estimate_cpu_time(str(input), wildcards, attempt),
        ref = lambda wildcards, attempt, input: get_ref(wildcards, attempt),
        logdir = LOG_DIR,
        queue = "short",
    shell:
        "/usr/bin/time -v {SAFETY_PATH} -f {input.fasta} --threads {threads} --alpha {ALPHA_SAFE} --delta {DELTA_SAFE} -r \"{resources.ref}\" > {output.out} 2> {output.benchmark}" if BENCHMARK else "{SAFETY_PATH} -f {input.fasta} --threads {threads} --alpha {ALPHA_SAFE} --delta {DELTA_SAFE} -r \"{resources.ref}\" > {output.out}"


rule identity:
    input:
        ancient(expand(os.path.join(WORK_DIR, "id", "{id}.out"), id=CLUSTER_IDS))

rule _identity:
    input:
        msa = ancient(os.path.join(WORK_DIR, "msa", "{id}.fasta"))
    output:
        os.path.join(WORK_DIR, "id", "{id}.out")
    threads: workflow.cores
    resources:
        logdir = LOG_DIR,
        mem_mb = 512,
        time = "00:05:00",
        queue = "short"
    shell:
        "{HMMER_PATH}/easel/miniapps/esl-alipid {input.msa} > {output}"

# search sequences agaianst a pfam domain database
rule hmmscan:
    input:
        ancient(expand(os.path.join(WORK_DIR, "hmmscan", "{id}.out"), id=CLUSTER_IDS))

rule _hmmscan:
    input:
        FAMILY_PROFILE + ".h3f",
        FAMILY_PROFILE + ".h3i",
        FAMILY_PROFILE + ".h3m",
        FAMILY_PROFILE + ".h3p",
        profile = FAMILY_PROFILE,
        fasta = ancient(os.path.join(WORK_DIR, "clean", "{id}.clean.fasta"))
    output:
        os.path.join(WORK_DIR, "hmmscan", "{id}.out")
    threads: workflow.cores
    resources:
        logdir = LOG_DIR,
        mem_mb = lambda wildcards, attempt, input: 512 * attempt,
        time = lambda wildcards, attempt, input: f"00:{15*attempt:02}:00",
        queue = "short"
    shell:
        "{HMMER_PATH}/src/hmmscan --tblout {output} --noali {input.profile} {input.fasta} > /dev/null"

# def get_pdbs():
#     clusters, accessions = glob_wildcards(os.path.join(WORK_DIR, "pdb", "{cluster}", "{accession}.pdb"))
#     return accessions

def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def parse_fasta(protein_fasta):
    id = protein_fasta.split()[0]
    sequence = ''.join(protein_fasta.split("\n")[1:])
    return (id, sequence)
    
def get_sequences():
    clusters, = glob_wildcards(os.path.join(WORK_DIR, "fasta", "{id}.fasta"))
    a = []
    not_found = 0
    for cluster in clusters:
        with open(os.path.join(WORK_DIR, "fasta", f"{cluster}.fasta"), "r") as f:
            db_fasta = ("\n" + f.read()).split("\n>")[1:]
            for protein in db_fasta:
                s = parse_fasta(protein)[0]
                if ":" in s:
                    accession = s.split(":")[1]
                else:
                    accession = s.split("|")[1]
                path = os.path.join(PDB_DIR, f"AF-{accession}-F1-model_v3.pdb")
                if os.path.isfile(path):
                    # print("found")
                    a.append(os.path.join(cluster, accession))
                else:
                    # print(path)
                    not_found += 1

    if len(a) + not_found > 0:
        print(f"Structures found: {len(a)}, not found: {not_found} ({not_found * 100.0/(len(a) + not_found):.2f}%)")
    return a


rule stride:
    input:
        expand(os.path.join(WORK_DIR, "stride", "{sequence}.out"), sequence=get_sequences())

rule _stride:
    input:
        lambda wildcards: os.path.join(PDB_DIR, "AF-" + wildcards.sequence.split("/")[-1] + "-F1-model_v3.pdb")
    output:
        os.path.join(WORK_DIR, "stride", "{sequence}.out")
    shell:
        "{STRIDE} {input} > {output}"

def get_pdbs():
    clusters, = glob_wildcards(os.path.join(WORK_DIR, "fasta", "{id}.fasta"))
    a = []
    for cluster in clusters:
        with open(os.path.join(WORK_DIR, "fasta", f"{cluster}.fasta"), "r") as f:
            db_fasta = ("\n" + f.read()).split("\n>")[1:]
            for protein in db_fasta:
                s = parse_fasta(protein)[0]
                if ":" in s:
                    accession = s.split(":")[1]
                else:
                    accession = s.split("|")[1]

                a.append(accession)

    return a

rule pdbs:
    input:
        expand(os.path.join(WORK_DIR, "pdb", "{accession}.pdb"), accession=get_pdbs())

rule _pdbs:
    output:
        os.path.join(WORK_DIR, "pdb", "{accession}.pdb")
    shell:
        "wget https://alphafold.ebi.ac.uk/files/AF-{wildcards.accession}-F1-model_v3.pdb -O {output}"

# rule rep_dssp:
#     input:
#         expand(os.path.join(WORK_DIR, "dssp", "{accession}.out"), accession=get_rep_pdbs())

# rule _rep_dssp:
#     input:
#         os.path.join(WORK_DIR, "rep_pdb", "{accession}.pdb")
#     output:
#         os.path.join(WORK_DIR, "dssp", "{accession}.out") 
#     shell:
#         "mkdir -p $(dirname {input}) & {DSSP} --output-format dssp {input} {output}"

rule score_splits:
    input:
        expand(os.path.join(WORK_DIR, f"tmscore.{ALPHA_SAFE_S}.d{DELTA_SAFE}", "{id}"), id=CLUSTER_IDS)

rule _score_splits:
    input:
        directory(os.path.join(WORK_DIR, f"split.{ALPHA_SAFE_S}.d{DELTA_SAFE}", "{id}"))
    output:
        directory(os.path.join(WORK_DIR, f"tmscore.{ALPHA_SAFE_S}.d{DELTA_SAFE}", "{id}"))
    shell:
        """
        for dir in {input}/*/safe/; do
            LEN=$(expr $(ls $dir | wc -l) / 2)
            SUB=$(echo $dir | cut -d'/' -f5-)
            mkdir -p {output}/"$SUB"
            for ((i=0; i < $LEN; i++)); do
               {TMSCORE_PATH} "$dir"ref_"$i".pdb "$dir"seq_"$i".pdb -outfmt 2 > {output}/"$SUB"tmscore_"$i".out
            done
        done
        for dir in {input}/*/unsafe/; do
            LEN=$(expr $(ls $dir | wc -l) / 2)
            SUB=$(echo $dir | cut -d'/' -f5-)
            mkdir -p {output}/"$SUB"
            for ((i=0; i < $LEN; i++)); do
               {TMSCORE_PATH} "$dir"ref_"$i".pdb "$dir"seq_"$i".pdb -outfmt 2 > {output}/"$SUB"tmscore_"$i".out
            done
        done
        """    

rule pdb_splits:
    input:
        expand(os.path.join(WORK_DIR, f"split.{ALPHA_SAFE_S}.d{DELTA_SAFE}", "{id}"), id=CLUSTER_IDS)

rule _pdb_splits:
    input:
        pdb_dir = directory(os.path.join(WORK_DIR, "pdb")),
        safe = os.path.join(WORK_DIR, f"safety.{ALPHA_SAFE_S}.d{DELTA_SAFE}", "{id}.out")
    output:
        directory(os.path.join(WORK_DIR, f"split.{ALPHA_SAFE_S}.d{DELTA_SAFE}", "{id}"))
    shell:
        "python3 scripts/split_safes.py {input} {output}"


rule dssp:
    input:
        expand(os.path.join(WORK_DIR, "dssp", "{sequence}.out"), sequence=get_sequences())

rule _dssp:
    input:
        lambda wildcards: os.path.join(PDB_DIR, "AF-" + wildcards.sequence.split("/")[-1] + "-F1-model_v3.pdb.gz")
    output:
        os.path.join(WORK_DIR, "dssp", "{sequence}.out") 
    shell:
        "mkdir -p $(dirname {input}) & {DSSP} --output-format dssp {input} {output}"
        


rule sample_x:
    input:
        ancient(expand(os.path.join(WORK_DIR, f"X.{ALPHA_SAFE_S}.d{DELTA_SAFE}", "{id}.fasta"), id=CLUSTER_IDS))

rule _sample_x:
    input:
        ancient(os.path.join(WORK_DIR, f"safety.{ALPHA_SAFE_S}.d{DELTA_SAFE}", "{id}.out"))
    output:
        os.path.join(WORK_DIR, f"X.{ALPHA_SAFE_S}.d{DELTA_SAFE}", "{id}.fasta")
    threads: workflow.cores
    resources:
        logdir = LOG_DIR,
        mem_mb = 512,
        time = "00:05:00",
        queue = "short"
    shell:
        "python3 scripts/safety_sampler.py {input} {output} -x"

# search sequences agaianst a pfam domain database
rule hmmscan_align:
    input:
        ancient(expand(os.path.join(WORK_DIR, "hmmscan_align", "{id}.out"), id=CLUSTER_IDS))

rule _hmmscan_align:
    input:
        FAMILY_PROFILE + ".h3f",
        FAMILY_PROFILE + ".h3i",
        FAMILY_PROFILE + ".h3m",
        FAMILY_PROFILE + ".h3p",
        profile = FAMILY_PROFILE,
        fasta = ancient(os.path.join(WORK_DIR, "clean", "{id}.clean.fasta"))
    output:
        os.path.join(WORK_DIR, "hmmscan_align", "{id}.out")
    threads: workflow.cores
    resources:
        logdir = LOG_DIR,
        mem_mb = lambda wildcards, attempt, input: 512 * attempt,
        time = lambda wildcards, attempt, input: f"00:{10*attempt:02}:00",
        queue = "short"
    shell:
        "{HMMER_PATH}/src/hmmscan -o {output} {input.profile} {input.fasta}"

rule build_pfam:
    input:
        db = FAMILY_DB,
    output:
        profile = FAMILY_PROFILE
    threads: 8
    resources:
        logdir = LOG_DIR,
        cpus = 8,
        mem_mb = 1024 * 8,
        time = "01:00:00",
        queue = "short",
    shell:
        "{HMMER_PATH}/src/hmmbuild --cpu {threads} {output.profile} {input.db}"


rule profile_pfam:
    input:
        profile = FAMILY_PROFILE
    output:
        FAMILY_PROFILE + ".h3f",
        FAMILY_PROFILE + ".h3i",
        FAMILY_PROFILE + ".h3m",
        FAMILY_PROFILE + ".h3p"
    threads: workflow.cores
    resources:
        mem_mb = 1024 * 8,
        time = "01:00:00",
        queue = "short",
    shell:
        "{HMMER_PATH}/src/hmmpress {input.profile}"

rule jackhmmer:
    input:
        expand(os.path.join(WORK_DIR, "jackhmmer", "{id}.out"), id=CLUSTER_IDS)
    
rule _jackhmmer:
    input:
        fasta = ancient(os.path.join(WORK_DIR, "clean", "{id}.clean.fasta")),
        db = DB_FILE
    output:
        os.path.join(WORK_DIR, "jackhmmer", "{id}.out")
    threads: 1
    resources:
        logdir = LOG_DIR,
        mem_mb = 256,
        time = "00:15:00",
        cpus = 1,
        queue = "short"
    shell:
        "{HMMER_PATH}/src/jackhmmer --tblout {output} --noali {input.fasta} {input.db} > /dev/null"

# phmmer for all clusters
rule phmmer:
    input:
        expand(os.path.join(WORK_DIR, "phmmer", "{id}.out"), id=CLUSTER_IDS)

rule _phmmer:
    input:
        ref = ancient(os.path.join(WORK_DIR, "refs", "{id}.ref.fasta")),
        fasta = ancient(os.path.join(WORK_DIR, "clean", "{id}.clean.fasta")),
    output:
        os.path.join(WORK_DIR, "phmmer", "{id}.out")
    threads: workflow.cores
    resources:
        logdir = LOG_DIR,
        mem_mb = 256,
        time = "00:05:00",
        cpus = 1,
        queue = "short"
    shell:
        "{HMMER_PATH}/src/phmmer --tblout {output} --noali {input.ref} {input.fasta} > /dev/null"

# hmmsearch for all clusters
rule hmmsearch:
    input:
        ancient(expand(os.path.join(WORK_DIR, "hmmsearch", "{id}.out"), id=CLUSTER_IDS))

rule _hmmsearch:
    input:
        fasta = ancient(os.path.join(WORK_DIR, "clean", "{id}.clean.fasta")),
        hmm = ancient(os.path.join(WORK_DIR, "profiles", "{id}.hmm"))
    output:
        os.path.join(WORK_DIR, "hmmsearch", "{id}.out")
    threads: workflow.cores
    resources:
        logdir = LOG_DIR,
        mem_mb = 512,
        time = "00:10:00",
        cpus = 1,
        queue = "short"
    shell:
        "{HMMER_PATH}/src/hmmsearch --tblout {output} --noali {input.hmm} {input.fasta} > /dev/null"

# creates needed hmm-profiles for hmmsearch
rule hmmbuild:
    input:
        expand(os.path.join(WORK_DIR, "profiles", "{id}.hmm"), id=CLUSTER_IDS)

rule:
    input:
        ancient(os.path.join(WORK_DIR, "msa", "{id}.fasta"))
    output:
        temp(os.path.join(WORK_DIR, "profiles", "{id}.hmm"))
    threads: workflow.cores
    resources:
        logdir = LOG_DIR,
        mem_mb = 512,
        time = "00:10:00",
        cpus = 1,
        queue = "short"
    shell:
        "{HMMER_PATH}/src/hmmbuild {output} {input}"

# MSA for all clusters with Muscle
rule msa:
    input:
        ancient(expand(os.path.join(WORK_DIR, "msa", "{id}.fasta"), id=CLUSTER_IDS))

rule _msa:
    input:
        ancient(os.path.join(WORK_DIR, "clean", "{id}.clean.fasta"))
    output:
        os.path.join(WORK_DIR, "msa", "{id}.fasta")
    threads: workflow.cores
    resources:
        logdir = LOG_DIR,
        mem_mb = 4024,
        time = "03:30:00",
        cpus = 1,
        queue = "short"
    shell:
        "{MUSCLE_PATH} -in {input} -out {output} -quiet"

# Separates all clusters (that fit the threshold) into separate fasta-files
checkpoint separate_clusters:
    input:
        f"{DB_FILE}",
        os.path.join(WORK_DIR, f"{DB_FILENAME}.{CLUSTER_ALGO}")
    output:
        os.path.join(WORK_DIR, "info.txt"),
        directory(os.path.join(WORK_DIR, "clean")),
        directory(os.path.join(WORK_DIR, "fasta")),
        directory(os.path.join(WORK_DIR, "refs"))
    params:
        n = config["cluster_number"],
        min_size = config["cluster_min_size"],
        max_size = config["cluster_max_size"]
    threads: workflow.cores
    resources:
        logdir = LOG_DIR,
        mem_mb = 1024,
        time = "00:30:00",
        queue = "short"
    shell:
        "python3 scripts/clusteread.py a {input} --n {params.n} --min {params.min_size} --max {params.max_size} {UNIFORM}"

CLUSTER_REQS = ""
if CLUSTER_ALGO == "mcl" or CLUSTER_ALGO == "multi-step":
    CLUSTER_REQS = os.path.join(DATA_DIR, f"{DB_FILENAME}.dmnd")
elif CLUSTER_ALGO == "mmseqs":
    CLUSTER_REQS = os.path.join(DATA_DIR, "mmseqs", f"{DB_FILENAME}")

# Generates diamond binary database for clustering
rule cluster_reqs:
    input:
        f"{DB_FILE}"
    output:
        CLUSTER_REQS
    resources:
        logdir = LOG_DIR,
        mem_mb = 1024,
        time = "00:05:00",
        cpus = 1,
        queue = "short"
    run:
        if CLUSTER_ALGO == "mcl" or CLUSTER_ALGO == "multi-step":
            shell("{DIAMOND_PATH} makedb --in {input} -d {output} --threads {resources.cpus}")
        elif CLUSTER_ALGO == "mmseqs":
            shell("{MMSEQS_PATH} createdb {input} {output}")

# Clustering
rule cluster:
    input:
        CLUSTER_REQS
    output:
        clu = os.path.join(WORK_DIR, f"{DB_FILENAME}.{CLUSTER_ALGO}"),
    params:
        identity=MIN_IDENTITY,
        algo=CLUSTER_ALGO,
        sens=CLUSTER_SENSITIVITY,
    threads: 56
    resources:
        logdir = LOG_DIR,
        mem_mb = 18500,
        time = "02:00:00",
        cpus = 56,
        queue = "short"
    shell:
        "{DIAMOND_PATH} cluster --cluster-algo {params.algo} --id {params.identity} {params.sens} -d {input} -o {output} --threads {resources.cpus} --tmpdir {TEMP_DIR} -v"
