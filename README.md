## Dependencies
- Diamond   https://github.com/bbuchfink/diamond
- Muscle    http://www.drive5.com/muscle/
- ete3      https://github.com/etetoolkit/ete
- Hmmer     https://github.com/EddyRivasLab/hmmer

Compile Diamond from source (version 2.0.8 does not have issues with clustering):
```
wget https://github.com/bbuchfink/diamond/archive/refs/tags/v2.0.8.tar.gz
tar -xvzf v2.0.8.tar.gz
cd diamond-2.0.8
mkdir bin
cd bin
module load GCC CMake       # Only on computer cluster
cmake -DEXTRA=ON ..
make
```
Compile Hmmer and easel from source:
```
git clone https://github.com/EddyRivasLab/hmmer
cd hmmer
git clone https://github.com/EddyRivasLab/easel
module load Autoconf        # Only on computer cluster
autoconf
./configure
make
make check                  # optional: run automated tests
cd easel
make
```
Download muscle:
<http://www.drive5.com/muscle/downloads.htm>
wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
tar -xvzf v2.0.8.tar.gz

ete3 and snakemake via Conda environment
```
cd alignment-safety/pipeline
conda env create -f environment.yaml
conda activate pipeline
```
(optional)
Compile Raxml-ng from source:
```
git clone https://github.com/amkozlov/raxml-ng
cd raxml-ng
mkdir build
cd build
cmake ..
make
```
(optional)
Follow instruction for MMseqs2 installation (compilation from source recommended for better performance):
<https://github.com/soedinglab/MMseqs2#installation>

---

## How to run the pipeline?

#### 1. Edit `parameters.yaml`:
- `db_file`<br/>
    Path to the protein sequence database (fasta format), e.g. `uniprot_sprot.fasta`.
- `family_db`<br/>
    Path to the Pfam protein domain database: `Pfam-A.seed`.
- `diamond`<br/>
    Path to DIAMOND, e.g. `./../diamond-2.08`.
- `raxml`<br/>
    Path to RaxML-ng, e.g.  `./../raxml-ng`.
- `muscle`<br/>
    Path to Muscle - multiple sequence alignment tool, e.g. `./../muscle`.
- `hmmer`<br/>
    Path to the hmmer programs, e.g. `./../hmmer-3.3.2`.
- `safety`<br/>
    Parent folder of the compiled safety-window program, `./../safety-windows`
- `min_identity`<br/>
    Identity threshhold-% minimum to report and alignment. Used in db clustering process. Tested 20%-90%.
- `sensitivity`<br/>
    `auto`, or any of the DIAMOND sensitivity options: <https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#sensitivity-modes>
- `clustering_algorithm`<br/>
    `mcl` (DIAMOND), `multi-step` (DIAMOND) or `mmseqs` (MMSeqs2).
- `clustering_min_size`<br/>
    Treshhold of the minimum cluster size to include.<br/>
    Warning: `< 20`, will produce large amount of files.
- `clustering_max_size`<br/>
    Treshhold of the maximum cluster size to include.
- `cluster_number`<br/>
    If less than available clusters, `cluster_number` amount of clusters will be chosen randomly to include.<br/>
    Speeds up debugging/testing.
- `ref_criterion`<br/>
    * `--clustering` - default, depends on clustering: mcl or multi-step
    * `--identity` - Highest mean pair-wise identity score
    * `--highlow` - Highest lowest pair-wise identity score
    * `--similarity` - Highest hmmsearch score (i.e. most similar sequence to the MSA of the cluster)
    * `--taxonomy` - Highest node in taxonomic tree
#### 2. Run `separate_clusters` snakemake rule:
    snakemake -j n separate_clusters
- `n`<br/>
    Number of parallel processes.
    Clusters the database and separates clusters satisfying the treshholds to `WORK_DIR/`
- `WORK_DIR/fasta/`<br/>
    Each fasta-file corresponds to one cluster. Name of the file as well as the first sequence in the file is the reference sequence of that cluster.
- `WORK_DIR/clean/`<br/>
    Same as `fasta/`, but fasta sequences are cleaned with no additional information, such as taxonomic id. Needed for some programs, such as Muscle.
- `WORK_DIR/refs/`<br/>
    One file containing each cluster's reference sequence in fasta-format. Needed for `phmmer`.
#### 3. Run all rules or some specific rule:
    snakemake -j n rule
- `all`<br/>
    Runs rules: `safe`, `identity`, `hmmsearch`, `hmmscan`, `phmmer` and their prerequisites.
- `safe`<br/>
    Runs Safety-Window-program on all clusters. Outputs to `WORK_DIR/safety/`.
- `identity`<br/>
    Runs `esl-alipid`-program on all clusters to calculate pairwise identities. Outputs to `WORK_DIR/id/`.
- `hmmsearch`<br/>
    Runs `hmmsearch` on all clusters. Outputs to `WORK_DIR/hmmsearch/`.
- `hmmscan`<br/>
    Runs `hmmscan` on all clusters. Outputs to `WORK_DIR/hmmscan/`.
- `phmmer`<br/>
    Runs `phmmer` on all clusters. Outputs to `WORK_DIR/phmmer/`.

## How to configure and run the pipeline on cluster?

#### 1. Install Mamba (Conda):
-   `wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh`
-   `bash Mambaforge-$(uname)-$(uname -m).sh`
-   install into `/proj/<username>/mamba`
#### 2. Download repository and dependencies somewhere into $PROJ on Turso

#### 3. Download Database e.g. swissprot into $WRKDIR

#### 4. Edit `turso/parameters.yaml`:
-   Dependencies should be located in `/proj/<username>/` ($PROJ)
-   Data such as Swissprotein and Pfam DB should be located somewhere in `/wrk-vakka/users/<username>` ($WRKDIR)
-    Work (wrkdir) and temporary (tempdir) directory should be located in `/wrk-vakka/users/<username>` ($WRKDIR)

#### 5. Run the Snakemake pipeline via shell script:
-   `turso/run.sh rule -j <num_of_maximum_parallel_processes>`
-    For clustering -j 56 should conclude in 30-60min
-    Separating and changing reference -j 1-4

#### 6. Follow the progress:
-   To stream the progress `less +F logs/latest/progress.log`
-   `ctrl + c` then `q` exits the stream. Pipeline is still running in background.
-   Logs are found in `logs/<datetime>/<job_id>.out`

#### 7. Useful slurm commands:
-   `slurm w q` to see running jobs
-   `slurm w qq` to see running jobs with resource usage
-   `scancel -M ukko2 <job_id>` to cancel job
-   `seff -M ukko2 <job_id>` to see used resources and run time of job
-   `squeue -o '%A %.28R %j' -u <username>` to see if you have any jobs running

## Issues:

#### It's not possible to run rules, such as, `all`, `msa`, `safe`, `identity`, `hmmsearch`, `hmmscan`, `phmmer` before separating clusters `separate_clusters`
-   This is due to that Snakemake will look for files in `fasta/` and `clean/` to execute these rules. Before separating clusters there are no files in those folders.

#### On turso, it's not possible to use `--taxonomy` as cluster reference out of the box
-   This is due to disk quota limits on Turso. Ete3 tries to download taxonomy database into `~` which is not meant for data storage. This exceeds disk quota limit and is interrupted.
-   Workaround solution to fix this is to run Ete3 on your local machine and copy contents of `~/.etetoolkit` on your local machine into $WKRDIR (e.g. `/wrk-vakka/users/<username>/ncbi`) and make symbolic link `ln -s /wrk-vakka/users/<username>/ncbi ~/.etetoolkit`

#### Problems with Diamond clustering for version > 2.0.8.
-   Use Diamond v2.0.8 for now

## Additional downloads:

#### Swiss-protein database:
-   If working on cluster download and extract somewhere in $WRKDIR
-   `wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz`
-   `gzip -d uniprot_sprot.fasta.gz`
-   Add path to parameters.yaml and/or turso/parameters.yaml

#### Pfam database (for hmmscan):
-   If working on cluster download and extract somewhere in $WRKDIR
-   `wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.seed.gz`
-   `gzip -d Pfam-A.seed.gz`
-   Add path to parameters.yaml and/or turso/parameters.yaml
---
-   University of Helsinki HPC user guide: <https://wiki.helsinki.fi/display/it4sci/HPC+Environment+User+Guide>
-   Hmmer documentation: <http://eddylab.org/software/hmmer/Userguide.pdf>
-   Guide for slurm and other useful documentation (everything might not be applicable for Turso) <https://scicomp.aalto.fi/>