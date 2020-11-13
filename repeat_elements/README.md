De novo discovery of repeat elements using Dfam TE Tools v1.2 (Docker image)
```
curl -sSLO https://github.com/Dfam-consortium/TETools/raw/master/dfam-tetools.sh
chmod +x dfam-tetools.sh
singularity run -B /Volumes/userdata/student_users/chunshenlim/,/Volumes/scratch/brownlab/cslim/fungi/ref/ensembl/FFP/annotations/ docker://dfam/tetools:1.2
BuildDatabase -name filename filename.fna
RepeatModeler -database filename -LTRStruct -srand 12345
```

Build profile HMM models and search for repeat elements in individual genomes using dfamscan.pl
```
ls *.stk | sed 's/-families\.stk//' \
| parallel -j16 \
"nl {}-families.stk | sed 's/ \+//;s/[0-9]\+\t#/#/;s/\t/_/;s|[0-9]\+_//|//|' > {}-families.sto;
hmmbuild {}-families.hmm {}-families.sto;
hmmpress {}-families.hmm;
perl ~/fungi/script/dfamscan.pl \
-T 10 --cpu 12 \
-fastafile {}.fna \
-hmmfile {}-families.hmm \
-dfam_outfile {}.rm.out"
```
