# sc-processing-kickstarter
I often get asked to assist with processing single cell data, starting from raw fastq files. The processing sequence is typically the same. In broad strokes this includes: 
- sequence alignment with 10X cellranger ➡️ align_libraries.py
- pre-processing and QC ➡️ process_libraries.py
- integration of all libraries ➡️ integrate_libraries.py
- analysis of the integrated data (only hints here) ➡️ process_all.ipynb

This repo contains code to get started with these tasks.
Here is what you need to do:
1. 👀 Skim over the code. There isn't a lot of it and most of it is self-explanatory. This gives you an idea of what is happening and if you need certain steps that this kickstarter does not include (some examples below).
2. 🔎 Search for "TODO" in the codebase. They tell you what you might need to do to customize the scripts.
3. ⚙️ Follow the steps in process_all.ipynb. Essentially you need to run the scripts in the order defined above.

Note that this isn't a fire and forget kind of pipeline where you just push a button and everything works. This is some code to help get you on the right track but there might be more to fill in and/or update depending on your single cell experiment. Some things that this kickstarter does not do are:
- handling multiplexed samples
- handling  samples that are sequenced as part of multiple libraries (different from multiplexing)
- handling any other single cell modality other than rna (so no handling for proteome, metabolome, immune repertoire data, etc.)
- filtering out ambient rna

Acknowledgements: Credit for some of this code goes to [@E-R](https://github.com/E-R).

Have fun! 😺