#!/bin/bash
ONT_PYTHON=/Applications/MinKNOW.app/Contents/Resources/ont-python/bin/python
${ONT_PYTHON} ./src/signal_based_analysis.py --reference ./data/Escherichia_coli_str_k_12_substr_w3110.ASM1024v1.dna.chromosome.Chromosome.fa \
	      --run_time 900 \
	      --model ./data/template_r9.4.csv\
	      --param ./data/parameters.csv\
	      --verbose 2>> ./result/2018-06-28-ecoli-900s.log 1>> ./result/2018-06-28-ecoli-900s.out