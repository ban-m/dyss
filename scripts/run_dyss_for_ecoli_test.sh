#!/bin/bash
ONT_PYTHON=/Applications/MinKNOW.app/Contents/Resources/ont-python/bin/python
$ONT_PYTHON ./src/dyss_debug.py --min_chunk_size 4000 \
	    --reference ./data/Escherichia_coli_str_k_12_substr_w3110.ASM1024v1.dna.chromosome.Chromosome.fa \
	    --model ./data/template_r9.2.csv \
	    --param ./data/mock_parameters.csv \
	    --power 20 \
	    --test ~/work/test_reads/temp_r9.2/
