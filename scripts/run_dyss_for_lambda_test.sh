#!/bin/bash
ONT_PYTHON=/Applications/MinKNOW.app/Contents/Resources/ont-python/bin/python
$ONT_PYTHON ./src/dyss_debug.py --min_chunk_size 4000 \
	    --reference ./data/lambda.fa \
	    --model ./data/template_r9.4.csv \
	    --param ./data/mock_parameters.csv \
	    --test ~/work/test_reads/temp/ 2>&1 | grep score | sort -k 2 -t',' -n > ./result/lambda_test.csv
