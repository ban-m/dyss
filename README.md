# Dyss --- a tiny package for selective sequencing on ONT MinION sequencer

## NOTES

If you are interested in the log files of our experiments, please see [this repo](https://bitbucket.org/ban-m/readuntillogfiles). These include IDs of reads marked as "rejected"/"accepted."

To benchmark my implementation by pre-recorded data,please see [this repo](https://bitbucket.org/ban-m/long_reference_speedcheck) for speed checking and [this one](https://bitbucket.org/ban-m/score_calculate) for precision/sensitivity. I recommend using these repositories rather than using read until API with a bulk-fast5 file.

## What this software can/can not? (DISCLAIMER)
- Dyss can select reads by using their first 250-500 "squiggles" rather than sequence comparison after base calling. These squiggles are computed from the raw nanopore signal by segmentation.
- Dyss can take reference up to its size 200Kbp(,which is very short, we admit).
- We intended to run Dyss on relatively "powerless" computers such as laptops/NUC without any additional computational resources(e.g., GPU).
- It is designed for reads(DNA) of low relative abundance. We also note that the total throughput may not be increased. For more details, please see our manuscript(published soon).
- If you are interested in "ordinary" selective sequencing, I think an on-the-fly base-calling approach would perform better. In fact, such workflow has been already bundled with the ReadUntilAPI itself.
- We don't take any responsibility for what you get/lose by using Dyss because this software modifies the behavior of the sequencer irreversibly.
- Currently, we confirm this package only on MacBookPro 13inch(2017) and Intel NUC7i7BNH(Ubuntu 17.10).
- We have tested only on R9.4/2 version flowcells(in other words, 1D protocol). I'm not sure whether or not our scripts runs on R9.5/R9.5.1 flowcell as well. In our opinion, the barcode attached to the head of the reads should be treated correctly.

## How I can use Dyss?
Please install Read Until API first. It is distributed by ONT(presumably on your requests).
During installation of Read Until API, you should modify [the line 408 of base.py](read_until_api/read_until/base.py https://github.com/nanoporetech/read_until_api/blob/release/read_until/base.py#L408)
as follows:
```python
                raw_data_type=self.msgs.GetLiveReadsRequest.UNCALIBRATED,
```
This edit is crucial step because Dyss package use Rust language to convert raw signals into events and Rust is very strongly typed.

Install [mercurial](https://www.mercurial-scm.org/) if it is not installed yet. For example,
```bash
sudo apt install mercurial
```
for Ubuntu and
```bash
brew install mercurial
```
for OSX with homebrew.
  
Install [rust language](https://www.rust-lang.org) if it is not installed yet.
```bash
curl https://sh.rustup.rs -sSf | sh
```

Install python3 if it is not installed yet.
After that, clone this repository and build it.
```bash
hg clone https://ban-m@bitbucket.org/ban-m/dyss
bash setup.sh
```
This doesn't require any root privilege. After execution of the code above,
please check there are compiled library named libdyss.dylib or libdyss.so in ./target/release/ directly.
This is the library to be called by Python during selective sequencing.

To confirm the library was properly compiled, test the library(in dyss directly). Type command below and you should see some output in your screen, such as score/result.
```bash
python3 ./src/dyss_debug.py --reference ./data/lambda.fa --model ./kmer_models/r9.4_180mv_450bps_6mer/template_median68pA.model --param ./data/parameters.csv --test ./data/test_reads/
```

Prep your sample, load it in your MinION sequencer, connect MinION to your PC, start MinKNOW GUI and fill in a few blanks in your MinKNOW application and click execution buttom.

After the Mux process finishes and the sequencing starts, in ./dyss/ directly, do
```bash
ONTPYTHON ./src/signal_based_analysis.py --reference [ref] --model [model] --param [param] --run-time [time(sec)]
```
where ONTPYTHON is the python binary bundled with MinKNOW( for Mac, /Applications/MinKNOW.app/Contents/Resources/ont-python/bin/python).

It will start selective sequencing targeting [ref]:0-200K region for [time] seconds to the pore with odd number. If you would like to apply selective sequencing to all the pores, add `--control-group 1024`. Please make sure that the script will not be interrupted during execution.
[ref] should be a fasta file. Model should be a parameter for k-mer hidden Markov Model(for example, ./data/template_r9.4.csv).
The parameter should be ./data/parameters.csv. If you want to test the selective sequencing, accept every single read, then this parameter should be ./data/mock_parameters.csv.

## Contact
- Mail: banmasutani@gmail.com

