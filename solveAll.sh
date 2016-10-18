#!/bin/bash
for filename in ./input/*.txt; do
	./Coiote_heuristic -i $filename -o ./output/summary.csv
done