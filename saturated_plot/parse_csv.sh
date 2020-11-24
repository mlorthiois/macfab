reports=$(find ./*/*/results/SQANTI3/* -name *_classification.txt)

echo "sample,reads,tool,isoforms,genes" > results_bench.csv

for report in $reports; do
    sample_id=$(echo $report | cut -f2 -d '/')
    reads=$(echo $report | cut -f3 -d '/' | head -c1)
    tool=$(echo $report | cut -f6 -d '/')
    isoforms=$(cat $report | wc -l)
    genes=$(cat $report | cut -d $'\t' -f7 | sort | uniq -c | wc -l)
    echo "$sample_id,$reads,$tool,$isoforms,$genes" >> results_bench.csv
done

samples=(CML10 popsi twiny Bear)
tools=(bambu stringtie talon flair)
for sample in ${samples[@]};do
    for name in ${tools[@]}; do
        echo "$sample,0,$name,0,0" >> results_bench.csv 
    done
done

echo "CSV file created"