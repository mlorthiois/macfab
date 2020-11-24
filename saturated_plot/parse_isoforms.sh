reports=$(find ./*/*/results/SQANTI3/* -name *_classification.txt)

echo "sample,reads,tool,isoforms" > results_bench.csv
for report in $reports; do
    sample_id=$(echo $report | cut -f2 -d '/')
    reads=$(echo $report | cut -f3 -d '/' | head -c1)
    tool=$(echo $report | cut -f6 -d '/')
    isoforms=$(cat $report | wc -l)
    echo "$sample_id,$reads,$tool,$isoforms" >> results_bench.csv
done

echo "CSV file created"