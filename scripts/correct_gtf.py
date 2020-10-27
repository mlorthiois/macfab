#! python
import glob
import csv

list_transcripts = {}
new_csv = ""

with open(snakemake.input[0]) as csvfile:
    reader = csv.reader(csvfile, delimiter="\t")
    for row in reader:
        line = ""
        for tag in row[8].split("; "):
            if tag.startswith("transcript_id"):
                transcript_id = tag.split(" ")[1]
                if transcript_id not in list_transcripts:
                    list_transcripts[transcript_id] = row[0] + row[6]
                    line = "\t".join(row) + "\n"
                else:
                    if list_transcripts[transcript_id] == row[0] + row[6]:
                        line += "\t".join(row) + "\n"
                    else:
                        print(f"Remove duplicated line : {row}")
        new_csv += line

with open(snakemake.output[0], "w") as f:
    f.write(new_csv)
