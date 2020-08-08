import glob
import re

filenames=" "
filenames=filenames.join(list(snakemake.input))

print('Searching files matching arg "' + filenames + '"...')
for type_value in "Values", "Sensitivity":
    results_file=open(type_value + ".gffparse.tsv", "w") # path + type + extension
    total_line=dict() # dict of [annot]=type,value;type2,value2;...
    
    if type_value=="Sensitivity":
        results_file.write("Source\tAnnotation\tMetric\tFeature\tValue\n")
        regex="(\w+\s\w+\s*\w*):\s+(\d+\.\d+)\s+\|\s+(\d+\.\d+)" # word:space:word:possible space:possible word:":":spaces:num:.:num:|:spaces:num:
    elif type_value=="Values":
        results_file.write("Source\tAnnotation\tType\tFeature\tValue\n")
        regex="(\w+\s\w+):\s+(\d+)\/(\d+)" # regex : word:space:word:":":spaces:num:/:num
    else:
        print("Wrong keyword. Valid keywords: Sensitivity, Values")

    for stats_filename in list(snakemake.input):
        stats_file=open(stats_filename)
        soft_ref=stats_filename.split("/")[-1].split(".") # 0 = soft 1 = annot
        for group in re.findall(regex, stats_file.read()): # return tuples of label,found,total or label,sensi,precision
            if type_value=="Sensitivity":
                results_file.write(soft_ref[0] + "\t" + soft_ref[1] + "\tSensitivity\t" + group[0] + "\t" + group[1] + "\n")
                results_file.write(soft_ref[0] + "\t" + soft_ref[1] + "\tPrecision\t" + group[0] + "\t" + group[2] + "\n")
            elif type_value=="Values":
                results_file.write(soft_ref[0] + "\t" + soft_ref[1] + "\tRaw\t" + group[0] + "\t" + group[1] + "\n") # file:label:value
                if not soft_ref[1] in total_line.keys():
                    total_line[soft_ref[1]]=set()
                if group[0].split(" ")[0] == "Novel": # first part of the label
                    total_line[soft_ref[1]].add((soft_ref[0],group[0],group[2]))
                else:
                    total_line[soft_ref[1]].add(("Reference",group[0],group[2]))
    if type_value=="Values":
        for key in total_line.keys():
            for element in total_line[key]:
                results_file.write(element[0] + "\t" + key + "\tTotal\t" + element[1] + "\t" + element[2] + "\n") # file:label:value
        
        filenames=" "
        filenames=filenames.join(list(snakemake.output))
        print("Data summarised and saved as " + filenames)
