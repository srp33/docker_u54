import sys

vcf_file_path = sys.argv[1]
dict_file_path = sys.argv[2]
out_vcf_file_path = sys.argv[3]

with open(dict_file_path) as dict_file:
    sequences = set()
    for line in dict_file:
        if not "SN:" in line:
            continue

        line_items = line.rstrip("\n").split("\t")
        sequences.add(line_items[1].replace("SN:", ""))

with open(vcf_file_path) as vcf_file:
    with open(out_vcf_file_path, 'w') as out_vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                if not line.startswith("##contig="):
                    out_vcf_file.write(line)
                continue

            if line.split("\t")[0] in sequences:
                out_vcf_file.write(line)
