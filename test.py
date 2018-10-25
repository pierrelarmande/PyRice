# #!/usr/bin/env python3

import query

print("Oryzabase Query: ")
for iter in query.query("oryzabase", ["Os03g0149100"]):
    print(iter)
# print("\nRapDB Query")
# for iter in query.query("rapdb", ["Os01g0102700"]):
#     print(iter)
# print("\nGramene Query")
# for iter in query.query("Gramene", ["Os03g0149100"]):
#     print(iter)
# # print("\nic4r Query")
# # for iter in query.query("ic4r", ["Os03g0149100"]):
# #     print(iter)
# print("\nplntfdb Query")
# for iter in query.query("plntfdb", ["321718"]):
#     print(iter)
# # print("\nSNP-Seek Query")
# # for iter in query.query("snpseek", ["chr01", "1", "43270923", "rap"]):
# #     print(iter)
# print("\nfunricegene Query")
# for iter in query.query("funricegene_genekeywords", ["","LOC_Os07g39750"]):
#     print(iter)
# print("\nfunricegene Query")
# for iter in query.query("funricegene_geneinfo", ["","LOC_Os07g39750"]):
#     print(iter)
# print("\nfunricegene Query")
# for iter in query.query("funricegene_faminfo", ["","LOC_Os07g39750"]):
#     print(iter)
# print("\nMSU Query")
# for iter in query.query("msu", ["LOC_Os10g01006"]):
#     print(iter)


 # import os
 # import sys
 # sys.path.append('inst/python')
 # import helper
 # import query
 # import csv
 # from multiprocessing import Pool
 # from multiprocessing.dummy import Pool as ThreadPool
 #
 # # Configuration:
 # num_of_thread = 2
 # path = '/storage/bio/coffee/input/'
 # set = ['G1_clustered_genes_annotations.tab','G2_clustered_genes_annotations.tab','G5_clustered_genes_annotations.tab','G6_ clustered_genes_annotations.tab']
 # #set = ['Control-1-NT_25628_clustered_genes_annotations.tab', 'Curcuma-longa-N0-L_25630_clustered_genes_annotations.tab',' Curcuma-longa-N350-H_25629_clustered_genes_annotations.tab']
 #
 # # Support functions:
 # def get_kegg_id(key):
 #     try:
 #        return {'KEGG_id': key, 'Definition':query.query('kegg', [key])[0]['Definition']}
 #     except:
 #             print(key)
 #             return None
 # # Get a list of unique gene
 # unique_gene_list = []
 # for data in set:
 #     with open(path + data, newline='') as input_file:
 #         dictReader = csv.DictReader(input_file, delimiter = "\t")
 #         for line in dictReader:
 #             for kegg_id in line['KEGG_ko'].split('|'):
 #                 if not(kegg_id in unique_gene_list): unique_gene_list.append(kegg_id)
 #             for kegg_id in line['KEGG_module'].split('|'):
 #                 if not(kegg_id in unique_gene_list): unique_gene_list.append(kegg_id)
 #             for kegg_id in line['KEGG_pathway'].split('|'):
 #                 if not(kegg_id in unique_gene_list): unique_gene_list.append(kegg_id)
 #
 # print('Querying')
 # # Query from unique list
 # # Using multiple (4) threads
 # pool = ThreadPool(num_of_thread)
 # unique_return = pool.map(get_kegg_id, unique_gene_list)