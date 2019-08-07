# -*- coding: utf-8 -*-
from MultiQuery import MultiQuery,search_text
import time
import re
import csv
from helper import *

if __name__ == "__main__":
    #Search gene and query loc
    # t = time.time()
    # test = MultiQuery()
    # t = time.time()
    # file_id = test.search_gene(chro="chr01", start_pos="1",
    #                     end_pos="30000",dbs=["iric"],save_path="./result/")
    # print("Time for search gene:",time.time() - t)
    # print("Ouput file",file_id)

    #Query with chromosome
    # t = time.time()
    # test = MultiQuery()
    # db = test.query_iric(chro="chr01", start_pos="1",
    #                     end_pos="100000",dbs="all")
    # test.save_file(db,save_path="./result/",format=["csv","html","json","pkl"],hyper_link=False)
    # print("Time for query:", time.time() - t)
    # print("Output database",db)

    ##Query with ids, locs and irics
    idents = []
    locs = []
    # with open('/Users/mac/Downloads/Test4Stef.csv') as mapping:
    #     for line in mapping:
    #         loc, annot = re.split(';', line)
    #         if loc != "":
    #             locs.append(loc)
    # print(len(locs),locs)
    # locs = locs[:0]
    # test = MultiQuery()
    # t = time.time()
    # db = test.query_ids_locs(idents=["Os08g0164400", "Os07g0586200","Os01g0100900","Os01g0311400"],locs=locs,irics=[],dbs=[
    #     "oryzabase", "Gramene", "funricegene_genekeywords", "funricegene_faminfo", "rapdb", "funricegene_geneinfo"])
    # print(time.time()-t)
    # test.save_file(db,'./result4/',format=["csv","html","json"])
    # print("Output database",db)

    #Query with ids, locs and irics
    # test = MultiQuery()
    # t = time.time()
    # db = test.query_ids_locs(idents=["Os04g0557800"]
    #                          ,locs=[],irics=[],dbs=[ "Gramene"])
    # print(time.time()-t)
    # test.save_file(db,'./fix/',format=["csv","html","json","pkl"],hyper_link=False)
    # print("Output database",db)

    test = MultiQuery()
    t = time.time()
    db = test.new_query(atts=['TRAES3BF001000010CFD'], dbs=['urgi'])


# chro="chr01", start_pos="1", end_pos="43270923"
# chro="chr02", start_pos="1", end_pos="35937250"
# chro="chr03", start_pos="1", end_pos="36413819"
# chro="chr04", start_pos="1", end_pos="35502694"
# chro="chr05", start_pos="1", end_pos="29958434"
# chro="chr06", start_pos="1", end_pos="31248787"
# chro="chr07", start_pos="1", end_pos="29697621"
# chro="chr08", start_pos="1", end_pos="28443022"
# chro="chr09", start_pos="1", end_pos="23012720"
# chro="chr10", start_pos="1", end_pos="23207287"
# chro="chr11", start_pos="1", end_pos="29021106"
# chro="chr12", start_pos="1", end_pos="27531856"
