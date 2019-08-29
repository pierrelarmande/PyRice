# -*- coding: utf-8 -*-
# from MultiQuery import MultiQuery
# from MultiQuery2 import MultiQuery
from MultiQuery3 import MultiQuery
import time
import re
import csv
from helper import *
from multiprocessing import cpu_count
from multiprocessing_logging import install_mp_handler
import logging
from argparse import ArgumentParser

parser = ArgumentParser("PyRice", conflict_handler='resolve')

parser.add_argument("--number_process", type=int)
parser.add_argument("--case",type=int)
parser.add_argument("--multi_processing",type=bool)
parser.add_argument("--multi_threading",type=bool)
parser.add_argument("--number_gene",type = float)
args = parser.parse_args()

if __name__ == "__main__":
    # Search gene and query loc
    # t = time.time()
    # test = MultiQuery()
    # t = time.time()
    # file_id = test.search_gene(chro="chr01", start_pos="10000",
    #                     end_pos="50000",number_process=cpu_count(),dbs=["iric"],save_path="./result/")
    # print("Time for search gene:",time.time() - t)
    # print("Ouput file",file_id)

    #Query with chromosome
    t = time.time()
    test = MultiQuery()
    db = test.query_iric(chro="chr01", start_pos="1",
                         end_pos=str(int(args.number_gene)), number_process=args.number_process, multi_processing=args.multi_processing, multi_threading=args.multi_threading,
                         dbs="all")
    test_time = time.time() - t
    print("Time for query:", test_time)

    # t = time.time()
    # test = MultiQuery()
    # db = test.query_iric(chro="chr01", start_pos="1",
    #                     end_pos="100000",number_process = cpu_count(),multi_processing=True,multi_threading=True,dbs="all")
    # test_time = time.time() - t
    # print("Time for query:", test_time)
    # test.save_file(db,save_path="./result4_"+str(round(test_time,2))+"/",format=["csv","html","json","pkl"],hyper_link=False)
    # print("Output database",db)

    ##Query with ids, locs and irics
    # idents = []
    # locs = []
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
    #                          ,locs=[],irics=[],number_process = cpu_count()*2,dbs=["oryzabase","rapdb","gramene"])
    # print(time.time()-t)
    # test.save_file(db,'./sdf/',format=["csv","html","json","pkl"],hyper_link=False)
    # print("Output database",db)

    # test = MultiQuery()
    # t = time.time(),
    # db = test.new_query(atts=['TRAES3BF001000010CFD'], dbs=['urgi'])


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
