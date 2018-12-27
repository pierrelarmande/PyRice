# -*- coding: utf-8 -*-
from MultiQuery import MultiQuery
import time

if __name__ == "__main__":
    #Search gene and query
    t = time.time()
    test = MultiQuery()
    file_id = test.search_gene(chro="chr01", start_pos="1",
                            end_pos="300000",dbs="all",save_path="./result/")
    print("Time for search gene:",time.time() - t)
    print("Ouput file",file_id)
    #Query iric name
    t = time.time()
    db = test.query_iric(file_id,dbs="all",save_path="./result/")
    print("Time for query:", time.time() - t)
    print("Output database",db)



