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
