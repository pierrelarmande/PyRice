# -*- coding: utf-8 -*-
from MultiQuery import MultiQuery
import time

if __name__ == "__main__":
    #Search gene and queryloc
    # t = time.time()
    # test = MultiQuery()
    # file_id = test.search_gene(chro="chr01", start_pos="1",
    #                         end_pos="300000",dbs=["iric"],save_path="./result/")
    # print("Time for search gene:",time.time() - t)
    # print("Ouput file",file_id)
    # #Query iric name
    # t = time.time()
    # db = test.query_iric(file_id,dbs="all",save_path="./result/")
    # print("Time for query:", time.time() - t)
    # print("Output database",db)


    # a,b,c,d = multi_query.check_gene(idents=["Os08g0164400", "Os07g0586200","Os01g0100900"]
    #                                          ,locs=["LOC_Os10g01006", "LOC_Os07g39750","LOC_Os10g13914","LOC_Os01g01019"])
    # print("a: ",a)
    # print("b: ",b)
    # print("c: ",c)
    # print("d: ",d)

    #Query with only id and
    multi_query = MultiQuery()
    test = multi_query.query_ids_locs(idents=["Os08g0164400", "Os07g0586200","Os01g0100900"]
                                             ,locs=["LOC_Os10g01006", "LOC_Os07g39750","LOC_Os10g13914","LOC_Os01g01019"],save_path="./result2/")
    print(test)
    #test = multi_query.query_ids_locs(idents=[],locs=["LOC_Os01g01019"],dbs=["Gramene","msu"])




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
