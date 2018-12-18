# -*- coding: utf-8 -*-
from MultiQuery import MultiQuery
import time

if __name__ == "__main__":
    #Check exist of gene
    # set_ids, set_locs, idents, locs = MultiQuery().check_gene(idents=["Os08g0164400", "Os07g0586200"],locs=["LOC_Os10g01006", "LOC_Os07g39750","LOC_Os10g13914"])
    # print("Set of ids {} \nSet of locs {}".format(set_ids,set_locs))

    #Query all databse
    # t = time.time()
    # multi_query = MultiQuery()
    # test = multi_query.query_all(idents=["Os08g0164400", "Os07g0586200","Os01g0100900"]
    #                                          ,locs=["LOC_Os10g01006", "LOC_Os07g39750","LOC_Os10g13914","LOC_Os01g01019"])
    # print(time.time() - t)
    #save_file
    #multi_query.save_file("./result/")

    # for i in test.keys():
    #     print("Gene:", i)
    #     for k,v in test[i].items():
    #         print(k,v)
    #     print("\n------------------------")

    ##Search gene and query
    t = time.time()
    test = MultiQuery().search_gene(chro="chr01", start_pos="1",
                                    end_pos="300000",dbs="all")
    print("Time for search gene:",time.time() - t)
    for i in test.keys():
        print("Database:", i)
        for k,v in test[i].items():
            print(k,v)
        print("\n------------------------")

    ##Get id,loc from snpseek database
    id1 = set()
    loc1 = set()
    for key,value in test.items():
        for db,list_gene in value.items():
            if db == 'msu7':
                id_att = "raprepName"
                loc_att = "msu7Name"
            elif db == 'rap':
                id_att = "uniquename"
                loc_att = "msu7Name"
            else:
                id_att = "raprepName"
                loc_att = "msu7Name"
            for gene in list_gene:
                if gene[id_att] !=None:
                #    print(gene[id_att],type(gene[id_att]))
                    for g in gene[id_att].split(','):
                        id1.add(g)
                if gene[loc_att] != None:
                #    print(gene[loc_att],type(gene[loc_att]))
                    for l in gene[loc_att].split(','):
                        loc1.add(l)
    print("List_id {} \nList_loc {}".format(id1,loc1))
    ##Query with id and loc
    t = time.time()
    multi_query = MultiQuery()
    test = multi_query.query_all(idents=id1 ,locs=loc1)
    multi_query.save_file("./result/")
    print("Time for query:",time.time() - t)



