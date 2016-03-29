75271874    5b  64  61  74  61  5d  0a  6e  0e  74  02  75  30  02  19  00
           [   d   a   t   a   ]  nl   n  so   t stx   u   0 stx  em nul

75271874 + 7 = 75271881


75271881 + 1000 * 17

split -b 75288881 tweets



75271881 + 100000 * 17
split -b 76971881 tweets_27M.nci




75271881 + 1000000 * 17
split -b 92271881 tweets_27M.nci



# 
    time besort -i tweets_27M.nci -s 17 -f 9 -k 4 -h 75271881 -t 4 -b 128m -m 2 -o tweets_27M_sorted.nci

    besort  v.0.1rc1
    chunks in level 0: 4
    chunks in level 1: 2
    chunks in level 2: 1
    sorting leaf chunk: /tmp/__besort_chunk_l0_i0.txt
    sorting leaf chunk: /tmp/__besort_chunk_l0_i1.txt
    sorting leaf chunk: /tmp/__besort_chunk_l0_i2.txt
    sorting leaf chunk: /tmp/__besort_chunk_l0_i3.txt
    finished sorting: /tmp/__besort_chunk_l0_i3.txt
    finished sorting: /tmp/__besort_chunk_l0_i0.txt
    finished sorting: /tmp/__besort_chunk_l0_i2.txt
    merging sorted chunks: /tmp/__besort_chunk_l1_i2.txt
    finished sorting: /tmp/__besort_chunk_l0_i1.txt
    merging sorted chunks: /tmp/__besort_chunk_l1_i0.txt
    finished merging: /tmp/__besort_chunk_l1_i2.txt
    finished merging: /tmp/__besort_chunk_l1_i0.txt
    merging sorted chunks: /tmp/__besort_chunk_l2_i0.txt
    finished merging: /tmp/__besort_chunk_l2_i0.txt
    besort -i tweets_27M.nci -s 17 -f 9 -k 4 -h 75271881 -t 4 -b 128m -m 2 -o   26.81s user 18.66s system 130% cpu 34.771 total


    nanocube-topk -d tweets_27M_sorted.nci -o tweets_27M_sorted.nc -f 100000
    nanocube-topk -d tweets_27M.nci -o tweets_27M.nc -f 100000

    data git:(master) time nanocube-topk -d tweets_27M.nci -o tweets_27M.nc -f 100000
    Nanocube Top-K version: 4.0.1r8
    (Info) loading header...
    (stdin     ) count:     100000 mem. res:        330MB. time(s):         10
    (stdin     ) count:     200000 mem. res:        439MB. time(s):         25
    (stdin     ) count:     300000 mem. res:        486MB. time(s):         38
    (stdin     ) count:     400000 mem. res:        593MB. time(s):         53
    (stdin     ) count:     500000 mem. res:        720MB. time(s):         71
    (stdin     ) count:     600000 mem. res:        874MB. time(s):         91
    (stdin     ) count:     700000 mem. res:       1002MB. time(s):        112
    (stdin     ) count:     800000 mem. res:       1133MB. time(s):        131
    (stdin     ) count:     900000 mem. res:       1273MB. time(s):        147
    (stdin     ) count:    1000000 mem. res:       1294MB. time(s):        167
    (stdin     ) count:    1100000 mem. res:       1094MB. time(s):        190

    (stdin     ) count:   27490000 mem. res:      47021MB. time(s):       4624
    (stdin:done) count:   27494686 mem. res:      47035MB. time(s):       4625 <--------- much faster!!!
    Nanocube with 27494686 records was saved into tweets_27M_sorted.nc

# Flickr sort by ID

    llins@lion5:/n/hadoop_llins/itopk$ grep -b "\[data\]" flickr_84M.nci
    35667607:[data]
    llins@lion5:/n/hadoop_llins/itopk$ grep -n "\[data\]" flickr_84M.nci
    1565120:[data]
    llins@lion5:/n/hadoop_llins/itopk$ od -A d -t x1 -t a -j 35667607 -N 100 flickr_84M.nci
    35667607  5b  64  61  74  61  5d  0a  2d  f5  50  03  ff  03  02  00  00
              [   d   a   t   a   ]  nl   -   u   P etx del etx stx nul nul
    35667623  c8  02  00  00  01  00  00  00  2d  f5  50  03  ff  03  02  00
              H stx nul nul soh nul nul nul   -   u   P etx del etx stx nul
    35667639  00  51  10  00  00  01  00  00  00  2d  f5  50  03  ff  03  02

    # add 7 to the number
    35667607 + 7 = 35667614

    time besort -i flickr_84M.nci -s 17 -f 9 -k 4 -h 35667614 -t 10 -b 128m -m 2 -o flickr_84M_sorted.nci

# Wikipedia

    grep -b "\[data\]" wikipedia-2005-2014_112M.nci
    85147577:[data]

    85147577 + 7 = 85147584

    time besort -i wikipedia-2005-2014_112M.nci -s 17 -f 9 -k 4 -h 85147584 -t 10 -b 128m -m 2 -o wikipedia_112M_sorted.nci


    finished merging: /tmp/__besort_chunk_l3_i0.txt
    merging sorted chunks: /tmp/__besort_chunk_l4_i0.txt
    finished merging: /tmp/__besort_chunk_l4_i0.txt
    besort -i wikipedia-2005-2014_112M.nci -s 17 -f 9 -k 4 -h 85147584 -t 10 -b    85.89s user 188.97s system 109% cpu 4:10.09 total




    # sorted!
    llins@nano2:~/data/topk$ screen -r tw
    (stdin     ) count:   27190000 mem. res:      46205MB. time(s):       3856
    (stdin     ) count:   27200000 mem. res:      46229MB. time(s):       3858
    (stdin     ) count:   27210000 mem. res:      46253MB. time(s):       3859

    # unsorted!
    (stdin     ) count:   27480000 mem. res:      46305MB. time(s):      10776
    (stdin     ) count:   27490000 mem. res:      46319MB. time(s):      10785
    (stdin:done) count:   27494686 mem. res:      46328MB. time(s):      10790
    Nanocube with 27494686 records was saved into tweets_27M.nc


## github

    (stdin     ) count:   18310000 mem. res:      13804MB. time(s):       8447
    (stdin     ) count:   18320000 mem. res:      13810MB. time(s):       8452
    (stdin     ) count:   18330000 mem. res:      13936MB. time(s):       8458
    (stdin:done) count:   18335154 mem. res:      13950MB. time(s):       8462
    Nanocube with 18335154 records was saved into github_84M.nc