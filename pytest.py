from phymmr_dedup import dedup_lines
from time import time

with open("SRR6824218.fa", "r") as fr:
    lines = fr.readlines()

    len_lines = len(lines)



headers = len([l for l in lines if l.strip()[0] == ">"])

print("File loaded. Counted headers.")

t1 = time()

dep = dedup_lines(lines, 90)

t2 = time()



print(f"SRR6824218.fa: There are {headers} headers. There are {len(lines)} lines. In {t2 - t1} returned {len(dep)} lines.")



