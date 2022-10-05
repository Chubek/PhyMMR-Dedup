from phymmr_dedup import dedup_lines
from time import time

t1 = time()

with open("test_file.fa", "r") as fr:
    lines = fr.readlines()

    len_lines = len(lines)

    dedup_lines(lines, 90, 5)


t2 = time()


print(f"Deduped {len_lines} lines in {t2 - t1} secs")