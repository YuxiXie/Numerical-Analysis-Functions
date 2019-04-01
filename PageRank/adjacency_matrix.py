import warnings
warnings.filterwarnings(action='ignore', category=UserWarning, module='gensim')
import random
import numpy as np
from gensim.corpora.wikicorpus import extract_pages, find_interlinks


def build_dict(N):
    tuple = extract_pages("enwiki-20181220-pages-articles-multistream.xml")
    page_dict = {}
    elect_id = random.randint(1, 500)
    id = 0
    cnt = 1
    for t in tuple:
        if (cnt > N):
            break
        id += 1
        if (id == elect_id):
            title = t[0]
            interlinks = find_interlinks(str(t))
            outlink_num = len(interlinks)
            page_dict[title] = [cnt, outlink_num, list(interlinks.keys())]
            cnt += 1
            elect_id += random.randint(1, 150)
    return page_dict


def obtain_A(N):
    A = np.zeros((N, N))
    page_dict = build_dict(N)
    page_name = list(page_dict.keys())
    for i in range(N):
        title = page_name[i]
        outlink_num = page_dict[title][1]
        for outlink in range(outlink_num):
            link_name = page_dict[title][2][outlink]
            if (page_dict.get(link_name) != None):
                j = page_dict[link_name][0]
                A[i, j] = 1
    return A


N = 30000
A = obtain_A(N)
np.save('adjacency_matrix.npy', A)
