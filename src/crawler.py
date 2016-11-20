#!/usr/local/bin/python3
import openpyxl
from bs4 import BeautifulSoup
import urllib
import urllib3
import os
import json
import re
import sys
from data import *
from url import *
from soup import *
from excel import *

http = urllib3.PoolManager()


def load_fusion(data):
    j = fetch_json(fusion_url, init_param(fusion_params, data.gen, data.gen_id, data.gen_seq))
    for arr in j_data(j):
        link = parse(arr[2]).a
        mutation = link.text
        data.fusions.append(Fusion(mutation, arr[3], int(arr[4]), link['href']))


def load_tissue_hist(data, uat, hist_type, variation_type, positive_name, negative_name):
    hist = uat.find('td', attrs={'data-hist-type': hist_type})
    pos, neg = hist['data-hist-muts'].split(',')
    total = hist['data-hist-tested']
    if total == '-' or total == '':
        return
    if pos != '-' and pos != '':
        data.tissues.append(Tissue(variation_type, positive_name, int(pos), int(total), True))
    if neg != '-' and neg != '':
        data.tissues.append(Tissue(variation_type, negative_name, int(neg), int(total), False))


def load_tissue(data):
    tissue_header = headers.copy()
    tissue_header['Accept'] = 'text/plain, */*; q=0.01'
    h = fetch_html(tissue_url, init_param(tissue_params, data.gen, data.gen_id, data.gen_seq), h=tissue_header)
    uat = h.find('td', attrs={'class': 't_hist', 'data-hist-title': 'Upper aerodigestive tract'}).parent
    load_tissue_hist(data, uat, 'cnv', 'CNV', 'Gain', 'Loss')
    load_tissue_hist(data, uat, 'expr', 'Expression', 'Over', 'Under')
    load_tissue_hist(data, uat, 'meth', 'Methylation', 'Hypo', 'Hyper')


def check_fs_somatic_status(mut_id):
    params = {'id': mut_id}
    h = fetch_html(mut_overview_url, params)
    overview = get_id(h, 'div', 'overview')
    pattern = re.compile('Ever confirmed somatic:')
    confirmed = overview.find('dt', text=pattern).find_next('dd').text
    return confirmed == 'Yes'


def get_pathogenic_score(mut_id):
    params = {'id': mut_id}
    h = fetch_html(mut_overview_url, params)
    overview = get_id(h, 'div', 'overview')
    pattern = re.compile('FATHMM prediction:')
    fathmm_str = overview.find('dt', text=pattern).find_next('dd').text
    return float(re.findall('\d+.\d+', fathmm_str)[0])


def get_snv_conseq(mut_id):
    params = {'id': mut_id}
    h = fetch_html(mut_overview_url, params)
    overview = get_id(h, 'div', 'overview')
    pattern = re.compile('AA Mutation:')
    conseq = overview.find('dt', text=pattern).find_next('dd').text
    if 'Missense' in conseq:
        return 'Missense'
    else:
        return 'Nonsense'


def load_fs_data(data, arr, variation_type, frameshift):
    mut_id = get_param(parse(arr[5]).a['href'], 'id')[0]
    transcript = parse(arr[1]).a.text
    aa_mut = parse(arr[5]).a.text
    cds_mut = parse(arr[6]).a.text
    link = parse(arr[5]).a['href']

    if frameshift:
        if check_fs_somatic_status(mut_id):
            if mut_id in data.frameshifts:
                data.frameshifts[mut_id].incr()
            else:
                references = get_reference(data, mut_id)
                data.frameshifts[mut_id] = Frameshift(mut_id, transcript, variation_type, 'Frameshift',
                                                      aa_mut, cds_mut, link, references)
    else:
        if mut_id in data.snvs:
            data.snvs[mut_id].incr()
        else:
            references = get_reference(data, mut_id)
            data.snvs[mut_id] = SNV(mut_id, transcript, get_snv_conseq(mut_id),
                                    aa_mut, cds_mut, get_pathogenic_score(mut_id), link, references)


def load_frameshift(data, mut, variation_type):
    params = init_param(frameshift_params, data.gen, data.gen_id, data.gen_seq)
    params['mut'] = mut
    j = fetch_json(frameshift_url, params)
    for arr in j_data(j):
        load_fs_data(data, arr, variation_type, True)


def load_snv(data):
    j = fetch_json(snv_url + '?' + urlencode(init_param(snv_params, data.gen, data.gen_id, data.gen_seq)))
    for arr in j_data(j):
        load_fs_data(data, arr, 'SNV', False)


def filter_reference(title):
    return any(x in title.lower() for x in reference_filter)


def get_reference(data, mut_id):
    j = fetch_json(reference_url, init_param(reference_params, None, mut_id, None))
    references = []
    for arr in j_data(j):
        title = parse(arr[0]).span['title']
        if filter_reference(title):
            references.append(title)
    return references


def load_gen_data(data):
    h = fetch_html('http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln=' + data.gen)
    data.gen_id = h.find('input', attrs={'name': 'id'})['value']
    data.gen_seq = h.find('input', attrs={'name': 'end'})['value']


def load_gen(gen):
    print('Loading', gen)
    data = Data(gen)
    load_gen_data(data)
    load_fusion(data)
    load_tissue(data)
    load_frameshift(data, 'deletion_frameshift', 'Deletion')
    load_frameshift(data, 'insertion_frameshift', 'Insertion')
    load_snv(data)
    output(data, '/Users/mu/' + data.gen + '.xlsx')


def filter_gen(arr, gen):
    gens = []
    # exact match
    if any(parse(x[0]).a.text == gen for x in arr):
        return [gen]
    pattern = re.compile(gen + '_ENST.*')
    for g in arr:
        gen_name = parse(g[0]).a.text
        if pattern.match(gen_name):
            gens.append(gen_name)
    return gens


def main(gen):
    search_param['ln'] = gen
    j = fetch_json(search_url, init_param(search_param, gen, None, None))
    if no_result(j):
        print('No result for ', gen)
    else:
        for gen in filter_gen(j_data(j), gen):
            load_gen(gen)


if __name__ == '__main__':
    main(sys.argv[1])